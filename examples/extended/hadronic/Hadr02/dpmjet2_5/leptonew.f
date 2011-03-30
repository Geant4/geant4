
*CMZ :          16/04/97  11.35.54  by  Unknown
*CMZ :  1.02/12 15/01/97  23.18.19  by  P. Zucchelli
*CMZ :  1.02/10 14/01/97  16.20.11  by  P. Zucchelli
*CMZ :  1.02/09 14/01/97  16.02.14  by  P. Zucchelli
*CMZ :  1.02/06 13/01/97  17.29.17  by  P. Zucchelli
*CMZ :  1.02/04 13/01/97  14.53.16  by  P. Zucchelli
*CMZ :  1.02/01 12/01/97  16.45.28  by  J. Brunner
*CMZ :  1.02/00 12/01/97  16.24.18  by  J. Brunner
*CMZ :  1.01/51 04/07/96  10.22.19  by  Piero Zucchelli
*CMZ :  1.01/50 23/05/96  12.34.50  by  Piero Zucchelli
*CMZ :  1.01/49 29/01/96  15.55.13  by  Piero Zucchelli
*CMZ :  1.01/48 15/01/96  15.30.50  by  Piero Zucchelli
*CMZ :  1.01/47 11/01/96  09.41.56  by  Piero Zucchelli
*CMZ :  1.01/46 09/01/96  11.43.51  by  Piero Zucchelli
*CMZ :  1.01/45 08/01/96  14.21.12  by  Piero Zucchelli
*CMZ :  1.01/44 05/01/96  18.02.33  by  Piero Zucchelli
*CMZ :  1.01/42 15/12/95 BA 14.58.49  by  Piero Zucchelli
*CMZ :  1.01/41 15/12/95  09.27.14  by  Piero Zucchelli
*CMZ :  1.01/40 12/12/95  15.36.35  by  Piero Zucchelli
*CMZ :  1.01/39 02/11/95  18.49.42  by  Piero Zucchelli
*CMZ :  1.01/38 18/10/95  18.46.18  by  Piero Zucchelli
*CMZ :  1.01/37 18/10/95  18.17.36  BY  PIERO ZUCCHELLI
*CMZ :  1.01/36 31/07/95  18.02.18  BY  PIERO ZUCCHELLI
*CMZ :  1.01/35 26/07/95  14.56.50  BY  PIERO ZUCCHELLI
*CMZ :  1.01/34 25/07/95  11.29.30  BY  PIERO ZUCCHELLI
*CMZ :  1.01/33 12/07/95  09.40.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/32 02/06/95  20.27.59  BY  PIERO ZUCCHELLI
*CMZ :  1.01/31 02/06/95  20.17.58  BY  PIERO ZUCCHELLI
*CMZ :  1.01/29 02/06/95  19.47.43  BY  PIERO ZUCCHELLI
*CMZ :  1.01/27 02/06/95  15.00.58  BY  PIERO ZUCCHELLI
*CMZ :  1.01/24 29/05/95  15.39.50  BY  PIERO ZUCCHELLI
*CMZ :  1.01/23 29/05/95  15.31.35  BY  PIERO ZUCCHELLI
*CMZ :  1.01/22 27/05/95  16.17.50  BY  PIERO ZUCCHELLI
*CMZ :  1.01/21 27/05/95  15.46.09  BY  PIERO ZUCCHELLI
*CMZ :  1.01/20 27/05/95  15.12.35  BY  PIERO ZUCCHELLI
*CMZ :  1.01/18 14/05/95  12.43.38  BY  PIERO ZUCCHELLI
*CMZ :  1.01/15 14/05/95  11.39.34  BY  PIERO ZUCCHELLI
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/13 14/05/95  11.21.51  BY  PIERO ZUCCHELLI
*CMZ :  1.01/11 13/05/95  18.30.15  BY  PIERO ZUCCHELLI
*CMZ :  1.01/10 28/04/95  14.19.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  19.33.11  BY  PIERO ZUCCHELLI
*CMZ :  1.01/06 05/03/95  10.42.02  BY  PIERO ZUCCHELLI
*CMZ :  1.01/05 05/03/95  01.15.37  BY  PIERO ZUCCHELLI
*CMZ :  1.01/04 05/03/95  01.00.25  BY  PIERO ZUCCHELLI
*CMZ :  1.01/03 05/03/95  00.09.33  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 23/09/94  13.01.24  BY  PIERO ZUCCHELLI
*CMZ :  1.01/00 08/09/94  09.46.35  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 22/08/94  14.08.45  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 20/07/94  12.28.44  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C     PROGRAM MAIN
      SUBROUTINE JETTARUN
*     IMPLICIT NONE
      PARAMETER(MEMOR=8000000)
      PARAMETER (NMXHEP=2000)
      PARAMETER (LUNBEAM=10)

*        INTEGER*4 MEMOR
      INTEGER*4 IBADEV,LEPIN,INTERACTION,MAXIEVT
      REAL*4 PPXYZ(3),CI,VDECY(4),BB(3),BRF(3),XM(3),PM(3),PH(3)
      COMMON /TAUPOS / NPA,NPB
      COMMON/BERI/JALLY,JEIN
      COMMON /SLATE/ISL(40)
      COMMON/RUNCOM/IMODE
      COMMON/SBEAM/ PNUMBER,NEUTYPE,VECT(3),GKIN(3),MESTYPE,G4MES(4)
      INTEGER JALLY(30),JARRY(30),IRAWHEAD(11),KEYLIST(50)
      DOUBLE PRECISION HH(4)
      INTEGER*4 DALUAEF(200),IPROT
      COMMON/DALUA/DALUAEF
*KEEP,LUDAT3.
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      SAVE /LUDAT3/
*KEEP,POLAR.
C--
      COMMON /POLARIZ/POL(4000,3)
      REAL POLARX(4)
*KEEP,LUDAT1.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,JETTA.
C--
         PARAMETER (ICENTO=100)
         PARAMETER (ISIZ=93)
         PARAMETER (IOF1=32)
         PARAMETER (IOF2=83)
         PARAMETER (LUX_LEVEL=4)
         INTEGER*4 JTAU(100),JPRI(100),JSTRO(100)
         REAL*4 FTUPLE(ISIZ)
         COMMON/JETTAGL/JTAU,JPRI,JSTRO
         COMMON/NTUPLA/FTUPLE,ISFIRST
         COMMON/BEAM/SPEC(ICENTO)
         COMMON /MAXSPEC/RMAXSPEC,RINTSPEC
         COMMON/SAV/XMINSAV(ICENTO),XMAXSAV(ICENTO),YMINSAV(ICENTO),
     &   YMAXSAV(ICENTO),Q2MINSAV(ICENTO),Q2MAXSAV(ICENTO),
     &   W2MINSAV(ICENTO),W2MAXSAV(ICENTO),PARIMAX(ICENTO),
     &   PPSAVE(ICENTO,3,4,5),PARICOR(ICENTO),INDEX,SIGMASAV(ICENTO),
     &   XMSIGMA,XSECT

*KEEP,FOREFI.
C--
	INTEGER*4 IEVT
         COMMON/FOREFICASS/IEVT


*KEEP,HEPEVT.
      DOUBLE PRECISION PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      SAVE /HEPEVT/

*KEEP,KEYS.
        COMMON/CFREAD/SPACE(5000)
 	COMMON/MYKEYS/IEVAR,IF5CC,INEUT,IINTE,IFERM,IFLAT,ICOUN,
     &  REFIX,IMUDO,IKAT1,IKAT2,IKAT3,IKAT4,IKAT5,IKAT6,
     &  INEVT,IJAK1,IJAK2,IITDK,RPTAU,RXK0D,NTGR,IDIMUON,IGLU,
     &  IQDEN,LOME(2),IFILES,ICCHA,ISEED,IDSUBS,IONEM,EHAC,IPSEL,
     &  IHIST


*KEEP,PERROR.
*-- AUTHOR :    PIERO ZUCCHELLI   01/09/94
        PARAMETER (CHARMSENS=10000)
 	COMMON/MYERR/ICRACK


*KEEP,zebra.

      PARAMETER (NNQ=1000000)
*
      DIMENSION LQ(NNQ),IQ(NNQ),Q(NNQ)
      EQUIVALENCE (Q(1),IQ(1),LQ(9),JSTRUC(8))
      COMMON /QUEST/IQUEST(100)
      COMMON /XQSTOR/IXEVT,IFENCE(16),JGEEV,JSTRUC(99),JREFER(100),
     +DIV12(NNQ)
      COMMON /FZLUN/LUNFZ
      COMMON/MZIOALL/IOGENF

*KEND.



      COMMON/WLIST/WW1,WW2,WW3,WW5
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)

      REAL*4 VNPALIFE(1000)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON/PAWC/H(MEMOR)
      CHARACTER *8 FTAGS(ISIZ)


      CALL READFFKY


* GRV HO/LO SET
      LST(15)=IQDEN
      LST(17)=IEVAR
      LST(32)=IF5CC
      INTERACTION=IINTE
      LST(33)=IFERM
      LST(37)=IFLAT
      MAXIEVT=ICOUN
      ENERGY_FIX=REFIX
      LST(38)=IMUDO
      LEPIN =INEUT
      LST(10)=NTGR
      LST(3)=3

* random number initialization
      CALL RLUXGO(LUX_LEVEL,ISEED,0,0)

* muon generation by charm and nu_e/anu_e
      IF (ABS(LEPIN).EQ.14) THEN
        IMUREQ=1
      ELSE
        IMUREQ=0
      ENDIF

      IF (IDSUBS.EQ.0) THEN
        IF (IDIMUON.EQ.2) THEN
* Dplus
          R1=0
          R1M=0
          DO IO=367,377
            MDME(IO,1)=0
            R1=R1+BRAT(IO)
          ENDDO
          DO IO=378,388
            MDME(IO,1)=1
            R1M=R1M+BRAT(IO)
          ENDDO
          DO IO=389,429
            MDME(IO,1)=0
            R1=R1+BRAT(IO)
          ENDDO

* D0
          R2=0
          R2M=0
          DO IO=430,437
            MDME(IO,1)=0
            R2=R2+BRAT(IO)
          ENDDO
          DO IO=438,445
            MDME(IO,1)=1
            R2M=R2M+BRAT(IO)
          ENDDO
          DO IO=446,490
            MDME(IO,1)=0
            R2=R2+BRAT(IO)
          ENDDO

* D+_s
          R3=0
          R3M=0
          DO IO=491,496
            R3=R3+BRAT(IO)
            MDME(IO,1)=0
          ENDDO
          DO IO=497,501
            R3M=R3M+BRAT(IO)
            MDME(IO,1)=1
          ENDDO
          DO IO=502,523
            R3=R3+BRAT(IO)
            MDME(IO,1)=0
          ENDDO

* Lambda_c
          R4=0
          R4M=0
          DO IO=997,1003
            R4=R4+BRAT(IO)
            MDME(IO,1)=0
          ENDDO
          DO IO=1004,1010
            R4M=R4M+BRAT(IO)
            MDME(IO,1)=1
          ENDDO
          DO IO=1011,1072
            R4=R4+BRAT(IO)
            MDME(IO,1)=0
          ENDDO
* now R* contain the branching ratios

          RA=R1M/(R1M+R1)
          RB=R2M/(R2M+R2)
          RC=R3M/(R3M+R3)
          RD=R4M/(R4M+R4)

          RMA=MAX(RA,RB,RC,RD)

          RA=RA/RMA
          RB=RB/RMA
          RC=RC/RMA
          RD=RD/RMA

* restore Dplus
          BRAT(377)=1-RA
          MDME(377,1)=1
          DO IO=378,388
            BRAT(IO)=BRAT(IO)*RA/R1M
          ENDDO
* restore D0
          BRAT(437)=1-RB
          MDME(437,1)=1
          DO IO=438,445
            BRAT(IO)=BRAT(IO)*RB/R2M
          ENDDO
* restore Ds
          BRAT(496)=1-RC
          MDME(496,1)=1
          DO IO=497,501
            BRAT(IO)=BRAT(IO)*RC/R3M
          ENDDO
* restore lambda
          BRAT(1003)=1-RD
          MDME(1003,1)=1
          DO IO=1004,1010
            BRAT(IO)=BRAT(IO)*RD/R4M
          ENDDO

        ENDIF
      ELSE
* D+_s forced to tau and then...
        MDME(491,1)=1
        DO IO=492,523
          MDME(IO,1)=0
        ENDDO
* so tau to chosen channel ...
        DO IO=83,123
          MDME(IO,1)=0
        ENDDO
        IF (IDSUBS.EQ.1) THEN
          MDME(84,1)=1
        ELSEIF (IDSUBS.EQ.2) THEN
          MDME(85,1)=1
        ENDIF


      ENDIF

********TAUOLA SECTION START

      CALL TAUINIT

*******TAUOLA SECTION END


* PARTIAL DECAY RULES:

      MSTJ(22)=2
* 2 MM
      PARJ(71)=2.0
*******INITIALIZE ZEBRA*****
      CALL FZINI
      CALL MZINI

      CALL HLIMIT(-MEMOR)
      CALL HROPEN(1,'histos','jetta.paw','N',1024,ISTAT)
      CALL HCDIR('//histos',' ')
      IF (ISTAT.NE.0) THEN
        WRITE(*,*)'+++MAIN: WARNING, HROPEN FAILED'
      ENDIF

* Differential cross section dsigma/dx/dy/de in pb
      FTAGS(1) ='SIGMA'
* Bjorken X
      FTAGS(2) ='X'
* Bjorken Y
      FTAGS(3) ='Y'
* Interacting neutrino energy (GeV)
      FTAGS(4) ='E'
* W**2 distribution (GeV**2)
      FTAGS(5) ='W2'
* Q**2 distribution (GeV**2)
      FTAGS(6) ='Q2'
* \nu distribution (GeV)
      FTAGS(7) ='U'
* Fermi motion of nucleon (X-Y-Z) cohordinates
      FTAGS(8) ='PFERX'
      FTAGS(9) ='PFERY'
      FTAGS(10)='PFERZ'
* Tau polarization vector
      FTAGS(11)='POLX'
      FTAGS(12)='POLY'
      FTAGS(13)='POLZ'
* Primary Charged lepton LAB momentum (i.e. Tau for nutau, mu for numu)
      FTAGS(14)='PLEPX'
      FTAGS(15)='PLEPY'
      FTAGS(16)='PLEPZ'
* Hadronic string momentum
      FTAGS(17)='PSTRX'
      FTAGS(18)='PSTRY'
      FTAGS(19)='PSTRZ'
* Number of Charged particles in tau decay
      FTAGS(20)='ICHGTAU'
* Number of Neutral particles in tau decay
      FTAGS(21)='INEUTAU'
* Number of Charged particles in Hadronic shower
      FTAGS(22)='ICHGPRI'
* Number of Neutral particles in Hadronic shower
      FTAGS(23)='INEUPRI'
* Energy of Charged particles in tau decay*
      FTAGS(24)='ECHGTAU'
* Energy of Neutral particles in tau decay
      FTAGS(25)='ENEUTAU'
* Energy of Charged particles in Hadronic shower
      FTAGS(26)='ECHGPRI'
* Energy of Neutral particles in Hadronic shower
      FTAGS(27)='ENEUPRI'
* 2*P*k relativistic invariant
      FTAGS(28)='TWOPK'
* 2*P*Q relativistic invariant
      FTAGS(29)='TWOPQ'
* Type of nucleon hit
      FTAGS(30)='NUCL'
* String energy
      FTAGS(31)='ESTR'
* Invariant mass of the string
      FTAGS(32)='WSTR'
* Primary Charged lepton Lab Energy (i.e. Tau for nutau, mu for numu)
      FTAGS(IOF1+1)='ELEP'
* Longitudinal polarization of tau
      FTAGS(IOF1+2)='LPOL'
* Momentum in tau rest frame of the Outgoing Charged particles (mu, pi-)
* summed over in 3 prong decays
      FTAGS(IOF1+3)='PVXRF'
      FTAGS(IOF1+4)='PVYRF'
      FTAGS(IOF1+5)='PVZRF'
* Event number
      FTAGS(IOF1+6)='EVT'
* Eventual Charm particle flag
      FTAGS(IOF1+7)='CHARM'
* mu minus momentum
      FTAGS(IOF1+8)='PMUM'
* mu plus momentum
      FTAGS(IOF1+9)='PMUP'
* charm decay length (lab frame)
      FTAGS(IOF1+10)='CHARLEN'
* charm type (Lund ID)
      FTAGS(IOF1+11)='CHARTYP'
* charm momentum, energy
      FTAGS(IOF1+12)='CHARP'
      FTAGS(IOF1+13)='CHARE'
* proton number (event-by-event beamfiles)
      FTAGS(IOF1+14)='PNUMB'
* charm three-momentum
      FTAGS(IOF1+15)='CHARPX'
      FTAGS(IOF1+16)='CHARPY'
      FTAGS(IOF1+17)='CHARPZ'

* Energy in tau rest frame of the Outgoing Charged particles (mu, pi-)
* Polarization Intensity
      FTAGS(51)='MPOL'
* Center of Mass Primary Lepton momentum, plus modulo
      FTAGS(52)='PLEPXCMS'
      FTAGS(53)='PLEPYCMS'
      FTAGS(54)='PLEPZCMS'
      FTAGS(55)='PLEPCMS'
* Center of mass Cos(theta) of the primary lepton
      FTAGS(56)='CSLEPCMS'
* Helicity, i.e. projection of Polarization vector on tau momentum
      FTAGS(57)='HELIX'
* visible momentum/energy from tau decay (neutrino excluded)
      FTAGS(58)='PVXTAU'
      FTAGS(59)='PVYTAU'
      FTAGS(60)='PVZTAU'
      FTAGS(61)='EVTAU'
* muon from tau decay: momentum
      FTAGS(62)='PMUX'
      FTAGS(63)='PMUY'
      FTAGS(64)='PMUZ'
* Total visible energy
      FTAGS(65)='ECAL'
* Total Electromagnetic energy
      FTAGS(66)='EEM'
* Total hadronic energy
      FTAGS(67)='EHAD'
* Total track multiplicity
      FTAGS(68)='NE'
* Total electromagnetic track multiplicity
      FTAGS(69)='NEEM'
* Total hadronic track multiplicity (ch+neutrals)
      FTAGS(70)='NEHAD'
* Total Neutral multiplicity
      FTAGS(71)='NEN'
* Total Neutral electromagnetic multiplicity
      FTAGS(72)='NENEM'
* Total Neutral hadronic multiplicity
      FTAGS(73)='NENHAD'
* Muon impact parater (distance from vertex)
      FTAGS(74)='BMU'
* Pion impact parameter (distance from vertex)
      FTAGS(75)='BPI'
* nucleon daughters, momentum and energy
      FTAGS(76)='PXNUC'
      FTAGS(77)='PYNUC'
      FTAGS(78)='PZNUC'
      FTAGS(79)='ENUC'
* tau decay pathlength (LAB frame)
      FTAGS(80)='TAULEN'
* electron in tau decay:momentum
      FTAGS(81)='PELECX'
      FTAGS(82)='PELECY'
      FTAGS(83)='PELECZ'
* random transverse position of events
      FTAGS(IOF2+1)='XEMU'
      FTAGS(IOF2+2)='YEMU'
* muon momentum OUT of shower plane
      FTAGS(IOF2+3)='PTMOP'
* shower momentum OUT of muon plane
      FTAGS(IOF2+4)='PTHOP'
* muon angle w.r.t. hadron plane
      FTAGS(IOF2+5)='THEMOP'
* hadron angle w.r.t. muon plane
      FTAGS(IOF2+6)='THEHOP'
* integrated cross section dsigma/de
      FTAGS(IOF2+7)='SIGMAPB'
* transverse position of events as from input-neutrino files
      FTAGS(IOF2+8)='EMUX'
      FTAGS(IOF2+9)='EMUY'




      CALL HBOOKN(11,'X-sect',ISIZ,'//histos',50000,FTAGS)


C...................... READ DATACARD AND INITIALIZE LEPTO
* INACTIVE LEPTON EVOLUTION IN ORDER TO STUDY THE PRIMARY HADRON SHOWER
*      LST(4)=0
* TEMPORARY: NO PHI ROTATION
*     LST(6)=0
* HADRONIZATION AND GENERATION OF SHOWER
      LST(8)=IGLU
*      LST(8)=2
* W**2 SCALE AS SUGGESTED FOR NOT SO LOW X INTERACTIONS
*      LST(9)=2
* ISOSCALAR/EMULSIONS TARGET
*      PARL(1)=2.
*      PARL(2)=1.
      PARL(1)=30.4
      PARL(2)=15.2
*INITIAL ENERGIES FROM P IN LUJETS
*      LST(17)=1

* INCREASE NUMBER OF ERRORS TO AVOID STOP
      MSTU(21)=1
      MSTU(22)=10
      CALL CATS

* TARGET REMNANTS ....WHO CARE? BUT WE PREFER A BARYON....
      LST(14)=1
* GRID SUITABLE FOR FIXED TARGET < 300 GEV
      LST(19)=1
* ADDED TAU LEPTON MASS EFFECTS
* THIS IS TAU CC
*      LST(32)=1
*      LEPIN=16
* THIS IS MU CC
*      LST(32)=0
*      LEPIN=14
* ADD FERMI MOMENTUM SMEARING ON THE NUCLEON
*      LST(33)=1
* LST(34)=1 INHIBITS THE TAU DECAY
      LST(34)=0
* LST(35) CONTAINS THE LOWER W2 LIMIT FOR PARTON SHOWER EVOLUTION
*~ 0.938+1 **2
      LST(35)=2.0
*      LST(35)=0.01
* TAUOLA TAU DECAY
      LST(36)=1
* FLAT ENERGY DISTRIBUTION
*      LST(37)=1
* GENERAL RUN PARAMETERS
*      INTERACTION=2
*      MAXIEVT=10000
* TO SET A MUON TAU SELECTION: 1=ONLY MU,2=NO MU,0=EVERYTHING
*       LST(38)=0
      IF (LST(36).EQ.1) LST(34)=1


      CALL GENTABLE (0,LEPIN,ENERGY_FIX,0.,INTERACTION)


      IF (LOME(1).LT.LOME(2)) THEN
        CALL LULIST(12)
      ENDIF
      IPROT=0
C...................... START OF EVENT

      ISTATUS=2
   10 CONTINUE
      CALL VZERO(FTUPLE,ISIZ)
      CALL MZWIPE(IXEVT)
      IF (LEPIN.EQ.16) NEUFORCE=51
      IF (LEPIN.EQ.14) NEUFORCE=51
      IF (LEPIN.EQ.-14) NEUFORCE=52
      IF (LEPIN.EQ.12) NEUFORCE=49
      IF (LEPIN.EQ.-12) NEUFORCE=50
C........................... CALL LEPTO

      IBADEV=0
   20 CONTINUE
      ICRACK=0

* take next neutrino

      IF (LST(17).GT.0) THEN

   39   CONTINUE
        IF (IHIST.EQ.0) THEN
*single neutrino operations
          CALL GETNEU(IPNUMBER,NEUTYPE,VECT,GKIN, MESTYPE,G4MES,
     +    NEUFORCE,ISTATUS)
          IF (ISTATUS.NE.4) THEN
            PNUMBER=IPNUMBER
          ELSE
            WRITE(*,*)' END OF NEUTRINO BEAMDATA:RELOOP'
            ISTATUS=2
            GOTO 39
          ENDIF
        ELSE
*simulate a real neutrino using histogram file

          CALL GETHNEU(IPNUMBER,NEUTYPE,VECT,GKIN, MESTYPE,G4MES,
     +    NEUFORCE,ISTATUS)
        ENDIF

C       PNUMBER = current proton number
C       NEUTYPE = type of neutrino;
C       51= muon         neutrino
C       52= muon     antineutrino
C       49= electron     neutrino
C       50= electron antineutrino
C       VECT(3) = coordinates of neutrino origin
C       GKIN(3) = momentum components of the neutrino


* to be changed when nu_E oscillations become of interest
        IF (LEPIN.EQ.16.AND.NEUTYPE.NE.51) GOTO 39
        IF (LEPIN.EQ.14.AND.NEUTYPE.NE.51) GOTO 39
        IF (LEPIN.EQ.-14.AND.NEUTYPE.NE.52) GOTO 39
        IF (LEPIN.EQ.12.AND.NEUTYPE.NE.49) GOTO 39
        IF (LEPIN.EQ.-12.AND.NEUTYPE.NE.50) GOTO 39

        PTEST=SQRT(GKIN(1)**2+GKIN(2)**2+GKIN(3)**2)
        IF (PTEST.LE.3.5.AND.LEPIN.EQ.16) GOTO 39
        PFIN=DISTRR(PTEST)

        IF (PFIN.EQ.0) GOTO 39

        CALL BTOCHO2(VECT,GKIN,EMUX,EMUY)

        FTUPLE(IOF2+8)=EMUX
        FTUPLE(IOf2+9)=EMUY

        P(1,1)=0.
        P(1,2)=0.
        P(1,3)=PFIN
*        P(1,3)=100.
        P(1,4)=P(1,3)
        P(1,5)=0.
      ENDIF

      CALL VZERO(PPXYZ,3)

c      IF (LST(33).EQ.1.AND.LST(17).NE.0) CALL FERMII(PPXYZ)
       IF( LST(33).EQ.1) CALL FERMII(PPXYZ)
      DO I=1,3
        FTUPLE(7+I)=PPXYZ(I)
      END DO

      P(2,1)=PPXYZ(1)
      P(2,2)=PPXYZ(2)
      P(2,3)=PPXYZ(3)
      RNUCKIN2=P(2,1)**2+P(2,2)**2+P(2,3)**2
      P(2,4)=SQRT( RNUCKIN2 +P(2,4)**2)
      P(2,5)=0.938


      IBADEV=IBADEV+1
*      WRITE(*,*)' NU_TAU ENERGY=',P(1,4)
*      WRITE(*,*)' NUCLEON ENERGY=',P(2,4)


      CALL CATS
      CALL LEPTO
*      CALL LULIST(3)

      IF (ICRACK.NE.0) THEN
        WRITE(*,*)'+++MAIN: ICRACK SAFETY TRAPPED'
        GOTO 20
      ENDIF

      IF (LST(21).NE.0) THEN
        IF (LST(21).NE.3131) WRITE (*,10100) LST(21),IBADEV
10100 FORMAT (/,10X,'!!!!!! LST(21)=',I10,' AFTER ',I2,' CALL TO LEPTO')
        GOTO 20
      ENDIF
*

      CALL PARUPD

      IF (ABS(P(2,3)+P(1,3)).LT.0.01) THEN
        WRITE(*,*)'ERROR CMS FRAME!! 28,29=',LST(28),LST(29)
        GOTO 20
      ENDIF

      IF (LST(32).EQ.1) THEN
* TRANSFORM TO THE LEPTON CMS TO CALCULATE ITS POLARIZATION
        DO I=1,3
          BB(I)=-P(4,I)/P(4,4)
        END DO
*       WRITE(*,*)'BB=',BB
        CALL LUROBO(0.0,0.0,BB(1),BB(2),BB(3))
* NOW PARTICLES ARE IN THE SCATTERED LEPTON CMS


        RML=P(4,5)
        RMM=P(2,5)
        EE=PARL(21)/2./RMM
        QQ2=+Q2
        QM2=QQ2+RML**2

        NU=U*RMM

        FRAC=QM2*WW1 + (2.*EE*(EE-U) - 0.5*QM2)*WW2 - 0.5/RMM**2*(2.*
     +  RMM*EE*QQ2 - NU*QM2)*WW3 - RML**2/RMM*EE*WW5

        FACTK=2.*WW1 -WW2 - EE/RMM*WW3 +(EE-U)/RMM*WW5
        FACTP=2.*EE/RMM*WW2 - QM2/2./RMM**2*(WW3+WW5)

        DO I=1,3
          POL(4,I)=RML*(FACTK*P(1,I)+FACTP*P(2,I))/FRAC
          POLARX(I)=POL(4,I)
        END DO

        PMODUL=0.
        DO I=1,3
          PMODUL=PMODUL+POL(4,I)**2
        END DO
        IF(PMODUL.GT.1.05) WRITE(*,*)'PMODUL>1 ',SQRT(PMODUL)

*      CALL LUROBO(0.,0.,-BB(1),-BB(2),-BB(3))
        CALL LUDBRB(1,4,0.,0.,DBLE(-BB(1)),DBLE(-BB(2)),DBLE(-BB(3)))

      ENDIF

      IF (LST(32).EQ.2) THEN
        DO I=1,3
          POLARX(I)=0.
        END DO
      ENDIF


C DECAY....
      IF (LST(36).EQ.1.AND.LEPIN.EQ.16) THEN
* REMEMEBER, NEEDS TAU POSITION AND TAU POLARIZATION:NPA,NPB
        NPB=0
        NPA=0
        DO I=1,N
          IF(ABS(K(I,2)).EQ.15) THEN
            NPA=I
            GOTO 30
          ENDIF
        END DO

   30   CONTINUE

        IF (NPA.NE.0) THEN
C...CHOOSE LIFETIME AND DETERMINE DECAY VERTEX.
          KFA=IABS(K(NPA,2))
          KC=LUCOMP(KFA)
          IF(K(NPA,1).EQ.5) THEN
            V(NPA,5)=0.
          ELSEIF(K(NPA,1).NE.4) THEN
            V(NPA,5)=-PMAS(KC,4)*LOG(RLU(0))
          ENDIF
          DO 40  J=1,4
            VDECY(J)=V(NPA,J)+V(NPA,5)*P(NPA,J)/P(NPA,5)
   40     CONTINUE

          NBAK=N


*        CALL LULIST(3)
          DO I=1,3
            BB(I)=-P(4,I)/P(4,4)
          END DO
*       WRITE(*,*)'BB=',BB
*        CALL LUROBO(0.0,0.0,BB(1),BB(2),BB(3))
          CALL LUDBRB(1,4,0.,0.,DBLE(BB(1)),DBLE(BB(2)),DBLE(BB(3)))

* PATCH TO RECOVER FROM LUHEPC BUG
          DO II=1,N
            VNPALIFE(II)=V(II,5)
          END DO
          CALL LUHEPC(1)
          CALL DEXAY(2,POLARX)
*       CALL DEKAY(2,HH)
*       CALL DEKAY(12,HH)
          CALL LUHEPC(2)
*
          DO IK=1,4
            FTUPLE(IOF1+2+IK)=0.
          END DO
          DO II=NBAK+1,N
            IF(ABS(K(II,2)).NE.16.AND.ABS(K(II,2)).NE.14.AND. ABS(K(II,
     +      2)).NE.12) THEN
              DO IK=1,4
                FTUPLE(IOF1+2+IK)=FTUPLE(IOF1+2+IK)+P(II,IK)
              END DO
            ENDIF
          END DO


*        CALL LULIST(3)
* APPARENT BUG IN LUHEPC
          CALL LUROBO(0.0,0.0,-BB(1),-BB(2),-BB(3))

          FTUPLE(58)=0.
          FTUPLE(59)=0.
          FTUPLE(60)=0.
          FTUPLE(61)=0.


          DO II=1,NBAK
            V(II,5)=VNPALIFE(II)
          END DO
          DO II=NBAK+1,N
            DO JJ=1,4
              V(II,JJ)=VDECY(JJ)
            END DO
          END DO
*        CALL LULIST(3)
* NEW STRATEGY: NO DECAYS EXCEPT TAU AND CHARM
          CALL LUEXEC
          DO II=NBAK+1,N
            IF (K(II,2).EQ.13) THEN
              DO IK=1,3
                FTUPLE(61+IK)=P(II,IK)
              END DO
            ELSEIF (K(II,2).EQ.11) THEN
              DO IK=1,3
                FTUPLE(80+IK)=P(II,IK)
              END DO
            ENDIF

            IF(ABS(K(II,2)).NE.16.AND.ABS(K(II,2)).NE.14.AND. ABS(K(II,
     +      2)).NE.12) THEN
              DO IK=1,4
                FTUPLE(57+IK)=FTUPLE(57+IK)+P(II,IK)
              END DO
            ENDIF
          END DO

        ENDIF
      ENDIF

      IMU=0
      ICHM=0
      IPSELFOUND=0
      DO II=1,N
        IAKI=ABS(K(II,2))
        IF (IAKI.EQ.92.OR.IAKI.EQ.91) THEN
          ESTR=P(II,4)
          IF (EHAC.NE.0) THEN
            IF (ESTR.GT.EHAC) THEN
              IEHAC=IEHAC+1
              GOTO 10
            ENDIF
          ENDIF
          IF (ESTR.EQ.0) THEN
            WRITE(*,*)'+++MAIN: STRING ENERGY=0,EVT=',IEVT
            GOTO 10
          ENDIF
        ENDIF

        IF (IAKI.GE.400
     +  .AND.IAKI.LE.500
     +  .OR.IAKI.EQ.4122) THEN
          IF (IDSUBS.EQ.0) THEN
            IF (ICCHA.EQ.0) THEN
              ICHM=ICHM+1
            ELSE
              IF (IAKI.EQ.411.OR.IAKI.EQ.431.OR.IAKI.EQ.4122) THEN
                ICHM=ICHM+1
              ENDIF
            ENDIF
          ELSE
            IF (IAKI.EQ.431) ICHM=ICHM+1
          ENDIF
        ENDIF

        IF (IAKI.EQ.ABS(IPSEL)) THEN
          IPSELFOUND=1
        ENDIF


        IF (IAKI.EQ.13) THEN
          IMU=IMU+1
          IF (K(II,2).EQ.13) THEN
            IHM=II
          ELSE
            IHMP=II
          ENDIF
        ENDIF
      END DO
   50 CONTINUE


      IPI=0
      DO II=NBAK+1,N
        IF (K(II,2).EQ.-211) THEN
          IPI=1
          IHP=II
          GOTO 60
        ENDIF
      END DO
   60 CONTINUE

      IF (IPSELFOUND.EQ.0..AND.IPSEL.NE.0) THEN
        IPSELTHROW=IPSELTHROW+1
        GOTO 20
      ENDIF

      IF (IDIMUON.GE.1.AND.ICHM.LT.1) THEN
        GOTO 20
      ENDIF

      IF (IDIMUON.GE.2.AND.IMU.LE.IMUREQ) THEN
*        WRITE(*,*)' TOO FEW MUONS DETECTED W.R.T. REQUIRED'
*        CALL LULIST(3)
        MUMISS=MUMISS+1
        GOTO 20
      ENDIF


      IF ((IMU.GE.1.AND.LST(38).EQ.2).OR.
     + (IMU.EQ.0.AND.LST(38).EQ.1)) THEN
        GOTO 20
      ENDIF


* FUNDAMENTAL HELICITY CROSS_CHECK


      DO JIJ=1,3
        BRF(JIJ)=(P(1,JIJ)+P(2,JIJ))/(P(1,4)+P(2,4))
      END DO

      CALL LUDBRB(1,4,0.,0.,DBLE(-BRF(1)),DBLE(-BRF(2)),DBLE(-BRF(3)))
*        CALL LUROBO(0.0,0.0,-BRF(1),-BRF(2),-BRF(3))
*        CALL LULIST(1)

      FTUPLE(55)=0.
      DO JIJ=1,3
        FTUPLE(51+JIJ)=P(4,JIJ)
        FTUPLE(55)=FTUPLE(55)+P(4,JIJ)**2
      END DO
      FTUPLE(55)=SQRT(FTUPLE(55))
      FTUPLE(51)=SQRT(POLARX(1)**2+POLARX(2)**2+POLARX(3)**2)

      FTUPLE(56)=(P(1,1)*P(4,1)+P(1,2)*P(4,2)+P(1,3)*P(4,3)) /P(1,4)/
     +FTUPLE(55)

      IF (IF5CC.EQ.1) THEN
        FTUPLE(57)=(POLARX(1)*FTUPLE(52)+POLARX(2)*FTUPLE(53)+
     +  POLARX(3)* FTUPLE(54))/FTUPLE(55)/FTUPLE(51)
      ENDIF

      CALL LUDBRB(1,4,0.,0.,DBLE(BRF(1)),DBLE(BRF(2)),DBLE(BRF(3)))
*        CALL LUROBO(0.0,0.0,BRF(1),BRF(2),BRF(3))

      IF (IMU.EQ.1) THEN
        POPI=P(IHM,1)**2+P(IHM,2)**2+P(IHM,3)**2
        TVAR=-(V(IHM,1)*P(IHM,1)+V(IHM,2)*P(IHM,2)+
     +  V(IHM,3)*P(IHM,3))/POPI

        FTUPLE(74)=SQRT(
     +  (V(IHM,1)+P(IHM,1)*TVAR)**2+
     +  (V(IHM,2)+P(IHM,2)*TVAR)**2+
     +  (V(IHM,3)+P(IHM,3)*TVAR)**2  )
      ENDIF

      IF (IMU.EQ.2) THEN
        FTUPLE(IOF1+8)=SQRT(P(IHM,1)**2+P(IHM,2)**2+P(IHM,3)**2)
        FTUPLE(IOF1+9)=SQRT(P(IHMP,1)**2+P(IHMP,2)**2+P(IHMP,3)**2)
      ENDIF


      IF (IPI.EQ.1) THEN
        POPI=P(IHP,1)**2+P(IHP,2)**2+P(IHP,3)**2
        TVAR=-(V(IHP,1)*P(IHP,1)+V(IHP,2)*P(IHP,2)+
     +  V(IHP,3)*P(IHP,3))/POPI

        FTUPLE(75)=SQRT(
     +  (V(IHP,1)+P(IHP,1)*TVAR)**2+
     +  (V(IHP,2)+P(IHP,2)*TVAR)**2+
     +  (V(IHP,3)+P(IHP,3)*TVAR)**2  )
      ENDIF


      TEMP=SQRT(P(4,1)**2+P(4,2)**2+P(4,3)**2)
      FTUPLE(IOF1+1)=P(4,4)
      FTUPLE(IOF1+2)=(POLARX(1)*P(4,1)+POLARX(2)*P(4,2)+
     +POLARX(3)*P(4,3))
     +/TEMP
      FTUPLE(4)=P(1,4)
      CALL HFILL(1000,P(1,4),0.,1.)
      CALL HFILL(1001,SQRT(RNUCKIN2),0.,1.)
      CALL HFILL(1002,SQRT(W2),0.,1.)
      CALL HFILL(1003,SQRT(Q2),0.,1.)
      CALL HFILL(1006,P(1,4),0.,FTUPLE(1))
      FTUPLE(5)=W2
      FTUPLE(6)=Q2
      FTUPLE(7)=U
      FTUPLE(28)=PARL(21)
      FTUPLE(29)=PARL(22)
      FTUPLE(30)=K(2,2)
      DO I=1,3
        FTUPLE(13+I)=P(4,I)
      END DO
      DO I=1,3
        FTUPLE(I+10)=POL(4,I)
      END DO

      CHK1=0
      CHK2=0
      RFACT=1
      DO I=1,N

        IF ((I.EQ.4.AND.V(I,5).EQ.0).OR.(K(I,1).LE.10.AND.I.GT.4)) THEN
          CHK1=CHK1+P(I,1)
          CHK2=CHK2+P(I,2)
        ENDIF

        IF (K(I,3).EQ.2) THEN
* NUCLEON REMNANTS
          DO IJ=1,4
            FTUPLE(75+IJ)=P(I,IJ)
          END DO
        ENDIF

        IF (K(I,2).EQ.91) THEN
* this takes into account, at the checking level,
* that clustering requires some "forced" momentum conservation
* arising from the initial transverse momenta of the quarks
* inside the nucleon
          RFACT=10
        ENDIF

        IF (K(I,2).EQ.92.OR.K(I,2).EQ.91) THEN
          DO IJ=1,3
            FTUPLE(16+IJ)=P(I,IJ)
            PH(IJ)=P(I,IJ)
          END DO
        ENDIF

        IF (K(I,2).EQ.13) THEN
          DO IJ=1,3
            PM(IJ)=P(I,IJ)
          END DO
        ENDIF


        IF (K(I,2).EQ.15) THEN
* THIS IS THE TAU
          CALL HFILL(1004,P(I,4),0.,1.)
        ENDIF

        IF ((K(I,2).LT.10).AND.K(I,3).EQ.0) THEN
* THIS IS THE FIRST PRIMARY PARTON TAKING AWAY ALL THE ENERGY
          CALL HFILL(1005,P(I,4),0.,1.)
        ENDIF

      END DO


      CHK1=CHK1-P(2,1)
      CHK2=CHK2-P(2,2)
      RELERR=0.03*RFACT
      IF (ABS(CHK1).GT.RELERR.OR.ABS(CHK2).GT.RELERR) THEN
        CALL LULIST(3)
        WRITE(*,*)'ERROR IN CHK',CHK1,CHK2
        GOTO 10
      ENDIF

   70 CONTINUE



      IEVT=IEVT+1

      CALL JETMC
* NOW WE HAVE A MONTECARLO SEPRATION OF THE TWO JETS...
      CALL JETTA
*     WRITE(*,*)'W2,U,Q2,X,Y=',W2,U,Q2,X,Y


      IF(IEVT.GE.LOME(1).AND.IEVT.LE.LOME(2)) THEN
        WRITE(*,*)'**********************************'
        WRITE(*,*)'********LOOK AT ME****************'
        WRITE(*,*)'EVENT NUMBER:',IEVT
        WRITE(*,*)'**********************************'
        CALL LULIST(3)
      ENDIF

C...................... END OF EVENT
      IF (MOD(IEVT,100).EQ.0) THEN
        WRITE(*,*)'EVENTS=',IEVT,'SIGMA=',PARL(24)
      ENDIF


      TEST=P(1,4)+P(2,4)-FTUPLE(31)-P(4,4)

      IF(TEST.LT.-0.2) THEN
        WRITE(*,*)'...HERE BOZZO:EVENT/LST(24)',IEVT,LST(24)
*        WRITE(*,*)'BEGIN BOZZO EVENT,EVT=',IEVT,' ISFIRST=',ISFIRST
*        CALL LULIST(3)
*        WRITE(*,*)'END BOZZO EVENT,EVT=',IEVT
      ELSE
*        WRITE(*,*)'GOOD EVENT,EVT=',IEVT,' ISFIRST=',ISFIRST
* this IS an event, book banks
        IF (JGEEV.NE.0) CALL MZDROP(IXSTOR,JGEEV,'.')
        CALL MZBOOK(IXEVT,JGEEV,JGEEV,2,'GEEV',4,4,0,2,0)
        CALL EVTINFO
        CALL JETTOUT
        CALL MZBOOK(IXEVT,JGELU,JGEEV,-2,'GELU',N,N,0,2,0)
        DO IAS=1,N
          JGELU=LQ(JGEEV-2)
          CALL MZBOOK(IXEVT,JGELN,JGELU,-IAS,'GELN',0,0,16,3,0)
          DO IAT=1,5
            Q(JGELN+IAT)=K(IAS,IAT)
            Q(JGELN+IAT+5)=P(IAS,IAT)
            Q(JGELN+IAT+10)=V(IAS,IAT)
            Q(JGELN+16)=DALUAEF(IAS)
          ENDDO
        ENDDO
        CALL ZVERIF (IXEVT,IFLRTN,'Verification')
        KEYLIST(1)=121
        KEYLIST(2)=1990
        KEYLIST(3)=-IEVT
        KEYLIST(4)=-99999
*        CALL DZSHOW('GEEV Bank:',0,JGEEV,'BLV',0,0,0,0)
        CALL FZOUT(LUNFZ,IXEVT,0,1,'Z',2,4,KEYLIST)
        CALL FZOUT(LUNFZ,IXEVT,JGEEV,0,' ',2,1,1990)
        CALL FZOUT(LUNFZ,IXEVT,0,0,'Z',2,1,-1)
        IF(IEVT.GE.LOME(1).AND.IEVT.LE.LOME(2)) THEN
          WRITE(*,*)'**********************************'
        ENDIF
      ENDIF


* OUT OF PLANE ANALYSIS

      PHL=SQRT(PH(1)**2+PH(2)**2+PH(3)**2)
      PML=SQRT(PM(1)**2+PM(2)**2+PM(3)**2)

      CALL ORTH(PO,PM,PH)
      IF (PHL.GT.0.AND.PML.GT.0) THEN
        FTUPLE(IOF2+3)=PO
      ENDIF
      IF (PML.GT.0) THEN
        FTUPLE(IOF2+5)=ACOS(PM(3)/PML)
      ENDIF

      CALL ORTH(PO,PH,PM)
      FTUPLE(IOF2+4)=PO
      IF (PHL.GT.0) THEN
        FTUPLE(IOF2+6)=ACOS(PH(3)/PHL)
      ENDIF

      FTUPLE(IOF2+7)=XSECT
      FTUPLE(IOF1+14)=PNUMBER

      IF(IEVT.LE.10000) THEN
        FTUPLE(IOF1+6)=IEVT
c        WRITE(*,5341) ievt
c        WRITE(*,5342) ineut,interaction,lst(22),ftuple(4)
        WRITE(88,5341) ievt,n
        WRITE(88,5342) ineut,interaction,lst(22),ftuple(4)
         do ms=1,n
           write(88,5343) ms,k(ms,1),k(ms,2),k(ms,3),k(ms,4),
     +k(ms,5),p(ms,1),p(ms,2),p(ms,3),p(ms,4)
         end do
5341      FORMAT(1x,i6,1x,i6)
5342      FORMAT(1x,3i6,e15.7)
5343      format(1x,i3,5i6,4e15.7)
        CALL HFN(11,FTUPLE)
      ENDIF
*       CALL LUEDIT(13)
*       CALL LULIST(3)
      JARRY(JEIN)=JARRY(JEIN)+1
      IF (IEVT.LT.MAXIEVT) GOTO 10

 3131 CONTINUE
      IF (LST(36).EQ.1) THEN
        CALL DEXAY(100,POLARX)
        DO IJI=1,30
          IJEJE=IJEJE+JARRY(IJI)
          IJAJA=IJAJA+JALLY(IJI)
        END DO
        DO IJI=1,30
          WRITE(*,*)'JARRY-JALLY JAKER SAYS',FLOAT(JARRY(IJI))/IJEJE,
     +    FLOAT(JALLY(IJI))/IJAJA,'FOR JAK',IJI
        END DO
      ENDIF
      CALL FZRUN(LUNFZ,-99999,0,0)
      CALL FZCLOS
      CALL HCDIR('//histos',' ')
      CALL HLDIR('//histos','T')
      CALL HROUT(11,ICYCLE,' ')
      CALL HREND('histos')
      WRITE(*,*)' MISSED MUONS in DIMUON GENERATION:',MUMISS,' IN ',
     + IEVT,' EVENTS (',FLOAT(MUMISS)/FLOAT(IEVT),'%)'
      WRITE(*,*)' REJECTED EVENTS FOR EHAD CUT:',IEHAC,' IN ',
     + IEVT,' EVENTS (',FLOAT(IEHAC)/FLOAT(IEVT),'%)'
      WRITE(*,*)' REJECTED EVENTS FOR PARTICLE ID=',IPSEL,
     +':KEPT ',
     + IEVT,' REJECTED',IPSELTHROW

      STOP
      END
*CMZ :          04/03/97  13.03.17  by  Unknown
*CMZ :  1.01/51 24/05/96  11.17.19  by  Piero Zucchelli
*CMZ :  1.01/40 11/12/95  19.22.51  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  19.28.31  BY  PIERO ZUCCHELLI
*CMZ :  1.01/04 05/03/95  00.45.11  BY  PIERO ZUCCHELLI
*CMZ :  1.01/03 05/03/95  00.22.25  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 15/08/94  03.58.02  BY  UNKNOWN
*CMZ :  1.00/00 20/07/94  12.12.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LEPTO



C...ADMINISTER THE GENERATION OF AN EVENT.
C...NOTE: IF ERROR FLAG LST(21) IS NON-ZERO, NO PROPER EVENT GENERATED.

      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
      COMMON /ARDAT1/ PARA(40),MSTA(40)
      COMMON/WLIST/WW1,WW2,WW3,WW5
*KEEP,POLAR.
C--
      COMMON /POLARIZ/POL(4000,3)
      REAL POLARX(4)
*KEEP,JETTA.
C--
         PARAMETER (ICENTO=100)
         PARAMETER (ISIZ=93)
         PARAMETER (IOF1=32)
         PARAMETER (IOF2=83)
         PARAMETER (LUX_LEVEL=4)
         INTEGER*4 JTAU(100),JPRI(100),JSTRO(100)
         REAL*4 FTUPLE(ISIZ)
         COMMON/JETTAGL/JTAU,JPRI,JSTRO
         COMMON/NTUPLA/FTUPLE,ISFIRST
         COMMON/BEAM/SPEC(ICENTO)
         COMMON /MAXSPEC/RMAXSPEC,RINTSPEC
         COMMON/SAV/XMINSAV(ICENTO),XMAXSAV(ICENTO),YMINSAV(ICENTO),
     &   YMAXSAV(ICENTO),Q2MINSAV(ICENTO),Q2MAXSAV(ICENTO),
     &   W2MINSAV(ICENTO),W2MAXSAV(ICENTO),PARIMAX(ICENTO),
     &   PPSAVE(ICENTO,3,4,5),PARICOR(ICENTO),INDEX,SIGMASAV(ICENTO),
     &   XMSIGMA,XSECT

*KEND.
      DOUBLE PRECISION DBETATMP(3)
      DOUBLE PRECISION DTHETA,DPHI,DBETA,DETOT,DARI29,DARI30
      REAL*4 BB(3)
      DIMENSION SPQ(17)
      DATA NUMMIS,NWARN/0,10/,DARI29,DARI30/2*0.D0/

      ISFIRST=1

   10 LST(21)=0
      DO 20 I=1,10
        DO 20 J=1,5
          K(I,J)=0
   20 V(I,J)=0.
      DO 30 I=1,4
        K(I,1)=21
   30 K(I,2)=KSAVE(I)
      K(4,1)=1
      N=2
      IF(LST(17).NE.0.AND.LST(2).GT.0) THEN
C...LEPTON AND/OR NUCLEON ENERGY MAY VARY FROM EVENT TO EVENT,
C...MOMENTUM VECTORS TAKEN FROM P(I,1), P(I,2) AND P(I,3), I=1,2
        DO 50 I=1,2
          DO 40 J=1,5
            IF(ISFIRST.EQ.1) THEN
              PSAVE(3,I,J)=P(I,J)
            ELSE
              P(I,J)=PSAVE(3,I,J)
            ENDIF
   40     CONTINUE
          P(I,5)=ULMASS(K(I,2))
          P(I,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2)
C...Momentum vectors from PSAVE if new try, i.e. jump back to 1
          DO 25 II=1,2
            DO 25 J=1,5
   25     P(II,J)=PSAVE(3,II,J)
   50  CONTINUE 
C...TRANSFORM TO CMS OF INCOMING PARTICLES, LEPTON ALONG +Z AXIS.
        DO 60 J=1,3
   60   DBETA(1,J)=(DBLE(P(1,J))+DBLE(P(2,J)))/
     +             (DBLE(P(1,4))+DBLE(P(2,4)))
        IF (ISFIRST.EQ.1) THEN
          ISFIRST=0
          DO 70 J=1,3
   70     DBETATMP(J)=DBETA(1,J)
        ENDIF

        CALL LUDBRB(0,0,0.,0.,-DBETA(1,1),-DBETA(1,2),-DBETA(1,3))
        SPHI(1)=ULANGL(P(1,1),P(1,2))
        CALL LUDBRB(0,0,0.,-SPHI(1),0.D0,0.D0,0.D0)
        STHETA(1)=ULANGL(P(1,3),P(1,1))
        CALL LUDBRB(0,0,-STHETA(1),0.,0.D0,0.D0,0.D0)
        LST(28)=2
        PARL(21)=2.*(P(1,4)*P(2,4)-P(1,3)*P(2,3))
      ELSE
C...INITIAL STATE MOMENTA FIXED FROM LINIT CALL.
        DO 90 I=1,2
          DO 80 J=1,5
   80     P(I,J)=PSAVE(3,I,J)
   90   IF(PSAVE(3,1,3).LT.0.) P(I,3)=-PSAVE(3,I,3)
        LST(28)=3
      ENDIF

      CALL LEPTOX
      IF (W2.LT.LST(35)) THEN
*      CALL RESONANCES REGION
*      ...FOR THE MOMENT,
        LST(21)=3131
*      WRITE(*,*)' ERROR 21- GOTO 1=',LST(21)
        IERR31=IERR31+1
        IF(MOD(IERR31,5000).EQ.0) WRITE(*,*)'STATUS  21: CASE 3131',IERR31
*      GOTO 1
      ENDIF

C...RETURN IF ERROR OR IF NO EVENT TO BE GENERATED.
c         print*,'lst(21),lst(2),lst(7)'
c        write(*,*) lst(21),lst(2),lst(7)
      IF(LST(21).NE.0.OR.LST(2).LE.0.OR.LST(7).EQ.-1) THEN
        RETURN
      ENDIF

      IF(PARI(29).LT.0.5) THEN
C...FOR FIRST CALL, RESET DOUBLE PRECISION COUNTERS.
        DARI29=0.D0
        DARI30=0.D0
      ENDIF
      DARI29=DARI29+1.D0
      PARI(29)=DARI29


C     CALL GULIST(-3,2)
C...SCATTERED LEPTON AND EXCHANGED BOSON ADDED TO EVENT RECORD IN LKINEM
C...TRANSFORM TO LEPTON-NUCLEON CMS IF NOT MADE EARLIER
      IF(LST(17).EQ.0) THEN
        DO 110 I=3,4
          DO 100 J=1,5
  100     PSAVE(3,I,J)=P(I,J)
  110   IF(PSAVE(3,1,3).LT.0.) PSAVE(3,I,3)=-P(I,3)
        CALL LUDBRB(0,0,0.,0.,0.D0,0.D0,-DBETA(1,3))
        LST(28)=2
      ENDIF
      DO 120 I=1,4
        DO 120 J=1,5
  120 PSAVE(2,I,J)=P(I,J)
C     CALL GULIST(-2,2)

C...PREPARE FOR PARTON CASCADE.
      IF(LST(8).GE.2.AND.MOD(LST(8),10).NE.9) CALL LSHOWR(0)

C...TRANSFORM TO HADRONIC CMS, BOOST PARAMETERS IN DOUBLE PRECISION.
      DETOT=DBLE(P(1,4))-DBLE(P(4,4))+DBLE(P(2,4))
      DBETA(2,1)=-DBLE(P(4,1))/DETOT
      DBETA(2,2)=-DBLE(P(4,2))/DETOT
      DBETA(2,3)=(DBLE(P(1,3))-DBLE(P(4,3))+DBLE(P(2,3)))/DETOT
      CALL LUDBRB(0,0,0.,0.,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3))
      SPHI(2)=0.
      STHETA(2)=ULANGL(P(3,3),P(3,1))
      CALL LUDBRB(0,0,-STHETA(2),0.,0.D0,0.D0,0.D0)
      LST(28)=1
      DO 130 I=1,4
        DO 130 J=1,5
  130 PSAVE(1,I,J)=P(I,J)
C...SAVE MOMENTUM OF EXCHANGED BOSON (USED IN SUBROUTINE LFRAME).
      DO 140 J=1,5
  140 PB(J)=P(3,J)
*     WRITE(*,*)'HADRONIC CMS....'
*      CALL LULIST(1)




  150 N=4
      MSTU(1)=N+1
      LST(26)=N+1
      LST(27)=0
      PARL(25)=ULALPS(Q2)
      IF(LST(8).EQ.1.OR.LST(8)/10.EQ.1.OR.MOD(LST(8),10).EQ.9) THEN
C...PROBABILITIES FOR HARD, FIRST ORDER QCD EVENTS.
        CALL LQCDPR(QG,QQB)
        DO 160 I=1,17
  160   SPQ(I)=PQ(I)
  170   SRLU=RLU(0)
        IF(SRLU.GT.QQB+QG) THEN
          DO 180 I=1,17
  180     PQ(I)=SPQ(I)
          CALL LQEV
        ELSEIF(SRLU.GT.QQB) THEN
          IF(LST(8).EQ.9) THEN
            DO 190 I=1,17
  190       PQ(I)=SPQ(I)
            CALL LQEV
          ELSE
            CALL LQGEV
          ENDIF
        ELSE
          CALL LQQBEV
          IF(LST(8).EQ.9.AND.LST(21).EQ.0) THEN
            IF(PLU(5,11).LT.Q2*PARA(20)) THEN
              DO 200 I=1,17
  200         PQ(I)=SPQ(I)
              CALL LQEVAR(K(5,2),K(7,2))
            ENDIF
          ENDIF
        ENDIF
        IF(LST(21).NE.0) GOTO 170
      ELSE
C...QPM MODEL WITHOUT QCD CORRECTIONS (CASCADE APPLIED LATER).
  210   CALL LQEV
*     WRITE(*,*)'AFTER QPM CALL TO LQEV'
        IF(LST(21).NE.0) GOTO 210
      ENDIF

      NS=MSTU(1)
      MSTU(1)=0
*     WRITE(*,*)' AFTER QPM:'
*      CALL LULIST(1)


      IF(LST(8).LE.1.OR.MOD(LST(8),10).EQ.9) THEN
C...NO PARTON CASCADE, INTRODUCE PRIMORDIAL KT.
        IF(PARL(3).GT.1.E-03) THEN
          CALL LPRIKT(PARL(3),PT,PHI)
          CALL LUDBRB(NS,N,0.,-PHI,0.D0,0.D0,0.D0)
          CALL LUDBRB(NS,N,ATAN(2.*PT/SQRT(W2)),PHI,0.D0,0.D0,0.D0)
        ENDIF
C...CHECK SYSTEM AGAINST FRAGMENTATION CUTS.
        MSTU(24)=0
        CALL LUPREP(0)
        IF(MSTU(24).NE.0) THEN
          IF(LST(3).GE.1) WRITE(6,*) ' LUPREP ERROR MSTU(24)= ',MSTU(24)
          LST(21)=11
        ENDIF
      ELSEIF(LST(24).EQ.1) THEN
C...INCLUDE PARTON CASCADES (+ REMNANT & KT) ON Q-EVENT
*       WRITE(*,*)'AND NOW GO INTO Q SHOWERING!'
        IF(LST(21).NE.0) THEN
*      WRITE(*,*)' ERROR 21- PRE LSHOWR(1)'
        ENDIF

* CUT ON SHOWER SIMULATION
        CALL LSHOWR(1)

*       WRITE(*,*)' AFTER SHOWER THE STATUS IS:'
*       CALL LULIST(1)
      ELSE
C...INCLUDE PARTON CASCADES (+ REMNANT & KT) ON QG- OR QQBAR-EVENT
*       WRITE(*,*)'AND NOW GO INTO QG-QQBAR SHOWERING!'
        CALL LMEPS
      ENDIF

  220 CONTINUE




      IF(LST(21).NE.0) THEN
*      WRITE(*,*)' ERROR 21- GOTO 1=',LST(21)
        IERR21=IERR21+1
        IF(MOD(IERR21,100).EQ.0) WRITE(*,*)'ERROR 21:',IERR21
        GOTO 10
      ENDIF

      DO 230 I=1,N
C...CORRECT ENERGY-MOMENTUM-MASS MISMATCH FOR REAL PARTICLE
        IF(P(I,5).LT.0.) GOTO 230
        ENERGY=SQRT(DBLE(P(I,5))**2+DBLE(P(I,1))**2+DBLE(P(I,2))**2+
     +  DBLE(P(I,3))**2)
        P2=DBLE(P(I,4))**2-DBLE(P(I,1))**2-DBLE(P(I,2))**2-DBLE(P(I,3))
     +  **2
        IF(ABS(ENERGY-P(I,4))/(PSAVE(3,1,4)+PSAVE(3,2,4)).GT.PARU(11))
     +  THEN
          NUMMIS=NUMMIS+1
C...FOR TESTING PURPOSES
C       IF(LST(3).GE.1.AND.NUMMIS.LE.NWARN) THEN
C         WRITE(6,1000) I,(K(I,J),J=1,2),(P(I,J),J=1,5),
C    &    SIGN(SQRT(ABS(P2)),P2),ENERGY,INT(DARI29),NWARN
C         IF(ABS(P2-P(I,5)**2).GT.400.) CALL LULIST(2)
C       ENDIF
          GOTO 150
        ENDIF
        P(I,4)=ENERGY
  230 CONTINUE

      DARI30=DARI30+1.D0
      PARI(30)=DARI30
      IF(LST(23).EQ.2) PARL(24)=PARL(24)*DARI30/DARI29

      DO 240 I=1,N
        DO 240 J=1,5
  240 V(I,J)=0.
      IF(LST(7).EQ.1) THEN
        IF(LST(34).EQ.1) THEN
          K(4,1)=21
        ENDIF
* NEW PHYILOSOPHY:ONLY TAU AND CHARM DECAYS
        CALL LUEXEC
*        CALL LUSTRF(IP)
        IF(MSTU(24).NE.0) THEN
          WRITE(*,*) ' ERROR FROM JETSET, NEW EVENT MADE'
          GOTO 150
        ENDIF
      ENDIF

C     CALL GULIST(-1,2)
C...TRANSFORM TO DESIRED FRAME
C     LST(28)=1
      LST(29)=0
      PHIR=6.2832*RLU(0)
*      WRITE(*,*)'SYNC...',PHIR
      IF(LST(17).EQ.0) THEN
        IF(LST(5).GE.2) CALL LFRAME(LST(5),0)
C...RESTORE MOMENTA (E,P,BOSON,L) DUE TO NUMERICAL ERRORS FROM BOOSTS
        DO 250 I=1,4
          DO 250 J=1,5
  250   P(I,J)=PSAVE(LST(28),I,J)
        IF(LST(6).EQ.1.AND.LST(28).GE.2) THEN
C...RANDOM ROTATION IN AZIMUTHAL ANGLE
          CALL LUDBRB(0,0,0.,PHIR,0.D0,0.D0,0.D0)
          LST(29)=1
        ENDIF
      ELSE
        IF(LST(5).GE.2) THEN
*        WRITE(*,*)'DBETA,28,5',DBETA,LST(28),LST(5)
          IF (DBETA(1,3).NE.0) THEN
            CALL LFRAME(LST(5),LST(6))
          ELSE
            WRITE(*,*)'0 DBETA!!!'
            DO J=1,3
              DBETA(1,J)=DBETATMP(J)
            END DO
            CALL LFRAME(LST(5),LST(6))
            WRITE(*,*)'1ST ATTEMPT RECOVERY',PHIR
            CALL LULIST(1)
            WRITE(*,*)'1ST ATTEMPT RECOVERY ENDED',PHIR
          ENDIF
*         WRITE(*,*)'CALLED LFRAME',LST(5),LST(6),PHIR
          IF (ABS(P(2,3)+P(1,3)).LT.0.2) THEN
            WRITE(*,*)'LEPTO ERROR CMS FRAME!!',PHIR
*          CALL LULIST(1)
            DO J=1,3
              DBETA(1,J)=DBETATMP(J)
            END DO
            WRITE(*,*)' TRY RECOVERY...',DBETATMP
            LST(28)=2
            CALL LFRAME(LST(5),LST(6))
*          CALL LULIST(1)
          ENDIF
        ENDIF
      ENDIF
C...DEACTIVATE SCATTERED LEPTON
      IF(MOD(LST(4),10).EQ.0) K(4,1)=21
C     CALL GULIST(0,2)

      RETURN
10000 FORMAT(' WARNING: TOO LARGE NUMERICAL MISMATCH IN ',
     +'PARTICLE ENERGY-MOMENTUM-MASS',
     +/,3X,'I K(I,1) ..2)  P(I,1)  P(I,2)  P(I,3)',
     +'  P(I,4)  P(I,5)    MASS  ENERGY',/,I4,2I6,7F8.3,/,
     +' EVENT NO.',I8,' REGENERATED. ONLY FIRST',I5,' WARNINGS PRINTED')
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION AMAS4(PP)
C     ******************
C ----------------------------------------------------------------------
C CALCULATES MASS OF PP
C
C     USED BY :
C ----------------------------------------------------------------------
      REAL  PP(4)
      AAA=PP(4)**2-PP(3)**2-PP(2)**2-PP(1)**2
      IF(AAA.NE.0.0) AAA=AAA/SQRT(ABS(AAA))
      AMAS4=AAA
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION AMAST(PP)
C ----------------------------------------------------------------------
C CALCULATES MASS OF PP (DOUBLE PRECISION)
C
C     USED BY : RADKOR
C ----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  PP(4)
      AAA=PP(4)**2-PP(3)**2-PP(2)**2-PP(1)**2
C
      IF(AAA.NE.0.0) AAA=AAA/SQRT(ABS(AAA))
      AMAST=AAA
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION ANGFI(X,Y)
C ----------------------------------------------------------------------
* CALCULATES ANGLE IN (0,2*PI) RANGE OUT OF X-Y
C
C     USED BY : KORALZ RADKOR
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA PI /3.141592653589793238462643D0/
C
      IF(ABS(Y).LT.ABS(X)) THEN
        THE=ATAN(ABS(Y/X))
        IF(X.LE.0D0) THE=PI-THE
      ELSE
        THE=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      IF(Y.LT.0D0) THE=2D0*PI-THE
      ANGFI=THE
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION ANGXY(X,Y)
C ----------------------------------------------------------------------
C
C     USED BY : KORALZ RADKOR
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA PI /3.141592653589793238462643D0/
C
      IF(ABS(Y).LT.ABS(X)) THEN
        THE=ATAN(ABS(Y/X))
        IF(X.LE.0D0) THE=PI-THE
      ELSE
        THE=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      ANGXY=THE
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE BOSTD3(EXE,PVEC,QVEC)
C ----------------------------------------------------------------------
C BOOST ALONG Z AXIS, EXE=EXP(ETA), ETA= HIPERBOLIC VELOCITY.
C
C     USED BY : KORALZ RADKOR
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PVEC(4),QVEC(4),RVEC(4)
C
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      RPL=RVEC(4)+RVEC(3)
      RMI=RVEC(4)-RVEC(3)
      QPL=RPL*EXE
      QMI=RMI/EXE
      QVEC(1)=RVEC(1)
      QVEC(2)=RVEC(2)
      QVEC(3)=(QPL-QMI)/2
      QVEC(4)=(QPL+QMI)/2
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE BOSTR3(EXE,PVEC,QVEC)
C ----------------------------------------------------------------------
C BOOST ALONG Z AXIS, EXE=EXP(ETA), ETA= HIPERBOLIC VELOCITY.
C
C     USED BY : TAUOLA KORALZ (?)
C ----------------------------------------------------------------------
      REAL*4 PVEC(4),QVEC(4),RVEC(4)
C
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      RPL=RVEC(4)+RVEC(3)
      RMI=RVEC(4)-RVEC(3)
      QPL=RPL*EXE
      QMI=RMI/EXE
      QVEC(1)=RVEC(1)
      QVEC(2)=RVEC(2)
      QVEC(3)=(QPL-QMI)/2
      QVEC(4)=(QPL+QMI)/2
      END
*CMZ :  1.01/39 03/11/95  16.36.14  by  Piero Zucchelli
*CMZ :  1.01/04 15/08/95  16.55.04  by  Stefania RICCIARDI
*CMZ :  1.01/02 02/08/95  12.14.00  by  Stefania RICCIARDI
*CMZ :  1.00/03 01/07/95  16.21.55  by  Stefania RICCIARDI
*-- Author :    Stefania RICCIARDI   26/06/95
      SUBROUTINE BTOCHO2(VIN,PIN,PTX,PTY)

      REAL VIN(3),PIN(3)
      DATA ZBEAM/82342/

      DZ = ZBEAM-VIN(3)
      PTY = VIN(2)+PIN(2)/PIN(3)*DZ
      PTX = VIN(1)+PIN(1)/PIN(3)*DZ
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      COMPLEX FUNCTION BWIG(S,M,G)
C **********************************************************
C     P-WAVE BREIT-WIGNER  FOR RHO
C **********************************************************
      REAL S,M,G
      REAL PI,PIM,QS,QM,W,GS
      DATA INIT /0/
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0) THEN
        INIT=1
        PI=3.141592654
        PIM=.139
C -------  BREIT-WIGNER -----------------------
      ENDIF
      IF (S.GT.4.*PIM**2) THEN
        QS=SQRT(ABS(ABS(S/4.-PIM**2)+(S/4.-PIM**2))/2.0)
        QM=SQRT(M**2/4.-PIM**2)
        W=SQRT(S)
        GS=G*(M/W)*(QS/QM)**3
      ELSE
        GS=0.0
      ENDIF
      BWIG=M**2/CMPLX(M**2-S,-M*GS)
      RETURN
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      COMPLEX FUNCTION BWIGM(S,M,G,XM1,XM2)
C **********************************************************
C     P-WAVE BREIT-WIGNER  FOR RHO
C **********************************************************
      REAL S,M,G,XM1,XM2
      REAL PI,QS,QM,W,GS
      DATA INIT /0/
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0) THEN
        INIT=1
        PI=3.141592654
C -------  BREIT-WIGNER -----------------------
      ENDIF
      IF (S.GT.(XM1+XM2)**2) THEN
        QS=SQRT(ABS((S -(XM1+XM2)**2)*(S -(XM1-XM2)**2)))/SQRT(S)
        QM=SQRT(ABS((M**2-(XM1+XM2)**2)*(M**2-(XM1-XM2)**2)))/M
        W=SQRT(S)
        GS=G*(M/W)**2*(QS/QM)**3
      ELSE
        GS=0.0
      ENDIF
      BWIGM=M**2/CMPLX(M**2-S,-SQRT(S)*GS)
      RETURN
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      COMPLEX FUNCTION BWIGS(S,M,G)
C **********************************************************
C     P-WAVE BREIT-WIGNER  FOR K*
C **********************************************************
      REAL S,M,G
      REAL PI,PIM,QS,QM,W,GS,MK
      DATA INIT /0/
      P(A,B,C)=SQRT(ABS(ABS(((A+B-C)**2-4.*A*B)/4./A)
     +                    +(((A+B-C)**2-4.*A*B)/4./A))/2.0)
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0) THEN
        INIT=1
        PI=3.141592654
        PIM=.139
        MK=.493667
C -------  BREIT-WIGNER -----------------------
      ENDIF
      QS=P(S,PIM**2,MK**2)
      QM=P(M**2,PIM**2,MK**2)
      W=SQRT(S)
      GS=G*(M/W)*(QS/QM)**3
      BWIGS=M**2/CMPLX(M**2-S,-M*GS)
      RETURN
      END
*CMZ :  1.01/24 29/05/95  15.40.09  BY  PIERO ZUCCHELLI
*CMZ :  1.01/22 27/05/95  16.39.05  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.40.52  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 20/09/94  14.44.05  BY  PIERO ZUCCHELLI
*-- AUTHOR :    PIERO ZUCCHELLI   20/09/94

      SUBROUTINE CATS
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LMINUI/ XKIN(4),UKIN(4),WKIN(4),AIN(4),BIN(4), MAXFIN,
     +RELUP,RELERR,RELER2,FCNMAX
* DO ESPLICITELY WWHAT IS DONE IN LEPTOD AND USER CUTS
      COMMON/LINPATCH/NCALLS,NCALL
      REAL*4 SAVECUT(14)

      IF (IMYFIRST.EQ.0) THEN
        IMYFIRST=1

        DO I=1,14
          SAVECUT(I)=CUT(I)
        END DO
      ENDIF

      DO I=1,14
        CUT(I)=SAVECUT(I)
      END DO
      DO I=1,4
        XKIN(I)=I
        UKIN(I)=0
        WKIN(I)=0
        AIN(I)=0
        BIN(I)=0
      END DO

      FCNMAX=0
      NCALLS=0
      NCALL=0

* USER CATS

      CUT(7)=0
      CUT(5)=0

      RETURN
      END
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 09/08/94  17.43.59  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE CHOICE(MNUM,RR,ICHAN,PROB1,PROB2,PROB3,
     +            AMRX,GAMRX,AMRA,GAMRA,AMRB,GAMRB)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      AMROP=1.1
      GAMROP=0.36
      AMOM=.782
      GAMOM=0.0084
C     XXXXA CORRESPOND TO S2 CHANNEL !
      IF(MNUM.EQ.0) THEN
        PROB1=0.5
        PROB2=0.5
        AMRX =AMA1
        GAMRX=GAMA1
        AMRA =AMRO
        GAMRA=GAMRO
        AMRB =AMRO
        GAMRB=GAMRO
      ELSEIF(MNUM.EQ.1) THEN
        PROB1=0.5
        PROB2=0.5
        AMRX =1.57
        GAMRX=0.9
        AMRB =AMKST
        GAMRB=GAMKST
        AMRA =AMRO
        GAMRA=GAMRO
      ELSEIF(MNUM.EQ.2) THEN
        PROB1=0.5
        PROB2=0.5
        AMRX =1.57
        GAMRX=0.9
        AMRB =AMKST
        GAMRB=GAMKST
        AMRA =AMRO
        GAMRA=GAMRO
      ELSEIF(MNUM.EQ.3) THEN
        PROB1=0.5
        PROB2=0.5
        AMRX =1.27
        GAMRX=0.3
        AMRA =AMKST
        GAMRA=GAMKST
        AMRB =AMKST
        GAMRB=GAMKST
      ELSEIF(MNUM.EQ.4) THEN
        PROB1=0.5
        PROB2=0.5
        AMRX =1.27
        GAMRX=0.3
        AMRA =AMKST
        GAMRA=GAMKST
        AMRB =AMKST
        GAMRB=GAMKST
      ELSEIF(MNUM.EQ.5) THEN
        PROB1=0.5
        PROB2=0.5
        AMRX =1.27
        GAMRX=0.3
        AMRA =AMKST
        GAMRA=GAMKST
        AMRB =AMRO
        GAMRB=GAMRO
      ELSEIF(MNUM.EQ.6) THEN
        PROB1=0.4
        PROB2=0.4
        AMRX =1.27
        GAMRX=0.3
        AMRA =AMRO
        GAMRA=GAMRO
        AMRB =AMKST
        GAMRB=GAMKST
      ELSEIF(MNUM.EQ.7) THEN
        PROB1=0.0
        PROB2=1.0
        AMRX =1.27
        GAMRX=0.9
        AMRA =AMRO
        GAMRA=GAMRO
        AMRB =AMRO
        GAMRB=GAMRO
      ELSEIF(MNUM.EQ.8) THEN
        PROB1=0.0
        PROB2=1.0
        AMRX =AMROP
        GAMRX=GAMROP
        AMRB =AMOM
        GAMRB=GAMOM
        AMRA =AMRO
        GAMRA=GAMRO
      ELSEIF(MNUM.EQ.101) THEN
        PROB1=.35
        PROB2=.35
        AMRX =1.2
        GAMRX=.46
        AMRB =AMOM
        GAMRB=GAMOM
        AMRA =AMOM
        GAMRA=GAMOM
      ELSEIF(MNUM.EQ.102) THEN
        PROB1=0.0
        PROB2=0.0
        AMRX =1.4
        GAMRX=.6
        AMRB =AMOM
        GAMRB=GAMOM
        AMRA =AMOM
        GAMRA=GAMOM
      ELSE
        PROB1=0.0
        PROB2=0.0
        AMRX =AMA1
        GAMRX=GAMA1
        AMRA =AMRO
        GAMRA=GAMRO
        AMRB =AMRO
        GAMRB=GAMRO
      ENDIF
C
      IF    (RR.LE.PROB1) THEN
        ICHAN=1
      ELSEIF(RR.LE.(PROB1+PROB2)) THEN
        ICHAN=2
        AX   =AMRA
        GX   =GAMRA
        AMRA =AMRB
        GAMRA=GAMRB
        AMRB =AX
        GAMRB=GX
        PX   =PROB1
        PROB1=PROB2
        PROB2=PX
      ELSE
        ICHAN=3
      ENDIF
C
      PROB3=1.0-PROB1-PROB2
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE CLAXI(HJ,PN,PIA)
C ----------------------------------------------------------------------
* CALCULATES THE "AXIAL TYPE"  PI-VECTOR  PIA
* NOTE THAT THE NEUTRINO MOM. PN IS ASSUMED TO BE ALONG Z-AXIS
C SIGN IS CHOSEN +/- FOR DECAY OF TAU +/- RESPECTIVELY
C     CALLED BY : DAMPAA, CLNUT
C ----------------------------------------------------------------------
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / IDFC  / IDFF
      REAL PIA(4),PN(4)
      COMPLEX HJ(4),HJC(4)
C     DET2(I,J)=AIMAG(HJ(I)*HJC(J)-HJ(J)*HJC(I))
C -- HERE WAS AN ERROR (ZW, 21.11.1991)
      DET2(I,J)=AIMAG(HJC(I)*HJ(J)-HJC(J)*HJ(I))
C -- IT WAS AFFECTING SIGN OF A_LR ASYMMETRY IN A1 DECAY.
C -- NOTE ALSO COLLISION OF NOTATION OF GAMMA_VA AS DEFINED IN
C -- TAUOLA PAPER AND J.H. KUHN AND SANTAMARIA Z. PHYS C 48 (1990) 445
* -----------------------------------
      IF     (KTOM.EQ.1.OR.KTOM.EQ.-1) THEN
        SIGN= IDFF/ABS(IDFF)
      ELSEIF (KTOM.EQ.2) THEN
        SIGN=-IDFF/ABS(IDFF)
      ELSE
        PRINT *, 'STOP IN CLAXI: KTOM=',KTOM
        STOP
      ENDIF
C
      DO 10 I=1,4
 10   HJC(I)=CONJG(HJ(I))
      PIA(1)= -2.*PN(3)*DET2(2,4)+2.*PN(4)*DET2(2,3)
      PIA(2)= -2.*PN(4)*DET2(1,3)+2.*PN(3)*DET2(1,4)
      PIA(3)=  2.*PN(4)*DET2(1,2)
      PIA(4)=  2.*PN(3)*DET2(1,2)
C ALL FOUR INDICES ARE UP SO  PIA(3) AND PIA(4) HAVE SAME SIGN
      DO 20 I=1,4
  20  PIA(I)=PIA(I)*SIGN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE CLNUT(HJ,B,HV)
C ----------------------------------------------------------------------
* CALCULATES THE CONTRIBUTION BY NEUTRINO MASS
* NOTE THE TAU IS ASSUMED TO BE AT REST
C
C     CALLED BY : DAMPAA
C ----------------------------------------------------------------------
      COMPLEX HJ(4)
      REAL HV(4),P(4)
      DATA P /3*0.,1.0/
C
      CALL CLAXI(HJ,P,HV)
      B=REAL( HJ(4)*AIMAG(HJ(4)) - HJ(3)*AIMAG(HJ(3))
     &      - HJ(2)*AIMAG(HJ(2)) - HJ(1)*AIMAG(HJ(1))  )
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE CLVEC(HJ,PN,PIV)
C ----------------------------------------------------------------------
* CALCULATES THE "VECTOR TYPE"  PI-VECTOR  PIV
* NOTE THAT THE NEUTRINO MOM. PN IS ASSUMED TO BE ALONG Z-AXIS
C
C     CALLED BY : DAMPAA
C ----------------------------------------------------------------------
      REAL PIV(4),PN(4)
      COMPLEX HJ(4),HN
C
      HN= HJ(4)*CMPLX(PN(4))-HJ(3)*CMPLX(PN(3))
      HH= REAL(HJ(4)*CONJG(HJ(4))-HJ(3)*CONJG(HJ(3))
     $        -HJ(2)*CONJG(HJ(2))-HJ(1)*CONJG(HJ(1)))
      DO 10 I=1,4
   10 PIV(I)=4.*REAL(HN*CONJG(HJ(I)))-2.*HH*PN(I)
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE CURR(MNUM,PIM1,PIM2,PIM3,PIM4,HADCUR)
C     ==================================================================
C     HADRONIC CURRENT FOR 4 PI FINAL STATE
C     R. FISHER, J. WESS AND F. WAGNER Z. PHYS C3 (1980) 313
C     R. DECKER Z. PHYS C36 (1987) 487.
C     M. GELL-MANN, D. SHARP, W. WAGNER PHYS. REV. LETT 8 (1962) 261.
C     ==================================================================

      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
C ARBITRARY FIXING OF THE FOUR PI X-SECTION NORMALIZATION
      COMMON /ARBIT/ ARFLAT,AROMEG
      REAL  PIM1(4),PIM2(4),PIM3(4),PIM4(4),PAA(4)
      COMPLEX HADCUR(4),FORM1,FORM2,FORM3,FPIKM
      COMPLEX BWIGN
      REAL PA(4),PB(4)
      REAL AA(4,4),PP(4,4)
      DATA PI /3.141592653589793238462643/
      DATA  FPI /93.3E-3/
      BWIGN(A,XM,XG)=1.0/CMPLX(A-XM**2,XM*XG)
C
C --- MASSES AND CONSTANTS
      G1=12.924
      G2=1475.98
      G =G1*G2
      ELPHA=-.1
      AMROP=1.7
      GAMROP=0.26
      AMOM=.782
      GAMOM=0.0085
      ARFLAT=1.0
      AROMEG=1.0
C
      FRO=0.266*AMRO**2
      COEF1=2.0*SQRT(3.0)/FPI**2*ARFLAT
      COEF2=FRO*G*AROMEG
C --- INITIALIZATION OF FOUR VECTORS
      DO 20 K=1,4
        DO 10 L=1,4
   10   AA(K,L)=0.0
        HADCUR(K)=CMPLX(0.0)
        PAA(K)=PIM1(K)+PIM2(K)+PIM3(K)+PIM4(K)
        PP(1,K)=PIM1(K)
        PP(2,K)=PIM2(K)
        PP(3,K)=PIM3(K)
   20 PP(4,K)=PIM4(K)
C
      IF (MNUM.EQ.1) THEN
C ===================================================================
C PI- PI- P0 PI+ CASE                                            ====
C ===================================================================
        QQ=PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2
C --- LOOP OVER THRE CONTRIBUTION OF THE NON-OMEGA CURRENT
        DO 80  K=1,3
          SK=(PP(K,4)+PIM4(4))**2-(PP(K,3)+PIM4(3))**2 -(PP(K,2)+
     +    PIM4(2))**2-(PP(K,1)+PIM4(1))**2
C -- DEFINITION OF AA MATRIX
C -- CRONECKER DELTA
          DO 40  I=1,4
            DO 30  J=1,4
   30       AA(I,J)=0.0
   40     AA(I,I)=1.0
C ... AND THE REST ...
          DO 60  L=1,3
            IF (L.NE.K) THEN
              DENOM=(PAA(4)-PP(L,4))**2-(PAA(3)-PP(L,3))**2 -(PAA(2)-
     +        PP(L,2))**2-(PAA(1)-PP(L,1))**2
              DO 50  I=1,4
                DO 50  J=1,4
                  SIG= 1.0
                  IF(J.NE.4) SIG=-SIG
                  AA(I,J)=AA(I,J) -SIG*(PAA(I)-2.0*PP(L,I))*(PAA(J)-
     +            PP(L,J))/DENOM
   50         CONTINUE
            ENDIF
   60     CONTINUE
C --- LET'S ADD SOMETHING TO HADCURR
          FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKM(SQRT(QQ),AMPI,AMPI)
C       FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKMD(SQRT(QQ),AMPI,AMPI)
CCCCCCCCCCCCCCCCC       FORM1=WIGFOR(SK,AMRO,GAMRO)      (TESTS)
C
          FIX=1.0
          IF (K.EQ.3) FIX=-2.0
          DO 70  I=1,4
            DO 70  J=1,4
              HADCUR(I)= HADCUR(I)+CMPLX(FIX*COEF1)*FORM1*AA(I,J)*
     +        (PP(K,J)-PP(4,J))
   70     CONTINUE
C --- END OF THE NON OMEGA CURRENT (3 POSSIBILITIES)
   80   CONTINUE
C
C
C --- THERE ARE TWO POSSIBILITIES FOR OMEGA CURRENT
C --- PA PB ARE CORRESPONDING FIRST AND SECOND PI-'S
        DO 120 KK=1,2
          DO 90  I=1,4
            PA(I)=PP(KK,I)
            PB(I)=PP(3-KK,I)
   90     CONTINUE
C --- LORENTZ INVARIANTS
          QQA=0.0
          SS23=0.0
          SS24=0.0
          SS34=0.0
          QP1P2=0.0
          QP1P3=0.0
          QP1P4=0.0
          P1P2 =0.0
          P1P3 =0.0
          P1P4 =0.0
          DO 100 K=1,4
            SIGN=-1.0
            IF (K.EQ.4) SIGN= 1.0
            QQA=QQA+SIGN*(PAA(K)-PA(K))**2
            SS23=SS23+SIGN*(PB(K) +PIM3(K))**2
            SS24=SS24+SIGN*(PB(K) +PIM4(K))**2
            SS34=SS34+SIGN*(PIM3(K)+PIM4(K))**2
            QP1P2=QP1P2+SIGN*(PAA(K)-PA(K))*PB(K)
            QP1P3=QP1P3+SIGN*(PAA(K)-PA(K))*PIM3(K)
            QP1P4=QP1P4+SIGN*(PAA(K)-PA(K))*PIM4(K)
            P1P2=P1P2+SIGN*PA(K)*PB(K)
            P1P3=P1P3+SIGN*PA(K)*PIM3(K)
            P1P4=P1P4+SIGN*PA(K)*PIM4(K)
  100     CONTINUE
C
          FORM2=COEF2*(BWIGN(QQ,AMRO,GAMRO)+ELPHA*BWIGN(QQ,AMROP,
     +    GAMROP))
C        FORM3=BWIGN(QQA,AMOM,GAMOM)*(BWIGN(SS23,AMRO,GAMRO)+
C     $        BWIGN(SS24,AMRO,GAMRO)+BWIGN(SS34,AMRO,GAMRO))
          FORM3=BWIGN(QQA,AMOM,GAMOM)
C
          DO 110 K=1,4
            HADCUR(K)=HADCUR(K)+FORM2*FORM3*( PB (K)*(QP1P3*P1P4-QP1P4*
     +      P1P3) +PIM3(K)*(QP1P4*P1P2-QP1P2*P1P4) +PIM4(K)*(QP1P2*
     +      P1P3-QP1P3*P1P2) )
  110     CONTINUE
  120   CONTINUE
C
      ELSE
C ===================================================================
C PI0 PI0 P0 PI- CASE                                            ====
C ===================================================================
        QQ=PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2
        DO 180 K=1,3
C --- LOOP OVER THRE CONTRIBUTION OF THE NON-OMEGA CURRENT
          SK=(PP(K,4)+PIM4(4))**2-(PP(K,3)+PIM4(3))**2 -(PP(K,2)+
     +    PIM4(2))**2-(PP(K,1)+PIM4(1))**2
C -- DEFINITION OF AA MATRIX
C -- CRONECKER DELTA
          DO 140 I=1,4
            DO 130 J=1,4
  130       AA(I,J)=0.0
  140     AA(I,I)=1.0
C
C ... AND THE REST ...
          DO 160 L=1,3
            IF (L.NE.K) THEN
              DENOM=(PAA(4)-PP(L,4))**2-(PAA(3)-PP(L,3))**2 -(PAA(2)-
     +        PP(L,2))**2-(PAA(1)-PP(L,1))**2
              DO 150 I=1,4
                DO 150 J=1,4
                  SIG=1.0
                  IF(J.NE.4) SIG=-SIG
                  AA(I,J)=AA(I,J) -SIG*(PAA(I)-2.0*PP(L,I))*(PAA(J)-
     +            PP(L,J))/DENOM
  150         CONTINUE
            ENDIF
  160     CONTINUE
C --- LET'S ADD SOMETHING TO HADCURR
          FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKM(SQRT(QQ),AMPI,AMPI)
C       FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKMD(SQRT(QQ),AMPI,AMPI)
CCCCCCCCCCCCC       FORM1=WIGFOR(SK,AMRO,GAMRO)        (TESTS)
          DO 170 I=1,4
            DO 170 J=1,4
              HADCUR(I)= HADCUR(I)+CMPLX(COEF1)*FORM1*AA(I,J)*(PP(K,J)-
     +        PP(4,J))
  170     CONTINUE
C --- END OF THE NON OMEGA CURRENT (3 POSSIBILITIES)
  180   CONTINUE
      ENDIF
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DADMAA(MODE,ISGN,HHV,PNU,PAA,PIM1,PIM2,PIPL,JAA)
C ----------------------------------------------------------------------
* A1 DECAY UNWEIGHTED EVENTS
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      REAL  HHV(4)
      REAL  HV(4),PAA(4),PNU(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4),PDUM5(4)
      REAL*4 RRR(3)
      REAL*8 SWT, SSWT
      DATA PI /3.141592653589793238462643/
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        NEVRAW=0
        NEVACC=0
        NEVOVR=0
        SWT=0
        SSWT=0
        WTMAX=1E-20
        DO 10 I=1,500
          CALL DPHSAA(WT,HV,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5,JAA)
          IF(WT.GT.WTMAX/1.2) WTMAX=WT*1.2
   10   CONTINUE
CC      CALL HBOOK1(801,'WEIGHT DISTRIBUTION  DADMAA    $',100,0,2)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
   20   CONTINUE
        IF(IWARM.EQ.0) GOTO 40
        CALL DPHSAA(WT,HV,PNU,PAA,PIM1,PIM2,PIPL,JAA)
CC      CALL HFILL(801,WT/WTMAX)
        NEVRAW=NEVRAW+1
        SWT=SWT+WT
        SSWT=SSWT+WT**2
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX) NEVOVR=NEVOVR+1
        IF(RN*WTMAX.GT.WT) GOTO 20
C ROTATIONS TO BASIC TAU REST FRAME
        COSTHE=-1.+2.*RRR(2)
        THET=ACOS(COSTHE)
        PHI =2*PI*RRR(3)
        CALL ROTPOL(THET,PHI,PNU)
        CALL ROTPOL(THET,PHI,PAA)
        CALL ROTPOL(THET,PHI,PIM1)
        CALL ROTPOL(THET,PHI,PIM2)
        CALL ROTPOL(THET,PHI,PIPL)
        CALL ROTPOL(THET,PHI,HV)
        DO 30 I=1,3
   30   HHV(I)=-ISGN*HV(I)
        NEVACC=NEVACC+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVRAW.EQ.0) RETURN
        PARGAM=SWT/FLOAT(NEVRAW+1)
        ERROR=0
        IF(NEVRAW.NE.0) ERROR=SQRT(SSWT/SWT**2-1./FLOAT(NEVRAW))
        RAT=PARGAM/GAMEL
        WRITE(IOUT, 10100) NEVRAW,NEVACC,NEVOVR,PARGAM,RAT,ERROR
CC      CALL HPRINT(801)
        GAMPMC(5)=RAT
        GAMPER(5)=ERROR
CAM     NEVDEC(5)=NEVACC
      ENDIF
C     =====
      RETURN
10000 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'******** DADMAA INITIALISATION ********',9X,1H*
     + /,' *',E20.5,5X,'WTMAX  = MAXIMUM WEIGHT                ',9X,1H*
     +  /,1X,15(5H*****)/)
10100 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'******** DADMAA FINAL REPORT  ******** ',9X,1H*
     + /,' *',I20  ,5X,'NEVRAW = NO. OF A1  DECAYS TOTAL       ',9X,1H*
     + /,' *',I20  ,5X,'NEVACC = NO. OF A1   DECS. ACCEPTED    ',9X,1H*
     + /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
     + /,' *',E20.5,5X,'PARTIAL WTDTH (A1  DECAY) IN GEV UNITS ',9X,1H*
     + /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
     + /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
     +  /,1X,15(5H*****)/)
   40 WRITE(IOUT, 10200)
10200 FORMAT(' ----- DADMAA: LACK OF INITIALISATION')
      STOP
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DADMEL(MODE,ISGN,HHV,PNU,PWB,Q1,Q2,PHX)
C ----------------------------------------------------------------------
C
C     CALLED BY : DEXEL,(DEKAY,DEKAY1)
C ----------------------------------------------------------------------
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      REAL*4         PHX(4)
      COMMON / INOUT / INUT,IOUT
      REAL  HHV(4),HV(4),PWB(4),PNU(4),Q1(4),Q2(4)
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4),PDUM5(4)
      REAL*4 RRR(3)
      REAL*8 SWT, SSWT
      DATA PI /3.141592653589793238462643/
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        NEVRAW=0
        NEVACC=0
        NEVOVR=0
        SWT=0
        SSWT=0
        WTMAX=1E-20
        DO 10 I=1,500
          CALL DPHSEL(WT,HV,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          IF(WT.GT.WTMAX/1.2) WTMAX=WT*1.2
   10   CONTINUE
CC      CALL HBOOK1(803,'WEIGHT DISTRIBUTION  DADMEL    $',100,0,2)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
   20   CONTINUE
        IF(IWARM.EQ.0) GOTO 40
        NEVRAW=NEVRAW+1
        CALL DPHSEL(WT,HV,PNU,PWB,Q1,Q2,PHX)
CC      CALL HFILL(803,WT/WTMAX)
        SWT=SWT+WT
        SSWT=SSWT+WT**2
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX) NEVOVR=NEVOVR+1
        IF(RN*WTMAX.GT.WT) GOTO 20
C ROTATIONS TO BASIC TAU REST FRAME
        RR2=RRR(2)
        COSTHE=-1.+2.*RR2
        THET=ACOS(COSTHE)
        RR3=RRR(3)
        PHI =2*PI*RR3
        CALL ROTOR2(THET,PNU,PNU)
        CALL ROTOR3( PHI,PNU,PNU)
        CALL ROTOR2(THET,PWB,PWB)
        CALL ROTOR3( PHI,PWB,PWB)
        CALL ROTOR2(THET,Q1,Q1)
        CALL ROTOR3( PHI,Q1,Q1)
        CALL ROTOR2(THET,Q2,Q2)
        CALL ROTOR3( PHI,Q2,Q2)
        CALL ROTOR2(THET,HV,HV)
        CALL ROTOR3( PHI,HV,HV)
        CALL ROTOR2(THET,PHX,PHX)
        CALL ROTOR3( PHI,PHX,PHX)
        DO 30 ,I=1,3
   30   HHV(I)=-ISGN*HV(I)
        NEVACC=NEVACC+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVRAW.EQ.0) RETURN
        PARGAM=SWT/FLOAT(NEVRAW+1)
        ERROR=0
        IF(NEVRAW.NE.0) ERROR=SQRT(SSWT/SWT**2-1./FLOAT(NEVRAW))
        RAT=PARGAM/GAMEL
        WRITE(IOUT, 10000) NEVRAW,NEVACC,NEVOVR,PARGAM,RAT,ERROR
CC      CALL HPRINT(803)
        GAMPMC(1)=RAT
        GAMPER(1)=ERROR
CAM     NEVDEC(1)=NEVACC
      ENDIF
C     =====
      RETURN
10000 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'******** DADMEL FINAL REPORT  ******** ',9X,1H*
     + /,' *',I20  ,5X,'NEVRAW = NO. OF EL  DECAYS TOTAL       ',9X,1H*
     + /,' *',I20  ,5X,'NEVACC = NO. OF EL   DECS. ACCEPTED    ',9X,1H*
     + /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
     + /,' *',E20.5,5X,'PARTIAL WTDTH ( ELECTRON) IN GEV UNITS ',9X,1H*
     + /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
     + /,' *',F20.9,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
     + /,' *',25X,     'COMPLETE QED CORRECTIONS INCLUDED      ',9X,1H*
     + /,' *',25X,     'BUT ONLY V-A CUPLINGS                  ',9X,1H*
     +  /,1X,15(5H*****)/)
   40 WRITE(IOUT, 10100)
10100 FORMAT(' ----- DADMEL: LACK OF INITIALISATION')
      STOP
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DADMKK(MODE,ISGN,HV,PKK,PNU)
C ----------------------------------------------------------------------
C FZ
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      REAL  PKK(4),PNU(4),HV(4)
      DATA PI /3.141592653589793238462643/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        NEVTOT=0
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
        NEVTOT=NEVTOT+1
        EKK= (AMTAU**2+AMK**2-AMNUTA**2)/(2*AMTAU)
        ENU= (AMTAU**2-AMK**2+AMNUTA**2)/(2*AMTAU)
        XKK= SQRT(EKK**2-AMK**2)
C K MOMENTUM
        CALL SPHERA(XKK,PKK)
        PKK(4)=EKK
C TAU-NEUTRINO MOMENTUM
        DO 10 I=1,3
   10   PNU(I)=-PKK(I)
        PNU(4)=ENU
        PXQ=AMTAU*EKK
        PXN=AMTAU*ENU
        QXN=PKK(4)*PNU(4)-PKK(1)*PNU(1)-PKK(2)*PNU(2)-PKK(3)*PNU(3)
        BRAK=(GV**2+GA**2)*(2*PXQ*QXN-AMK**2*PXN)
     &      +(GV**2-GA**2)*AMTAU*AMNUTA*AMK**2
        DO 20 I=1,3
   20   HV(I)=-ISGN*2*GA*GV*AMTAU*(2*PKK(I)*QXN-PNU(I)*AMK**2)/BRAK
        HV(4)=1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVTOT.EQ.0) RETURN
        FKK=0.0354
CFZ THERE WAS BRAK/AMTAU**4 BEFORE
C        GAMM=(GFERMI*FKK)**2/(16.*PI)*AMTAU**3*
C     *       (BRAK/AMTAU**4)**2
CZW 7.02.93 HERE WAS AN ERROR AFFECTING NON STANDARD MODEL
C       CONFIGURATIONS ONLY
        GAMM=(GFERMI*FKK)**2/(16.*PI)*AMTAU**3*
     $       (BRAK/AMTAU**4)*
     $       SQRT((AMTAU**2-AMK**2-AMNUTA**2)**2
     $            -4*AMK**2*AMNUTA**2           )/AMTAU**2
        ERROR=0

        ERROR=0
        RAT=GAMM/GAMEL
        WRITE(IOUT, 10000) NEVTOT,GAMM,RAT,ERROR
        GAMPMC(6)=RAT
        GAMPER(6)=ERROR
CAM     NEVDEC(6)=NEVTOT
      ENDIF
C     =====
      RETURN
10000 FORMAT(///1X,15(5H*****)
     $ /,' *',     25X,'******** DADMKK FINAL REPORT   ********',9X,1H*
     $ /,' *',I20  ,5X,'NEVTOT = NO. OF K  DECAYS TOTAL        ',9X,1H*,
     $ /,' *',E20.5,5X,'PARTIAL WTDTH ( K DECAY) IN GEV UNITS  ',9X,1H*,
     $ /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
     $ /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH (STAT.)',9X,1H*
     $  /,1X,15(5H*****)/)
      END
*CMZ :  1.01/50 19/04/96  09.49.04  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DADMKS(MODE,ISGN,HHV,PNU,PKS,PKK,PPI,JKST)
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      COMMON / INOUT / INUT,IOUT
      REAL  HHV(4)
      REAL  HV(4),PKS(4),PNU(4),PKK(4),PPI(4)
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4)
      REAL*4 RRR(3),RTEMP(1)
      REAL*8 SWT, SSWT
      DATA PI /3.141592653589793238462643/
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        NEVRAW=0
        NEVACC=0
        NEVOVR=0
        SWT=0
        SSWT=0
        WTMAX=1E-20
        DO 10 I=1,500
C THE INITIALISATION IS DONE WITH THE 66.7% MODE
          JKST=10
          CALL DPHSKS(WT,HV,PDUM1,PDUM2,PDUM3,PDUM4,JKST)
          IF(WT.GT.WTMAX/1.2) WTMAX=WT*1.2
   10   CONTINUE
CC      CALL HBOOK1(801,'WEIGHT DISTRIBUTION  DADMKS    $',100,0,2)
CC      PRINT 7003,WTMAX
CC      CALL HBOOK1(112,'-------- K* MASS -------- $',100,0.,2.)
      ELSEIF(MODE.EQ. 0) THEN
C     =====================================
        IF(IWARM.EQ.0) GOTO 40
C  HERE WE CHOOSE RANDOMLY BETWEEN K0 PI+_ (66.7%)
C  AND K+_ PI0 (33.3%)
        DEC1=BRKS
   20   CONTINUE
        RTEMP(1)=RMOD
        CALL RANMAR(RTEMP,1)
        RMOD=RTEMP(1)
        IF(RMOD.LT.DEC1) THEN
          JKST=10
        ELSE
          JKST=20
        ENDIF
        CALL DPHSKS(WT,HV,PNU,PKS,PKK,PPI,JKST)
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX) NEVOVR=NEVOVR+1
        NEVRAW=NEVRAW+1
        SWT=SWT+WT
        SSWT=SSWT+WT**2
        IF(RN*WTMAX.GT.WT) GOTO 20
C ROTATIONS TO BASIC TAU REST FRAME
        COSTHE=-1.+2.*RRR(2)
        THET=ACOS(COSTHE)
        PHI =2*PI*RRR(3)
        CALL ROTOR2(THET,PNU,PNU)
        CALL ROTOR3( PHI,PNU,PNU)
        CALL ROTOR2(THET,PKS,PKS)
        CALL ROTOR3( PHI,PKS,PKS)
        CALL ROTOR2(THET,PKK,PKK)
        CALL ROTOR3(PHI,PKK,PKK)
        CALL ROTOR2(THET,PPI,PPI)
        CALL ROTOR3( PHI,PPI,PPI)
        CALL ROTOR2(THET,HV,HV)
        CALL ROTOR3( PHI,HV,HV)
        DO 30 I=1,3
   30   HHV(I)=-ISGN*HV(I)
        NEVACC=NEVACC+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVRAW.EQ.0) RETURN
        PARGAM=SWT/FLOAT(NEVRAW+1)
        ERROR=0
        IF(NEVRAW.NE.0) ERROR=SQRT(SSWT/SWT**2-1./FLOAT(NEVRAW))
        RAT=PARGAM/GAMEL
        WRITE(IOUT, 10100) NEVRAW,NEVACC,NEVOVR,PARGAM,RAT,ERROR
CC      CALL HPRINT(801)
        GAMPMC(7)=RAT
        GAMPER(7)=ERROR
CAM     NEVDEC(7)=NEVACC
      ENDIF
C     =====
      RETURN
10000 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'******** DADMKS INITIALISATION ********',9X,1H*
     + /,' *',E20.5,5X,'WTMAX  = MAXIMUM WEIGHT                ',9X,1H*
     +  /,1X,15(5H*****)/)
10100 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'******** DADMKS FINAL REPORT   ********',9X,1H*
     + /,' *',I20  ,5X,'NEVRAW = NO. OF K* DECAYS TOTAL        ',9X,1H*,
     + /,' *',I20  ,5X,'NEVACC = NO. OF K*  DECS. ACCEPTED     ',9X,1H*,
     + /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
     + /,' *',E20.5,5X,'PARTIAL WTDTH (K* DECAY) IN GEV UNITS  ',9X,1H*,
     + /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
     + /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
     +  /,1X,15(5H*****)/)
   40 WRITE(IOUT, 10200)
10200 FORMAT(' ----- DADMKS: LACK OF INITIALISATION')
      STOP
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DADMMU(MODE,ISGN,HHV,PNU,PWB,Q1,Q2,PHX)
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      REAL*4         PHX(4)
      REAL  HHV(4),HV(4),PNU(4),PWB(4),Q1(4),Q2(4)
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4),PDUM5(4)
      REAL*4 RRR(3)
      REAL*8 SWT, SSWT
      DATA PI /3.141592653589793238462643/
      DATA IWARM /0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        NEVRAW=0
        NEVACC=0
        NEVOVR=0
        SWT=0
        SSWT=0
        WTMAX=1E-20
        DO 10 I=1,500
          CALL DPHSMU(WT,HV,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          IF(WT.GT.WTMAX/1.2) WTMAX=WT*1.2
   10   CONTINUE
CC      CALL HBOOK1(802,'WEIGHT DISTRIBUTION  DADMMU    $',100,0,2)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
   20   CONTINUE
        IF(IWARM.EQ.0) GOTO 40
        NEVRAW=NEVRAW+1
        CALL DPHSMU(WT,HV,PNU,PWB,Q1,Q2,PHX)
CC      CALL HFILL(802,WT/WTMAX)
        SWT=SWT+WT
        SSWT=SSWT+WT**2
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX) NEVOVR=NEVOVR+1
        IF(RN*WTMAX.GT.WT) GOTO 20
C ROTATIONS TO BASIC TAU REST FRAME
        COSTHE=-1.+2.*RRR(2)
        THET=ACOS(COSTHE)
        PHI =2*PI*RRR(3)
        CALL ROTOR2(THET,PNU,PNU)
        CALL ROTOR3( PHI,PNU,PNU)
        CALL ROTOR2(THET,PWB,PWB)
        CALL ROTOR3( PHI,PWB,PWB)
        CALL ROTOR2(THET,Q1,Q1)
        CALL ROTOR3( PHI,Q1,Q1)
        CALL ROTOR2(THET,Q2,Q2)
        CALL ROTOR3( PHI,Q2,Q2)
        CALL ROTOR2(THET,HV,HV)
        CALL ROTOR3( PHI,HV,HV)
        CALL ROTOR2(THET,PHX,PHX)
        CALL ROTOR3( PHI,PHX,PHX)
        DO 30 ,I=1,3
   30   HHV(I)=-ISGN*HV(I)
        NEVACC=NEVACC+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVRAW.EQ.0) RETURN
        PARGAM=SWT/FLOAT(NEVRAW+1)
        ERROR=0
        IF(NEVRAW.NE.0) ERROR=SQRT(SSWT/SWT**2-1./FLOAT(NEVRAW))
        RAT=PARGAM/GAMEL
        WRITE(IOUT, 10000) NEVRAW,NEVACC,NEVOVR,PARGAM,RAT,ERROR
CC      CALL HPRINT(802)
        GAMPMC(2)=RAT
        GAMPER(2)=ERROR
CAM     NEVDEC(2)=NEVACC
      ENDIF
C     =====
      RETURN
10000 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'******** DADMMU FINAL REPORT  ******** ',9X,1H*
     + /,' *',I20  ,5X,'NEVRAW = NO. OF MU  DECAYS TOTAL       ',9X,1H*
     + /,' *',I20  ,5X,'NEVACC = NO. OF MU   DECS. ACCEPTED    ',9X,1H*
     + /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
     + /,' *',E20.5,5X,'PARTIAL WTDTH (MU  DECAY) IN GEV UNITS ',9X,1H*
     + /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
     + /,' *',F20.9,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
     + /,' *',25X,     'COMPLETE QED CORRECTIONS INCLUDED      ',9X,1H*
     + /,' *',25X,     'BUT ONLY V-A CUPLINGS                  ',9X,1H*
     +  /,1X,15(5H*****)/)
   40 WRITE(IOUT, 10100)
10100 FORMAT(' ----- DADMMU: LACK OF INITIALISATION')
      STOP
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DADMPI(MODE,ISGN,HV,PPI,PNU)
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      REAL  PPI(4),PNU(4),HV(4)
      DATA PI /3.141592653589793238462643/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        NEVTOT=0
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
        NEVTOT=NEVTOT+1
        EPI= (AMTAU**2+AMPI**2-AMNUTA**2)/(2*AMTAU)
        ENU= (AMTAU**2-AMPI**2+AMNUTA**2)/(2*AMTAU)
        XPI= SQRT(EPI**2-AMPI**2)
C PI MOMENTUM
        CALL SPHERA(XPI,PPI)
        PPI(4)=EPI
C TAU-NEUTRINO MOMENTUM
        DO 10 I=1,3
   10   PNU(I)=-PPI(I)
        PNU(4)=ENU
        PXQ=AMTAU*EPI
        PXN=AMTAU*ENU
        QXN=PPI(4)*PNU(4)-PPI(1)*PNU(1)-PPI(2)*PNU(2)-PPI(3)*PNU(3)
        BRAK=(GV**2+GA**2)*(2*PXQ*QXN-AMPI**2*PXN)
     &      +(GV**2-GA**2)*AMTAU*AMNUTA*AMPI**2
        DO 20 I=1,3
   20   HV(I)=-ISGN*2*GA*GV*AMTAU*(2*PPI(I)*QXN-PNU(I)*AMPI**2)/BRAK
        HV(4)=1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVTOT.EQ.0) RETURN
        FPI=0.1284
C        GAMM=(GFERMI*FPI)**2/(16.*PI)*AMTAU**3*
C     *       (BRAK/AMTAU**4)**2
CZW 7.02.93 HERE WAS AN ERROR AFFECTING NON STANDARD MODEL
C       CONFIGURATIONS ONLY
        GAMM=(GFERMI*FPI)**2/(16.*PI)*AMTAU**3*
     $       (BRAK/AMTAU**4)*
     $       SQRT((AMTAU**2-AMPI**2-AMNUTA**2)**2
     $            -4*AMPI**2*AMNUTA**2           )/AMTAU**2
        ERROR=0
        RAT=GAMM/GAMEL
        WRITE(IOUT, 10000) NEVTOT,GAMM,RAT,ERROR
        GAMPMC(3)=RAT
        GAMPER(3)=ERROR
CAM     NEVDEC(3)=NEVTOT
      ENDIF
C     =====
      RETURN
10000 FORMAT(///1X,15(5H*****)
     $ /,' *',     25X,'******** DADMPI FINAL REPORT  ******** ',9X,1H*
     $ /,' *',I20  ,5X,'NEVTOT = NO. OF PI  DECAYS TOTAL       ',9X,1H*
     $ /,' *',E20.5,5X,'PARTIAL WTDTH ( PI DECAY) IN GEV UNITS ',9X,1H*
     $ /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
     $ /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH (STAT.)',9X,1H*
     $  /,1X,15(5H*****)/)
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DADMRO(MODE,ISGN,HHV,PNU,PRO,PIC,PIZ)
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      REAL  HHV(4)
      REAL  HV(4),PRO(4),PNU(4),PIC(4),PIZ(4)
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4)
      REAL*4 RRR(3)
      REAL*8 SWT, SSWT
      DATA PI /3.141592653589793238462643/
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        NEVRAW=0
        NEVACC=0
        NEVOVR=0
        SWT=0
        SSWT=0
        WTMAX=1E-20
        DO 10 I=1,500
          CALL DPHSRO(WT,HV,PDUM1,PDUM2,PDUM3,PDUM4)
          IF(WT.GT.WTMAX/1.2) WTMAX=WT*1.2
   10   CONTINUE
CC      CALL HBOOK1(801,'WEIGHT DISTRIBUTION  DADMRO    $',100,0,2)
CC      PRINT 7003,WTMAX
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
   20   CONTINUE
        IF(IWARM.EQ.0) GOTO 40
        CALL DPHSRO(WT,HV,PNU,PRO,PIC,PIZ)
CC      CALL HFILL(801,WT/WTMAX)
        NEVRAW=NEVRAW+1
        SWT=SWT+WT
        SSWT=SSWT+WT**2
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX) NEVOVR=NEVOVR+1
        IF(RN*WTMAX.GT.WT) GOTO 20
C ROTATIONS TO BASIC TAU REST FRAME
        COSTHE=-1.+2.*RRR(2)
        THET=ACOS(COSTHE)
        PHI =2*PI*RRR(3)
        CALL ROTOR2(THET,PNU,PNU)
        CALL ROTOR3( PHI,PNU,PNU)
        CALL ROTOR2(THET,PRO,PRO)
        CALL ROTOR3( PHI,PRO,PRO)
        CALL ROTOR2(THET,PIC,PIC)
        CALL ROTOR3( PHI,PIC,PIC)
        CALL ROTOR2(THET,PIZ,PIZ)
        CALL ROTOR3( PHI,PIZ,PIZ)
        CALL ROTOR2(THET,HV,HV)
        CALL ROTOR3( PHI,HV,HV)
        DO 30 I=1,3
   30   HHV(I)=-ISGN*HV(I)
        NEVACC=NEVACC+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVRAW.EQ.0) RETURN
        PARGAM=SWT/FLOAT(NEVRAW+1)
        ERROR=0
        IF(NEVRAW.NE.0) ERROR=SQRT(SSWT/SWT**2-1./FLOAT(NEVRAW))
        RAT=PARGAM/GAMEL
        WRITE(IOUT, 10100) NEVRAW,NEVACC,NEVOVR,PARGAM,RAT,ERROR
CC      CALL HPRINT(801)
        GAMPMC(4)=RAT
        GAMPER(4)=ERROR
CAM     NEVDEC(4)=NEVACC
      ENDIF
C     =====
      RETURN
10000 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'******** DADMRO INITIALISATION ********',9X,1H*
     + /,' *',E20.5,5X,'WTMAX  = MAXIMUM WEIGHT                ',9X,1H*
     +  /,1X,15(5H*****)/)
10100 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'******** DADMRO FINAL REPORT  ******** ',9X,1H*
     + /,' *',I20  ,5X,'NEVRAW = NO. OF RHO DECAYS TOTAL       ',9X,1H*
     + /,' *',I20  ,5X,'NEVACC = NO. OF RHO  DECS. ACCEPTED    ',9X,1H*
     + /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
     + /,' *',E20.5,5X,'PARTIAL WTDTH (RHO DECAY) IN GEV UNITS ',9X,1H*
     + /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
     + /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
     +  /,1X,15(5H*****)/)
   40 WRITE(IOUT, 10200)
10200 FORMAT(' ----- DADMRO: LACK OF INITIALISATION')
      STOP
      END
*CMZ :  1.01/50 22/05/96  18.06.08  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DADNEW(MODE,ISGN,HV,PNU,PWB,PNPI,JNPI)
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / MDECOMP /IDFFIN(9,NMODE),MULPIK(NMODE)
     +                ,NAMES
      CHARACTER NAMES(NMODE)*31

      REAL*4 PNU(4),PWB(4),PNPI(4,9),HV(4),HHV(4)
      REAL*4 PDUM1(4),PDUM2(4),PDUMI(4,9)
      REAL*4 RRR(3)
      REAL*4 WTMAX(NMODE)
      REAL*8              SWT(NMODE),SSWT(NMODE)
      DIMENSION NEVRAW(NMODE),NEVOVR(NMODE),NEVACC(NMODE)
C
      DATA PI /3.141592653589793238462643/
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
C -- AT THE MOMENT ONLY TWO DECAY MODES OF MULTIPIONS HAVE M. ELEM
        NMOD=NMODE
        IWARM=1
C       PRINT 7003
        DO 10 JNPI=1,NMOD
          NEVRAW(JNPI)=0
          NEVACC(JNPI)=0
          NEVOVR(JNPI)=0
          SWT(JNPI)=0
          SSWT(JNPI)=0
          WTMAX(JNPI)=-1.
          DO I=1,500
            IF (JNPI.LE.0) THEN
              GOTO 60
            ELSEIF(JNPI.LE.NM4) THEN
              CALL DPH4PI(WT,HV,PDUM1,PDUM2,PDUMI,JNPI)
            ELSEIF(JNPI.LE.NM4+NM5) THEN
              CALL DPH5PI(WT,HV,PDUM1,PDUM2,PDUMI,JNPI)
            ELSEIF(JNPI.LE.NM4+NM5+NM6) THEN
              CALL DPHNPI(WT,HV,PDUM1,PDUM2,PDUMI,JNPI)
            ELSEIF(JNPI.LE.NM4+NM5+NM6+NM3) THEN
              INUM=JNPI-NM4-NM5-NM6
              CALL DPHSPK(WT,HV,PDUM1,PDUM2,PDUMI,INUM)
            ELSEIF(JNPI.LE.NM4+NM5+NM6+NM3+NM2) THEN
              INUM=JNPI-NM4-NM5-NM6-NM3
              CALL DPHSRK(WT,HV,PDUM1,PDUM2,PDUMI,INUM)
            ELSE
              GOTO 60
            ENDIF
            IF(WT.GT.WTMAX(JNPI)/1.2) WTMAX(JNPI)=WT*1.2
          ENDDO
C       CALL HBOOK1(801,'WEIGHT DISTRIBUTION  DADNPI    $',100,0.,2.,.0)
C       PRINT 7004,WTMAX(JNPI)
   10   CONTINUE
        WRITE(IOUT,10200)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
        IF(IWARM.EQ.0) GOTO 50
C
   20   CONTINUE
        IF (JNPI.LE.0) THEN
          GOTO 60
        ELSEIF(JNPI.LE.NM4) THEN
          CALL DPH4PI(WT,HHV,PNU,PWB,PNPI,JNPI)
        ELSEIF(JNPI.LE.NM4+NM5) THEN
          CALL DPH5PI(WT,HHV,PNU,PWB,PNPI,JNPI)
        ELSEIF(JNPI.LE.NM4+NM5+NM6) THEN
          CALL DPHNPI(WT,HHV,PNU,PWB,PNPI,JNPI)
        ELSEIF(JNPI.LE.NM4+NM5+NM6+NM3) THEN
          INUM=JNPI-NM4-NM5-NM6
          CALL DPHSPK(WT,HHV,PNU,PWB,PNPI,INUM)
        ELSEIF(JNPI.LE.NM4+NM5+NM6+NM3+NM2) THEN
          INUM=JNPI-NM4-NM5-NM6-NM3
          CALL DPHSRK(WT,HHV,PNU,PWB,PNPI,INUM)
        ELSE
          GOTO 60
        ENDIF
        DO I=1,4
          HV(I)=-ISGN*HHV(I)
        ENDDO
C       CALL HFILL(801,WT/WTMAX(JNPI))
        NEVRAW(JNPI)=NEVRAW(JNPI)+1
        SWT(JNPI)=SWT(JNPI)+WT
        SSWT(JNPI)=SSWT(JNPI)+WT**2
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX(JNPI)) NEVOVR(JNPI)=NEVOVR(JNPI)+1
        IF(RN*WTMAX(JNPI).GT.WT) GOTO 20
C ROTATIONS TO BASIC TAU REST FRAME
        COSTHE=-1.+2.*RRR(2)
        THET=ACOS(COSTHE)
        PHI =2*PI*RRR(3)
        CALL ROTOR2(THET,PNU,PNU)
        CALL ROTOR3( PHI,PNU,PNU)
        CALL ROTOR2(THET,PWB,PWB)
        CALL ROTOR3( PHI,PWB,PWB)
        CALL ROTOR2(THET,HV,HV)
        CALL ROTOR3( PHI,HV,HV)
        ND=MULPIK(JNPI)
        DO 30  I=1,ND
          CALL ROTOR2(THET,PNPI(1,I),PNPI(1,I))
          CALL ROTOR3( PHI,PNPI(1,I),PNPI(1,I))
   30   CONTINUE
        NEVACC(JNPI)=NEVACC(JNPI)+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        DO 40  JNPI=1,NMOD
          IF(NEVRAW(JNPI).EQ.0) GOTO 40
          PARGAM=SWT(JNPI)/FLOAT(NEVRAW(JNPI)+1)
          ERROR=0
          IF(NEVRAW(JNPI).NE.0)
     +    ERROR=SQRT(SSWT(JNPI)/SWT(JNPI)**2-1./FLOAT(NEVRAW(JNPI)))
          RAT=PARGAM/GAMEL
          WRITE(IOUT, 10300) NAMES(JNPI), NEVRAW(JNPI),NEVACC(JNPI),
     +    NEVOVR(JNPI),PARGAM,RAT,ERROR
CC        CALL HPRINT(801)
          GAMPMC(8+JNPI-1)=RAT
          GAMPER(8+JNPI-1)=ERROR
CAM       NEVDEC(8+JNPI-1)=NEVACC(JNPI)
   40   CONTINUE
      ENDIF
C     =====
      RETURN
10000 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'******** DADNEW INITIALISATION ********',9X,1H*
     + )
10100 FORMAT(' *',E20.5,5X,'WTMAX  = MAXIMUM WEIGHT  ',9X,1H*/)
10200 FORMAT(
     +  /,1X,15(5H*****)/)
10300 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'******** DADNEW FINAL REPORT  ******** ',9X,1H*
     + /,' *',     25X,'CHANNEL:',A31                           ,9X,1H*
     + /,' *',I20  ,5X,'NEVRAW = NO. OF DECAYS TOTAL           ',9X,1H*
     + /,' *',I20  ,5X,'NEVACC = NO. OF DECAYS ACCEPTED        ',9X,1H*
     + /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
     + /,' *',E20.5,5X,'PARTIAL WTDTH IN GEV UNITS             ',9X,1H*
     + /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
     + /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
     +  /,1X,15(5H*****)/)
   50 WRITE(IOUT, 10400)
10400 FORMAT(' ----- DADNEW: LACK OF INITIALISATION')
      STOP
   60 WRITE(IOUT, 10500) JNPI,MODE
10500 FORMAT(' ----- DADNEW: WRONG JNPI',2I5)
      STOP
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DAM4PI(MNUM,PT,PN,PIM1,PIM2,PIM3,PIM4,AMPLIT,HV)
C ----------------------------------------------------------------------
* CALCULATES DIFFERENTIAL CROSS SECTION AND POLARIMETER VECTOR
* FOR TAU DECAY INTO 4 PI MODES
* ALL SPIN EFFECTS IN THE FULL DECAY CHAIN ARE TAKEN INTO ACCOUNT.
* CALCULATIONS DONE IN TAU REST FRAME WITH Z-AXIS ALONG NEUTRINO MOMENT
C MNUM DECAY MODE IDENTIFIER.
C
C     CALLED BY : DPHSAA
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL  HV(4),PT(4),PN(4),PIM1(4),PIM2(4),PIM3(4),PIM4(4)
      REAL  PIVEC(4),PIAKS(4),HVM(4)
      COMPLEX HADCUR(4),FORM1,FORM2,FORM3,FORM4,FORM5
      EXTERNAL FORM1,FORM2,FORM3,FORM4,FORM5
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
C
      CALL CURR(MNUM,PIM1,PIM2,PIM3,PIM4,HADCUR)
C
* CALCULATE PI-VECTORS: VECTOR AND AXIAL
      CALL CLVEC(HADCUR,PN,PIVEC)
      CALL CLAXI(HADCUR,PN,PIAKS)
      CALL CLNUT(HADCUR,BRAKM,HVM)
* SPIN INDEPENDENT PART OF DECAY DIFF-CROSS-SECT. IN TAU REST  FRAME
      BRAK= (GV**2+GA**2)*PT(4)*PIVEC(4) +2.*GV*GA*PT(4)*PIAKS(4)
     +     +2.*(GV**2-GA**2)*AMNUTA*AMTAU*BRAKM
      AMPLIT=(CCABIB*GFERMI)**2*BRAK/2.
C POLARIMETER VECTOR IN TAU REST FRAME
      DO 10 I=1,3
        HV(I)=-(AMTAU*((GV**2+GA**2)*PIAKS(I)+2.*GV*GA*PIVEC(I)))
     +  +(GV**2-GA**2)*AMNUTA*AMTAU*HVM(I)
C HV IS DEFINED FOR TAU-    WITH GAMMA=B+HV*POL
        IF (BRAK.NE.0.0) HV(I)=-HV(I)/BRAK
   10 CONTINUE
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DAMPAA(PT,PN,PIM1,PIM2,PIPL,AMPLIT,HV)
C ----------------------------------------------------------------------
* CALCULATES DIFFERENTIAL CROSS SECTION AND POLARIMETER VECTOR
* FOR TAU DECAY INTO A1, A1 DECAYS NEXT INTO RHO+PI AND RHO INTO PI+PI.
* ALL SPIN EFFECTS IN THE FULL DECAY CHAIN ARE TAKEN INTO ACCOUNT.
* CALCULATIONS DONE IN TAU REST FRAME WITH Z-AXIS ALONG NEUTRINO MOMENT
* THE ROUTINE IS WRITEN FOR ZERO NEUTRINO MASS.
C
C     CALLED BY : DPHSAA
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON /TESTA1/ KEYA1
      REAL  HV(4),PT(4),PN(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PAA(4),VEC1(4),VEC2(4)
      REAL  PIVEC(4),PIAKS(4),HVM(4)
      COMPLEX BWIGN,HADCUR(4),FPIK
      DATA ICONT /1/
C
* F CONSTANTS FOR A1, A1-RHO-PI, AND RHO-PI-PI
*
      DATA  FPI /93.3E-3/
* THIS INLINE FUNCT. CALCULATES THE SCALAR PART OF THE PROPAGATOR
      BWIGN(XM,AM,GAMMA)=1./CMPLX(XM**2-AM**2,GAMMA*AM)
C
* FOUR MOMENTUM OF A1
      DO 10 I=1,4
   10 PAA(I)=PIM1(I)+PIM2(I)+PIPL(I)
* MASSES OF A1, AND OF TWO PI-PAIRS WHICH MAY FORM RHO
      XMAA   =SQRT(ABS(PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2))
      XMRO1  =SQRT(ABS((PIPL(4)+PIM1(4))**2-(PIPL(1)+PIM1(1))**2
     +                -(PIPL(2)+PIM1(2))**2-(PIPL(3)+PIM1(3))**2))
      XMRO2  =SQRT(ABS((PIPL(4)+PIM2(4))**2-(PIPL(1)+PIM2(1))**2
     +                -(PIPL(2)+PIM2(2))**2-(PIPL(3)+PIM2(3))**2))
* ELEMENTS OF HADRON CURRENT
      PROD1  =PAA(4)*(PIM1(4)-PIPL(4))-PAA(1)*(PIM1(1)-PIPL(1))
     +       -PAA(2)*(PIM1(2)-PIPL(2))-PAA(3)*(PIM1(3)-PIPL(3))
      PROD2  =PAA(4)*(PIM2(4)-PIPL(4))-PAA(1)*(PIM2(1)-PIPL(1))
     +       -PAA(2)*(PIM2(2)-PIPL(2))-PAA(3)*(PIM2(3)-PIPL(3))
      DO 20 I=1,4
        VEC1(I)= PIM1(I)-PIPL(I) -PAA(I)*PROD1/XMAA**2
   20 VEC2(I)= PIM2(I)-PIPL(I) -PAA(I)*PROD2/XMAA**2
* HADRON CURRENT SATURATED WITH A1 AND RHO RESONANCES
      IF (KEYA1.EQ.1) THEN
        FA1=9.87
        FAROPI=1.0
        FRO2PI=1.0
        FNORM=FA1/SQRT(2.)*FAROPI*FRO2PI
        DO 30 I=1,4
          HADCUR(I)= CMPLX(FNORM) *AMA1**2*BWIGN(XMAA,AMA1,GAMA1)
     +    *(CMPLX(VEC1(I))*AMRO**2*BWIGN(XMRO1,AMRO,GAMRO) +
     +    CMPLX(VEC2(I))*AMRO**2*BWIGN(XMRO2,AMRO,GAMRO))
   30   CONTINUE
      ELSE
        FNORM=2.0*SQRT(2.)/3.0/FPI
        GAMAX=GAMA1*GFUN(XMAA**2)/GFUN(AMA1**2)
        DO 40 I=1,4
          HADCUR(I)= CMPLX(FNORM) *AMA1**2*BWIGN(XMAA,AMA1,GAMAX)
     +    *(CMPLX(VEC1(I))*FPIK(XMRO1) +CMPLX(VEC2(I))*FPIK(XMRO2))
   40   CONTINUE
      ENDIF
C
* CALCULATE PI-VECTORS: VECTOR AND AXIAL
      CALL CLVEC(HADCUR,PN,PIVEC)
      CALL CLAXI(HADCUR,PN,PIAKS)
      CALL CLNUT(HADCUR,BRAKM,HVM)
* SPIN INDEPENDENT PART OF DECAY DIFF-CROSS-SECT. IN TAU REST  FRAME
      BRAK= (GV**2+GA**2)*PT(4)*PIVEC(4) +2.*GV*GA*PT(4)*PIAKS(4)
     +     +2.*(GV**2-GA**2)*AMNUTA*AMTAU*BRAKM
      AMPLIT=(GFERMI*CCABIB)**2*BRAK/2.
C THE STATISTICAL FACTOR FOR IDENTICAL PI'S WAS CANCELLED WITH
C TWO, FOR TWO MODES OF A1 DECAY NAMELLY PI+PI-PI- AND PI-PI0PI0
C POLARIMETER VECTOR IN TAU REST FRAME
      DO 50 I=1,3
        HV(I)=-(AMTAU*((GV**2+GA**2)*PIAKS(I)+2.*GV*GA*PIVEC(I)))
     +  +(GV**2-GA**2)*AMNUTA*AMTAU*HVM(I)
C HV IS DEFINED FOR TAU-    WITH GAMMA=B+HV*POL
        HV(I)=-HV(I)/BRAK
   50 CONTINUE
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DAMPOG(PT,PN,PIM1,PIM2,PIPL,AMPLIT,HV)
C ----------------------------------------------------------------------
* CALCULATES DIFFERENTIAL CROSS SECTION AND POLARIMETER VECTOR
* FOR TAU DECAY INTO A1, A1 DECAYS NEXT INTO RHO+PI AND RHO INTO PI+PI.
* ALL SPIN EFFECTS IN THE FULL DECAY CHAIN ARE TAKEN INTO ACCOUNT.
* CALCULATIONS DONE IN TAU REST FRAME WITH Z-AXIS ALONG NEUTRINO MOMENT
* THE ROUTINE IS WRITEN FOR ZERO NEUTRINO MASS.
C
C     CALLED BY : DPHSAA
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON /TESTA1/ KEYA1
      REAL  HV(4),PT(4),PN(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PAA(4),VEC1(4),VEC2(4)
      REAL  PIVEC(4),PIAKS(4),HVM(4)
      COMPLEX BWIGN,HADCUR(4),FNORM,FORMOM
      DATA ICONT /1/
* THIS INLINE FUNCT. CALCULATES THE SCALAR PART OF THE PROPAGATOR
      BWIGN(XM,AM,GAMMA)=1./CMPLX(XM**2-AM**2,GAMMA*AM)
C
* FOUR MOMENTUM OF A1
      DO 10 I=1,4
        VEC1(I)=0.0
        VEC2(I)=0.0
        HV(I) =0.0
   10 PAA(I)=PIM1(I)+PIM2(I)+PIPL(I)
      VEC1(1)=1.0
* MASSES OF A1, AND OF TWO PI-PAIRS WHICH MAY FORM RHO
      XMAA   =SQRT(ABS(PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2))
      XMOM   =SQRT(ABS( (PIM2(4)+PIPL(4))**2-(PIM2(3)+PIPL(3))**2
     +                 -(PIM2(2)+PIPL(2))**2-(PIM2(1)+PIPL(1))**2   ))
      XMRO2  =(PIPL(1))**2 +(PIPL(2))**2 +(PIPL(3))**2
* ELEMENTS OF HADRON CURRENT
      PROD1  =VEC1(1)*PIPL(1)
      PROD2  =VEC2(2)*PIPL(2)
      P12    =PIM1(4)*PIM2(4)-PIM1(1)*PIM2(1)
     +       -PIM1(2)*PIM2(2)-PIM1(3)*PIM2(3)
      P1PL   =PIM1(4)*PIPL(4)-PIM1(1)*PIPL(1)
     +       -PIM1(2)*PIPL(2)-PIM1(3)*PIPL(3)
      P2PL   =PIPL(4)*PIM2(4)-PIPL(1)*PIM2(1)
     +       -PIPL(2)*PIM2(2)-PIPL(3)*PIM2(3)
      DO 20 I=1,3
        VEC1(I)= (VEC1(I)-PROD1/XMRO2*PIPL(I))
   20 CONTINUE
      GNORM=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
      DO 30 I=1,3
        VEC1(I)= VEC1(I)/GNORM
   30 CONTINUE
      VEC2(1)=(VEC1(2)*PIPL(3)-VEC1(3)*PIPL(2))/SQRT(XMRO2)
      VEC2(2)=(VEC1(3)*PIPL(1)-VEC1(1)*PIPL(3))/SQRT(XMRO2)
      VEC2(3)=(VEC1(1)*PIPL(2)-VEC1(2)*PIPL(1))/SQRT(XMRO2)
      P1VEC1   =PIM1(4)*VEC1(4)-PIM1(1)*VEC1(1)
     +         -PIM1(2)*VEC1(2)-PIM1(3)*VEC1(3)
      P2VEC1   =VEC1(4)*PIM2(4)-VEC1(1)*PIM2(1)
     +         -VEC1(2)*PIM2(2)-VEC1(3)*PIM2(3)
      P1VEC2   =PIM1(4)*VEC2(4)-PIM1(1)*VEC2(1)
     +         -PIM1(2)*VEC2(2)-PIM1(3)*VEC2(3)
      P2VEC2   =VEC2(4)*PIM2(4)-VEC2(1)*PIM2(1)
     +         -VEC2(2)*PIM2(2)-VEC2(3)*PIM2(3)
* HADRON CURRENT
      FNORM=FORMOM(XMAA,XMOM)
      BRAK=0.0
      DO 60  JJ=1,2
        DO 40 I=1,4
          IF (JJ.EQ.1) THEN
            HADCUR(I) = FNORM *( VEC1(I)*(AMPI**2*P1PL-P2PL*(P12-P1PL))
     +      -PIM2(I)*(P2VEC1*P1PL-P1VEC1*P2PL) +PIPL(I)*(P2VEC1*P12 -
     +      P1VEC1*(AMPI**2+P2PL)) )
          ELSE
            HADCUR(I) = FNORM *( VEC2(I)*(AMPI**2*P1PL-P2PL*(P12-P1PL))
     +      -PIM2(I)*(P2VEC2*P1PL-P1VEC2*P2PL) +PIPL(I)*(P2VEC2*P12 -
     +      P1VEC2*(AMPI**2+P2PL)) )
          ENDIF
   40   CONTINUE
C
* CALCULATE PI-VECTORS: VECTOR AND AXIAL
        CALL CLVEC(HADCUR,PN,PIVEC)
        CALL CLAXI(HADCUR,PN,PIAKS)
        CALL CLNUT(HADCUR,BRAKM,HVM)
* SPIN INDEPENDENT PART OF DECAY DIFF-CROSS-SECT. IN TAU REST  FRAME
        BRAK=BRAK+(GV**2+GA**2)*PT(4)*PIVEC(4) +2.*GV*GA*PT(4)*PIAKS(4)
     +  +2.*(GV**2-GA**2)*AMNUTA*AMTAU*BRAKM
        DO 50 I=1,3
          HV(I)=HV(I)-(AMTAU*((GV**2+GA**2)*PIAKS(I)+2.*GV*GA*PIVEC(I))
     +    ) +(GV**2-GA**2)*AMNUTA*AMTAU*HVM(I)
   50   CONTINUE
C HV IS DEFINED FOR TAU-    WITH GAMMA=B+HV*POL
   60 CONTINUE
      AMPLIT=(GFERMI*CCABIB)**2*BRAK/2.
C THE STATISTICAL FACTOR FOR IDENTICAL PI'S WAS CANCELLED WITH
C TWO, FOR TWO MODES OF A1 DECAY NAMELLY PI+PI-PI- AND PI-PI0PI0
C POLARIMETER VECTOR IN TAU REST FRAME
      DO 70 I=1,3
        HV(I)=-HV(I)/BRAK
   70 CONTINUE

      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DAMPPK(MNUM,PT,PN,PIM1,PIM2,PIM3,AMPLIT,HV)
C ----------------------------------------------------------------------
* CALCULATES DIFFERENTIAL CROSS SECTION AND POLARIMETER VECTOR
* FOR TAU DECAY INTO K K PI, K PI PI.
* ALL SPIN EFFECTS IN THE FULL DECAY CHAIN ARE TAKEN INTO ACCOUNT.
* CALCULATIONS DONE IN TAU REST FRAME WITH Z-AXIS ALONG NEUTRINO MOMENT
C MNUM DECAY MODE IDENTIFIER.
C
C     CALLED BY : DPHSAA
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL  HV(4),PT(4),PN(4),PIM1(4),PIM2(4),PIM3(4)
      REAL  PAA(4),VEC1(4),VEC2(4),VEC3(4),VEC4(4),VEC5(4)
      REAL  PIVEC(4),PIAKS(4),HVM(4)
      REAL FNORM(0:7),COEF(1:5,0:7)
      COMPLEX HADCUR(4),FORM1,FORM2,FORM3,FORM4,FORM5,UROJ
      EXTERNAL FORM1,FORM2,FORM3,FORM4,FORM5
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
C
      DATA  FPI /93.3E-3/
      IF (ICONT.EQ.0) THEN
        ICONT=1
        UROJ=CMPLX(0.0,1.0)
        DWAPI0=SQRT(2.0)
        FNORM(0)=CCABIB/FPI
        FNORM(1)=CCABIB/FPI
        FNORM(2)=CCABIB/FPI
        FNORM(3)=CCABIB/FPI
        FNORM(4)=SCABIB/FPI/DWAPI0
        FNORM(5)=SCABIB/FPI
        FNORM(6)=SCABIB/FPI
        FNORM(7)=CCABIB/FPI
C
        COEF(1,0)= 2.0*SQRT(2.)/3.0
        COEF(2,0)=-2.0*SQRT(2.)/3.0
        COEF(3,0)= 0.0
        COEF(4,0)= FPI
        COEF(5,0)= 0.0
C
        COEF(1,1)=-SQRT(2.)/3.0
        COEF(2,1)= SQRT(2.)/3.0
        COEF(3,1)= 0.0
        COEF(4,1)= FPI
        COEF(5,1)= SQRT(2.)
C
        COEF(1,2)=-SQRT(2.)/3.0
        COEF(2,2)= SQRT(2.)/3.0
        COEF(3,2)= 0.0
        COEF(4,2)= 0.0
        COEF(5,2)=-SQRT(2.)
C
        COEF(1,3)= 0.0
        COEF(2,3)=-1.0
        COEF(3,3)= 0.0
        COEF(4,3)= 0.0
        COEF(5,3)= 0.0
C
        COEF(1,4)= 1.0/SQRT(2.)/3.0
        COEF(2,4)=-1.0/SQRT(2.)/3.0
        COEF(3,4)= 0.0
        COEF(4,4)= 0.0
        COEF(5,4)= 0.0
C
        COEF(1,5)=-SQRT(2.)/3.0
        COEF(2,5)= SQRT(2.)/3.0
        COEF(3,5)= 0.0
        COEF(4,5)= 0.0
        COEF(5,5)=-SQRT(2.)
C
        COEF(1,6)= 0.0
        COEF(2,6)=-1.0
        COEF(3,6)= 0.0
        COEF(4,6)= 0.0
        COEF(5,6)=-2.0
C
        COEF(1,7)= 0.0
        COEF(2,7)= 0.0
        COEF(3,7)= 0.0
        COEF(4,7)= 0.0
        COEF(5,7)=-SQRT(2.0/3.0)
C
      ENDIF
C
      DO 10 I=1,4
   10 PAA(I)=PIM1(I)+PIM2(I)+PIM3(I)
      XMAA   =SQRT(ABS(PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2))
      XMRO1  =SQRT(ABS((PIM3(4)+PIM2(4))**2-(PIM3(1)+PIM2(1))**2
     +                -(PIM3(2)+PIM2(2))**2-(PIM3(3)+PIM2(3))**2))
      XMRO2  =SQRT(ABS((PIM3(4)+PIM1(4))**2-(PIM3(1)+PIM1(1))**2
     +                -(PIM3(2)+PIM1(2))**2-(PIM3(3)+PIM1(3))**2))
      XMRO3  =SQRT(ABS((PIM1(4)+PIM2(4))**2-(PIM1(1)+PIM2(1))**2
     +                -(PIM1(2)+PIM2(2))**2-(PIM1(3)+PIM2(3))**2))
* ELEMENTS OF HADRON CURRENT
      PROD1  =PAA(4)*(PIM2(4)-PIM3(4))-PAA(1)*(PIM2(1)-PIM3(1))
     +       -PAA(2)*(PIM2(2)-PIM3(2))-PAA(3)*(PIM2(3)-PIM3(3))
      PROD2  =PAA(4)*(PIM3(4)-PIM1(4))-PAA(1)*(PIM3(1)-PIM1(1))
     +       -PAA(2)*(PIM3(2)-PIM1(2))-PAA(3)*(PIM3(3)-PIM1(3))
      PROD3  =PAA(4)*(PIM1(4)-PIM2(4))-PAA(1)*(PIM1(1)-PIM2(1))
     +       -PAA(2)*(PIM1(2)-PIM2(2))-PAA(3)*(PIM1(3)-PIM2(3))
      DO 20 I=1,4
        VEC1(I)= PIM2(I)-PIM3(I) -PAA(I)*PROD1/XMAA**2
        VEC2(I)= PIM3(I)-PIM1(I) -PAA(I)*PROD2/XMAA**2
        VEC3(I)= PIM1(I)-PIM2(I) -PAA(I)*PROD3/XMAA**2
   20 VEC4(I)= PIM1(I)+PIM2(I)+PIM3(I)
      CALL PROD5(PIM1,PIM2,PIM3,VEC5)
* HADRON CURRENT
C BE AWARE THAT SIGN OF VEC2 IS OPPOSITE TO SIGN OF VEC1 IN A1 CASE
      DO 30 I=1,4
        HADCUR(I)= CMPLX(FNORM(MNUM)) * ( CMPLX(VEC1(I)*COEF(1,MNUM))*
     +  FORM1(MNUM,XMAA**2,XMRO1**2,XMRO2**2)+CMPLX(VEC2(I)*COEF(2,
     +  MNUM))*FORM2(MNUM,XMAA**2,XMRO2**2,XMRO1**2)+CMPLX(VEC3(I)*
     +  COEF(3,MNUM))*FORM3(MNUM,XMAA**2,XMRO3**2,XMRO1**2)+(-1.0*UROJ)
     +  * CMPLX(VEC4(I)*COEF(4,MNUM))*FORM4(MNUM,XMAA**2,XMRO1**2,
     +  XMRO2**2,XMRO3**2) +(-1.0)*UROJ/4.0/PI**2/FPI**2* CMPLX(VEC5(I)
     +  *COEF(5,MNUM))*FORM5(MNUM,XMAA**2,XMRO1**2,XMRO2**2))
   30 CONTINUE
C
* CALCULATE PI-VECTORS: VECTOR AND AXIAL
      CALL CLVEC(HADCUR,PN,PIVEC)
      CALL CLAXI(HADCUR,PN,PIAKS)
      CALL CLNUT(HADCUR,BRAKM,HVM)
* SPIN INDEPENDENT PART OF DECAY DIFF-CROSS-SECT. IN TAU REST  FRAME
      BRAK= (GV**2+GA**2)*PT(4)*PIVEC(4) +2.*GV*GA*PT(4)*PIAKS(4)
     +     +2.*(GV**2-GA**2)*AMNUTA*AMTAU*BRAKM
      AMPLIT=(GFERMI)**2*BRAK/2.
      IF (MNUM.GE.9) THEN
        PRINT *, 'MNUM=',MNUM
        ZNAK=-1.0
        XM1=0.0
        XM2=0.0
        XM3=0.0
        DO 40 K=1,4
          IF (K.EQ.4) ZNAK=1.0
          XM1=ZNAK*PIM1(K)**2+XM1
          XM2=ZNAK*PIM2(K)**2+XM2
          XM3=ZNAK*PIM3(K)**2+XM3
   40   PRINT *, 'PIM1=',PIM1(K),'PIM2=',PIM2(K),'PIM3=',PIM3(K)
        PRINT *, 'XM1=',SQRT(XM1),'XM2=',SQRT(XM2),'XM3=',SQRT(XM3)
        PRINT *, '************************************************'
      ENDIF
C POLARIMETER VECTOR IN TAU REST FRAME
      DO 50 I=1,3
        HV(I)=-(AMTAU*((GV**2+GA**2)*PIAKS(I)+2.*GV*GA*PIVEC(I)))
     +  +(GV**2-GA**2)*AMNUTA*AMTAU*HVM(I)
C HV IS DEFINED FOR TAU-    WITH GAMMA=B+HV*POL
        HV(I)=-HV(I)/BRAK
   50 CONTINUE
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DAMPRY(ITDKRC,XK0DEC,XK,XA,QP,XN,AMPLIT,HV)
      IMPLICIT REAL*8 (A-H,O-Z)
C ----------------------------------------------------------------------
C IT CALCULATES MATRIX ELEMENT FOR THE
C TAU --> MU(E) NU NUBAR DECAY MODE
C INCLUDING COMPLETE ORDER ALPHA QED CORRECTIONS.
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      REAL*8  HV(4),QP(4),XN(4),XA(4),XK(4)
C
      HV(4)=1.D0
      AK0=XK0DEC*AMTAU
      IF(XK(4).LT.0.1D0*AK0) THEN
        AMPLIT=THB(ITDKRC,QP,XN,XA,AK0,HV)
      ELSE
        AMPLIT=SQM2(ITDKRC,QP,XN,XA,XK,AK0,HV)
      ENDIF
      RETURN
      END
*CMZ :  1.00/00 09/08/94  17.43.59  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION DCDMAS(IDENT)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      IF      (IDENT.EQ. 1) THEN
        APKMAS=AMPI
      ELSEIF  (IDENT.EQ.-1) THEN
        APKMAS=AMPI
      ELSEIF  (IDENT.EQ. 2) THEN
        APKMAS=AMPIZ
      ELSEIF  (IDENT.EQ.-2) THEN
        APKMAS=AMPIZ
      ELSEIF  (IDENT.EQ. 3) THEN
        APKMAS=AMK
      ELSEIF  (IDENT.EQ.-3) THEN
        APKMAS=AMK
      ELSEIF  (IDENT.EQ. 4) THEN
        APKMAS=AMKZ
      ELSEIF  (IDENT.EQ.-4) THEN
        APKMAS=AMKZ
      ELSEIF  (IDENT.EQ. 8) THEN
        APKMAS=0.0001
      ELSEIF  (IDENT.EQ.-8) THEN
        APKMAS=0.0001
      ELSEIF  (IDENT.EQ. 9) THEN
        APKMAS=0.5488
      ELSEIF  (IDENT.EQ.-9) THEN
        APKMAS=0.5488
      ELSE
        PRINT *, 'STOP IN APKMAS, WRONG IDENT=',IDENT
        STOP
      ENDIF
      DCDMAS=APKMAS
      END
*CMZ :          04/03/97  14.19.02  by  Unknown
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 14/07/94  16.40.21  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      FUNCTION DCROSS(V1,V2)

C...DIFFERENTIAL CROSS-SECTION DSIGMA/DV1DV2; V1=X, V2=Q2 OR Y OR W2.
C...USED FOR NUMERICAL INTEGRATION ETC.
C...NOTE, NON-ZERO RESULT ONLY FOR REGION DEFINED BY CUTS THROUGH CUT.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON /LINTEG/ NTOT,NPASS
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)


      DCROSS=0.
      NTOT=NTOT+1
C...VARIABLE V1 IS X, VARIABLE V2 IS EITHER Q**2, Y OR W**2
      X=V1
      IF(X.LT.XMIN.OR.X.GT.XMAX) RETURN
      S=PARL(21)

c      write(*,*) s
      PM2=PSAVE(3,2,5)**2
      IF(LST(31).EQ.1) THEN
        Q2=V2
        Y=Q2/(PARL(21)*X)
        W2=(1.-X)*Y*PARL(21)+PSAVE(3,2,5)**2
      ELSEIF(LST(31).EQ.2) THEN
        Y=V2
        Q2=Y*X*PARL(21)
        W2=(1.-X)*Y*PARL(21)+PSAVE(3,2,5)**2
      ELSEIF(LST(31).EQ.3) THEN
        W2=V2
        Y=(W2-PSAVE(3,2,5)**2)/((1.-X)*PARL(21))
        Q2=X*Y*PARL(21)
      ENDIF
      Q2LOW=MAX(Q2MIN,X*YMIN*S,(W2MIN-PM2)*X/(1.-X))
      Q2UPP=MIN(Q2MAX,X*YMAX*S,(W2MAX-PM2)*X/(1.-X))
      YLOW=MAX(YMIN,Q2MIN/(S*X),(W2MIN-PM2)/(S*(1.-X)))
      YUPP=MIN(YMAX,Q2MAX/(S*X),(W2MAX-PM2)/(S*(1.-X)))
      W2LOW=MAX(W2MIN,(1.-X)*YMIN*S+PM2,Q2MIN*(1.-X)/X+PM2)
      W2UPP=MIN(W2MAX,(1.-X)*YMAX*S+PM2,Q2MAX*(1.-X)/X+PM2)
      IF(Q2.LT.Q2LOW.OR.Q2.GT.Q2UPP) RETURN
      IF(Y.LT.YLOW.OR.Y.GT.YUPP) RETURN
      IF(W2.LT.W2LOW.OR.W2.GT.W2UPP) RETURN
      LST2=LST(2)
      LST(2)=-2
c      print*,' calling lepto'
      CALL LEPTO
      LST(2)=LST2
c      write(*,*) lst(21)
      IF(LST(21).NE.0) RETURN
      NPASS=NPASS+1
      DCROSS=PARI(31)*PQ(17)*COMFAC
c      write(*,*) dcross
      RETURN
      END
*CMZ :  1.01/50 20/03/96  12.40.51  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/00 01/09/94  17.30.38  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 09/08/94  18.46.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DECTES(KTORY)
C     ************************
*KEEP,KEYS.
        COMMON/CFREAD/SPACE(5000)
 	COMMON/MYKEYS/IEVAR,IF5CC,INEUT,IINTE,IFERM,IFLAT,ICOUN,
     &  REFIX,IMUDO,IKAT1,IKAT2,IKAT3,IKAT4,IKAT5,IKAT6,
     &  INEVT,IJAK1,IJAK2,IITDK,RPTAU,RXK0D,NTGR,IDIMUON,IGLU,
     &  IQDEN,LOME(2),IFILES,ICCHA,ISEED,IDSUBS,IONEM,EHAC,IPSEL,
     &  IHIST


*KEND.
      REAL POL(4)
      DOUBLE PRECISION HH(4)
C SWITCHES FOR TAUOLA;
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / IDFC  / IDFF
C I/O UNITS  NUMBERS
      COMMON / INOUT /  INUT,IOUT
C LUND TYPE IDENTIFIER FOR A1
      COMMON / IDPART / IA1
C /PTAU/ IS USED IN ROUTINE TRALO4
      COMMON /PTAU/ PTAU
      COMMON / TAURAD / XK0DEC,ITDKRC
      REAL*8            XK0DEC
      COMMON /TESTA1/ KEYA1
C SPECIAL SWITCH FOR TESTS OF DGAMMA/DQ**2 IN A1 DECAY
C KEYA1=1 CONSTANT WIDTH OF A1 AND RHO
C KEYA1=2 FREE CHOICE OF RHO PROPAGATOR (DEFINED IN FUNCTION FPIK)
C         AND FREE CHOICE OF A1 MASS AND WIDTH. FUNCTION G(Q**2)
C         (SEE FORMULA 3.48 IN COMP. PHYS. COMM. 64 (1991) 275)
C         HARD CODED BOTH IN MONTE CARLO AND IN TESTING DISTRIBUTION.
C KEYA1=3 FUNCTION G(Q**2) HARDCODED IN THE MONTE CARLO
C         (IT IS TIMY TO CALCULATE!), BUT APPROPRIATELY ADJUSTED IN
C         TESTING DISTRIBUTION.
C-----------------------------------------------------------------------
C          INITIALIZATION
C-----------------------------------------------------------------------
C======================================
      NINP=INUT
      NOUT=IOUT
10000 FORMAT(A80)
10100 FORMAT(8I2)
10200 FORMAT(I10)
10300 FORMAT(F10.0)
      IF (KTORY.EQ.1) THEN
*      READ( NINP,3000) TESTIT
*      WRITE(NOUT,3000) TESTIT
*      READ( NINP,3001) KAT1,KAT2,KAT3,KAT4,KAT5,KAT6
*      READ( NINP,3002) NEVT,JAK1,JAK2,ITDKRC
*      READ( NINP,3003) PTAU,XK0DEC

        KAT1=IKAT1
        KAT2=IKAT2
        KAT3=IKAT3
        KAT4=IKAT4
        KAT5=IKAT5
        KAT6=IKAT6
        NEVT=INEVT
        JAK1=IJAK1
        JAK2=IJAK2
        ITDKRC=IITDK
        PTAU=RPTAU
        XK0DEC=RXK0D




      ENDIF
C======================================
C CONTROL OUTPUT
      WRITE(NOUT,'(6A6/6I6)')
     + 'KAT1','KAT2','KAT3','KAT4','KAT5','KAT6',
     +  KAT1 , KAT2 , KAT3 , KAT4 , KAT5 , KAT6
      WRITE(NOUT,'(4A12/4I12)')
     +  'NEVT','JAK1','JAK2','ITDKRC',
     +   NEVT,  JAK1 , JAK2 , ITDKRC
      WRITE(NOUT,'(2A12/2F12.6)')
     + 'PTAU','XK0DEC',
     +  PTAU , XK0DEC
C======================================
      JAK=0
C      JAK1=5
C      JAK2=5
C LUND IDENTIFIER (FOR TAU+) -15
      IF (KTORY.EQ.1) THEN
        IDFF=-15
      ELSE
        IDFF= 15
      ENDIF
C KTO=1 DENOTES TAU DEFINED BY IDFF (I.E. TAU+)
C KTO=2 DENOTES THE OPPOSITE        (I.E. TAU-)
*PZ PATCH
      IDFF=-15

      KTO=2
      IF (KTO.NE.2) THEN
        PRINT *, 'FOR THE SAKE OF THESE TESTS KTO HAS TO BE 2'
        PRINT *, 'TO CHANGE TAU- TO TAU+ CHANGE IDFF FROM -15 TO 15'
        STOP
      ENDIF
C TAU POLARIZATION IN ITS RESTFRAME;
      POL(1)=0.
      POL(2)=0.
      POL(3)=.9
C TAU MOMENTUM IN GEV;
C      PTAU=CMSENE/2.D0
C NUMBER OF EVENTS TO BE GENERATED;
      NEVTES=0
C      NEVTES=NEVT
      PRINT *, 'NEVTES= ',NEVTES
      WRITE(IOUT,10800) KEYA1
C
      IF (KTORY.EQ.1) THEN
        WRITE(IOUT,10400) JAK,IDFF,POL(3),PTAU
      ELSE
        WRITE(IOUT,10700) JAK,IDFF,POL(3),PTAU
      ENDIF
C INITIALISATION OF TAU DECAY PACKAGE TAUOLA
C ******************************************
      CALL INIMAS
      CALL INITDK


      CALL INIPHY(0.1D0)
      IF (KTORY.EQ.1) THEN
        CALL DEXAY(-1,POL)
      ELSE
        CALL DEKAY(-1,HH)
      ENDIF

      RETURN
10400 FORMAT(//4(/1X,15(5H=====))
     + /,' ',     19X,'  TEST OF RAD. CORR IN ELECTRON DECAY   ',9X,1H ,
     + /,' ',     19X,'    TESTS OF TAU DECAY ROUTINES         ',9X,1H ,
     + /,' ',     19X,'    INTERFACE OF THE KORAL-Z TYPE       ',9X,1H ,
     +  2(/,1X,15(5H=====)),
     + /,5X ,'JAK   =',I7  ,'  KEY DEFINING DECAY TYPE         ',9X,1H ,
     + /,5X ,'IDFF  =',I7  ,'  LUND IDENTIFIER FOR FIRST TAU   ',9X,1H ,
     + /,5X ,'POL(3)=',F7.2,'  THIRD COMPONENT OF TAU POLARIZ. ',9X,1H ,
     + /,5X ,'PTAU  =',F7.2,'  THIRD COMPONENT OF TAU MOM. GEV ',9X,1H ,
     +  2(/,1X,15(5H=====))/)
10500 FORMAT(///1X, '===== EVENT NO.',I4,1X,5H=====)
10600 FORMAT(5X,'POLARIMETRIC VECTOR: ',
     +       7X,'HH(1)',7X,'HH(2)',7X,'HH(3)',7X,'HH(4)',
     + /,    5X,'                     ', 4(1X,F11.8)   )
10700 FORMAT(//4(/1X,15(5H=====))
     + /,'  ',     19X,'  TEST OF RAD. CORR IN ELECTRON DECAY  ',9X,1H ,
     + /,'  ',     19X,'    TESTS OF TAU DECAY ROUTINES        ',9X,1H ,
     + /,'  ',     19X,'    INTERFACE OF THE KORAL-B TYPE      ',9X,1H ,
     +  2(/,1X,15(5H=====)),
     + /,5X ,'JAK   =',I7  ,'  KEY DEFINING DECAY TYPE         ',9X,1H ,
     + /,5X ,'IDFF  =',I7  ,'  LUND IDENTIFIER FOR FIRST TAU   ',9X,1H ,
     + /,5X ,'POL(3)=',F7.2,'  THIRD COMPONENT OF TAU POLARIZ. ',9X,1H ,
     + /,5X ,'PTAU  =',F7.2,'  THIRD COMPONENT OF TAU MOM. GEV ',9X,1H ,
     +  2(/,1X,15(5H=====))/)
10800 FORMAT(///1X, '===== TYPE OF CURRENT',I4,1X,5H=====)
      END
*CMZ :  1.01/50 22/05/96  18.06.08  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.39  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEKAY(KTO,HX)
C     ***********************
C THIS DEKAY IS IN SPIRIT OF THE 'DECAY' WHICH
C WAS INCLUDED IN KORAL-B PROGRAM, COMP. PHYS. COMMUN.
C VOL. 36 (1985) 191, SEE COMMENTS  ON GENERAL PHILOSOPHY THERE.
C KTO=0 INITIALISATION (OBLIGATORY)
C KTO=1,11 DENOTES TAU+ AND KTO=2,12 TAU-
C DEKAY(1,H) AND DEKAY(2,H) IS CALLED INTERNALLY BY MC GENERATOR.
C H DENOTES THE POLARIMETRIC VECTOR, USED BY THE HOST PROGRAM FOR
C CALCULATION OF THE SPIN WEIGHT.
C USER MAY OPTIONALLY CALL DEKAY(11,H) DEKAY(12,H) IN ORDER
C TO TRANSFORM DECAY PRODUCTS TO CMS AND WRITE LUND RECORD IN /LUJETS/.
C KTO=100, PRINT FINAL REPORT  (OPTIONAL).
C DECAY MODES:
C JAK=1 ELECTRON DECAY
C JAK=2 MU  DECAY
C JAK=3 PI  DECAY
C JAK=4 RHO DECAY
C JAK=5 A1  DECAY
C JAK=6 K   DECAY
C JAK=7 K*  DECAY
C JAK=8 NPI DECAY
C JAK=0 INCLUSIVE:  JAK=1,2,3,4,5,6,7,8
      REAL  H(4)
      REAL*8 HX(4)
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / IDFC  / IDF
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / MDECOMP /IDFFIN(9,NMODE),MULPIK(NMODE)
     +                ,NAMES
      CHARACTER NAMES(NMODE)*31
      COMMON / INOUT / INUT,IOUT
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4),PDUM5(4),HDUM(4)
      REAL  PDUMX(4,9)
      DATA IWARM/0/
      KTOM=KTO
      IF(KTO.EQ.-1) THEN
C     ==================
C       INITIALISATION OR REINITIALISATION
        KTOM=1
        IF (IWARM.EQ.1) X=5/(IWARM-1)
        IWARM=1
        WRITE(IOUT,10000) JAK1,JAK2
        NEVTOT=0
        NEV1=0
        NEV2=0
        IF(JAK1.NE.-1.OR.JAK2.NE.-1) THEN
          CALL DADMEL(-1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DADMMU(-1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DADMPI(-1,IDUM,PDUM,PDUM1,PDUM2)
          CALL DADMRO(-1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4)
          CALL DADMAA(-1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5,JDUM)
          CALL DADMKK(-1,IDUM,PDUM,PDUM1,PDUM2)
          CALL DADMKS(-1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,JDUM)
          CALL DADNEW(-1,IDUM,HDUM,PDUM1,PDUM2,PDUMX,JDUM)
        ENDIF
        DO 10 I=1,30
          NEVDEC(I)=0
          GAMPMC(I)=0
   10   GAMPER(I)=0
      ELSEIF(KTO.EQ.1) THEN
C     =====================
C DECAY OF TAU+ IN THE TAU REST FRAME
        NEVTOT=NEVTOT+1
        IF(IWARM.EQ.0) GOTO 30
        ISGN= IDF/IABS(IDF)
        CALL DEKAY1(0,H,ISGN)
      ELSEIF(KTO.EQ.2) THEN
C     =================================
C DECAY OF TAU- IN THE TAU REST FRAME
        NEVTOT=NEVTOT+1
        IF(IWARM.EQ.0) GOTO 30
        ISGN=-IDF/IABS(IDF)
        CALL DEKAY2(0,H,ISGN)
      ELSEIF(KTO.EQ.11) THEN
C     ======================
C REST OF DECAY PROCEDURE FOR ACCEPTED TAU+ DECAY
        NEV1=NEV1+1
        ISGN= IDF/IABS(IDF)
        CALL DEKAY1(1,H,ISGN)
      ELSEIF(KTO.EQ.12) THEN
C     ======================
C REST OF DECAY PROCEDURE FOR ACCEPTED TAU- DECAY
        NEV2=NEV2+1
        ISGN=-IDF/IABS(IDF)
        CALL DEKAY2(1,H,ISGN)
      ELSEIF(KTO.EQ.100) THEN
C     =======================
        IF(JAK1.NE.-1.OR.JAK2.NE.-1) THEN
          CALL DADMEL( 1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DADMMU( 1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DADMPI( 1,IDUM,PDUM,PDUM1,PDUM2)
          CALL DADMRO( 1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4)
          CALL DADMAA( 1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5,JDUM)
          CALL DADMKK( 1,IDUM,PDUM,PDUM1,PDUM2)
          CALL DADMKS( 1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,JDUM)
          CALL DADNEW( 1,IDUM,HDUM,PDUM1,PDUM2,PDUMX,JDUM)
          WRITE(IOUT,10100) NEV1,NEV2,NEVTOT
          WRITE(IOUT,10200) (NEVDEC(I),GAMPMC(I),GAMPER(I),I= 1,7)
          WRITE(IOUT,10300) (NEVDEC(I),GAMPMC(I),GAMPER(I),NAMES(I-7),
     +    I=8,7+NMODE)
          WRITE(IOUT,10400)
        ENDIF
      ELSE
C     ====
        GOTO 40
      ENDIF
C     =====
      DO 20 K=1,4
   20 HX(K)=H(K)
      RETURN

10000 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'*****TAUOLA LIBRARY: VERSION 2.5 ******',9X,1H*,
     + /,' *',     25X,'***********JUNE     1994***************',9X,1H*,
     + /,' *',     25X,'**AUTHORS: S.JADACH, Z.WAS*************',9X,1H*,
     + /,' *',     25X,'**R. DECKER, M. JEZABEK, J.H.KUEHN*****',9X,1H*,
     + /,' *',     25X,'**AVAILABLE FROM: WASM AT CERNVM ******',9X,1H*,
     + /,' *',     25X,'***** PUBLISHED IN COMP. PHYS. COMM.***',9X,1H*,
     + /,' *',     25X,'*******CERN-TH-5856 SEPTEMBER 1990*****',9X,1H*,
     + /,' *',     25X,'*******CERN-TH-6195 SEPTEMBER 1991*****',9X,1H*,
     + /,' *',     25X,'*******CERN TH-6793 NOVEMBER  1992*****',9X,1H*,
     + /,' *',     25X,'**5 OR MORE PI DEC.: PRECISION LIMITED ',9X,1H*,
     + /,' *',     25X,'****DEKAY ROUTINE: INITIALIZATION******',9X,1H*,
     + /,' *',I20  ,5X,'JAK1   = DECAY MODE TAU+               ',9X,1H*,
     + /,' *',I20  ,5X,'JAK2   = DECAY MODE TAU-               ',9X,1H*,
     +  /,1X,15(5H*****)/)
10100 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'*****TAUOLA LIBRARY: VERSION 2.5 ******',9X,1H*,
     + /,' *',     25X,'***********JUNE     1994***************',9X,1H*,
     + /,' *',     25X,'**AUTHORS: S.JADACH, Z.WAS*************',9X,1H*,
     + /,' *',     25X,'**R. DECKER, M. JEZABEK, J.H.KUEHN*****',9X,1H*,
     + /,' *',     25X,'**AVAILABLE FROM: WASM AT CERNVM ******',9X,1H*,
     + /,' *',     25X,'***** PUBLISHED IN COMP. PHYS. COMM.***',9X,1H*,
     + /,' *',     25X,'*******CERN-TH-5856 SEPTEMBER 1990*****',9X,1H*,
     + /,' *',     25X,'*******CERN-TH-6195 SEPTEMBER 1991*****',9X,1H*,
     + /,' *',     25X,'*******CERN TH-6793 NOVEMBER  1992*****',9X,1H*,
     + /,' *',     25X,'*****DEKAY ROUTINE: FINAL REPORT*******',9X,1H*,
     + /,' *',I20  ,5X,'NEV1   = NO. OF TAU+ DECS. ACCEPTED    ',9X,1H*,
     + /,' *',I20  ,5X,'NEV2   = NO. OF TAU- DECS. ACCEPTED    ',9X,1H*,
     + /,' *',I20  ,5X,'NEVTOT = SUM                           ',9X,1H*,
     + /,' *','    NOEVTS ',
     +   ' PART.WIDTH     ERROR       ROUTINE    DECAY MODE    ',9X,1H*)
10200 FORMAT(1X,'*'
     +       ,I10,2F12.7       ,'     DADMEL     ELECTRON      ',9X,1H*
     + /,' *',I10,2F12.7       ,'     DADMMU     MUON          ',9X,1H*
     + /,' *',I10,2F12.7       ,'     DADMPI     PION          ',9X,1H*
     + /,' *',I10,2F12.7,       '     DADMRO     RHO (->2PI)   ',9X,1H*
     + /,' *',I10,2F12.7,       '     DADMAA     A1  (->3PI)   ',9X,1H*
     + /,' *',I10,2F12.7,       '     DADMKK     KAON          ',9X,1H*
     + /,' *',I10,2F12.7,       '     DADMKS     K*            ',9X,1H*)
10300 FORMAT(1X,'*'
     +       ,I10,2F12.7,A31                                    ,8X,1H*)
10400 FORMAT(1X,'*'
     +       ,20X,'THE ERROR IS RELATIVE AND  PART.WIDTH      ',10X,1H*
     + /,' *',20X,'IN UNITS GFERMI**2*MASS**5/192/PI**3       ',10X,1H*
     +  /,1X,15(5H*****)/)
   30 PRINT 10500
10500 FORMAT(' ----- DEKAY: LACK OF INITIALISATION')
      STOP
   40 PRINT 10600
10600 FORMAT(' ----- DEKAY: WRONG VALUE OF KTO ')
      STOP
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.39  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEKAY1(IMOD,HH,ISGN)
C     *******************************
C THIS ROUTINE  SIMULATES TAU+  DECAY
      COMMON / DECP4 / PP1(4),PP2(4),KF1,KF2
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      REAL  HH(4)
      REAL  HV(4),PNU(4),PPI(4)
      REAL  PWB(4),PMU(4),PNM(4)
      REAL  PRHO(4),PIC(4),PIZ(4)
      REAL  PAA(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PKK(4),PKS(4)
      REAL  PNPI(4,9)
      REAL  PHOT(4)
      REAL  PDUM(4)
      DATA NEV,NPRIN/0,10/
      KTO=1
      IF(JAK1.EQ.-1) RETURN
      IMD=IMOD
      IF(IMD.EQ.0) THEN
C     =================
        JAK=JAK1
        IF(JAK1.EQ.0) CALL JAKER(JAK)
        IF(JAK.EQ.1) THEN
          CALL DADMEL(0, ISGN,HV,PNU,PWB,PMU,PNM,PHOT)
        ELSEIF(JAK.EQ.2) THEN
          CALL DADMMU(0, ISGN,HV,PNU,PWB,PMU,PNM,PHOT)
        ELSEIF(JAK.EQ.3) THEN
          CALL DADMPI(0, ISGN,HV,PPI,PNU)
        ELSEIF(JAK.EQ.4) THEN
          CALL DADMRO(0, ISGN,HV,PNU,PRHO,PIC,PIZ)
        ELSEIF(JAK.EQ.5) THEN
          CALL DADMAA(0, ISGN,HV,PNU,PAA,PIM1,PIM2,PIPL,JAA)
        ELSEIF(JAK.EQ.6) THEN
          CALL DADMKK(0, ISGN,HV,PKK,PNU)
        ELSEIF(JAK.EQ.7) THEN
          CALL DADMKS(0, ISGN,HV,PNU,PKS ,PKK,PPI,JKST)
        ELSE
          CALL DADNEW(0, ISGN,HV,PNU,PWB,PNPI,JAK-7)
        ENDIF
        DO 10 I=1,3
   10   HH(I)=HV(I)
        HH(4)=1.0

      ELSEIF(IMD.EQ.1) THEN
C     =====================
        NEV=NEV+1
        IF (JAK.LT.31) THEN
          NEVDEC(JAK)=NEVDEC(JAK)+1
        ENDIF
        DO 20 I=1,4
   20   PDUM(I)=.0
        IF(JAK.EQ.1) THEN
          CALL DWLUEL(1,ISGN,PNU,PWB,PMU,PNM)
          CALL DWRPH(KTOM,PHOT)
          DO 30 I=1,4
   30     PP1(I)=PMU(I)

        ELSEIF(JAK.EQ.2) THEN
          CALL DWLUMU(1,ISGN,PNU,PWB,PMU,PNM)
          CALL DWRPH(KTOM,PHOT)
          DO 40 I=1,4
   40     PP1(I)=PMU(I)

        ELSEIF(JAK.EQ.3) THEN
          CALL DWLUPI(1,ISGN,PPI,PNU)
          DO 50 I=1,4
   50     PP1(I)=PPI(I)

        ELSEIF(JAK.EQ.4) THEN
          CALL DWLURO(1,ISGN,PNU,PRHO,PIC,PIZ)
          DO 60 I=1,4
   60     PP1(I)=PRHO(I)

        ELSEIF(JAK.EQ.5) THEN
          CALL DWLUAA(1,ISGN,PNU,PAA,PIM1,PIM2,PIPL,JAA)
          DO 70 I=1,4
   70     PP1(I)=PAA(I)
        ELSEIF(JAK.EQ.6) THEN
          CALL DWLUKK(1,ISGN,PKK,PNU)
          DO 80 I=1,4
   80     PP1(I)=PKK(I)
        ELSEIF(JAK.EQ.7) THEN
          CALL DWLUKS(1,ISGN,PNU,PKS,PKK,PPI,JKST)
          DO 90 I=1,4
   90     PP1(I)=PKS(I)
        ELSE
CAM     MULTIPION DECAY
          CALL DWLNEW(1,ISGN,PNU,PWB,PNPI,JAK)
          DO 100 I=1,4
  100     PP1(I)=PWB(I)
        ENDIF

      ENDIF
C     =====
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.39  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEKAY2(IMOD,HH,ISGN)
C     *******************************
C THIS ROUTINE  SIMULATES TAU-  DECAY
      COMMON / DECP4 / PP1(4),PP2(4),KF1,KF2
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      REAL  HH(4)
      REAL  HV(4),PNU(4),PPI(4)
      REAL  PWB(4),PMU(4),PNM(4)
      REAL  PRHO(4),PIC(4),PIZ(4)
      REAL  PAA(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PKK(4),PKS(4)
      REAL  PNPI(4,9)
      REAL  PHOT(4)
      REAL  PDUM(4)
      DATA NEV,NPRIN/0,10/
      KTO=2
      IF(JAK2.EQ.-1) RETURN
      IMD=IMOD
      IF(IMD.EQ.0) THEN
C     =================
        JAK=JAK2
        IF(JAK2.EQ.0) CALL JAKER(JAK)
        IF(JAK.EQ.1) THEN
          CALL DADMEL(0, ISGN,HV,PNU,PWB,PMU,PNM,PHOT)
        ELSEIF(JAK.EQ.2) THEN
          CALL DADMMU(0, ISGN,HV,PNU,PWB,PMU,PNM,PHOT)
        ELSEIF(JAK.EQ.3) THEN
          CALL DADMPI(0, ISGN,HV,PPI,PNU)
        ELSEIF(JAK.EQ.4) THEN
          CALL DADMRO(0, ISGN,HV,PNU,PRHO,PIC,PIZ)
        ELSEIF(JAK.EQ.5) THEN
          CALL DADMAA(0, ISGN,HV,PNU,PAA,PIM1,PIM2,PIPL,JAA)
        ELSEIF(JAK.EQ.6) THEN
          CALL DADMKK(0, ISGN,HV,PKK,PNU)
        ELSEIF(JAK.EQ.7) THEN
          CALL DADMKS(0, ISGN,HV,PNU,PKS ,PKK,PPI,JKST)
        ELSE
          CALL DADNEW(0, ISGN,HV,PNU,PWB,PNPI,JAK-7)
        ENDIF
        DO 10 I=1,3
   10   HH(I)=HV(I)
        HH(4)=1.0
      ELSEIF(IMD.EQ.1) THEN
C     =====================
        NEV=NEV+1
        IF (JAK.LT.31) THEN
          NEVDEC(JAK)=NEVDEC(JAK)+1
        ENDIF
        DO 20 I=1,4
   20   PDUM(I)=.0
        IF(JAK.EQ.1) THEN
          CALL DWLUEL(2,ISGN,PNU,PWB,PMU,PNM)
          CALL DWRPH(KTOM,PHOT)
          DO 30 I=1,4
   30     PP2(I)=PMU(I)

        ELSEIF(JAK.EQ.2) THEN
          CALL DWLUMU(2,ISGN,PNU,PWB,PMU,PNM)
          CALL DWRPH(KTOM,PHOT)
          DO 40 I=1,4
   40     PP2(I)=PMU(I)

        ELSEIF(JAK.EQ.3) THEN
          CALL DWLUPI(2,ISGN,PPI,PNU)
          DO 50 I=1,4
   50     PP2(I)=PPI(I)

        ELSEIF(JAK.EQ.4) THEN
          CALL DWLURO(2,ISGN,PNU,PRHO,PIC,PIZ)
          DO 60 I=1,4
   60     PP2(I)=PRHO(I)

        ELSEIF(JAK.EQ.5) THEN
          CALL DWLUAA(2,ISGN,PNU,PAA,PIM1,PIM2,PIPL,JAA)
          DO 70 I=1,4
   70     PP2(I)=PAA(I)
        ELSEIF(JAK.EQ.6) THEN
          CALL DWLUKK(2,ISGN,PKK,PNU)
          DO 80 I=1,4
   80     PP1(I)=PKK(I)
        ELSEIF(JAK.EQ.7) THEN
          CALL DWLUKS(2,ISGN,PNU,PKS,PKK,PPI,JKST)
          DO 90 I=1,4
   90     PP1(I)=PKS(I)
        ELSE
CAM     MULTIPION DECAY
          CALL DWLNEW(2,ISGN,PNU,PWB,PNPI,JAK)
          DO 100 I=1,4
  100     PP1(I)=PWB(I)
        ENDIF
C
      ENDIF
C     =====
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEXAA(MODE,ISGN,POL,PNU,PAA,PIM1,PIM2,PIPL,JAA)
C ----------------------------------------------------------------------
* THIS SIMULATES TAU DECAY IN TAU REST FRAME
* INTO NU A1, NEXT A1 DECAYS INTO RHO PI AND FINALLY RHO INTO PI PI.
* OUTPUT FOUR MOMENTA: PNU   TAUNEUTRINO,
*                      PAA   A1
*                      PIM1  PION MINUS (OR PI0) 1      (FOR TAU MINUS)
*                      PIM2  PION MINUS (OR PI0) 2
*                      PIPL  PION PLUS  (OR PI-)
*                      (PIPL,PIM1) FORM A RHO
C ----------------------------------------------------------------------
      COMMON / INOUT / INUT,IOUT
      REAL  POL(4),HV(4),PAA(4),PNU(4),PIM1(4),PIM2(4),PIPL(4)
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        CALL DADMAA( -1,ISGN,HV,PNU,PAA,PIM1,PIM2,PIPL,JAA)
CC      CALL HBOOK1(816,'WEIGHT DISTRIBUTION  DEXAA    $',100,-2.,2.)
C
      ELSEIF(MODE.EQ. 0) THEN
*     =======================
   10   CONTINUE
        IF(IWARM.EQ.0) GOTO 20
        CALL DADMAA(  0,ISGN,HV,PNU,PAA,PIM1,PIM2,PIPL,JAA)
        WT=(1+POL(1)*HV(1)+POL(2)*HV(2)+POL(3)*HV(3))/2.
CC      CALL HFILL(816,WT)
        CALL RANMAR(RN,1)
        IF(RN.GT.WT) GOTO 10
C
      ELSEIF(MODE.EQ. 1) THEN
*     =======================
        CALL DADMAA(  1,ISGN,HV,PNU,PAA,PIM1,PIM2,PIPL,JAA)
CC      CALL HPRINT(816)
      ENDIF
C     =====
      RETURN
   20 WRITE(IOUT, 10000)
10000 FORMAT(' ----- DEXAA: LACK OF INITIALISATION')
      STOP
      END
*CMZ :  1.01/50 22/05/96  18.06.08  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.39  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEXAY(KTO,POL)
C ----------------------------------------------------------------------
C THIS 'DEXAY' IS A ROUTINE WHICH GENERATES DECAY OF THE SINGLE
C POLARIZED TAU,  POL IS A POLARIZATION VECTOR (NOT A POLARIMETER
C VECTOR AS IN DEKAY) OF THE TAU AND IT IS AN INPUT PARAMETER.
C KTO=0 INITIALISATION (OBLIGATORY)
C KTO=1 DENOTES TAU+ AND KTO=2 TAU-
C DEXAY(1,POL) AND DEXAY(2,POL) ARE CALLED INTERNALLY BY MC GENERATOR.
C DECAY PRODUCTS ARE TRANSFORMED READILY
C TO CMS AND WRITEN IN THE  LUND RECORD IN /LUJETS/
C KTO=100, PRINT FINAL REPORT (OPTIONAL).
C
C     CALLED BY : KORALZ
C ----------------------------------------------------------------------
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / IDFC  / IDFF
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / MDECOMP /IDFFIN(9,NMODE),MULPIK(NMODE)
     +                ,NAMES
      CHARACTER NAMES(NMODE)*31
      COMMON / INOUT / INUT,IOUT
      REAL  POL(4)
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4),PDUM5(4)
      REAL  PDUM(4)
      REAL  PDUMI(4,9)
      DATA IWARM/0/
      KTOM=KTO
C
      IF(KTO.EQ.-1) THEN
C     ==================
C       INITIALISATION OR REINITIALISATION
        IWARM=1
        WRITE(IOUT, 10000) JAK1,JAK2
        NEVTOT=0
        NEV1=0
        NEV2=0
        IF(JAK1.NE.-1.OR.JAK2.NE.-1) THEN
          CALL DEXEL(-1,IDUM,PDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DEXMU(-1,IDUM,PDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DEXPI(-1,IDUM,PDUM,PDUM1,PDUM2)
          CALL DEXRO(-1,IDUM,PDUM,PDUM1,PDUM2,PDUM3,PDUM4)
          CALL DEXAA(-1,IDUM,PDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5,IDUM)
          CALL DEXKK(-1,IDUM,PDUM,PDUM1,PDUM2)
          CALL DEXKS(-1,IDUM,PDUM,PDUM1,PDUM2,PDUM3,PDUM4,IDUM)
          CALL DEXNEW(-1,IDUM,PDUM,PDUM1,PDUM2,PDUMI,IDUM)
        ENDIF
        DO 10 I=1,30
          NEVDEC(I)=0
          GAMPMC(I)=0
   10   GAMPER(I)=0
      ELSEIF(KTO.EQ.1) THEN
C     =====================
C DECAY OF TAU+ IN THE TAU REST FRAME
        NEVTOT=NEVTOT+1
        NEV1=NEV1+1
        IF(IWARM.EQ.0) GOTO 20
        ISGN=IDFF/IABS(IDFF)
CAM     CALL DEXAY1(POL,ISGN)
        CALL DEXAY1(KTO,JAK1,JAKP,POL,ISGN)
      ELSEIF(KTO.EQ.2) THEN
C     =================================
C DECAY OF TAU- IN THE TAU REST FRAME
        NEVTOT=NEVTOT+1
        NEV2=NEV2+1
        IF(IWARM.EQ.0) GOTO 20
        ISGN=-IDFF/IABS(IDFF)
CAM     CALL DEXAY2(POL,ISGN)
        CALL DEXAY1(KTO,JAK2,JAKM,POL,ISGN)
      ELSEIF(KTO.EQ.100) THEN
C     =======================
        IF(JAK1.NE.-1.OR.JAK2.NE.-1) THEN
          CALL DEXEL( 1,IDUM,PDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DEXMU( 1,IDUM,PDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DEXPI( 1,IDUM,PDUM,PDUM1,PDUM2)
          CALL DEXRO( 1,IDUM,PDUM,PDUM1,PDUM2,PDUM3,PDUM4)
          CALL DEXAA( 1,IDUM,PDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5,IDUM)
          CALL DEXKK( 1,IDUM,PDUM,PDUM1,PDUM2)
          CALL DEXKS( 1,IDUM,PDUM,PDUM1,PDUM2,PDUM3,PDUM4,IDUM)
          CALL DEXNEW( 1,IDUM,PDUM,PDUM1,PDUM2,PDUMI,IDUM)
          WRITE(IOUT,10100) NEV1,NEV2,NEVTOT
          WRITE(IOUT,10200) (NEVDEC(I),GAMPMC(I),GAMPER(I),I= 1,7)
          WRITE(IOUT,10300) (NEVDEC(I),GAMPMC(I),GAMPER(I),NAMES(I-7),
     +    I=8,7+NMODE)
          WRITE(IOUT,10400)
        ENDIF
      ELSE
        GOTO 30
      ENDIF
      RETURN
10000 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'*****TAUOLA LIBRARY: VERSION 2.5 ******',9X,1H*,
     + /,' *',     25X,'***********JUNE     1994***************',9X,1H*,
     + /,' *',     25X,'**AUTHORS: S.JADACH, Z.WAS*************',9X,1H*,
     + /,' *',     25X,'**R. DECKER, M. JEZABEK, J.H.KUEHN*****',9X,1H*,
     + /,' *',     25X,'**AVAILABLE FROM: WASM AT CERNVM ******',9X,1H*,
     + /,' *',     25X,'***** PUBLISHED IN COMP. PHYS. COMM.***',9X,1H*,
     + /,' *',     25X,'*******CERN-TH-5856 SEPTEMBER 1990*****',9X,1H*,
     + /,' *',     25X,'*******CERN-TH-6195 SEPTEMBER 1991*****',9X,1H*,
     + /,' *',     25X,'*******CERN-TH-6793 NOVEMBER  1992*****',9X,1H*,
     + /,' *',     25X,'**5 OR MORE PI DEC.: PRECISION LIMITED ',9X,1H*,
     + /,' *',     25X,'******DEXAY ROUTINE: INITIALIZATION****',9X,1H*
     + /,' *',I20  ,5X,'JAK1   = DECAY MODE FERMION1 (TAU+)    ',9X,1H*
     + /,' *',I20  ,5X,'JAK2   = DECAY MODE FERMION2 (TAU-)    ',9X,1H*
     +  /,1X,15(5H*****)/)
CHBU  FORMAT 7010 HAD MORE THAN 19 CONTINUATION LINES
CHBU  SPLIT INTO TWO
10100 FORMAT(///1X,15(5H*****)
     + /,' *',     25X,'*****TAUOLA LIBRARY: VERSION 2.5 ******',9X,1H*,
     + /,' *',     25X,'***********JUNE     1994***************',9X,1H*,
     + /,' *',     25X,'**AUTHORS: S.JADACH, Z.WAS*************',9X,1H*,
     + /,' *',     25X,'**R. DECKER, M. JEZABEK, J.H.KUEHN*****',9X,1H*,
     + /,' *',     25X,'**AVAILABLE FROM: WASM AT CERNVM ******',9X,1H*,
     + /,' *',     25X,'***** PUBLISHED IN COMP. PHYS. COMM.***',9X,1H*,
     + /,' *',     25X,'*******CERN-TH-5856 SEPTEMBER 1990*****',9X,1H*,
     + /,' *',     25X,'*******CERN-TH-6195 SEPTEMBER 1991*****',9X,1H*,
     + /,' *',     25X,'*******CERN-TH-6793 NOVEMBER  1992*****',9X,1H*,
     + /,' *',     25X,'******DEXAY ROUTINE: FINAL REPORT******',9X,1H*
     + /,' *',I20  ,5X,'NEV1   = NO. OF TAU+ DECS. ACCEPTED    ',9X,1H*
     + /,' *',I20  ,5X,'NEV2   = NO. OF TAU- DECS. ACCEPTED    ',9X,1H*
     + /,' *',I20  ,5X,'NEVTOT = SUM                           ',9X,1H*
     + /,' *','    NOEVTS ',
     +   ' PART.WIDTH     ERROR       ROUTINE    DECAY MODE    ',9X,1H*)
10200 FORMAT(1X,'*'
     +       ,I10,2F12.7       ,'     DADMEL     ELECTRON      ',9X,1H*
     + /,' *',I10,2F12.7       ,'     DADMMU     MUON          ',9X,1H*
     + /,' *',I10,2F12.7       ,'     DADMPI     PION          ',9X,1H*
     + /,' *',I10,2F12.7,       '     DADMRO     RHO (->2PI)   ',9X,1H*
     + /,' *',I10,2F12.7,       '     DADMAA     A1  (->3PI)   ',9X,1H*
     + /,' *',I10,2F12.7,       '     DADMKK     KAON          ',9X,1H*
     + /,' *',I10,2F12.7,       '     DADMKS     K*            ',9X,1H*)
10300 FORMAT(1X,'*'
     +       ,I10,2F12.7,A31                                    ,8X,1H*)
10400 FORMAT(1X,'*'
     +       ,20X,'THE ERROR IS RELATIVE AND  PART.WIDTH      ',10X,1H*
     + /,' *',20X,'IN UNITS GFERMI**2*MASS**5/192/PI**3       ',10X,1H*
     +  /,1X,15(5H*****)/)
   20 WRITE(IOUT, 10500)
10500 FORMAT(' ----- DEXAY: LACK OF INITIALISATION')
      STOP
   30 WRITE(IOUT, 10600)
10600 FORMAT(' ----- DEXAY: WRONG VALUE OF KTO ')
      STOP
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.39  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEXAY1(KTO,JAKIN,JAK,POL,ISGN)
C ---------------------------------------------------------------------
C THIS ROUTINE  SIMULATES TAU+-  DECAY
C
C     CALLED BY : DEXAY
C ---------------------------------------------------------------------
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      REAL  POL(4),POLAR(4)
      REAL  PNU(4),PPI(4)
      REAL  PRHO(4),PIC(4),PIZ(4)
      REAL  PWB(4),PMU(4),PNM(4)
      REAL  PAA(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PKK(4),PKS(4)
      REAL  PNPI(4,9)
      REAL PHOT(4)
      REAL PDUM(4)
C
      IF(JAKIN.EQ.-1) RETURN
      DO 10 I=1,3
   10 POLAR(I)=POL(I)
      POLAR(4)=0.
      DO 20 I=1,4
   20 PDUM(I)=.0
      JAK=JAKIN
      IF(JAK.EQ.0) CALL JAKER(JAK)
CAM
      IF(JAK.EQ.1) THEN
        CALL DEXEL(0, ISGN,POLAR,PNU,PWB,PMU,PNM,PHOT)
        CALL DWLUEL(KTO,ISGN,PNU,PWB,PMU,PNM)
        CALL DWRPH(KTO,PHOT )
      ELSEIF(JAK.EQ.2) THEN
        CALL DEXMU(0, ISGN,POLAR,PNU,PWB,PMU,PNM,PHOT)
        CALL DWLUMU(KTO,ISGN,PNU,PWB,PMU,PNM)
        CALL DWRPH(KTO,PHOT )
      ELSEIF(JAK.EQ.3) THEN
        CALL DEXPI(0, ISGN,POLAR,PPI,PNU)
        CALL DWLUPI(KTO,ISGN,PPI,PNU)
      ELSEIF(JAK.EQ.4) THEN
        CALL DEXRO(0, ISGN,POLAR,PNU,PRHO,PIC,PIZ)
        CALL DWLURO(KTO,ISGN,PNU,PRHO,PIC,PIZ)
      ELSEIF(JAK.EQ.5) THEN
        CALL DEXAA(0, ISGN,POLAR,PNU,PAA,PIM1,PIM2,PIPL,JAA)
        CALL DWLUAA(KTO,ISGN,PNU,PAA,PIM1,PIM2,PIPL,JAA)
      ELSEIF(JAK.EQ.6) THEN
        CALL DEXKK(0, ISGN,POLAR,PKK,PNU)
        CALL DWLUKK(KTO,ISGN,PKK,PNU)
      ELSEIF(JAK.EQ.7) THEN
        CALL DEXKS(0, ISGN,POLAR,PNU,PKS,PKK,PPI,JKST)
        CALL DWLUKS(KTO,ISGN,PNU,PKS,PKK,PPI,JKST)
      ELSE
        JNPI=JAK-7
        CALL DEXNEW(0, ISGN,POLAR,PNU,PWB,PNPI,JNPI)
        CALL DWLNEW(KTO,ISGN,PNU,PWB,PNPI,JAK)
      ENDIF
      NEVDEC(JAK)=NEVDEC(JAK)+1
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.39  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEXEL(MODE,ISGN,POL,PNU,PWB,Q1,Q2,PH)
C ----------------------------------------------------------------------
C THIS SIMULATES TAU DECAY IN TAU REST FRAME
C INTO ELECTRON AND TWO NEUTRINOS
C
C     CALLED BY : DEXAY,DEXAY1
C ----------------------------------------------------------------------
      REAL  POL(4),HV(4),PWB(4),PNU(4),Q1(4),Q2(4),PH(4)
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        CALL DADMEL( -1,ISGN,HV,PNU,PWB,Q1,Q2,PH)
CC      CALL HBOOK1(813,'WEIGHT DISTRIBUTION  DEXEL    $',100,0,2)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
   10   CONTINUE
        IF(IWARM.EQ.0) GOTO 20
        CALL DADMEL(  0,ISGN,HV,PNU,PWB,Q1,Q2,PH)
        WT=(1+POL(1)*HV(1)+POL(2)*HV(2)+POL(3)*HV(3))/2.
CC      CALL HFILL(813,WT)
        CALL RANMAR(RN,1)
        IF(RN.GT.WT) GOTO 10
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        CALL DADMEL(  1,ISGN,HV,PNU,PWB,Q1,Q2,PH)
CC      CALL HPRINT(813)
      ENDIF
C     =====
      RETURN
   20 PRINT 10000
10000 FORMAT(' ----- DEXEL: LACK OF INITIALISATION')
      STOP
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEXKK(MODE,ISGN,POL,PKK,PNU)
C ----------------------------------------------------------------------
C TAU DECAY INTO KAON  AND TAU-NEUTRINO
C IN TAU REST FRAME
C OUTPUT FOUR MOMENTA: PNU   TAUNEUTRINO,
C                      PKK   KAON CHARGED
C ----------------------------------------------------------------------
      REAL  POL(4),HV(4),PNU(4),PKK(4)
C
      IF(MODE.EQ.-1) THEN
C     ===================
        CALL DADMKK(-1,ISGN,HV,PKK,PNU)
CC      CALL HBOOK1(815,'WEIGHT DISTRIBUTION  DEXPI    $',100,0,2)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
   10   CONTINUE
        CALL DADMKK( 0,ISGN,HV,PKK,PNU)
        WT=(1+POL(1)*HV(1)+POL(2)*HV(2)+POL(3)*HV(3))/2.
CC      CALL HFILL(815,WT)
        CALL RANMAR(RN,1)
        IF(RN.GT.WT) GOTO 10
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        CALL DADMKK( 1,ISGN,HV,PKK,PNU)
CC      CALL HPRINT(815)
      ENDIF
C     =====
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEXKS(MODE,ISGN,POL,PNU,PKS,PKK,PPI,JKST)
C ----------------------------------------------------------------------
C THIS SIMULATES TAU DECAY IN TAU REST FRAME
C INTO NU K*, THEN K* DECAYS INTO PI0,K+-(JKST=20)
C OR PI+-,K0(JKST=10).
C OUTPUT FOUR MOMENTA: PNU   TAUNEUTRINO,
C                      PKS   K* CHARGED
C                      PK0   K ZERO
C                      PKC   K CHARGED
C                      PIC   PION CHARGED
C                      PIZ   PION ZERO
C ----------------------------------------------------------------------
      COMMON / INOUT / INUT,IOUT
      REAL  POL(4),HV(4),PKS(4),PNU(4),PKK(4),PPI(4)
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
CFZ INITIALISATION DONE WITH THE GHARGED PION NEUTRAL KAON MODE(JKST=10
        CALL DADMKS( -1,ISGN,HV,PNU,PKS,PKK,PPI,JKST)
CC      CALL HBOOK1(816,'WEIGHT DISTRIBUTION  DEXKS    $',100,0,2)
CC      CALL HBOOK1(916,'ABS2 OF HV IN ROUTINE DEXKS   $',100,0,2)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
   10   CONTINUE
        IF(IWARM.EQ.0) GOTO 20
        CALL DADMKS(  0,ISGN,HV,PNU,PKS,PKK,PPI,JKST)
        WT=(1+POL(1)*HV(1)+POL(2)*HV(2)+POL(3)*HV(3))/2.
CC      CALL HFILL(816,WT)
CC      XHELP=HV(1)**2+HV(2)**2+HV(3)**2
CC      CALL HFILL(916,XHELP)
        CALL RANMAR(RN,1)
        IF(RN.GT.WT) GOTO 10
C
      ELSEIF(MODE.EQ. 1) THEN
C     ======================================
        CALL DADMKS( 1,ISGN,HV,PNU,PKS,PKK,PPI,JKST)
CC      CALL HPRINT(816)
CC      CALL HPRINT(916)
      ENDIF
C     =====
      RETURN
   20 WRITE(IOUT, 10000)
10000 FORMAT(' ----- DEXKS: LACK OF INITIALISATION')
      STOP
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.39  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEXMU(MODE,ISGN,POL,PNU,PWB,Q1,Q2,PH)
C ----------------------------------------------------------------------
C THIS SIMULATES TAU DECAY IN ITS REST FRAME
C INTO MUON AND TWO NEUTRINOS
C OUTPUT FOUR MOMENTA: PNU   TAUNEUTRINO,
C                      PWB   W-BOSON
C                      Q1    MUON
C                      Q2    MUON-NEUTRINO
C ----------------------------------------------------------------------
      COMMON / INOUT / INUT,IOUT
      REAL  POL(4),HV(4),PWB(4),PNU(4),Q1(4),Q2(4),PH(4)
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        CALL DADMMU( -1,ISGN,HV,PNU,PWB,Q1,Q2,PH)
CC      CALL HBOOK1(814,'WEIGHT DISTRIBUTION  DEXMU    $',100,0,2)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
   10   CONTINUE
        IF(IWARM.EQ.0) GOTO 20
        CALL DADMMU(  0,ISGN,HV,PNU,PWB,Q1,Q2,PH)
        WT=(1+POL(1)*HV(1)+POL(2)*HV(2)+POL(3)*HV(3))/2.
CC      CALL HFILL(814,WT)
        CALL RANMAR(RN,1)
        IF(RN.GT.WT) GOTO 10
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        CALL DADMMU(  1,ISGN,HV,PNU,PWB,Q1,Q2,PH)
CC      CALL HPRINT(814)
      ENDIF
C     =====
      RETURN
   20 WRITE(IOUT, 10000)
10000 FORMAT(' ----- DEXMU: LACK OF INITIALISATION')
      STOP
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEXNEW(MODE,ISGN,POL,PNU,PAA,PNPI,JNPI)
C ----------------------------------------------------------------------
* THIS SIMULATES TAU DECAY IN TAU REST FRAME
* INTO NU A1, NEXT A1 DECAYS INTO RHO PI AND FINALLY RHO INTO PI PI.
* OUTPUT FOUR MOMENTA: PNU   TAUNEUTRINO,
*                      PAA   A1
*                      PIM1  PION MINUS (OR PI0) 1      (FOR TAU MINUS)
*                      PIM2  PION MINUS (OR PI0) 2
*                      PIPL  PION PLUS  (OR PI-)
*                      (PIPL,PIM1) FORM A RHO
C ----------------------------------------------------------------------
      COMMON / INOUT / INUT,IOUT
      REAL  POL(4),HV(4),PAA(4),PNU(4),PNPI(4,9)
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        CALL DADNEW( -1,ISGN,HV,PNU,PAA,PNPI,JDUMM)
CC      CALL HBOOK1(816,'WEIGHT DISTRIBUTION  DEXAA    $',100,-2.,2.)
C
      ELSEIF(MODE.EQ. 0) THEN
*     =======================
   10   CONTINUE
        IF(IWARM.EQ.0) GOTO 20
        CALL DADNEW( 0,ISGN,HV,PNU,PAA,PNPI,JNPI)
        WT=(1+POL(1)*HV(1)+POL(2)*HV(2)+POL(3)*HV(3))/2.
CC      CALL HFILL(816,WT)
        CALL RANMAR(RN,1)
        IF(RN.GT.WT) GOTO 10
C
      ELSEIF(MODE.EQ. 1) THEN
*     =======================
        CALL DADNEW( 1,ISGN,HV,PNU,PAA,PNPI,JDUMM)
CC      CALL HPRINT(816)
      ENDIF
C     =====
      RETURN
   20 WRITE(IOUT, 10000)
10000 FORMAT(' ----- DEXNEW: LACK OF INITIALISATION')
      STOP
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEXPI(MODE,ISGN,POL,PPI,PNU)
C ----------------------------------------------------------------------
C TAU DECAY INTO PION AND TAU-NEUTRINO
C IN TAU REST FRAME
C OUTPUT FOUR MOMENTA: PNU   TAUNEUTRINO,
C                      PPI   PION CHARGED
C ----------------------------------------------------------------------
      REAL  POL(4),HV(4),PNU(4),PPI(4)
CC
      IF(MODE.EQ.-1) THEN
C     ===================
        CALL DADMPI(-1,ISGN,HV,PPI,PNU)
CC      CALL HBOOK1(815,'WEIGHT DISTRIBUTION  DEXPI    $',100,0,2)

      ELSEIF(MODE.EQ. 0) THEN
C     =======================
   10   CONTINUE
        CALL DADMPI( 0,ISGN,HV,PPI,PNU)
        WT=(1+POL(1)*HV(1)+POL(2)*HV(2)+POL(3)*HV(3))/2.
CC      CALL HFILL(815,WT)
        CALL RANMAR(RN,1)
        IF(RN.GT.WT) GOTO 10
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        CALL DADMPI( 1,ISGN,HV,PPI,PNU)
CC      CALL HPRINT(815)
      ENDIF
C     =====
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DEXRO(MODE,ISGN,POL,PNU,PRO,PIC,PIZ)
C ----------------------------------------------------------------------
C THIS SIMULATES TAU DECAY IN TAU REST FRAME
C INTO NU RHO, NEXT RHO DECAYS INTO PION PAIR.
C OUTPUT FOUR MOMENTA: PNU   TAUNEUTRINO,
C                      PRO   RHO
C                      PIC   PION CHARGED
C                      PIZ   PION ZERO
C ----------------------------------------------------------------------
      COMMON / INOUT / INUT,IOUT
      REAL  POL(4),HV(4),PRO(4),PNU(4),PIC(4),PIZ(4)
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        CALL DADMRO( -1,ISGN,HV,PNU,PRO,PIC,PIZ)
CC      CALL HBOOK1(816,'WEIGHT DISTRIBUTION  DEXRO    $',100,0,2)
CC      CALL HBOOK1(916,'ABS2 OF HV IN ROUTINE DEXRO   $',100,0,2)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
   10   CONTINUE
        IF(IWARM.EQ.0) GOTO 20
        CALL DADMRO(  0,ISGN,HV,PNU,PRO,PIC,PIZ)
        WT=(1+POL(1)*HV(1)+POL(2)*HV(2)+POL(3)*HV(3))/2.
CC      CALL HFILL(816,WT)
CC      XHELP=HV(1)**2+HV(2)**2+HV(3)**2
CC      CALL HFILL(916,XHELP)
        CALL RANMAR(RN,1)
        IF(RN.GT.WT) GOTO 10
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        CALL DADMRO(  1,ISGN,HV,PNU,PRO,PIC,PIZ)
CC      CALL HPRINT(816)
CC      CALL HPRINT(916)
      ENDIF
C     =====
      RETURN
   20 WRITE(IOUT, 10000)
10000 FORMAT(' ----- DEXRO: LACK OF INITIALISATION')
      STOP
      END
*CMZ :  1.02/00 12/01/97  15.00.23  by  J. Brunner
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      DOUBLE PRECISION FUNCTION DFUN(NDIM,X)
      INTEGER NDIM
      DOUBLE PRECISION X(NDIM)
      DFUN=RIWFUN(X)
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION DILOG(X)
C     *****************
      IMPLICIT REAL*8(A-H,O-Z)
CERN      C304      VERSION    29/07/71 DILOG        59                C
      Z=-1.64493406684822
      IF(X .LT.-1.0) GO TO 10
      IF(X .LE. 0.5) GO TO 20
      IF(X .EQ. 1.0) GO TO 30
      IF(X .LE. 2.0) GO TO 40
      Z=3.2898681336964
   10 T=1.0/X
      S=-0.5
      Z=Z-0.5* LOG(ABS(X))**2
      GO TO 50
   20 T=X
      S=0.5
      Z=0.
      GO TO 50
   30 DILOG=1.64493406684822
      RETURN
   40 T=1.0-X
      S=-0.5
      Z=1.64493406684822 - LOG(X)* LOG(ABS(T))
   50 Y=2.66666666666666 *T+0.66666666666666
      B=      0.00000 00000 00001
      A=Y*B  +0.00000 00000 00004
      B=Y*A-B+0.00000 00000 00011
      A=Y*B-A+0.00000 00000 00037
      B=Y*A-B+0.00000 00000 00121
      A=Y*B-A+0.00000 00000 00398
      B=Y*A-B+0.00000 00000 01312
      A=Y*B-A+0.00000 00000 04342
      B=Y*A-B+0.00000 00000 14437
      A=Y*B-A+0.00000 00000 48274
      B=Y*A-B+0.00000 00001 62421
      A=Y*B-A+0.00000 00005 50291
      B=Y*A-B+0.00000 00018 79117
      A=Y*B-A+0.00000 00064 74338
      B=Y*A-B+0.00000 00225 36705
      A=Y*B-A+0.00000 00793 87055
      B=Y*A-B+0.00000 02835 75385
      A=Y*B-A+0.00000 10299 04264
      B=Y*A-B+0.00000 38163 29463
      A=Y*B-A+0.00001 44963 00557
      B=Y*A-B+0.00005 68178 22718
      A=Y*B-A+0.00023 20021 96094
      B=Y*A-B+0.00100 16274 96164
      A=Y*B-A+0.00468 63619 59447
      B=Y*A-B+0.02487 93229 24228
      A=Y*B-A+0.16607 30329 27855
      A=Y*A-B+1.93506 43008 6996
      DILOG=S*T*(A-B)+Z
      RETURN
C=======================================================================
C===================END OF CPC PART ====================================
C=======================================================================
      END
*CMZ :  1.01/50 29/02/96  09.49.49  by  Piero Zucchelli
*CMZ :  1.01/43 15/12/95  18.02.33  by  Piero Zucchelli
*CMZ :  1.01/40 09/11/95  16.09.04  by  Piero Zucchelli
*CMZ :  1.01/39 20/10/95  18.27.32  by  Piero Zucchelli
*CMZ :  1.01/32 02/06/95  20.27.41  BY  PIERO ZUCCHELLI
*CMZ :  1.01/31 02/06/95  20.17.17  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.35.13  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 23/09/94  12.02.06  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 15/08/94  07.15.40  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 11/07/94  09.26.54  BY  PIERO ZUCCHELLI
*-- AUTHOR :    PIERO ZUCCHELLI   04/07/94

      REAL*4 FUNCTION DISTRR(DUMMY)
*KEEP,JETTA.
C--
         PARAMETER (ICENTO=100)
         PARAMETER (ISIZ=93)
         PARAMETER (IOF1=32)
         PARAMETER (IOF2=83)
         PARAMETER (LUX_LEVEL=4)
         INTEGER*4 JTAU(100),JPRI(100),JSTRO(100)
         REAL*4 FTUPLE(ISIZ)
         COMMON/JETTAGL/JTAU,JPRI,JSTRO
         COMMON/NTUPLA/FTUPLE,ISFIRST
         COMMON/BEAM/SPEC(ICENTO)
         COMMON /MAXSPEC/RMAXSPEC,RINTSPEC
         COMMON/SAV/XMINSAV(ICENTO),XMAXSAV(ICENTO),YMINSAV(ICENTO),
     &   YMAXSAV(ICENTO),Q2MINSAV(ICENTO),Q2MAXSAV(ICENTO),
     &   W2MINSAV(ICENTO),W2MAXSAV(ICENTO),PARIMAX(ICENTO),
     &   PPSAVE(ICENTO,3,4,5),PARICOR(ICENTO),INDEX,SIGMASAV(ICENTO),
     &   XMSIGMA,XSECT

*KEEP,KEYS.
        COMMON/CFREAD/SPACE(5000)
 	COMMON/MYKEYS/IEVAR,IF5CC,INEUT,IINTE,IFERM,IFLAT,ICOUN,
     &  REFIX,IMUDO,IKAT1,IKAT2,IKAT3,IKAT4,IKAT5,IKAT6,
     &  INEVT,IJAK1,IJAK2,IITDK,RPTAU,RXK0D,NTGR,IDIMUON,IGLU,
     &  IQDEN,LOME(2),IFILES,ICCHA,ISEED,IDSUBS,IONEM,EHAC,IPSEL,
     &  IHIST


*KEND.
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)

10    CONTINUE
      DISTR=DUMMY
* now we have to understand to what index corresponds
*     ENE=(I-1)*3. + 1.5

      ITEST=(DUMMY-1.5)/3.+1
* this is the lower bin, we have to determine its xsection
* and then accept/reject it ,
* if accepted we fill the numbers....
      DICE=RNDMM(ISEED)*XMSIGMA
      ENE1=(ITEST-1)*3. + 1.5
      ENE2=(ITEST)*3. + 1.5
      IF (ITEST.EQ.0) ENE1=0
      IF (ITEST.NE.0) THEN
        S1=SIGMASAV(ITEST)
      ELSE
        S1=0
      ENDIF
      S2=SIGMASAV(ITEST+1)
      XSECT=(DUMMY-ENE1)/(ENE2-ENE1)*(S2-S1)+S1


      IF (DICE.GT.XSECT) THEN
        DISTR=0.
        RETURN
      ENDIF

      INDEX=ITEST

      IF (INDEX.LE.0.OR.INDEX.GT.ICENTO) THEN
        WRITE(*,*)'+++DISTR WARNING: OUTOFRANGE',INDEX
        DISTR=0.
        RETURN
      ENDIF

      PARI(32)=PARICOR(INDEX)
      PARI(LST(23))=PARIMAX(INDEX)
      XMIN=XMINSAV(INDEX)
      XMAX=XMAXSAV(INDEX)
      YMIN=YMINSAV(INDEX)
      YMAX=YMAXSAV(INDEX)
      Q2MIN=Q2MINSAV(INDEX)
      Q2MAX=Q2MAXSAV(INDEX)
      W2MIN=W2MINSAV(INDEX)
      W2MAX=W2MAXSAV(INDEX)
      DO 20 IA=1,2
        DO 20 JA=1,5
   20 PSAVE(3,IA,JA)=PPSAVE(INDEX,3,IA,JA)

      RETURN
      END
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      FUNCTION DLOWER(V1)

C...LOWER LIMIT ON SECOND VARIABLE (Y, Q**2 OR W**2) DEPENDING ON FIRST
C...VARIABLE X=V1. USED FOR INTEGRATING DIFFERENTIAL CROSS-SECTION.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
C...CMS ENERGY SQUARED AND TARGET NUCLEON MASS.
      S=PARL(21)
      PM2=PSAVE(3,2,5)**2
      IF(LST(31).EQ.1) THEN
        DLOWER=MAX(Q2MIN,V1*YMIN*S,(W2MIN-PM2)*V1/MAX(1.-V1,1.E-22))
      ELSEIF(LST(31).EQ.2) THEN
        DLOWER=MAX(YMIN,Q2MIN/(S*V1),(W2MIN-PM2)/MAX(S*(1.-V1),1.E-22))
      ELSEIF(LST(31).EQ.3) THEN
        DLOWER=MAX(W2MIN,(1.-V1)*YMIN*S+PM2,
     &  Q2MIN*(1.-V1)/MAX(V1,1.E-22)+PM2)
      ENDIF
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DPH4PI(DGAMT,HV,PN,PAA,PMULT,JNPI)
C ----------------------------------------------------------------------
* IT SIMULATES A1  DECAY IN TAU REST FRAME WITH
* Z-AXIS ALONG A1  MOMENTUM
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL  HV(4),PT(4),PN(4),PAA(4),PIM1(4),PIM2(4),PIPL(4),PMULT(4,9)
      REAL  PR(4),PIZ(4)
      REAL*4 RRR(9)
      REAL*8 UU,FF,FF1,FF2,FF3,FF4,GG1,GG2,GG3,GG4,RR
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
      XLAM(X,Y,Z)=SQRT(ABS((X-Y-Z)**2-4.0*Y*Z))
C AMRO, GAMRO IS ONLY A PARAMETER FOR GETING HIGHT EFFICIENCY
C
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1./2**23/PI**11
      PHSP=1./2**5/PI**2
      IF (JNPI.EQ.1) THEN
        PREZ=0.7
        AMP1=AMPI
        AMP2=AMPI
        AMP3=AMPI
        AMP4=AMPIZ
        AMRX=0.782
        GAMRX=0.0084
        AMROP =1.2
        GAMROP=.46

      ELSE
        PREZ=0.0
        AMP1=AMPIZ
        AMP2=AMPIZ
        AMP3=AMPIZ
        AMP4=AMPI
        AMRX=1.4
        GAMRX=.6
        AMROP =AMRX
        GAMROP=GAMRX

      ENDIF
      RR=0.3
      CALL CHOICE(100+JNPI,RR,ICHAN,PROB1,PROB2,PROB3,
     +            AMROP,GAMROP,AMRX,GAMRX,AMRB,GAMRB)
      PREZ=PROB1+PROB2
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C
      CALL RANMAR(RRR,9)
C
* MASSES OF 4, 3 AND 2 PI SYSTEMS
C 3 PI WITH SAMPLING FOR RESONANCE
CAM
      RR1=RRR(6)
      AMS1=(AMP1+AMP2+AMP3+AMP4)**2
      AMS2=(AMTAU-AMNUTA)**2
      ALP1=ATAN((AMS1-AMROP**2)/AMROP/GAMROP)
      ALP2=ATAN((AMS2-AMROP**2)/AMROP/GAMROP)
      ALP=ALP1+RR1*(ALP2-ALP1)
      AM4SQ =AMROP**2+AMROP*GAMROP*TAN(ALP)
      AM4 =SQRT(AM4SQ)
      PHSPAC=PHSPAC* ((AM4SQ-AMROP**2)**2+(AMROP*GAMROP)**2)/(AMROP*
     +GAMROP)
      PHSPAC=PHSPAC*(ALP2-ALP1)

C
      RR1=RRR(1)
      AMS1=(AMP2+AMP3+AMP4)**2
      AMS2=(AM4-AMP1)**2
      IF (RRR(9).GT.PREZ) THEN
        AM3SQ=AMS1+ RR1*(AMS2-AMS1)
        AM3 =SQRT(AM3SQ)
C --- THIS PART OF JACOBIAN WILL BE RECOVERED LATER
        FF1=AMS2-AMS1
      ELSE
* PHASE SPACE WITH SAMPLING FOR OMEGA RESONANCE,
        ALP1=ATAN((AMS1-AMRX**2)/AMRX/GAMRX)
        ALP2=ATAN((AMS2-AMRX**2)/AMRX/GAMRX)
        ALP=ALP1+RR1*(ALP2-ALP1)
        AM3SQ =AMRX**2+AMRX*GAMRX*TAN(ALP)
        AM3 =SQRT(AM3SQ)
C --- THIS PART OF THE JACOBIAN WILL BE RECOVERED LATER ---------------
        FF1=((AM3SQ-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)
        FF1=FF1*(ALP2-ALP1)
      ENDIF
C MASS OF 2
      RR2=RRR(2)
      AMS1=(AMP3+AMP4)**2
      AMS2=(AM3-AMP2)**2
* FLAT PHASE SPACE;
      AM2SQ=AMS1+ RR2*(AMS2-AMS1)
      AM2 =SQRT(AM2SQ)
C --- THIS PART OF JACOBIAN WILL BE RECOVERED LATER
      FF2=(AMS2-AMS1)
*  2 RESTFRAME, DEFINE PIZ AND PIPL
      ENQ1=(AM2SQ-AMP3**2+AMP4**2)/(2*AM2)
      ENQ2=(AM2SQ+AMP3**2-AMP4**2)/(2*AM2)
      PPI= ENQ1**2-AMP4**2
      PPPI=SQRT(ABS(ENQ1**2-AMP4**2))
      PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM2)
* PIZ   MOMENTUM IN 2 REST FRAME
      CALL SPHERA(PPPI,PIZ)
      PIZ(4)=ENQ1
* PIPL  MOMENTUM IN 2 REST FRAME
      DO 10 I=1,3
   10 PIPL(I)=-PIZ(I)
      PIPL(4)=ENQ2
* 3 REST FRAME, DEFINE PIM1
*       PR   MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AM3)*(AM3**2+AM2**2-AMP2**2)
      PR(3)= SQRT(ABS(PR(4)**2-AM2**2))
      PPI = PR(4)**2-AM2**2
*       PIM1  MOMENTUM
      PIM1(1)=0
      PIM1(2)=0
      PIM1(4)=1./(2*AM3)*(AM3**2-AM2**2+AMP2**2)
      PIM1(3)=-PR(3)
C --- THIS PART OF JACOBIAN WILL BE RECOVERED LATER
      FF3=(4*PI)*(2*PR(3)/AM3)
* OLD PIONS BOOSTED FROM 2 REST FRAME TO 3 REST FRAME
      EXE=(PR(4)+PR(3))/AM2
      CALL BOSTR3(EXE,PIZ,PIZ)
      CALL BOSTR3(EXE,PIPL,PIPL)
      RR3=RRR(3)
      RR4=RRR(4)
      THET =ACOS(-1.+2*RR3)
      PHI = 2*PI*RR4
      CALL ROTPOL(THET,PHI,PIPL)
      CALL ROTPOL(THET,PHI,PIM1)
      CALL ROTPOL(THET,PHI,PIZ)
      CALL ROTPOL(THET,PHI,PR)
* 4  REST FRAME, DEFINE PIM2
*       PR   MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AM4)*(AM4**2+AM3**2-AMP1**2)
      PR(3)= SQRT(ABS(PR(4)**2-AM3**2))
      PPI = PR(4)**2-AM3**2
*       PIM2 MOMENTUM
      PIM2(1)=0
      PIM2(2)=0
      PIM2(4)=1./(2*AM4)*(AM4**2-AM3**2+AMP1**2)
      PIM2(3)=-PR(3)
C --- THIS PART OF JACOBIAN WILL BE RECOVERED LATER
      FF4=(4*PI)*(2*PR(3)/AM4)
* OLD PIONS BOOSTED FROM 3 REST FRAME TO 4 REST FRAME
      EXE=(PR(4)+PR(3))/AM3
      CALL BOSTR3(EXE,PIZ,PIZ)
      CALL BOSTR3(EXE,PIPL,PIPL)
      CALL BOSTR3(EXE,PIM1,PIM1)
      RR3=RRR(7)
      RR4=RRR(8)
      THET =ACOS(-1.+2*RR3)
      PHI = 2*PI*RR4
      CALL ROTPOL(THET,PHI,PIPL)
      CALL ROTPOL(THET,PHI,PIM1)
      CALL ROTPOL(THET,PHI,PIM2)
      CALL ROTPOL(THET,PHI,PIZ)
      CALL ROTPOL(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE PAA AND NEUTRINO MOMENTA
* PAA  MOMENTUM
      PAA(1)=0
      PAA(2)=0
      PAA(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AM4**2)
      PAA(3)= SQRT(ABS(PAA(4)**2-AM4**2))
      PPI   =          PAA(4)**2-AM4**2
      PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAU)
      PHSP=PHSP*(4*PI)*(2*PAA(3)/AMTAU)
* TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AM4**2)
      PN(3)=-PAA(3)
C WE INCLUDE REMAINING PART OF THE JACOBIAN
C --- FLAT CHANNEL
      AM3SQ=(PIM1(4)+PIZ(4)+PIPL(4))**2-(PIM1(3)+PIZ(3)+PIPL(3))**2
     +-(PIM1(2)+PIZ(2)+PIPL(2))**2-(PIM1(1)+PIZ(1)+PIPL(1))**2
      AMS2=(AM4-AMP2)**2
      AMS1=(AMP1+AMP3+AMP4)**2
      FF1=(AMS2-AMS1)
      AMS1=(AMP3+AMP4)**2
      AMS2=(SQRT(AM3SQ)-AMP1)**2
      FF2=AMS2-AMS1
      FF3=(4*PI)*(XLAM(AM2**2,AMP1**2,AM3SQ)/AM3SQ)
      FF4=(4*PI)*(XLAM(AM3SQ,AMP2**2,AM4**2)/AM4**2)
      UU=FF1*FF2*FF3*FF4
C --- FIRST CHANNEL
      AM3SQ=(PIM1(4)+PIZ(4)+PIPL(4))**2-(PIM1(3)+PIZ(3)+PIPL(3))**2
     +-(PIM1(2)+PIZ(2)+PIPL(2))**2-(PIM1(1)+PIZ(1)+PIPL(1))**2
      AMS2=(AM4-AMP2)**2
      AMS1=(AMP1+AMP3+AMP4)**2
      ALP1=ATAN((AMS1-AMRX**2)/AMRX/GAMRX)
      ALP2=ATAN((AMS2-AMRX**2)/AMRX/GAMRX)
      FF1=((AM3SQ-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)
      FF1=FF1*(ALP2-ALP1)
      AMS1=(AMP3+AMP4)**2
      AMS2=(SQRT(AM3SQ)-AMP1)**2
      FF2=AMS2-AMS1
      FF3=(4*PI)*(XLAM(AM2**2,AMP1**2,AM3SQ)/AM3SQ)
      FF4=(4*PI)*(XLAM(AM3SQ,AMP2**2,AM4**2)/AM4**2)
      FF=FF1*FF2*FF3*FF4
C --- SECOND CHANNEL
      AM3SQ=(PIM2(4)+PIZ(4)+PIPL(4))**2-(PIM2(3)+PIZ(3)+PIPL(3))**2
     +-(PIM2(2)+PIZ(2)+PIPL(2))**2-(PIM2(1)+PIZ(1)+PIPL(1))**2
      AMS2=(AM4-AMP1)**2
      AMS1=(AMP2+AMP3+AMP4)**2
      ALP1=ATAN((AMS1-AMRX**2)/AMRX/GAMRX)
      ALP2=ATAN((AMS2-AMRX**2)/AMRX/GAMRX)
      GG1=((AM3SQ-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)
      GG1=GG1*(ALP2-ALP1)
      AMS1=(AMP3+AMP4)**2
      AMS2=(SQRT(AM3SQ)-AMP2)**2
      GG2=AMS2-AMS1
      GG3=(4*PI)*(XLAM(AM2**2,AMP2**2,AM3SQ)/AM3SQ)
      GG4=(4*PI)*(XLAM(AM3SQ,AMP1**2,AM4**2)/AM4**2)
      GG=GG1*GG2*GG3*GG4
C --- JACOBIAN AVERAGED OVER THE TWO
      IF ( ( (FF+GG)*UU+FF*GG ).GT.0.0D0) THEN
        RR=FF*GG*UU/(0.5*PREZ*(FF+GG)*UU+(1.0-PREZ)*FF*GG)
        PHSPAC=PHSPAC*RR
      ELSE
        PHSPAC=0.0
      ENDIF
* MOMENTA OF THE TWO PI-MINUS ARE RANDOMLY SYMMETRISED
      IF (JNPI.EQ.1) THEN
        RR5= RRR(5)
        IF(RR5.LE.0.5) THEN
          DO 20 I=1,4
            X=PIM1(I)
            PIM1(I)=PIM2(I)
   20     PIM2(I)=X
        ENDIF
        PHSPAC=PHSPAC/2.
      ELSE
C MOMENTA OF PI0'S ARE GENERATED UNIFORMLY ONLY IF PREZ=0.0
        RR5= RRR(5)
        IF(RR5.LE.0.5) THEN
          DO 30 I=1,4
            X=PIM1(I)
            PIM1(I)=PIM2(I)
   30     PIM2(I)=X
        ENDIF
        PHSPAC=PHSPAC/6.
      ENDIF
* ALL PIONS BOOSTED FROM  4  REST FRAME TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO NEUTRINO MOMENTUM
      EXE=(PAA(4)+PAA(3))/AM4
      CALL BOSTR3(EXE,PIZ,PIZ)
      CALL BOSTR3(EXE,PIPL,PIPL)
      CALL BOSTR3(EXE,PIM1,PIM1)
      CALL BOSTR3(EXE,PIM2,PIM2)
      CALL BOSTR3(EXE,PR,PR)
C PARTIAL WIDTH CONSISTS OF PHASE SPACE AND AMPLITUDE
C CHECK ON CONSISTENCY WITH DADNPI, THEN, CODE BREAKES UNIFORM PION
C DISTRIBUTION IN HADRONIC SYSTEM
CAM     ASSUME NEUTRINO MASS=0. AND SUM OVER FINAL POLARISATION
C      AMX2=AM4**2
C      BRAK= 2*(AMTAU**2-AMX2) * (AMTAU**2+2.*AMX2)
C      AMPLIT=CCABIB**2*GFERMI**2/2. * BRAK * AMX2*SIGEE(AMX2,1)
      IF     (JNPI.EQ.1) THEN
        CALL DAM4PI(JNPI,PT,PN,PIM1,PIM2,PIZ,PIPL,AMPLIT,HV)
      ELSEIF (JNPI.EQ.2) THEN
        CALL DAM4PI(JNPI,PT,PN,PIM1,PIM2,PIPL,PIZ,AMPLIT,HV)
      ENDIF
      DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
C PHASE SPACE CHECK
C      DGAMT=PHSPAC
      DO 40 K=1,4
        PMULT(K,1)=PIM1(K)
        PMULT(K,2)=PIM2(K)
        PMULT(K,3)=PIPL(K)
        PMULT(K,4)=PIZ (K)
   40 CONTINUE
      END
*CMZ :  1.02/01 12/01/97  16.41.51  by  J. Brunner
*CMZ :  1.01/50 22/05/96  18.06.08  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DPH5PI(DGAMT,HV,PN,PAA,PMULT,JNPI)
C ----------------------------------------------------------------------
* IT SIMULATES 5PI DECAY IN TAU REST FRAME WITH
* Z-AXIS ALONG 5PI MOMENTUM
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST


      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / MDECOMP /IDFFIN(9,NMODE),MULPIK(NMODE)
     +                ,NAMES
      CHARACTER NAMES(NMODE)*31
      REAL  HV(4),PT(4),PN(4),PAA(4),PMULT(4,9)
      REAL*4 PR(4),PI1(4),PI2(4),PI3(4),PI4(4),PI5(4)
      REAL*8 AMP1,AMP2,AMP3,AMP4,AMP5,AMS1,AMS2,AMOM,GAMOM
      REAL*8 AM5SQ,AM4SQ,AM3SQ,AM2SQ,AM5,AM4,AM3
      REAL*4 RRR(10)
      REAL*8 GG1,GG2,GG3,FF1,FF2,FF3,FF4,ALP,ALP1,ALP2
      REAL*8 XM,AM,GAMMA
      COMPLEX BWIGN
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
      DATA FPI /93.3E-3/
C
C
      BWIGN(XM,AM,GAMMA)=XM**2/CMPLX(XM**2-AM**2,GAMMA*AM)
C
      AMOM=.782
      GAMOM=0.0085
C
C 6 BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1./2**29/PI**14
C     PHSPAC=1./2**5/PI**2
C INIT 5PI DECAY MODE (JNPI)
      AMP1=DCDMAS(IDFFIN(1,JNPI))
      AMP2=DCDMAS(IDFFIN(2,JNPI))
      AMP3=DCDMAS(IDFFIN(3,JNPI))
      AMP4=DCDMAS(IDFFIN(4,JNPI))
      AMP5=DCDMAS(IDFFIN(5,JNPI))
C
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C
      CALL RANMAR(RRR,10)
C
C MASSES OF 5, 4, 3 AND 2 PI SYSTEMS
C 3 PI WITH SAMPLING FOR OMEGA RESONANCE
CAM
C MASS OF 5   (12345)
      RR1=RRR(10)
      AMS1=(AMP1+AMP2+AMP3+AMP4+AMP5)**2
      AMS2=(AMTAU-AMNUTA)**2
      AM5SQ=AMS1+   RR1*(AMS2-AMS1)
      AM5 =SQRT(AM5SQ)
      PHSPAC=PHSPAC*(AMS2-AMS1)
C
C MASS OF 4   (2345)
C FLAT PHASE SPACE
      RR1=RRR(9)
      AMS1=(AMP2+AMP3+AMP4+AMP5)**2
      AMS2=(AM5-AMP1)**2
      AM4SQ=AMS1+   RR1*(AMS2-AMS1)
      AM4 =SQRT(AM4SQ)
      GG1=AMS2-AMS1
C
C MASS OF 3   (234)
C PHASE SPACE WITH SAMPLING FOR OMEGA RESONANCE
      RR1=RRR(1)
      AMS1=(AMP2+AMP3+AMP4)**2
      AMS2=(AM4-AMP5)**2
      ALP1=ATAN((AMS1-AMOM**2)/AMOM/GAMOM)
      ALP2=ATAN((AMS2-AMOM**2)/AMOM/GAMOM)
      ALP=ALP1+RR1*(ALP2-ALP1)
      AM3SQ =AMOM**2+AMOM*GAMOM*TAN(ALP)
      AM3 =SQRT(AM3SQ)
C --- THIS PART OF THE JACOBIAN WILL BE RECOVERED LATER ---------------
      GG2=((AM3SQ-AMOM**2)**2+(AMOM*GAMOM)**2)/(AMOM*GAMOM)
      GG2=GG2*(ALP2-ALP1)
C FLAT PHASE SPACE;
C      AM3SQ=AMS1+   RR1*(AMS2-AMS1)
C      AM3 =SQRT(AM3SQ)
C --- THIS PART OF JACOBIAN WILL BE RECOVERED LATER
C      GG2=AMS2-AMS1
C
C MASS OF 2  (34)
      RR2=RRR(2)
      AMS1=(AMP3+AMP4)**2
      AMS2=(AM3-AMP2)**2
C FLAT PHASE SPACE;
      AM2SQ=AMS1+   RR2*(AMS2-AMS1)
      AM2 =SQRT(AM2SQ)
C --- THIS PART OF JACOBIAN WILL BE RECOVERED LATER
      GG3=AMS2-AMS1
C
C (34) RESTFRAME, DEFINE PI3 AND PI4
      ENQ1=(AM2SQ+AMP3**2-AMP4**2)/(2*AM2)
      ENQ2=(AM2SQ-AMP3**2+AMP4**2)/(2*AM2)
      PPI=          ENQ1**2-AMP3**2
      PPPI=SQRT(ABS(ENQ1**2-AMP3**2))
      FF1=(4*PI)*(2*PPPI/AM2)
C PI3   MOMENTUM IN (34) REST FRAME
      CALL SPHERA(PPPI,PI3)
      PI3(4)=ENQ1
C PI4   MOMENTUM IN (34) REST FRAME
      DO 10 I=1,3
   10 PI4(I)=-PI3(I)
      PI4(4)=ENQ2
C
C (234) REST FRAME, DEFINE PI2
C PR   MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AM3)*(AM3**2+AM2**2-AMP2**2)
      PR(3)= SQRT(ABS(PR(4)**2-AM2**2))
      PPI  =          PR(4)**2-AM2**2
C PI2   MOMENTUM
      PI2(1)=0
      PI2(2)=0
      PI2(4)=1./(2*AM3)*(AM3**2-AM2**2+AMP2**2)
      PI2(3)=-PR(3)
C --- THIS PART OF JACOBIAN WILL BE RECOVERED LATER
      FF2=(4*PI)*(2*PR(3)/AM3)
C OLD PIONS BOOSTED FROM 2 REST FRAME TO 3 REST FRAME
      EXE=(PR(4)+PR(3))/AM2
      CALL BOSTR3(EXE,PI3,PI3)
      CALL BOSTR3(EXE,PI4,PI4)
      RR3=RRR(3)
      RR4=RRR(4)
      THET =ACOS(-1.+2*RR3)
      PHI = 2*PI*RR4
      CALL ROTPOL(THET,PHI,PI2)
      CALL ROTPOL(THET,PHI,PI3)
      CALL ROTPOL(THET,PHI,PI4)
C
C (2345)  REST FRAME, DEFINE PI5
C PR   MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AM4)*(AM4**2+AM3**2-AMP5**2)
      PR(3)= SQRT(ABS(PR(4)**2-AM3**2))
      PPI  =          PR(4)**2-AM3**2
C PI5  MOMENTUM
      PI5(1)=0
      PI5(2)=0
      PI5(4)=1./(2*AM4)*(AM4**2-AM3**2+AMP5**2)
      PI5(3)=-PR(3)
C --- THIS PART OF JACOBIAN WILL BE RECOVERED LATER
      FF3=(4*PI)*(2*PR(3)/AM4)
C OLD PIONS BOOSTED FROM 3 REST FRAME TO 4 REST FRAME
      EXE=(PR(4)+PR(3))/AM3
      CALL BOSTR3(EXE,PI2,PI2)
      CALL BOSTR3(EXE,PI3,PI3)
      CALL BOSTR3(EXE,PI4,PI4)
      RR3=RRR(5)
      RR4=RRR(6)
      THET =ACOS(-1.+2*RR3)
      PHI = 2*PI*RR4
      CALL ROTPOL(THET,PHI,PI2)
      CALL ROTPOL(THET,PHI,PI3)
      CALL ROTPOL(THET,PHI,PI4)
      CALL ROTPOL(THET,PHI,PI5)
C
C (12345)  REST FRAME, DEFINE PI1
C PR   MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AM5)*(AM5**2+AM4**2-AMP1**2)
      PR(3)= SQRT(ABS(PR(4)**2-AM4**2))
      PPI  =          PR(4)**2-AM4**2
C PI1  MOMENTUM
      PI1(1)=0
      PI1(2)=0
      PI1(4)=1./(2*AM5)*(AM5**2-AM4**2+AMP1**2)
      PI1(3)=-PR(3)
C --- THIS PART OF JACOBIAN WILL BE RECOVERED LATER
      FF4=(4*PI)*(2*PR(3)/AM5)
C OLD PIONS BOOSTED FROM 4 REST FRAME TO 5 REST FRAME
      EXE=(PR(4)+PR(3))/AM4
      CALL BOSTR3(EXE,PI2,PI2)
      CALL BOSTR3(EXE,PI3,PI3)
      CALL BOSTR3(EXE,PI4,PI4)
      CALL BOSTR3(EXE,PI5,PI5)
      RR3=RRR(7)
      RR4=RRR(8)
      THET =ACOS(-1.+2*RR3)
      PHI = 2*PI*RR4
      CALL ROTPOL(THET,PHI,PI1)
      CALL ROTPOL(THET,PHI,PI2)
      CALL ROTPOL(THET,PHI,PI3)
      CALL ROTPOL(THET,PHI,PI4)
      CALL ROTPOL(THET,PHI,PI5)
C
* NOW TO THE TAU REST FRAME, DEFINE PAA AND NEUTRINO MOMENTA
* PAA  MOMENTUM
      PAA(1)=0
      PAA(2)=0
C     PAA(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AM5**2)
C     PAA(3)= SQRT(ABS(PAA(4)**2-AM5**2))
C     PPI   =          PAA(4)**2-AM5**2
      PAA(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AM5SQ)
      PAA(3)= SQRT(ABS(PAA(4)**2-AM5SQ))
      PPI   =          PAA(4)**2-AM5SQ
      PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAU)
* TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AM5**2)
      PN(3)=-PAA(3)
C
      PHSPAC=PHSPAC * GG1*GG2*GG3*FF1*FF2*FF3*FF4
C
C ALL PIONS BOOSTED FROM  5  REST FRAME TO TAU REST FRAME
C Z-AXIS ANTIPARALLEL TO NEUTRINO MOMENTUM
      EXE=(PAA(4)+PAA(3))/AM5
      CALL BOSTR3(EXE,PI1,PI1)
      CALL BOSTR3(EXE,PI2,PI2)
      CALL BOSTR3(EXE,PI3,PI3)
      CALL BOSTR3(EXE,PI4,PI4)
      CALL BOSTR3(EXE,PI5,PI5)
C
C PARTIAL WIDTH CONSISTS OF PHASE SPACE AND AMPLITUDE
C AMPLITUDE  (CF YS.TSAI PHYS.REV.D4,2821(1971)
C    OR F.GILMAN SH.RHIE PHYS.REV.D31,1066(1985)
C
      PXQ=AMTAU*PAA(4)
      PXN=AMTAU*PN(4)
      QXN=PAA(4)*PN(4)-PAA(1)*PN(1)-PAA(2)*PN(2)-PAA(3)*PN(3)
      BRAK=2*(GV**2+GA**2)*(2*PXQ*QXN+AM5SQ*PXN)
     +    -6*(GV**2-GA**2)*AMTAU*AMNUTA*AM5SQ
      FOMPP = CABS(BWIGN(AM3,AMOM,GAMOM))**2
C NORMALISATION FACTOR (TO SOME NUMERICAL UNDIMENSIONED FACTOR;
C CF R.FISCHER ET AL ZPHYS C3, 313 (1980))
      FNORM = 1/FPI**6
C     AMPLIT=CCABIB**2*GFERMI**2/2. * BRAK * AM5SQ*SIGEE(AM5SQ,JNPI)
      AMPLIT=CCABIB**2*GFERMI**2/2. * BRAK
      AMPLIT = AMPLIT * FOMPP * FNORM
C PHASE SPACE TEST
C     AMPLIT = AMPLIT * FNORM
      DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
C IGNORE SPIN TERMS
      DO 20 I=1,3
   20 HV(I)=0.
C
      DO 30 K=1,4
        PMULT(K,1)=PI1(K)
        PMULT(K,2)=PI2(K)
        PMULT(K,3)=PI3(K)
        PMULT(K,4)=PI4(K)
        PMULT(K,5)=PI5(K)
   30 CONTINUE
      RETURN
C MISSING: TRANSPOSITION OF IDENTICAL PARTICLES, STARTISTICAL FACTORS
C FOR IDENTICAL MATRICES, POLARIMETRIC VECTOR. MATRIX ELEMENT RATHER NAIVE.
C FLAT PHASE SPACE IN PION SYSTEM + WITH BREIT WIGNER FOR OMEGA
C ANYWAY IT IS BETTER THAN NOTHING, AND CODE IS IMPROVABLE.
      END
*CMZ :  1.01/50 22/05/96  18.06.08  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DPHNPI(DGAMT,HVX,PNX,PRX,PPIX,JNPI)
C ----------------------------------------------------------------------
C IT SIMULATES MULTIPI DECAY IN TAU REST FRAME WITH
C Z-AXIS OPPOSITE TO NEUTRINO MOMENTUM
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / MDECOMP /IDFFIN(9,NMODE),MULPIK(NMODE)
     +                ,NAMES
      CHARACTER NAMES(NMODE)*31
      REAL*8 WETMAX(20)
C
      REAL*8  PN(4),PR(4),PPI(4,9),HV(4)
      REAL*4  PNX(4),PRX(4),PPIX(4,9),HVX(4)
      REAL*8  PV(5,9),PT(4),UE(3),BE(3)
      REAL*8  PAWT,AMX,AMS1,AMS2,PA,PHS,PHSMAX,PMIN,PMAX
      REAL*8  GAM,BEP,PHI,A,B,C
      REAL*8  AMPIK
      REAL*4 RRR(9),RRX(2),RTEMP(1)
C
      DATA PI /3.141592653589793238462643/
      DATA WETMAX /20*1D-15/
C
CC--      PAWT(A,B,C)=SQRT((A**2-(B+C)**2)*(A**2-(B-C)**2))/(2.*A)
C
      PAWT(A,B,C)=
     +  SQRT(MAX(0.D0,(A**2-(B+C)**2)*(A**2-(B-C)**2)))/(2.D0*A)
C
      AMPIK(I,J)=DCDMAS(IDFFIN(I,J))
C
C
      IF ((JNPI.LE.0).OR.JNPI.GT.20) THEN
        WRITE(6,*) 'JNPI OUTSIDE RANGE DEFINED BY WETMAX; JNPI=',JNPI
        STOP
      ENDIF

C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C
   10 CONTINUE
C MASS OF VIRTUAL W
      ND=MULPIK(JNPI)
      PS=0.
      PHSPAC = 1./2.**5 /PI**2
      DO 20 I=1,ND
   20 PS  =PS+AMPIK(I,JNPI)
      RTEMP(1)=RR1
      CALL RANMAR(RTEMP,1)
      RR1=RTEMP(1)
      AMS1=PS**2
      AMS2=(AMTAU-AMNUTA)**2
C
C
      AMX2=AMS1+   RR1*(AMS2-AMS1)
      AMX =SQRT(AMX2)
      AMW =AMX
      PHSPAC=PHSPAC * (AMS2-AMS1)
C
C TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AMX2)
      PN(3)=-SQRT((PN(4)-AMNUTA)*(PN(4)+AMNUTA))
C W MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AMX2)
      PR(3)=-PN(3)
      PHSPAC=PHSPAC * (4.*PI) * (2.*PR(3)/AMTAU)
C
C AMPLITUDE  (CF YS.TSAI PHYS.REV.D4,2821(1971)
C    OR F.GILMAN SH.RHIE PHYS.REV.D31,1066(1985)
C
      PXQ=AMTAU*PR(4)
      PXN=AMTAU*PN(4)
      QXN=PR(4)*PN(4)-PR(1)*PN(1)-PR(2)*PN(2)-PR(3)*PN(3)
C HERE WAS AN ERROR. 20.10.91 (ZW)
C       BRAK=2*(GV**2+GA**2)*(2*PXQ*PXN+AMX2*QXN)
      BRAK=2*(GV**2+GA**2)*(2*PXQ*QXN+AMX2*PXN) -6*(GV**2-GA**2)*AMTAU*
     +AMNUTA*AMX2
CAM     ASSUME NEUTRINO MASS=0. AND SUM OVER FINAL POLARISATION
C     BRAK= 2*(AMTAU**2-AMX2) * (AMTAU**2+2.*AMX2)
      AMPLIT=CCABIB**2*GFERMI**2/2. * BRAK * AMX2*SIGEE(AMX2,JNPI)
      DGAMT=1./(2.*AMTAU)*AMPLIT*PHSPAC
C
C   ISOTROPIC W DECAY IN W REST FRAME
      PHSMAX = 1.
      DO 30  I=1,4
   30 PV(I,1)=PR(I)
      PV(5,1)=AMW
      PV(5,ND)=AMPIK(ND,JNPI)
C    COMPUTE MAX. PHASE SPACE FACTOR
      PMAX=AMW-PS+AMPIK(ND,JNPI)
      PMIN=.0
      DO 40  IL=ND-1,1,-1
        PMAX=PMAX+AMPIK(IL,JNPI)
        PMIN=PMIN+AMPIK(IL+1,JNPI)
   40 PHSMAX=PHSMAX*PAWT(PMAX,PMIN,AMPIK(IL,JNPI))/PMAX

C --- 2.02.94 ZW  9 LINES
      AMX=AMW
      DO 60  IL=1,ND-2
        AMS1=.0
        DO 50  JL=IL+1,ND
   50   AMS1=AMS1+AMPIK(JL,JNPI)
        AMS1=AMS1**2
        AMX =(AMX-AMPIK(IL,JNPI))
        AMS2=(AMX)**2
        PHSMAX=PHSMAX * (AMS2-AMS1)
   60 CONTINUE
      NCONT=0
   70 CONTINUE
      NCONT=NCONT+1
CAM  GENERATE ND-2 EFFECTIVE MASSES
      PHS=1.D0
      PHSPAC = 1./2.**(6*ND-7) /PI**(3*ND-4)
      AMX=AMW
      CALL RANMAR(RRR,ND-2)
      DO 90  IL=1,ND-2
        AMS1=.0D0
        DO 80  JL=IL+1,ND
   80   AMS1=AMS1+AMPIK(JL,JNPI)
        AMS1=AMS1**2
        AMS2=(AMX-AMPIK(IL,JNPI))**2
        RR1=RRR(IL)
        AMX2=AMS1+ RR1*(AMS2-AMS1)
        AMX=SQRT(AMX2)
        PV(5,IL+1)=AMX
        PHSPAC=PHSPAC * (AMS2-AMS1)
C ---  2.02.94 ZW 1 LINE
        PHS=PHS* (AMS2-AMS1)
        PA=PAWT(PV(5,IL),PV(5,IL+1),AMPIK(IL,JNPI))
        PHS =PHS *PA/PV(5,IL)
   90 CONTINUE
      PA=PAWT(PV(5,ND-1),AMPIK(ND-1,JNPI),AMPIK(ND,JNPI))
      PHS   =PHS    *PA/PV(5,ND-1)
      RTEMP(1)=RN
      CALL RANMAR(RTEMP,1)
      RN=RTEMP(1)
      WETMAX(JNPI)=1.2D0*MAX(WETMAX(JNPI)/1.2D0,PHS/PHSMAX)
      IF (NCONT.EQ.500 000) THEN
        XNPI=0.0
        DO KK=1,ND
          XNPI=XNPI+AMPIK(KK,JNPI)
        ENDDO
        WRITE(6,*) 'ROUNDING INSTABILITY IN DPHNPI ?'
        WRITE(6,*) 'AMW=',AMW,'XNPI=',XNPI
        WRITE(6,*) 'IF =AMW= IS NEARLY EQUAL =XNPI= THAT IS IT'
        WRITE(6,*) 'PHS=',PHS,'PHSMAX=',PHSMAX
        GOTO 10
      ENDIF
      IF(RN*PHSMAX*WETMAX(JNPI).GT.PHS) GO TO 70
C...PERFORM SUCCESSIVE TWO-PARTICLE DECAYS IN RESPECTIVE CM FRAME
  100 DO 120 IL=1,ND-1
        PA=PAWT(PV(5,IL),PV(5,IL+1),AMPIK(IL,JNPI))
        CALL RANMAR(RRX,2)
        UE(3)=2.*RRX(1)-1.
        PHI=2.*PI*RRX(2)
        UE(1)=SQRT(1.D0-UE(3)**2)*COS(PHI)
        UE(2)=SQRT(1.D0-UE(3)**2)*SIN(PHI)
        DO 110 J=1,3
          PPI(J,IL)=PA*UE(J)
  110   PV(J,IL+1)=-PA*UE(J)
        PPI(4,IL)=SQRT(PA**2+AMPIK(IL,JNPI)**2)
        PV(4,IL+1)=SQRT(PA**2+PV(5,IL+1)**2)
        PHSPAC=PHSPAC *(4.*PI)*(2.*PA/PV(5,IL))
  120 CONTINUE
C...LORENTZ TRANSFORM DECAY PRODUCTS TO TAU FRAME
      DO 130 J=1,4
  130 PPI(J,ND)=PV(J,ND)
      DO 160 IL=ND-1,1,-1
        DO 140 J=1,3
  140   BE(J)=PV(J,IL)/PV(4,IL)
        GAM=PV(4,IL)/PV(5,IL)
        DO 160 I=IL,ND
          BEP=BE(1)*PPI(1,I)+BE(2)*PPI(2,I)+BE(3)*PPI(3,I)
          DO 150 J=1,3
  150     PPI(J,I)=PPI(J,I)+GAM*(GAM*BEP/(1.D0+GAM)+PPI(4,I))*BE(J)
          PPI(4,I)=GAM*(PPI(4,I)+BEP)
  160 CONTINUE
C
      HV(4)=1.
      HV(3)=0.
      HV(2)=0.
      HV(1)=0.
      DO K=1,4
        PNX(K)=PN(K)
        PRX(K)=PR(K)
        HVX(K)=HV(K)
        DO L=1,ND
          PPIX(K,L)=PPI(K,L)
        ENDDO
      ENDDO
      RETURN
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DPHSAA(DGAMT,HV,PN,PAA,PIM1,PIM2,PIPL,JAA)
C ----------------------------------------------------------------------
* IT SIMULATES A1  DECAY IN TAU REST FRAME WITH
* Z-AXIS ALONG A1  MOMENTUM
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      REAL  HV(4),PN(4),PAA(4),PIM1(4),PIM2(4),PIPL(4)


      REAL*4 RRR(1)
C MATRIX ELEMENT NUMBER:
      MNUM=0
C TYPE OF THE GENERATION:
      KEYT=1
      CALL RANMAR(RRR,1)
      RMOD=RRR(1)
      IF (RMOD.LT.BRA1) THEN
        JAA=1
        AMP1=AMPI
        AMP2=AMPI
        AMP3=AMPI
      ELSE
        JAA=2
        AMP1=AMPIZ
        AMP2=AMPIZ
        AMP3=AMPI
      ENDIF

      CALL
     +   DPHTRE(DGAMT,HV,PN,PAA,PIM1,AMP1,PIM2,AMP2,PIPL,AMP3,KEYT,MNUM)
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DPHSEL(DGAMX,HVX,XNX,PAAX,QPX,XAX,PHX)
C XNX,XNA WAS FLIPPED IN PARAMETERS OF DPHSEL AND DPHSMU
C *********************************************************************
C *   ELECTRON DECAY MODE                                             *
C *********************************************************************
      REAL*4         PHX(4)
      REAL*4  HVX(4),PAAX(4),XAX(4),QPX(4),XNX(4)
      REAL*8  HV(4),PH(4),PAA(4),XA(4),QP(4),XN(4)
      REAL*8  DGAMT
      IELMU=1
      CALL DRCMU(DGAMT,HV,PH,PAA,XA,QP,XN,IELMU)
      DO 10 K=1,4
        HVX(K)=HV(K)
        PHX(K)=PH(K)
        PAAX(K)=PAA(K)
        XAX(K)=XA(K)
        QPX(K)=QP(K)
        XNX(K)=XN(K)
   10 CONTINUE
      DGAMX=DGAMT
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DPHSKS(DGAMT,HV,PN,PKS,PKK,PPI,JKST)
C ----------------------------------------------------------------------
C IT SIMULATES KAON* DECAY IN TAU REST FRAME WITH
C Z-AXIS ALONG KAON* MOMENTUM
C     JKST=10 FOR K* --->K0 + PI+-
C     JKST=20 FOR K* --->K+- + PI0
C ----------------------------------------------------------------------
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      REAL  HV(4),PT(4),PN(4),PKS(4),PKK(4),PPI(4),QQ(4)
      COMPLEX BWIGS
      DATA PI /3.141592653589793238462643/
C
      DATA ICONT /0/
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
      PHSPAC=1./2**11/PI**5
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
      CALL RANMAR(RR1,1)
C HERE BEGIN THE K0,PI+_ DECAY
      IF(JKST.EQ.10)THEN
C     ==================
C MASS OF (REAL/VIRTUAL) K*
        AMS1=(AMPI+AMKZ)**2
        AMS2=(AMTAU-AMNUTA)**2
C FLAT PHASE SPACE
C       AMX2=AMS1+   RR1*(AMS2-AMS1)
C       AMX=SQRT(AMX2)
C       PHSPAC=PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH SAMPLING FOR K* RESONANCE
        ALP1=ATAN((AMS1-AMKST**2)/AMKST/GAMKST)
        ALP2=ATAN((AMS2-AMKST**2)/AMKST/GAMKST)
        ALP=ALP1+RR1*(ALP2-ALP1)
        AMX2=AMKST**2+AMKST*GAMKST*TAN(ALP)
        AMX=SQRT(AMX2)
        PHSPAC=PHSPAC*((AMX2-AMKST**2)**2+(AMKST*GAMKST)**2)
     &                /(AMKST*GAMKST)
        PHSPAC=PHSPAC*(ALP2-ALP1)
C
C TAU-NEUTRINO MOMENTUM
        PN(1)=0
        PN(2)=0
        PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AMX**2)
        PN(3)=-SQRT((PN(4)-AMNUTA)*(PN(4)+AMNUTA))
C
C K* MOMENTUM
        PKS(1)=0
        PKS(2)=0
        PKS(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AMX**2)
        PKS(3)=-PN(3)
        PHSPAC=PHSPAC*(4*PI)*(2*PKS(3)/AMTAU)
C
CAM
        ENPI=( AMX**2+AMPI**2-AMKZ**2 ) / ( 2*AMX )
        PPPI=SQRT((ENPI-AMPI)*(ENPI+AMPI))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AMX)
C CHARGED PI MOMENTUM IN KAON* REST FRAME
        CALL SPHERA(PPPI,PPI)
        PPI(4)=ENPI
C NEUTRAL KAON MOMENTUM IN K* REST FRAME
        DO 10 I=1,3
   10   PKK(I)=-PPI(I)
        PKK(4)=( AMX**2+AMKZ**2-AMPI**2 ) / ( 2*AMX )
        EXE=(PKS(4)+PKS(3))/AMX
C PION AND K  BOOSTED FROM K* REST FRAME TO TAU REST FRAME
        CALL BOSTR3(EXE,PPI,PPI)
        CALL BOSTR3(EXE,PKK,PKK)
        DO 20 I=1,4
   20   QQ(I)=PPI(I)-PKK(I)
C QQ TRANSVERSE TO PKS
        PKSD =PKS(4)*PKS(4)-PKS(3)*PKS(3)-PKS(2)*PKS(2)-PKS(1)*PKS(1)
        QQPKS=PKS(4)* QQ(4)-PKS(3)* QQ(3)-PKS(2)* QQ(2)-PKS(1)* QQ(1)
        DO 30 I=1,4
   30   QQ(I)=QQ(I)-PKS(I)*QQPKS/PKSD
C AMPLITUDE
        PRODPQ=PT(4)*QQ(4)
        PRODNQ=PN(4)*QQ(4)-PN(1)*QQ(1)-PN(2)*QQ(2)-PN(3)*QQ(3)
        PRODPN=PT(4)*PN(4)
        QQ2= QQ(4)**2-QQ(1)**2-QQ(2)**2-QQ(3)**2
        BRAK=(GV**2+GA**2)*(2*PRODPQ*PRODNQ-PRODPN*QQ2)
     &      +(GV**2-GA**2)*AMTAU*AMNUTA*QQ2
C A SIMPLE BREIT-WIGNER IS CHOSEN FOR K* RESONANCE
        FKS=CABS(BWIGS(AMX2,AMKST,GAMKST))**2
        AMPLIT=(GFERMI*SCABIB)**2*BRAK*2*FKS
        DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
        DO 40 I=1,3
   40   HV(I)=2*GV*GA*AMTAU*(2*PRODNQ*QQ(I)-QQ2*PN(I))/BRAK
C
C HERE BEGIN THE K+-,PI0 DECAY
      ELSEIF(JKST.EQ.20)THEN
C     ======================
C MASS OF (REAL/VIRTUAL) K*
        AMS1=(AMPIZ+AMK)**2
        AMS2=(AMTAU-AMNUTA)**2
C FLAT PHASE SPACE
C       AMX2=AMS1+   RR1*(AMS2-AMS1)
C       AMX=SQRT(AMX2)
C       PHSPAC=PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH SAMPLING FOR K* RESONANCE
        ALP1=ATAN((AMS1-AMKST**2)/AMKST/GAMKST)
        ALP2=ATAN((AMS2-AMKST**2)/AMKST/GAMKST)
        ALP=ALP1+RR1*(ALP2-ALP1)
        AMX2=AMKST**2+AMKST*GAMKST*TAN(ALP)
        AMX=SQRT(AMX2)
        PHSPAC=PHSPAC*((AMX2-AMKST**2)**2+(AMKST*GAMKST)**2)
     &                /(AMKST*GAMKST)
        PHSPAC=PHSPAC*(ALP2-ALP1)
C
C TAU-NEUTRINO MOMENTUM
        PN(1)=0
        PN(2)=0
        PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AMX**2)
        PN(3)=-SQRT((PN(4)-AMNUTA)*(PN(4)+AMNUTA))
C KAON* MOMENTUM
        PKS(1)=0
        PKS(2)=0
        PKS(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AMX**2)
        PKS(3)=-PN(3)
        PHSPAC=PHSPAC*(4*PI)*(2*PKS(3)/AMTAU)
C
CAM
        ENPI=( AMX**2+AMPIZ**2-AMK**2 ) / ( 2*AMX )
        PPPI=SQRT((ENPI-AMPIZ)*(ENPI+AMPIZ))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AMX)
C NEUTRAL PI MOMENTUM IN K* REST FRAME
        CALL SPHERA(PPPI,PPI)
        PPI(4)=ENPI
C CHARGED KAON MOMENTUM IN K* REST FRAME
        DO 50 I=1,3
   50   PKK(I)=-PPI(I)
        PKK(4)=( AMX**2+AMK**2-AMPIZ**2 ) / ( 2*AMX )
        EXE=(PKS(4)+PKS(3))/AMX
C PION AND K  BOOSTED FROM K* REST FRAME TO TAU REST FRAME
        CALL BOSTR3(EXE,PPI,PPI)
        CALL BOSTR3(EXE,PKK,PKK)
        DO 60 I=1,4
   60   QQ(I)=PKK(I)-PPI(I)
C QQ TRANSVERSE TO PKS
        PKSD =PKS(4)*PKS(4)-PKS(3)*PKS(3)-PKS(2)*PKS(2)-PKS(1)*PKS(1)
        QQPKS=PKS(4)* QQ(4)-PKS(3)* QQ(3)-PKS(2)* QQ(2)-PKS(1)* QQ(1)
        DO 70 I=1,4
   70   QQ(I)=QQ(I)-PKS(I)*QQPKS/PKSD
C AMPLITUDE
        PRODPQ=PT(4)*QQ(4)
        PRODNQ=PN(4)*QQ(4)-PN(1)*QQ(1)-PN(2)*QQ(2)-PN(3)*QQ(3)
        PRODPN=PT(4)*PN(4)
        QQ2= QQ(4)**2-QQ(1)**2-QQ(2)**2-QQ(3)**2
        BRAK=(GV**2+GA**2)*(2*PRODPQ*PRODNQ-PRODPN*QQ2)
     &      +(GV**2-GA**2)*AMTAU*AMNUTA*QQ2
C A SIMPLE BREIT-WIGNER IS CHOSEN FOR THE K* RESONANCE
        FKS=CABS(BWIGS(AMX2,AMKST,GAMKST))**2
        AMPLIT=(GFERMI*SCABIB)**2*BRAK*2*FKS
        DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
        DO 80 I=1,3
   80   HV(I)=2*GV*GA*AMTAU*(2*PRODNQ*QQ(I)-QQ2*PN(I))/BRAK
      ENDIF
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DPHSMU(DGAMX,HVX,XNX,PAAX,QPX,XAX,PHX)
C XNX,XNA WAS FLIPPED IN PARAMETERS OF DPHSEL AND DPHSMU
C *********************************************************************
C *   MUON     DECAY MODE                                             *
C *********************************************************************
      REAL*4         PHX(4)
      REAL*4  HVX(4),PAAX(4),XAX(4),QPX(4),XNX(4)
      REAL*8  HV(4),PH(4),PAA(4),XA(4),QP(4),XN(4)
      REAL*8  DGAMT
      IELMU=2
      CALL DRCMU(DGAMT,HV,PH,PAA,XA,QP,XN,IELMU)
      DO 10 K=1,4
        HVX(K)=HV(K)
        PHX(K)=PH(K)
        PAAX(K)=PAA(K)
        XAX(K)=XA(K)
        QPX(K)=QP(K)
        XNX(K)=XN(K)
   10 CONTINUE
      DGAMX=DGAMT
      END
*CMZ :  1.01/50 22/05/96  18.06.08  by  Piero Zucchelli
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DPHSPK(DGAMT,HV,PN,PAA,PNPI,JAA)
C ----------------------------------------------------------------------
* IT SIMULATES THREE PI (K) DECAY IN THE TAU REST FRAME
* Z-AXIS ALONG HADRONIC SYSTEM
C ----------------------------------------------------------------------
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / MDECOMP /IDFFIN(9,NMODE),MULPIK(NMODE)
     +                ,NAMES
      CHARACTER NAMES(NMODE)*31

      REAL  HV(4),PN(4),PAA(4),PIM1(4),PIM2(4),PIPL(4),PNPI(4,9)
C MATRIX ELEMENT NUMBER:
      MNUM=JAA
C TYPE OF THE GENERATION:
      KEYT=4
      IF(JAA.EQ.7) KEYT=3
C --- MASSES OF THE DECAY PRODUCTS
      AMP1=DCDMAS(IDFFIN(1,JAA+NM4+NM5+NM6))
      AMP2=DCDMAS(IDFFIN(2,JAA+NM4+NM5+NM6))
      AMP3=DCDMAS(IDFFIN(3,JAA+NM4+NM5+NM6))
      CALL
     +   DPHTRE(DGAMT,HV,PN,PAA,PIM1,AMP1,PIM2,AMP2,PIPL,AMP3,KEYT,MNUM)
      DO I=1,4
        PNPI(I,1)=PIM1(I)
        PNPI(I,2)=PIM2(I)
        PNPI(I,3)=PIPL(I)
      ENDDO
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DPHSRK(DGAMT,HV,PN,PR,PMULT,INUM)
C ----------------------------------------------------------------------
C IT SIMULATES RHO DECAY IN TAU REST FRAME WITH
C Z-AXIS ALONG RHO MOMENTUM
C RHO DECAYS TO K KBAR
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL  HV(4),PT(4),PN(4),PR(4),PKC(4),PKZ(4),QQ(4),PMULT(4,9)
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
C
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
      PHSPAC=1./2**11/PI**5
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C MASS OF (REAL/VIRTUAL) RHO
      AMS1=(AMK+AMKZ)**2
      AMS2=(AMTAU-AMNUTA)**2
C FLAT PHASE SPACE
      CALL RANMAR(RR1,1)
      AMX2=AMS1+   RR1*(AMS2-AMS1)
      AMX=SQRT(AMX2)
      PHSPAC=PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH SAMPLING FOR RHO RESONANCE
C     ALP1=ATAN((AMS1-AMRO**2)/AMRO/GAMRO)
C     ALP2=ATAN((AMS2-AMRO**2)/AMRO/GAMRO)
CAM
   10 CONTINUE
C     CALL RANMAR(RR1,1)
C     ALP=ALP1+RR1*(ALP2-ALP1)
C     AMX2=AMRO**2+AMRO*GAMRO*TAN(ALP)
C     AMX=SQRT(AMX2)
C     IF(AMX.LT.(AMK+AMKZ)) GO TO 100
CAM
C     PHSPAC=PHSPAC*((AMX2-AMRO**2)**2+(AMRO*GAMRO)**2)/(AMRO*GAMRO)
C     PHSPAC=PHSPAC*(ALP2-ALP1)
C
C TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AMX**2)
      PN(3)=-SQRT((PN(4)-AMNUTA)*(PN(4)+AMNUTA))
C RHO MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AMX**2)
      PR(3)=-PN(3)
      PHSPAC=PHSPAC*(4*PI)*(2*PR(3)/AMTAU)
C
CAM
      ENQ1=(AMX2+AMK**2-AMKZ**2)/(2.*AMX)
      ENQ2=(AMX2-AMK**2+AMKZ**2)/(2.*AMX)
      PPPI=SQRT((ENQ1-AMK)*(ENQ1+AMK))
      PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AMX)
C CHARGED PI MOMENTUM IN RHO REST FRAME
      CALL SPHERA(PPPI,PKC)
      PKC(4)=ENQ1
C NEUTRAL PI MOMENTUM IN RHO REST FRAME
      DO 20 I=1,3
   20 PKZ(I)=-PKC(I)
      PKZ(4)=ENQ2
      EXE=(PR(4)+PR(3))/AMX
C PIONS BOOSTED FROM RHO REST FRAME TO TAU REST FRAME
      CALL BOSTR3(EXE,PKC,PKC)
      CALL BOSTR3(EXE,PKZ,PKZ)
      DO 30 I=1,4
   30 QQ(I)=PKC(I)-PKZ(I)
C QQ TRANSVERSE TO PR
      PKSD =PR(4)*PR(4)-PR(3)*PR(3)-PR(2)*PR(2)-PR(1)*PR(1)
      QQPKS=PR(4)* QQ(4)-PR(3)* QQ(3)-PR(2)* QQ(2)-PR(1)* QQ(1)
      DO 40 I=1,4
   40 QQ(I)=QQ(I)-PR(I)*QQPKS/PKSD
C AMPLITUDE
      PRODPQ=PT(4)*QQ(4)
      PRODNQ=PN(4)*QQ(4)-PN(1)*QQ(1)-PN(2)*QQ(2)-PN(3)*QQ(3)
      PRODPN=PT(4)*PN(4)
      QQ2= QQ(4)**2-QQ(1)**2-QQ(2)**2-QQ(3)**2
      BRAK=(GV**2+GA**2)*(2*PRODPQ*PRODNQ-PRODPN*QQ2)
     +    +(GV**2-GA**2)*AMTAU*AMNUTA*QQ2
      AMPLIT=(GFERMI*CCABIB)**2*BRAK*2*FPIRK(AMX)
      DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
      DO 50 I=1,3
   50 HV(I)=2*GV*GA*AMTAU*(2*PRODNQ*QQ(I)-QQ2*PN(I))/BRAK
      DO 60 K=1,4
        PMULT(K,1)=PKC(K)
        PMULT(K,2)=PKZ(K)
   60 CONTINUE
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DPHSRO(DGAMT,HV,PN,PR,PIC,PIZ)
C ----------------------------------------------------------------------
C IT SIMULATES RHO DECAY IN TAU REST FRAME WITH
C Z-AXIS ALONG RHO MOMENTUM
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL  HV(4),PT(4),PN(4),PR(4),PIC(4),PIZ(4),QQ(4)
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
C
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
      PHSPAC=1./2**11/PI**5
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C MASS OF (REAL/VIRTUAL) RHO
      AMS1=(AMPI+AMPIZ)**2
      AMS2=(AMTAU-AMNUTA)**2
C FLAT PHASE SPACE
C     AMX2=AMS1+   RR1*(AMS2-AMS1)
C     AMX=SQRT(AMX2)
C     PHSPAC=PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH SAMPLING FOR RHO RESONANCE
      ALP1=ATAN((AMS1-AMRO**2)/AMRO/GAMRO)
      ALP2=ATAN((AMS2-AMRO**2)/AMRO/GAMRO)
CAM
   10 CONTINUE
      CALL RANMAR(RR1,1)
      ALP=ALP1+RR1*(ALP2-ALP1)
      AMX2=AMRO**2+AMRO*GAMRO*TAN(ALP)
      AMX=SQRT(AMX2)
      IF(AMX.LT.2.*AMPI) GO TO 10
CAM
      PHSPAC=PHSPAC*((AMX2-AMRO**2)**2+(AMRO*GAMRO)**2)/(AMRO*GAMRO)
      PHSPAC=PHSPAC*(ALP2-ALP1)
C
C TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AMX**2)
      PN(3)=-SQRT((PN(4)-AMNUTA)*(PN(4)+AMNUTA))
C RHO MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AMX**2)
      PR(3)=-PN(3)
      PHSPAC=PHSPAC*(4*PI)*(2*PR(3)/AMTAU)
C
CAM
      ENQ1=(AMX2+AMPI**2-AMPIZ**2)/(2.*AMX)
      ENQ2=(AMX2-AMPI**2+AMPIZ**2)/(2.*AMX)
      PPPI=SQRT((ENQ1-AMPI)*(ENQ1+AMPI))
      PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AMX)
C CHARGED PI MOMENTUM IN RHO REST FRAME
      CALL SPHERA(PPPI,PIC)
      PIC(4)=ENQ1
C NEUTRAL PI MOMENTUM IN RHO REST FRAME
      DO 20 I=1,3
   20 PIZ(I)=-PIC(I)
      PIZ(4)=ENQ2
      EXE=(PR(4)+PR(3))/AMX
C PIONS BOOSTED FROM RHO REST FRAME TO TAU REST FRAME
      CALL BOSTR3(EXE,PIC,PIC)
      CALL BOSTR3(EXE,PIZ,PIZ)
      DO 30 I=1,4
   30 QQ(I)=PIC(I)-PIZ(I)
C AMPLITUDE
      PRODPQ=PT(4)*QQ(4)
      PRODNQ=PN(4)*QQ(4)-PN(1)*QQ(1)-PN(2)*QQ(2)-PN(3)*QQ(3)
      PRODPN=PT(4)*PN(4)
      QQ2= QQ(4)**2-QQ(1)**2-QQ(2)**2-QQ(3)**2
      BRAK=(GV**2+GA**2)*(2*PRODPQ*PRODNQ-PRODPN*QQ2)
     &    +(GV**2-GA**2)*AMTAU*AMNUTA*QQ2
      AMPLIT=(GFERMI*CCABIB)**2*BRAK*2*FPIRHO(AMX)
      DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
      DO 40 I=1,3
   40 HV(I)=2*GV*GA*AMTAU*(2*PRODNQ*QQ(I)-QQ2*PN(I))/BRAK
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE
     +   DPHTRE(DGAMT,HV,PN,PAA,PIM1,AMPA,PIM2,AMPB,PIPL,AMP3,KEYT,MNUM)
C ----------------------------------------------------------------------
* IT SIMULATES A1  DECAY IN TAU REST FRAME WITH
* Z-AXIS ALONG A1  MOMENTUM
* IT CAN BE ALSO USED TO GENERATE K K PI AND K PI PI TAU DECAYS.
* INPUT PARAMETERS
* KEYT - ALGORITHM CONTROLLING SWITCH
*  2   - FLAT PHASE SPACE PIM1 PIM2 SYMMETRIZED STATISTICAL FACTOR 1/2
*  1   - LIKE 1 BUT PEAKED AROUND A1 AND RHO (TWO CHANNELS) MASSES.
*  3   - PEAKED AROUND OMEGA, ALL PARTICLES DIFFERENT
* OTHER- FLAT PHASE SPACE, ALL PARTICLES DIFFERENT
* AMP1 - MASS OF FIRST PI, ETC. (1-3)
* MNUM - MATRIX ELEMENT TYPE
*  0   - A1 MATRIX ELEMENT
* 1-6  - MATRIX ELEMENT FOR K PI PI, K K PI DECAY MODES
*  7   - PI- PI0 GAMMA MATRIX ELEMENT
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL  HV(4),PT(4),PN(4),PAA(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PR(4)
      REAL*4 RRR(5)
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
      XLAM(X,Y,Z)=SQRT(ABS((X-Y-Z)**2-4.0*Y*Z))
C AMRO, GAMRO IS ONLY A PARAMETER FOR GETING HIGHT EFFICIENCY
C
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1./2**17/PI**8
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C
      CALL RANMAR(RRR,5)
      RR=RRR(5)
C
      CALL CHOICE(MNUM,RR,ICHAN,PROB1,PROB2,PROB3,
     +            AMRX,GAMRX,AMRA,GAMRA,AMRB,GAMRB)
      IF     (ICHAN.EQ.1) THEN
        AMP1=AMPB
        AMP2=AMPA
      ELSEIF (ICHAN.EQ.2) THEN
        AMP1=AMPA
        AMP2=AMPB
      ELSE
        AMP1=AMPB
        AMP2=AMPA
      ENDIF
CAM
      RR1=RRR(1)
      AMS1=(AMP1+AMP2+AMP3)**2
      AMS2=(AMTAU-AMNUTA)**2
* PHASE SPACE WITH SAMPLING FOR A1  RESONANCE
      ALP1=ATAN((AMS1-AMRX**2)/AMRX/GAMRX)
      ALP2=ATAN((AMS2-AMRX**2)/AMRX/GAMRX)
      ALP=ALP1+RR1*(ALP2-ALP1)
      AM3SQ =AMRX**2+AMRX*GAMRX*TAN(ALP)
      AM3 =SQRT(AM3SQ)
      PHSPAC=PHSPAC*((AM3SQ-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)
      PHSPAC=PHSPAC*(ALP2-ALP1)
C MASS OF (REAL/VIRTUAL) RHO -
      RR2=RRR(2)
      AMS1=(AMP2+AMP3)**2
      AMS2=(AM3-AMP1)**2
      IF (ICHAN.LE.2) THEN
* PHASE SPACE WITH SAMPLING FOR RHO RESONANCE,
        ALP1=ATAN((AMS1-AMRA**2)/AMRA/GAMRA)
        ALP2=ATAN((AMS2-AMRA**2)/AMRA/GAMRA)
        ALP=ALP1+RR2*(ALP2-ALP1)
        AM2SQ =AMRA**2+AMRA*GAMRA*TAN(ALP)
        AM2 =SQRT(AM2SQ)
C --- THIS PART OF THE JACOBIAN WILL BE RECOVERED LATER ---------------
C     PHSPAC=PHSPAC*(ALP2-ALP1)
C     PHSPAC=PHSPAC*((AM2SQ-AMRA**2)**2+(AMRA*GAMRA)**2)/(AMRA*GAMRA)
C----------------------------------------------------------------------
      ELSE
* FLAT PHASE SPACE;
        AM2SQ=AMS1+   RR2*(AMS2-AMS1)
        AM2 =SQRT(AM2SQ)
        PHF0=(AMS2-AMS1)
      ENDIF
* RHO RESTFRAME, DEFINE PIPL AND PIM1
      ENQ1=(AM2SQ-AMP2**2+AMP3**2)/(2*AM2)
      ENQ2=(AM2SQ+AMP2**2-AMP3**2)/(2*AM2)
      PPI= ENQ1**2-AMP3**2
      PPPI=SQRT(ABS(ENQ1**2-AMP3**2))
C --- THIS PART OF JACOBIAN WILL BE RECOVERED LATER
      PHF1=(4*PI)*(2*PPPI/AM2)
* PI MINUS MOMENTUM IN RHO REST FRAME
      CALL SPHERA(PPPI,PIPL)
      PIPL(4)=ENQ1
* PI0 1 MOMENTUM IN RHO REST FRAME
      DO 10 I=1,3
   10 PIM1(I)=-PIPL(I)
      PIM1(4)=ENQ2
* A1 REST FRAME, DEFINE PIM2
*       RHO  MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AM3)*(AM3**2+AM2**2-AMP1**2)
      PR(3)= SQRT(ABS(PR(4)**2-AM2**2))
      PPI = PR(4)**2-AM2**2
*       PI0 2 MOMENTUM
      PIM2(1)=0
      PIM2(2)=0
      PIM2(4)=1./(2*AM3)*(AM3**2-AM2**2+AMP1**2)
      PIM2(3)=-PR(3)
      PHF2=(4*PI)*(2*PR(3)/AM3)
* OLD PIONS BOOSTED FROM RHO REST FRAME TO A1 REST FRAME
      EXE=(PR(4)+PR(3))/AM2
      CALL BOSTR3(EXE,PIPL,PIPL)
      CALL BOSTR3(EXE,PIM1,PIM1)
      RR3=RRR(3)
      RR4=RRR(4)
CAM   THET =PI*RR3
      THET =ACOS(-1.+2*RR3)
      PHI = 2*PI*RR4
      CALL ROTPOL(THET,PHI,PIPL)
      CALL ROTPOL(THET,PHI,PIM1)
      CALL ROTPOL(THET,PHI,PIM2)
      CALL ROTPOL(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE A1 AND NEUTRINO MOMENTA
* A1  MOMENTUM
      PAA(1)=0
      PAA(2)=0
      PAA(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AM3**2)
      PAA(3)= SQRT(ABS(PAA(4)**2-AM3**2))
      PPI   =          PAA(4)**2-AM3**2
      PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAU)
* TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AM3**2)
      PN(3)=-PAA(3)
C HERE WE CORRECT FOR THE JACOBIANS OF THE TWO CHAINS
C ---FIRST CHANNEL ------- PIM1+PIPL
      AMS1=(AMP2+AMP3)**2
      AMS2=(AM3-AMP1)**2
      ALP1=ATAN((AMS1-AMRA**2)/AMRA/GAMRA)
      ALP2=ATAN((AMS2-AMRA**2)/AMRA/GAMRA)
      XPRO = (PIM1(3)+PIPL(3))**2 +(PIM1(2)+PIPL(2))**2+(PIM1(1)+
     +PIPL(1))**2
      AM2SQ=-XPRO+(PIM1(4)+PIPL(4))**2
C JACOBIAN OF SPEEDING
      FF1 = ((AM2SQ-AMRA**2)**2+(AMRA*GAMRA)**2)/(AMRA*GAMRA)
      FF1 =FF1 *(ALP2-ALP1)
C LAMBDA OF RHO DECAY
      GG1 = (4*PI)*(XLAM(AM2SQ,AMP2**2,AMP3**2)/AM2SQ)
C LAMBDA OF A1 DECAY
      GG1 =GG1 *(4*PI)*SQRT(4*XPRO/AM3SQ)
      XJAJE=GG1*(AMS2-AMS1)
C ---SECOND CHANNEL ------ PIM2+PIPL
      AMS1=(AMP1+AMP3)**2
      AMS2=(AM3-AMP2)**2
      ALP1=ATAN((AMS1-AMRB**2)/AMRB/GAMRB)
      ALP2=ATAN((AMS2-AMRB**2)/AMRB/GAMRB)
      XPRO = (PIM2(3)+PIPL(3))**2 +(PIM2(2)+PIPL(2))**2+(PIM2(1)+
     +PIPL(1))**2
      AM2SQ=-XPRO+(PIM2(4)+PIPL(4))**2
      FF2 = ((AM2SQ-AMRB**2)**2+(AMRB*GAMRB)**2)/(AMRB*GAMRB)
      FF2 =FF2 *(ALP2-ALP1)
      GG2 = (4*PI)*(XLAM(AM2SQ,AMP1**2,AMP3**2)/AM2SQ)
      GG2 =GG2 *(4*PI)*SQRT(4*XPRO/AM3SQ)
      XJADW=GG2*(AMS2-AMS1)
C
      A1=0.0
      A2=0.0
      A3=0.0
      XJAC1=FF1*GG1
      XJAC2=FF2*GG2
      IF (ICHAN.EQ.2) THEN
        XJAC3=XJADW
      ELSE
        XJAC3=XJAJE
      ENDIF
      IF (XJAC1.NE.0.0) A1=PROB1/XJAC1
      IF (XJAC2.NE.0.0) A2=PROB2/XJAC2
      IF (XJAC3.NE.0.0) A3=PROB3/XJAC3
C
      IF (A1+A2+A3.NE.0.0) THEN
        PHSPAC=PHSPAC/(A1+A2+A3)
      ELSE
        PHSPAC=0.0
      ENDIF
      IF(ICHAN.EQ.2) THEN
        DO 20 I=1,4
          X=PIM1(I)
          PIM1(I)=PIM2(I)
   20   PIM2(I)=X
      ENDIF
* ALL PIONS BOOSTED FROM A1  REST FRAME TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO NEUTRINO MOMENTUM
      EXE=(PAA(4)+PAA(3))/AM3
      CALL BOSTR3(EXE,PIPL,PIPL)
      CALL BOSTR3(EXE,PIM1,PIM1)
      CALL BOSTR3(EXE,PIM2,PIM2)
      CALL BOSTR3(EXE,PR,PR)
C PARTIAL WIDTH CONSISTS OF PHASE SPACE AND AMPLITUDE
      IF (MNUM.EQ.8) THEN
        CALL DAMPOG(PT,PN,PIM1,PIM2,PIPL,AMPLIT,HV)
C      ELSEIF (MNUM.EQ.0) THEN
C        CALL DAMPAA(PT,PN,PIM1,PIM2,PIPL,AMPLIT,HV)
      ELSE
        CALL DAMPPK(MNUM,PT,PN,PIM1,PIM2,PIPL,AMPLIT,HV)
      ENDIF
      IF (KEYT.EQ.1.OR.KEYT.EQ.2) THEN
C THE STATISTICAL FACTOR FOR IDENTICAL PI'S IS CANCELLED WITH
C TWO, FOR TWO MODES OF A1 DECAY NAMELLY PI+PI-PI- AND PI-PI0PI0
        PHSPAC=PHSPAC*2.0
        PHSPAC=PHSPAC/2.
      ENDIF
      DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      FUNCTION DQCD(ICOSFI,IPART,IP,XP,ZP,Y)

C...FIRST ORDER QCD MATRIX ELEMENTS FROM R.D. PECCEI AND R. RUCKL:
C...NUCL. PHYS. B162 (1980) 125

C...CONSTANTS C1 TO C5 ARE RESP. 2/3/PI, 1/4/PI, 4/3/PI, 1/2/PI, 1/PI
      DATA C1,C2,C3,C4,C5/0.2122066,0.0795775,0.4244132,0.1591549,
     &     0.3183099/

      IF(ICOSFI.EQ.0) THEN
        IF(IPART.EQ.1) THEN
          IF(IP.EQ.1) THEN
            DQCD=C1*((ZP**2+XP**2)/(1.-XP)/(1.-ZP)+2.*(XP*ZP+1.))
          ELSEIF(IP.EQ.2) THEN
            DQCD=C1*4.*XP*ZP
          ELSEIF(IP.EQ.3) THEN
            DQCD=C1*((ZP**2+XP**2)/(1.-XP)/(1.-ZP)+2.*(XP+ZP))
          ELSE
            WRITE(6,10000) ICOSFI,IPART,IP
          ENDIF
        ELSEIF(IPART.EQ.2) THEN
          IF(IP.EQ.1) THEN
            DQCD=C2*(XP**2+(1.-XP)**2)*(ZP**2+(1.-ZP)**2)/(1.-ZP)/ZP
          ELSEIF(IP.EQ.2) THEN
            DQCD=C2*8.*XP*(1.-XP)
          ELSEIF(IP.EQ.3) THEN
            DQCD=C2*(XP**2+(1.-XP)**2)*(ZP-(1.-ZP))/(1.-ZP)/ZP
          ELSE
            WRITE(6,10000) ICOSFI,IPART,IP
          ENDIF
        ELSE
          WRITE(6,10000) ICOSFI,IPART,IP
        ENDIF

      ELSEIF(ICOSFI.EQ.1) THEN
        IF(IPART.EQ.1) THEN
          IF(IP.EQ.1) THEN
            DQCD=C3*Y*SQRT((1.-Y)*XP*ZP/(1.-XP)/(1.-ZP))*
     &      (1.-2./Y)*(1.-ZP-XP+2.*XP*ZP)
          ELSEIF(IP.EQ.3) THEN
            DQCD=C3*Y*SQRT((1.-Y)*XP*ZP/(1.-XP)/(1.-ZP))*
     &      (1.-XP-ZP)
          ELSE
            WRITE(6,10000) ICOSFI,IPART,IP
          ENDIF
        ELSEIF(IPART.EQ.2) THEN
          IF(IP.EQ.1) THEN
            DQCD=C4*Y*SQRT((1.-Y)*XP*(1.-XP)/ZP/(1.-ZP))*
     &      (1.-2./Y)*(1.-2.*ZP)*(1.-2.*XP)
          ELSEIF(IP.EQ.3) THEN
            DQCD=C4*Y*SQRT((1.-Y)*XP*(1.-XP)/ZP/(1.-ZP))*
     &      (1.-2.*XP)
          ELSE
            WRITE(6,10000) ICOSFI,IPART,IP
          ENDIF
        ENDIF

      ELSEIF(ICOSFI.EQ.2) THEN
        IF(IPART.EQ.1) THEN
          DQCD=C3*(1.-Y)*XP*ZP
        ELSEIF(IPART.EQ.2) THEN
          DQCD=C5*(1.-Y)*XP*(1.-XP)
        ELSE
          WRITE(6,10000) ICOSFI,IPART,IP
        ENDIF

      ELSE
        WRITE(6,10000) ICOSFI,IPART,IP
      ENDIF
      RETURN

10000 FORMAT(' ERROR IN ROUTINE DQCD     ',
     &' ICOSFI, IPART, IP = ',3I10)
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      FUNCTION DQCDI(IPART,IP,XP,ZPMIN,ZPMAX)

C...FIRST ORDER QCD MATRIX ELEMENTS AS IN FUNCTION DQCD BUT ANALYTICALLY
C...INTEGRATED OVER ZP FROM ZPMIN TO ZPMAX, ALSO A FACTOR 1/(1-XP) IS
C...FACTORED OUT (SINCE XP CHOSEN RANDOMLY ACCORDING TO 1/(1-XP) DISTR.)

      DATA C1,C2/0.2122066,0.0795775/

      IF(IPART.EQ.1) THEN
        IF(IP.EQ.1) THEN
          ZLOG=ALOG(ZPMAX/ZPMIN)
          DQCDI=C1*(XP**2*ZLOG+ZPMIN-ZPMAX+(ZPMIN**2-ZPMAX**2)/2.+ZLOG+
     &    XP*(1.-XP)*(ZPMAX**2-ZPMIN**2)+2.*(1.-XP)*(ZPMAX-ZPMIN))
        ELSEIF(IP.EQ.2) THEN
          DQCDI=C1*2.*XP*(1.-XP)*(ZPMAX**2-ZPMIN**2)
        ELSEIF(IP.EQ.3) THEN
          ZLOG=ALOG(ZPMAX/ZPMIN)
          DQCDI=C1*(XP**2*ZLOG+ZPMIN-ZPMAX+(ZPMIN**2-ZPMAX**2)/2.+ZLOG+
     &    2.*XP*(1.-XP)*(ZPMAX-ZPMIN)+(1.-XP)*(ZPMAX**2-ZPMIN**2))
        ELSE
          WRITE(6,10000) IPART,IP
        ENDIF

      ELSEIF(IPART.EQ.2) THEN
        IF(IP.EQ.1) THEN
          DQCDI=C2*(1.-XP)*(XP**2+(1.-XP)**2)*(2.*(ZPMIN-ZPMAX)+
     &    ALOG(ZPMAX**2/ZPMIN**2))
        ELSEIF(IP.EQ.2) THEN
          DQCDI=C2*8.*XP*(1.-XP)**2*(ZPMAX-ZPMIN)
        ELSEIF(IP.EQ.3) THEN
          DQCDI=0.
        ELSE
          WRITE(6,10000) IPART,IP
        ENDIF

      ELSE
        WRITE(6,10000) IPART,IP
      ENDIF
      RETURN

10000 FORMAT(' ERROR IN ROUTINE DQCDI     ',
     &' IPART, IP = ',2I10)
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DRCMU(DGAMT,HV,PH,PAA,XA,QP,XN,IELMU)
      IMPLICIT REAL*8 (A-H,O-Z)
C ----------------------------------------------------------------------
* IT SIMULATES E,MU CHANNELS OF TAU  DECAY IN ITS REST FRAME WITH
* QED ORDER ALPHA CORRECTIONS
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / INOUT / INUT,IOUT
      COMMON / TAURAD / XK0DEC,ITDKRC
      REAL*8            XK0DEC
      REAL*8  HV(4),PT(4),PH(4),PAA(4),XA(4),QP(4),XN(4)
      REAL*8  PR(4)
      REAL*4 RRR(6)
      LOGICAL IHARD
      DATA PI /3.141592653589793238462643D0/
      XLAM(X,Y,Z)=SQRT((X-Y-Z)**2-4.0*Y*Z)
C AMRO, GAMRO IS ONLY A PARAMETER FOR GETING HIGHT EFFICIENCY
C
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1./2**17/PI**8
      AMTAX=AMTAU
C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAX
C
      CALL RANMAR(RRR,6)
C
      IF (IELMU.EQ.1) THEN
        AMU=AMEL
      ELSE
        AMU=AMMU
      ENDIF
C
      PRHARD=0.30D0
      IF ( ITDKRC.EQ.0) PRHARD=0D0
      PRSOFT=1.-PRHARD
      IF(PRSOFT.LT.0.1) THEN
        PRINT *, 'ERROR IN DRCMU; PRSOFT=',PRSOFT
        STOP
      ENDIF
C
      RR5=RRR(5)
      IHARD=(RR5.GT.PRSOFT)
      IF (IHARD) THEN
C                     TAU DECAY TO `TAU+PHOTON'
        RR1=RRR(1)
        AMS1=(AMU+AMNUTA)**2
        AMS2=(AMTAX)**2
        XK1=1-AMS1/AMS2
        XL1=LOG(XK1/2/XK0DEC)
        XL0=LOG(2*XK0DEC)
        XK=EXP(XL1*RR1+XL0)
        AM3SQ=(1-XK)*AMS2
        AM3 =SQRT(AM3SQ)
        PHSPAC=PHSPAC*AMS2*XL1*XK
        PHSPAC=PHSPAC/PRHARD
      ELSE
        AM3=AMTAX
        PHSPAC=PHSPAC*2**6*PI**3
        PHSPAC=PHSPAC/PRSOFT
      ENDIF
C MASS OF NEUTRINA SYSTEM
      RR2=RRR(2)
      AMS1=(AMNUTA)**2
      AMS2=(AM3-AMU)**2
CAM
CAM
* FLAT PHASE SPACE;
      AM2SQ=AMS1+   RR2*(AMS2-AMS1)
      AM2 =SQRT(AM2SQ)
      PHSPAC=PHSPAC*(AMS2-AMS1)
* NEUTRINA REST FRAME, DEFINE XN AND XA
      ENQ1=(AM2SQ+AMNUTA**2)/(2*AM2)
      ENQ2=(AM2SQ-AMNUTA**2)/(2*AM2)
      PPI= ENQ1**2-AMNUTA**2
      PPPI=SQRT(ABS(ENQ1**2-AMNUTA**2))
      PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM2)
* NU TAU IN NUNU REST FRAME
      CALL SPHERD(PPPI,XN)
      XN(4)=ENQ1
* NU LIGHT IN NUNU REST FRAME
      DO 10 I=1,3
   10 XA(I)=-XN(I)
      XA(4)=ENQ2
* TAU' REST FRAME, DEFINE QP (MUON
*       NUNU  MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1.D0/(2*AM3)*(AM3**2+AM2**2-AMU**2)
      PR(3)= SQRT(ABS(PR(4)**2-AM2**2))
      PPI = PR(4)**2-AM2**2
*       MUON MOMENTUM
      QP(1)=0
      QP(2)=0
      QP(4)=1.D0/(2*AM3)*(AM3**2-AM2**2+AMU**2)
      QP(3)=-PR(3)
      PHSPAC=PHSPAC*(4*PI)*(2*PR(3)/AM3)
* NEUTRINA BOOSTED FROM THEIR FRAME TO TAU' REST FRAME
      EXE=(PR(4)+PR(3))/AM2
      CALL BOSTD3(EXE,XN,XN)
      CALL BOSTD3(EXE,XA,XA)
      RR3=RRR(3)
      RR4=RRR(4)
      IF (IHARD) THEN
        EPS=4*(AMU/AMTAX)**2
        XL1=LOG((2+EPS)/EPS)
        XL0=LOG(EPS)
        ETA  =EXP(XL1*RR3+XL0)
        CTHET=1+EPS-ETA
        THET =ACOS(CTHET)
        PHSPAC=PHSPAC*XL1/2*ETA
        PHI = 2*PI*RR4
        CALL ROTPOX(THET,PHI,XN)
        CALL ROTPOX(THET,PHI,XA)
        CALL ROTPOX(THET,PHI,QP)
        CALL ROTPOX(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE TAU' AND GAMMA MOMENTA
* TAU'  MOMENTUM
        PAA(1)=0
        PAA(2)=0
        PAA(4)=1/(2*AMTAX)*(AMTAX**2+AM3**2)
        PAA(3)= SQRT(ABS(PAA(4)**2-AM3**2))
        PPI   =          PAA(4)**2-AM3**2
        PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAX)
* GAMMA MOMENTUM
        PH(1)=0
        PH(2)=0
        PH(4)=PAA(3)
        PH(3)=-PAA(3)
* ALL MOMENTA BOOSTED FROM TAU' REST FRAME TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO PHOTON MOMENTUM
        EXE=(PAA(4)+PAA(3))/AM3
        CALL BOSTD3(EXE,XN,XN)
        CALL BOSTD3(EXE,XA,XA)
        CALL BOSTD3(EXE,QP,QP)
        CALL BOSTD3(EXE,PR,PR)
      ELSE
        THET =ACOS(-1.+2*RR3)
        PHI = 2*PI*RR4
        CALL ROTPOX(THET,PHI,XN)
        CALL ROTPOX(THET,PHI,XA)
        CALL ROTPOX(THET,PHI,QP)
        CALL ROTPOX(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE TAU' AND GAMMA MOMENTA
* TAU'  MOMENTUM
        PAA(1)=0
        PAA(2)=0
        PAA(4)=AMTAX
        PAA(3)=0
* GAMMA MOMENTUM
        PH(1)=0
        PH(2)=0
        PH(4)=0
        PH(3)=0
      ENDIF
C PARTIAL WIDTH CONSISTS OF PHASE SPACE AND AMPLITUDE
      CALL DAMPRY(ITDKRC,XK0DEC,PH,XA,QP,XN,AMPLIT,HV)
      DGAMT=1/(2.*AMTAX)*AMPLIT*PHSPAC
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 25/07/94  17.30.49  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 15/07/94  14.09.39  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      FUNCTION DSIGMA(XP)

C...DIFFERENTIAL CROSS SECTION FOR FIRST ORDER QCD PROCESSES.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      DIMENSION XPQ(-6:6),PQH(17,2)

      DSIGMA=0.
      DO 10 I=1,17
        PQH(I,1)=0.
        PQH(I,2)=0.
   10 PQ(I)=0.

      AMU=ULMASS(1)
      IF(LST(20).EQ.0.OR.LST(17).EQ.0) THEN
        IL=1
        IU=3
        IF(LST(23).EQ.1.OR.LST(24).EQ.3) IU=2
      ELSE
        IL=LST(20)
        IU=LST(20)
      ENDIF
      XI=X/XP
      ZPMIN=(1.-X)*XP/(XP-X)*PARL(27)
      IF(ZPMIN.GE.0.5) RETURN
      ZPMAX=1.D0-DBLE(ZPMIN)
      CALL LNSTRF(XI,Q2,XPQ)
      IF(LST(24).EQ.3) GOTO 80

C...GLUON BREMSSTRAHLUNG PROCESS, I.E. QG-EVENT.
   20 DO 60   IP=IL,IU
        SIG=DQCDI(1,IP,XP,ZPMIN,ZPMAX)
        SGN=SIGN(1.,5.-2.*IP)
        DO 50   IH=1,2
          IF(IH.EQ.1) THEN
            IF(PARL(6).GT.0.99) GOTO 50
            IF(LST(20).EQ.0.AND.LST(30).NE.-1) GOTO 50
          ELSEIF(IH.EQ.2) THEN
            IF(PARL(6).LT.-0.99) GOTO 50
            IF(LST(20).EQ.0.AND.LST(30).NE.1) GOTO 50
          ENDIF
          IF(LST(20).NE.0) LST(30)=SIGN(1.,IH-1.5)
          IF(LST(23).NE.2) THEN
            DO 30   I=1,LST(12)
              WQ=XPQ(I)*SIG*(EWQC(1,IH,I)+SGN*EWQC(2,IH,I))
              WQB=XPQ(-I)*SIG*SGN*(EWQC(1,IH,I)+SGN*EWQC(2,IH,I))
C...INCLUDE Y-DEPENDENCE.
              WQ=WQ*PARI(23+IP)
              WQB=WQB*PARI(23+IP)
              PQH(I,IH)=PQH(I,IH)+WQ
              PQH(I+LST(12),IH)=PQH(I+LST(12),IH)+WQB
              PQH(17,IH)=PQH(17,IH)+WQ+WQB
   30       CONTINUE
          ELSEIF(LST(23).EQ.2) THEN
C...ZERO CC CROSS-SECTION FOR ONE HELICITY STATE.
            IF(KSAVE(1).LT.0.AND.IH.EQ.1 .OR.KSAVE(1).GT.0.AND.IH.EQ.2)
     +      GOTO 50
            IF(IP.EQ.3) THEN
              TQ=-LST(30)
              TQB=-TQ
            ELSE
              TQ=1.
              TQB=1.
            ENDIF
            DO 40   I=1,LST(12)
              T1=-K(3,2)*QC(I)
              IF(T1.GT.0) THEN
                WQ=XPQ(I)*SIG*TQ
                WQB=0.
              ELSE
                WQB=XPQ(-I)*SIG*TQB
                WQ=0.
              ENDIF
C...INCLUDE Y-DEPENDENCE.
              WQ=WQ*PARI(23+IP)
              WQB=WQB*PARI(23+IP)
              PQH(I,IH)=PQH(I,IH)+WQ
              PQH(I+LST(12),IH)=PQH(I+LST(12),IH)+WQB
              PQH(17,IH)=PQH(17,IH)+WQ+WQB
   40       CONTINUE
          ENDIF
   50   CONTINUE
   60 CONTINUE
      DO 70   I=1,17
   70 PQ(I)=(1.-PARL(6))/2.*PQH(I,1)+(1.+PARL(6))/2.*PQH(I,2)
      IF(LST(20).EQ.0) THEN
C...SIMULATION: CROSS SECTION FOR CHOSEN HELICITY STATE ONLY.
        IH=1
        IF(LST(30).EQ.1) IH=2
        DSIGMA=PQH(17,IH)
*     WRITE(*,*)'DSIGMA_1=',DSIGMA
      ELSE
C...INTEGRATION: NORMALIZE AND INCLUDE ALPHA_S*1/(1-XP)
        DSIGMA=PQ(17)/PARI(20)*PARL(25)/(1.-XP)
*     WRITE(*,*)'DSIGMA_2=',DSIGMA
        IF(LST(17).EQ.0) THEN
C...FIXED BEAM ENERGY, MAX OF DSIGMA/DXP FOR L- AND R-HANDED LEPTON.
          IF(PQH(17,1).GT.PARI(15)) PARI(15)=PQH(17,1)
          IF(PQH(17,2).GT.PARI(16)) PARI(16)=PQH(17,2)
        ELSE
C...VARIABLE BEAM ENERGY, MAX OF DSIGMA/DXP FOR S, T, I CONTRIBUTIONS.
          IF(PQ(17)/PARI(23+LST(20)).GT.PARI(14+LST(20)))
     +    PARI(14+LST(20))=PQ(17)/PARI(23+LST(20))
        ENDIF
      ENDIF
      RETURN

C...BOSON-GLUON FUSION, I.E. QQ-EVENT.
   80 S13=Q2*(1.-XP)/XP
      IF(S13.LT.4.*AMU**2) RETURN
      DO 120  IP=IL,IU
        SIG=XPQ(0)*DQCDI(2,IP,XP,ZPMIN,ZPMAX)
        DO 110  IH=1,2
          IF(IH.EQ.1) THEN
            IF(PARL(6).GT.0.99) GOTO 110
            IF(LST(20).EQ.0.AND.LST(30).NE.-1) GOTO 110
          ELSEIF(IH.EQ.2) THEN
            IF(PARL(6).LT.-0.99) GOTO 110
            IF(LST(20).EQ.0.AND.LST(30).NE.1) GOTO 110
          ENDIF
          IF(LST(20).NE.0) LST(30)=SIGN(1.,IH-1.5)
          IF(LST(23).NE.2) THEN
            DO 90   I=1,LST(13)
              IF(S13.LT.4.*ULMASS(I)**2) GOTO 90
              WQ=SIG/2.*(EWQC(1,IH,I)+EWQC(2,IH,I))
              WQB=WQ
C...INCLUDE Y-DEPENDENCE.
              WQ=WQ*PARI(23+IP)
              WQB=WQB*PARI(23+IP)
              PQH(I,IH)=PQH(I,IH)+WQ
              PQH(I+LST(13),IH)=PQH(I+LST(13),IH)+WQB
              PQH(17,IH)=PQH(17,IH)+WQ+WQB
   90       CONTINUE
          ELSEIF(LST(23).EQ.2) THEN
C...ZERO CC CROSS-SECTION FOR ONE HELICITY STATE.
            IF(KSAVE(1).LT.0.AND.IH.EQ.1 .OR.KSAVE(1).GT.0.AND.IH.EQ.2)
     +      GOTO 110
            DO 100  I=1,LST(13)
              IF(S13.LT.(AMU+ULMASS(I))**2) GOTO 100
              IF(K(3,2)*QC(I).LT.0) THEN
                WQ=SIG
                WQB=0.
              ELSE
                WQB=SIG
                WQ=0.
              ENDIF
C...INCLUDE Y-DEPENDENCE.
              WQ=WQ*PARI(23+IP)
              WQB=WQB*PARI(23+IP)
              PQH(I,IH)=PQH(I,IH)+WQ
              PQH(I+LST(13),IH)=PQH(I+LST(13),IH)+WQB
              PQH(17,IH)=PQH(17,IH)+WQ+WQB
  100       CONTINUE
          ENDIF
  110   CONTINUE
  120 CONTINUE
      DO 130  I=1,17
  130 PQ(I)=(1.-PARL(6))/2.*PQH(I,1)+(1.+PARL(6))/2.*PQH(I,2)
      IF(LST(20).EQ.0) THEN
C...SIMULATION: CROSS SECTION FOR CHOSEN HELICITY STATE ONLY.
        IH=1
        IF(LST(30).EQ.1) IH=2
        DSIGMA=PQH(17,IH)
*     WRITE(*,*)'DSIGMA_3=',DSIGMA
      ELSE
C...INTEGRATION: NORMALIZE AND INCLUDE ALPHA_S*1/(1-XP)
        DSIGMA=PQ(17)/PARI(20)*PARL(25)/(1.-XP)
*     WRITE(*,*)'DSIGMA_4=',DSIGMA
        IF(LST(17).EQ.0) THEN
C...FIXED BEAM ENERGY, MAX OF DSIGMA/DXP FOR L- AND R-HANDED LEPTON.
          IF(PQH(17,1).GT.PARI(18)) PARI(18)=PQH(17,1)
          IF(PQH(17,2).GT.PARI(19)) PARI(19)=PQH(17,2)
        ELSE
C...VARIABLE BEAM ENERGY, MAX OF DSIGMA/DXP FOR S, T, I CONTRIBUTIONS.
          IF(PQ(17)/PARI(23+LST(20)).GT.PARI(17+LST(20)))
     +    PARI(17+LST(20))=PQ(17)/PARI(23+LST(20))
        ENDIF
      ENDIF
      RETURN
      END
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      FUNCTION DUPPER(V1)

C...UPPER LIMIT ON SECOND VARIABLE (Y, Q**2 OR W**2) DEPENDING ON FIRST
C...VARIABLE X=V1. USED FOR INTEGRATING DIFFERENTIAL CROSS-SECTION.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
C...CMS ENERGY SQUARED AND TARGET NUCLEON MASS.
      S=PARL(21)
      PM2=PSAVE(3,2,5)**2
      IF(LST(31).EQ.1) THEN
        DUPPER=MIN(Q2MAX,V1*YMAX*S,(W2MAX-PM2)*V1/MAX(1.-V1,1.E-22))
      ELSEIF(LST(31).EQ.2) THEN
        DUPPER=MIN(YMAX,Q2MAX/(S*V1),(W2MAX-PM2)/MAX(S*(1.-V1),1.E-22))
      ELSEIF(LST(31).EQ.3) THEN
        DUPPER=MIN(W2MAX,(1.-V1)*YMAX*S+PM2,
     &  Q2MAX*(1.-V1)/MAX(V1,1.E-22)+PM2)
      ENDIF
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE DVNOPT

C...CHANGE OF DEFAULT OPTIONS IN DIVONNE

      COMMON /PRINT/ IPRDIV
      COMMON /LPFLAG/ LST3
      IPRDIV=0
      IF(LST3.GE.2) IPRDIV=1000
      IF(LST3.GE.4) IPRDIV=10
      IF(LST3.GE.4) WRITE(6,10000) IPRDIV
      RETURN
10000 FORMAT(5X,'DIVON4 PRINT FLAG CHANGED: IPRDIV =',I5)
      END
*CMZ :  1.01/50 22/05/96  18.06.08  by  Piero Zucchelli
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DWLNEW(KTO,ISGN,PNU,PWB,PNPI,MODE)
C ----------------------------------------------------------------------
C LORENTZ TRANSFORMATION TO CMSYSTEM AND
C UPDATING OF HEPEVT RECORD
C
C ISGN = 1/-1 FOR TAU-/TAU+
C
C     CALLED BY : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / MDECOMP /IDFFIN(9,NMODE),MULPIK(NMODE)
     &                ,NAMES
      CHARACTER NAMES(NMODE)*31
      REAL  PNU(4),PWB(4),PNPI(4,9)
      REAL  PPI(4)
C
      JNPI=MODE-7
C POSITION OF DECAYING PARTICLE
      IF(KTO.EQ. 1) THEN
        NPS=3
      ELSE
        NPS=4
      ENDIF
C
C TAU NEUTRINO (NU_TAU IS 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C W BOSON (W+ IS 24)
      CALL TRALO4(KTO,PWB,PWB,AM)
      CALL FILHEP(0,1,-24*ISGN,NPS,NPS,0,0,PWB,AM,.TRUE.)
C
C MULTI PI MODE JNPI
C
C GET MULTIPLICITY OF MODE JNPI
      ND=MULPIK(JNPI)
      DO I=1,ND
        KFPI=LUNPIK(IDFFIN(I,JNPI),-ISGN)
C FOR CHARGED CONJUGATE CASE, CHANGE CHARGED PIONS ONLY
C        IF(KFPI.NE.111)KFPI=KFPI*ISGN
        DO J=1,4
          PPI(J)=PNPI(J,I)
        END DO
        CALL TRALO4(KTO,PPI,PPI,AM)
        CALL FILHEP(0,1,KFPI,-I,-I,0,0,PPI,AM,.TRUE.)
      END DO
C
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DWLUAA(KTO,ISGN,PNU,PAA,PIM1,PIM2,PIPL,JAA)
C ----------------------------------------------------------------------
C LORENTZ TRANSFORMATION TO CMSYSTEM AND
C UPDATING OF HEPEVT RECORD
C
C ISGN = 1/-1 FOR TAU-/TAU+
C JAA  = 1 (2) FOR A_1- DECAY TO PI+ 2PI- (PI- 2PI0)
C
C     CALLED BY : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PAA(4),PIM1(4),PIM2(4),PIPL(4)
C
C POSITION OF DECAYING PARTICLE:
      IF(KTO.EQ. 1) THEN
        NPS=3
      ELSE
        NPS=4
      ENDIF
C
C TAU NEUTRINO (NU_TAU IS 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C CHARGED A_1 MESON (A_1+ IS 20213)
      CALL TRALO4(KTO,PAA,PAA,AM)
      CALL FILHEP(0,1,-20213*ISGN,NPS,NPS,0,0,PAA,AM,.TRUE.)
C
C TWO POSSIBLE DECAYS OF THE CHARGED A1 MESON
      IF(JAA.EQ.1) THEN
C
C A1  --> PI+ PI-  PI- (OR CHARGED CONJUGATE)
C
C PI MINUS (OR C.C.) (PI+ IS 211)
        CALL TRALO4(KTO,PIM2,PIM2,AM)
        CALL FILHEP(0,1,-211*ISGN,-1,-1,0,0,PIM2,AM,.TRUE.)
C
C PI MINUS (OR C.C.) (PI+ IS 211)
        CALL TRALO4(KTO,PIM1,PIM1,AM)
        CALL FILHEP(0,1,-211*ISGN,-2,-2,0,0,PIM1,AM,.TRUE.)
C
C PI PLUS (OR C.C.) (PI+ IS 211)
        CALL TRALO4(KTO,PIPL,PIPL,AM)
        CALL FILHEP(0,1, 211*ISGN,-3,-3,0,0,PIPL,AM,.TRUE.)
C
      ELSE IF (JAA.EQ.2) THEN
C
C A1  --> PI- PI0  PI0 (OR CHARGED CONJUGATE)
C
C PI ZERO (PI0 IS 111)
        CALL TRALO4(KTO,PIM2,PIM2,AM)
        CALL FILHEP(0,1,111,-1,-1,0,0,PIM2,AM,.TRUE.)
C
C PI ZERO (PI0 IS 111)
        CALL TRALO4(KTO,PIM1,PIM1,AM)
        CALL FILHEP(0,1,111,-2,-2,0,0,PIM1,AM,.TRUE.)
C
C PI MINUS (OR C.C.) (PI+ IS 211)
        CALL TRALO4(KTO,PIPL,PIPL,AM)
        CALL FILHEP(0,1,-211*ISGN,-3,-3,0,0,PIPL,AM,.TRUE.)
C
      ENDIF
C
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DWLUEL(KTO,ISGN,PNU,PWB,PEL,PNE)
C ----------------------------------------------------------------------
C LORENTZ TRANSFORMATION TO CMSYSTEM AND
C UPDATING OF HEPEVT RECORD
C
C ISGN = 1/-1 FOR TAU-/TAU+
C
C     CALLED BY : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PWB(4),PEL(4),PNE(4)
C
C POSITION OF DECAYING PARTICLE:
      IF(KTO.EQ. 1) THEN
        NPS=3
      ELSE
        NPS=4
      ENDIF
C
C TAU NEUTRINO (NU_TAU IS 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C W BOSON (W+ IS 24)
      CALL TRALO4(KTO,PWB,PWB,AM)
C     CALL FILHEP(0,2,-24*ISGN,NPS,NPS,0,0,PWB,AM,.TRUE.)
C
C ELECTRON (E- IS 11)
      CALL TRALO4(KTO,PEL,PEL,AM)
      CALL FILHEP(0,1,11*ISGN,NPS,NPS,0,0,PEL,AM,.FALSE.)
C
C ANTI ELECTRON NEUTRINO (NU_E IS 12)
      CALL TRALO4(KTO,PNE,PNE,AM)
      CALL FILHEP(0,1,-12*ISGN,NPS,NPS,0,0,PNE,AM,.TRUE.)
C
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DWLUKK (KTO,ISGN,PKK,PNU)
C ----------------------------------------------------------------------
C LORENTZ TRANSFORMATION TO CMSYSTEM AND
C UPDATING OF HEPEVT RECORD
C
C ISGN = 1/-1 FOR TAU-/TAU+
C
C ----------------------------------------------------------------------
C
      REAL PKK(4),PNU(4)
C
C POSITION OF DECAYING PARTICLE
      IF (KTO.EQ.1) THEN
        NPS=3
      ELSE
        NPS=4
      ENDIF
C
C TAU NEUTRINO (NU_TAU IS 16)
      CALL TRALO4 (KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C K MESON (K+ IS 321)
      CALL TRALO4 (KTO,PKK,PKK,AM)
      CALL FILHEP(0,1,-321*ISGN,NPS,NPS,0,0,PKK,AM,.TRUE.)
C
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DWLUKS(KTO,ISGN,PNU,PKS,PKK,PPI,JKST)
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
C ----------------------------------------------------------------------
C LORENTZ TRANSFORMATION TO CMSYSTEM AND
C UPDATING OF HEPEVT RECORD
C
C ISGN = 1/-1 FOR TAU-/TAU+
C JKST=10 (20) CORRESPONDS TO K0B PI- (K- PI0) DECAY
C
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PKS(4),PKK(4),PPI(4)
C
C POSITION OF DECAYING PARTICLE
      IF(KTO.EQ. 1) THEN
        NPS=3
      ELSE
        NPS=4
      ENDIF
C
C TAU NEUTRINO (NU_TAU IS 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C CHARGED K* MESON (K*+ IS 323)
      CALL TRALO4(KTO,PKS,PKS,AM)
      CALL FILHEP(0,1,-323*ISGN,NPS,NPS,0,0,PKS,AM,.TRUE.)
C
C TWO POSSIBLE DECAY MODES OF CHARGED K*
      IF(JKST.EQ.10) THEN
C
C K*- --> PI- K0B (OR CHARGED CONJUGATE)
C
C CHARGED PI MESON  (PI+ IS 211)
        CALL TRALO4(KTO,PPI,PPI,AM)
        CALL FILHEP(0,1,-211*ISGN,-1,-1,0,0,PPI,AM,.TRUE.)
C
        BRAN=BRK0B
        IF (ISGN.EQ.-1) BRAN=BRK0
C K0 --> K0_LONG (IS 130) / K0_SHORT (IS 310) = 1/1
        CALL RANMAR(XIO,1)
        IF(XIO.GT.BRAN) THEN
          K0TYPE = 130
        ELSE
          K0TYPE = 310
        ENDIF
C
        CALL TRALO4(KTO,PKK,PKK,AM)
        CALL FILHEP(0,1,K0TYPE,-2,-2,0,0,PKK,AM,.TRUE.)
C
      ELSE IF(JKST.EQ.20) THEN
C
C K*- --> PI0 K-
C
C PI ZERO (PI0 IS 111)
        CALL TRALO4(KTO,PPI,PPI,AM)
        CALL FILHEP(0,1,111,-1,-1,0,0,PPI,AM,.TRUE.)
C
C CHARGED K MESON (K+ IS 321)
        CALL TRALO4(KTO,PKK,PKK,AM)
        CALL FILHEP(0,1,-321*ISGN,-2,-2,0,0,PKK,AM,.TRUE.)
C
      ENDIF
C
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DWLUMU(KTO,ISGN,PNU,PWB,PMU,PNM)
C ----------------------------------------------------------------------
C LORENTZ TRANSFORMATION TO CMSYSTEM AND
C UPDATING OF HEPEVT RECORD
C
C ISGN = 1/-1 FOR TAU-/TAU+
C
C     CALLED BY : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PWB(4),PMU(4),PNM(4)
C
C POSITION OF DECAYING PARTICLE:
      IF(KTO.EQ. 1) THEN
        NPS=3
      ELSE
        NPS=4
      ENDIF
C
C TAU NEUTRINO (NU_TAU IS 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C W BOSON (W+ IS 24)
      CALL TRALO4(KTO,PWB,PWB,AM)
C     CALL FILHEP(0,2,-24*ISGN,NPS,NPS,0,0,PWB,AM,.TRUE.)
C
C MUON (MU- IS 13)
      CALL TRALO4(KTO,PMU,PMU,AM)
      CALL FILHEP(0,1,13*ISGN,NPS,NPS,0,0,PMU,AM,.FALSE.)
C
C ANTI MUON NEUTRINO (NU_MU IS 14)
      CALL TRALO4(KTO,PNM,PNM,AM)
      CALL FILHEP(0,1,-14*ISGN,NPS,NPS,0,0,PNM,AM,.TRUE.)
C
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DWLUPH(KTO,PHOT)
C---------------------------------------------------------------------
C LORENTZ TRANSFORMATION TO CMSYSTEM AND
C UPDATING OF HEPEVT RECORD
C
C     CALLED BY : DEXAY1,(DEKAY1,DEKAY2)
C
C USED WHEN RADIATIVE CORRECTIONS IN DECAYS ARE GENERATED
C---------------------------------------------------------------------
C
      REAL  PHOT(4)
C
C CHECK ENERGY
      IF (PHOT(4).LE.0.0) RETURN
C
C POSITION OF DECAYING PARTICLE:
      IF((KTO.EQ. 1).OR.(KTO.EQ.11)) THEN
        NPS=3
      ELSE
        NPS=4
      ENDIF
C
      KTOS=KTO
      IF(KTOS.GT.10) KTOS=KTOS-10
C BOOST AND APPEND PHOTON (GAMMA IS 22)
      CALL TRALO4(KTOS,PHOT,PHOT,AM)
      CALL FILHEP(0,1,22,NPS,NPS,0,0,PHOT,0.0,.TRUE.)
C
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DWLUPI(KTO,ISGN,PPI,PNU)
C ----------------------------------------------------------------------
C LORENTZ TRANSFORMATION TO CMSYSTEM AND
C UPDATING OF HEPEVT RECORD
C
C ISGN = 1/-1 FOR TAU-/TAU+
C
C     CALLED BY : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PPI(4)
C
C POSITION OF DECAYING PARTICLE:
      IF(KTO.EQ. 1) THEN
        NPS=3
      ELSE
        NPS=4
      ENDIF
C
C TAU NEUTRINO (NU_TAU IS 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C CHARGED PI MESON (PI+ IS 211)
      CALL TRALO4(KTO,PPI,PPI,AM)
      CALL FILHEP(0,1,-211*ISGN,NPS,NPS,0,0,PPI,AM,.TRUE.)
C
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DWLURO(KTO,ISGN,PNU,PRHO,PIC,PIZ)
C ----------------------------------------------------------------------
C LORENTZ TRANSFORMATION TO CMSYSTEM AND
C UPDATING OF HEPEVT RECORD
C
C ISGN = 1/-1 FOR TAU-/TAU+
C
C     CALLED BY : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PRHO(4),PIC(4),PIZ(4)
C
C POSITION OF DECAYING PARTICLE:
      IF(KTO.EQ. 1) THEN
        NPS=3
      ELSE
        NPS=4
      ENDIF
C
C TAU NEUTRINO (NU_TAU IS 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C CHARGED RHO MESON (RHO+ IS 213)
      CALL TRALO4(KTO,PRHO,PRHO,AM)
      CALL FILHEP(0,2,-213*ISGN,NPS,NPS,0,0,PRHO,AM,.TRUE.)
C
C CHARGED PI MESON (PI+ IS 211)
      CALL TRALO4(KTO,PIC,PIC,AM)
      CALL FILHEP(0,1,-211*ISGN,-1,-1,0,0,PIC,AM,.TRUE.)
C
C PI0 MESON (PI0 IS 111)
      CALL TRALO4(KTO,PIZ,PIZ,AM)
      CALL FILHEP(0,1,111,-2,-2,0,0,PIZ,AM,.TRUE.)
C
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE DWRPH(KTO,PHX)
C
C -------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4         PHX(4)
      REAL*4 QHOT(4)
C
      DO 10 K=1,4
        QHOT(K) =0.0
   10 CONTINUE
C CASE OF TAU RADIATIVE DECAYS.
C FILLING OF THE LUND COMMON BLOCK.
      DO 20   I=1,4
   20 QHOT(I)=PHX(I)
      IF (QHOT(4).GT.1.E-5) CALL DWLUPH(KTO,QHOT)
      RETURN
      END
*CMZ :  1.02/04 13/01/97  14.47.40  by  P. Zucchelli
*CMZ :  1.02/02 12/01/97  17.44.36  by  P. Zucchelli
*CMZ :  1.01/50 19/04/96  11.22.32  by  Piero Zucchelli
*CMZ :  1.01/47 11/01/96  09.25.42  by  Piero Zucchelli
*CMZ :  1.01/45 08/01/96  14.21.55  by  Piero Zucchelli
*CMZ :  1.01/44 05/01/96  18.04.38  by  Piero Zucchelli
*CMZ :  1.01/40 09/11/95  16.12.24  by  Piero Zucchelli
*-- Author :    Piero Zucchelli   09/11/95
      SUBROUTINE EVTINFO

      CHARACTER*8 TITLE,VERSQQ
      CHARACTER*80 COOKIE
      INTEGER ITITLE(2)
*KEEP,zebra.

      PARAMETER (NNQ=1000000)
*
      DIMENSION LQ(NNQ),IQ(NNQ),Q(NNQ)
      EQUIVALENCE (Q(1),IQ(1),LQ(9),JSTRUC(8))
      COMMON /QUEST/IQUEST(100)
      COMMON /XQSTOR/IXEVT,IFENCE(16),JGEEV,JSTRUC(99),JREFER(100),
     +DIV12(NNQ)
      COMMON /FZLUN/LUNFZ
      COMMON/MZIOALL/IOGENF

*KEEP,info.
         COMMON/INFONEW/IRDATE,IRTIME
*KEEP,KEYS.
        COMMON/CFREAD/SPACE(5000)
 	COMMON/MYKEYS/IEVAR,IF5CC,INEUT,IINTE,IFERM,IFLAT,ICOUN,
     &  REFIX,IMUDO,IKAT1,IKAT2,IKAT3,IKAT4,IKAT5,IKAT6,
     &  INEVT,IJAK1,IJAK2,IITDK,RPTAU,RXK0D,NTGR,IDIMUON,IGLU,
     &  IQDEN,LOME(2),IFILES,ICCHA,ISEED,IDSUBS,IONEM,EHAC,IPSEL,
     &  IHIST


*KEEP,VIDQQ.
      CHARACTER*66 VIDQQ
      DATA VIDQQ/
     +'@(#)JETTA    1.02/13  15/01/97  23.59.31  C: 21/05/97  10.31.11'/
*KEND.
      IF (JGEEV.EQ.0) RETURN
      CALL MZBOOK(IXEVT,JGENF,JGEEV,-1,'GENF',1,1,7,IOGENF,0)
*KEEP,VERSQQ.
      VERSQQ = ' 1.02/13'
      IVERSQ =  10213
*KEEP,DATEQQ.
      IDATQQ = 970521
*KEEP,TIMEQQ.
      ITIMQQ =   1031
*KEND.
      IQ(JGENF+1)=IDATQQ
      IQ(JGENF+2)=ITIMQQ
      IQ(JGENF+3)=IVERSQ
      IQ(JGENF+4)=IRDATE
      IQ(JGENF+5)=IRTIME
      TITLE=
*KEEP,QFTITLCH,N= 8.
     +  'JETTA   '
*KEND.
      CALL UCTOH(TITLE,ITITLE,4,8)
      IQ(JGENF+6)=ITITLE(1)
      IQ(JGENF+7)=ITITLE(2)
      CALL MZBOOK(IXEVT,JGECR,JGENF,-1,'GECR',0,0,31,3,0)
      Q(JGECR+1)=ISEED
      Q(JGECR+2)=IGLU
      Q(JGECR+3)=IEVAR
      Q(JGECR+4)=IF5CC
      Q(JGECR+5)=INEUT
      Q(JGECR+6)=IINTE
      Q(JGECR+7)=IFERM
      Q(JGECR+8)=IFLAT
      Q(JGECR+9)=ICOUN
      Q(JGECR+10)=REFIX
      Q(JGECR+11)=IQDEN
      Q(JGECR+12)=IMUDO
      Q(JGECR+13)=NTGR
      Q(JGECR+14)=IDIMUON
      Q(JGECR+15)=ICCHA
      Q(JGECR+16)=LOME(1)
      Q(JGECR+17)=IFILES
      Q(JGECR+18)=IKAT1
      Q(JGECR+19)=IKAT2
      Q(JGECR+20)=IKAT3
      Q(JGECR+21)=IKAT4
      Q(JGECR+22)=IKAT5
      Q(JGECR+23)=IKAT6
      Q(JGECR+24)=INEVT
      Q(JGECR+25)=IJAK1
      Q(JGECR+26)=IJAK2
      Q(JGECR+27)=IITDK
      Q(JGECR+28)=RPTAU
      Q(JGECR+29)=RXKOD
      Q(JGECR+30)=LOME(2)
      Q(JGECR+31)=IHIST

      RETURN
      END
*CMZ :  1.01/44 15/12/95  18.07.48  by  Piero Zucchelli
*CMZ :  1.01/43 15/12/95  18.01.58  by  Piero Zucchelli
*CMZ :  1.01/22 27/05/95  16.23.58  BY  PIERO ZUCCHELLI
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.35.13  BY  PIERO ZUCCHELLI
*CMZ :  1.01/03 05/03/95  00.19.19  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 25/07/94  14.29.34  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 19/07/94  15.47.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE FERMII(F)
*KEEP,KEYS.
        COMMON/CFREAD/SPACE(5000)
 	COMMON/MYKEYS/IEVAR,IF5CC,INEUT,IINTE,IFERM,IFLAT,ICOUN,
     &  REFIX,IMUDO,IKAT1,IKAT2,IKAT3,IKAT4,IKAT5,IKAT6,
     &  INEVT,IJAK1,IJAK2,IITDK,RPTAU,RXK0D,NTGR,IDIMUON,IGLU,
     &  IQDEN,LOME(2),IFILES,ICCHA,ISEED,IDSUBS,IONEM,EHAC,IPSEL,
     &  IHIST


*KEND.
C.
C.
C.
      DIMENSION  F(3)
      REAL*4 F,ECI,PTR,CS,SN,PHI,FERM,FNUC,TMPX,TMPY,THE
C.
      DATA FERM/.027/
      DATA FNUC/.939/

   10 CONTINUE
      TMPX=RNDMM(ISEED)*FERM
      TMPY=RNDMM(ISEED)*SQRT(FERM)

      IF (SQRT(TMPX).GT.TMPY) THEN
        ECI=TMPX
      ELSE
        GOTO 10
      ENDIF




      PTR = SQRT(ECI*FNUC*2.)

      CS  = 2.*RNDMM(ISEED)-1.
      SN  = SQRT(1.-CS**2)
      PHI = 6.2832*RNDMM(ISEED)


      F(1) = PTR*SN*COS(PHI)
      F(2) = PTR*SN*SIN(PHI)
      F(3) = PTR*CS
      RETURN
      END
*CMZ :  1.02/03 13/01/97  13.40.31  by  P. Zucchelli
*CMZ :  1.01/50 23/05/96  10.19.16  by  Piero Zucchelli
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE FILHEP(N,IST,ID,JMO1,JMO2,JDA1,JDA2,P4,PINV,PHFLAG)
C ----------------------------------------------------------------------
C THIS SUBROUTINE FILLS ONE ENTRY INTO THE HEPEVT COMMON
C AND UPDATES THE INFORMATION FOR AFFECTED MOTHER ENTRIES
C
C WRITTEN BY MARTIN W. GRUENEWALD (91/01/28)
C
C     CALLED BY : ZTOHEP,BTOHEP,DWLUXY
C ----------------------------------------------------------------------
C
      PARAMETER (NMXHEP=2000)
*KEEP,HEPEVT.
      DOUBLE PRECISION PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      SAVE /HEPEVT/

*KEND.
      COMMON/PHOQED/QEDRAD(NMXHEP)
      LOGICAL QEDRAD
      SAVE /PHOQED/
      LOGICAL PHFLAG
C
      REAL*4  P4(4)
C
C CHECK ADDRESS MODE
      IF (N.EQ.0) THEN
C
C APPEND MODE
        IHEP=NHEP+1
      ELSE IF (N.GT.0) THEN
C
C ABSOLUTE POSITION
        IHEP=N
      ELSE
C
C RELATIVE POSITION
        IHEP=NHEP+N
      END IF
C
C CHECK ON IHEP
      IF ((IHEP.LE.0).OR.(IHEP.GT.NMXHEP)) RETURN
C
C ADD ENTRY
      NHEP=IHEP
      ISTHEP(IHEP)=IST
      IDHEP(IHEP)=ID
      JMOHEP(1,IHEP)=JMO1
      IF(JMO1.LT.0)JMOHEP(1,IHEP)=JMOHEP(1,IHEP)+IHEP
      JMOHEP(2,IHEP)=JMO2
      IF(JMO2.LT.0)JMOHEP(2,IHEP)=JMOHEP(2,IHEP)+IHEP
      JDAHEP(1,IHEP)=JDA1
      JDAHEP(2,IHEP)=JDA2
C
      DO I=1,4
        PHEP(I,IHEP)=P4(I)
C
C KORAL-B AND KORAL-Z DO NOT PROVIDE VERTEX AND/OR LIFETIME INFORMATIONS
        VHEP(I,IHEP)=0.0
      END DO
      PHEP(5,IHEP)=PINV
C FLAG FOR PHOTOS...
      QEDRAD(IHEP)=PHFLAG
C
C UPDATE PROCESS:
      DO IP=JMOHEP(1,IHEP),JMOHEP(2,IHEP)
        IF(IP.GT.0)THEN
C
C IF THERE IS A DAUGHTER AT IHEP, MOTHER ENTRY AT IP HAS DECAYED
          IF(ISTHEP(IP).EQ.1)ISTHEP(IP)=2
C
C AND DAUGHTER POINTERS OF MOTHER ENTRY MUST BE UPDATED
          IF(JDAHEP(1,IP).EQ.0)THEN
            JDAHEP(1,IP)=IHEP
            JDAHEP(2,IP)=IHEP
          ELSE
            JDAHEP(2,IP)=MAX(IHEP,JDAHEP(2,IP))
          END IF
        END IF
      END DO
C
      RETURN
      END
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      FUNCTION FLGINT(Z)

C...GLUON CONTRIBUTION INTEGRAND TO QCD LONGITUDINAL STRUCTURE FUNCTION.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTEG/ NTOT,NPASS
      DIMENSION XPQ(-6:6)
      DATA PI/3.14159/
      NTOT=NTOT+1
      CALL LNSTRF(Z,Q2,XPQ)
      FLGINT=20./9.*PARL(25)/PI*(X/Z)**2*(1.-X/Z)/Z*XPQ(0)
      NPASS=NPASS+1

      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE FLINTG(CFLQ,CFLG,CFLM)

C...EVENT-BY-EVENT CALCULATION OF CONTRIBUTION TO LONGITUDINAL
C...STRUCTURE FUNCTION FROM QCD AND TARGET MASS EFFECTS.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTEG/ NTOT,NPASS
      EXTERNAL FLQINT,FLGINT,FLTINT

      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      LHT=LST(11)/100
      PARL(25)=ULALPS(Q2)
      IF(LQCD.EQ.2) THEN
C...FL FROM QCD, QUARK AND GLUON CONTRIBUTIONS.
        ACCUR=PARL(11)
        IT=0
   10   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLQINT,EPS,CFLQ)
        IF(CFLQ.LT.1) THEN
          ACCUR=CFLQ*PARL(11)
          IF(IT.LT.2) GOTO 10
        ENDIF
        ACCUR=PARL(11)
        IT=0
   20   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLGINT,EPS,CFLG)
        IF(CFLG.LT.1.) THEN
          ACCUR=CFLG*PARL(11)
          IF(IT.LT.2) GOTO 20
        ENDIF
      ENDIF
      IF(LTM.EQ.2) THEN
        ACCUR=PARL(11)
        IT=0
   30   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLTINT,EPS,CFLM)
        IF(CFLM.LT.1.) THEN
          ACCUR=CFLM*PARL(11)
          IF(IT.LT.2) GOTO 30
        ENDIF
      ENDIF

      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE FLIPOL(FLQ,FLG,FLM)

C...QCD AND TARGET MASS CONTRIBUTIONS TO LONGITUDINAL STRUCTURE FUNCTION
C...FROM INTERPOLATION ON X,Q2 GRID.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /FLGRID/ NFX,NFQ,XR(2),QR(2),FLQT(41,16),FLGT(41,16),
     &FLMT(41,16)
      DATA NOUT/0/,NWARN/10/

      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      LHT=LST(11)/100
      XP=X
      Q2P=Q2
C...NOTE: TINY MISMATCH BETWEEN PRESENT X-VALUE AND THOSE ON GRID.
      QR(2)=X*PARL(21)
      IF(QR(1).GT.QR(2)) RETURN
      IF(X.LT.XR(1).OR.X.GT.XR(2).OR.
     +Q2.LT.QR(1).OR.Q2.GT.QR(2)) THEN
C...X AND/OR Q2 OUTSIDE GRID LIMITS, WRITE WARNING FOR FIRST NWARN CASES
        IF(LST(2).GE.0) THEN
          NOUT=NOUT+1
          IF(LST(3).GE.1.AND.NOUT.LE.NWARN) WRITE(6,10000) X,Q2,NWARN
        ENDIF
        IF(X.LT.XR(1)) XP=XR(1)
        IF(X.GT.XR(2)) XP=XR(2)
        IF(Q2.LT.QR(1)) Q2P=QR(1)
        IF(Q2.GT.QR(2)) Q2P=QR(2)
      ENDIF

      IX=(ALOG10(XP)-ALOG10(XR(1)))/
     &(ALOG10(XR(2))-ALOG10(XR(1)))*(NFX-1)+1
      IQ=(ALOG10(Q2P)-ALOG10(QR(1)))/
     &(ALOG10(QR(2))-ALOG10(QR(1)))*(NFQ-1)+1
      IX=MIN(IX,NFX-1)
      IQ=MIN(IQ,NFQ-1)
      Q2L=10**(ALOG10(QR(1))+(ALOG10(QR(2))-ALOG10(QR(1)))*
     &(IQ-1)/(NFQ-1))
      Q2H=10**(ALOG10(QR(1))+(ALOG10(QR(2))-ALOG10(QR(1)))*
     &(IQ  )/(NFQ-1))
      XL=10**(ALOG10(XR(1))+(ALOG10(XR(2))-ALOG10(XR(1)))*
     &(IX-1)/(NFX-1))
      XH=10**(ALOG10(XR(1))+(ALOG10(XR(2))-ALOG10(XR(1)))*
     &(IX  )/(NFX-1))
      QD=(Q2P-Q2L)/(Q2H-Q2L)
      XD=(XP-XL)/(XH-XL)

      IF(LQCD.EQ.1) THEN
        X1P=(FLQT(IX+1,IQ)-FLQT(IX,IQ))*XD+FLQT(IX,IQ)
        X2P=(FLQT(IX+1,IQ+1)-FLQT(IX,IQ+1))*XD+FLQT(IX,IQ+1)
        FLQ=(X2P-X1P)*QD+X1P
        X1P=(FLGT(IX+1,IQ)-FLGT(IX,IQ))*XD+FLGT(IX,IQ)
        X2P=(FLGT(IX+1,IQ+1)-FLGT(IX,IQ+1))*XD+FLGT(IX,IQ+1)
        FLG=(X2P-X1P)*QD+X1P
      ENDIF
      IF(LTM.EQ.1) THEN
        X1P=(FLMT(IX+1,IQ)-FLMT(IX,IQ))*XD+FLMT(IX,IQ)
        X2P=(FLMT(IX+1,IQ+1)-FLMT(IX,IQ+1))*XD+FLMT(IX,IQ+1)
        FLM=(X2P-X1P)*QD+X1P
      ENDIF

      RETURN
10000 FORMAT(' WARNING: X=',F7.4,' OR Q2=',F6.1,' OUTSIDE GRID,',
     &' FOR FL INTERPOLATION',/,10X,'VALUE ON GRID LIMIT USED.',
     &' ONLY FIRST',I5,' WARNINGS PRINTED.',/)
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      FUNCTION FLQINT(Z)

C...QUARK CONTRIBUTION INTEGRAND TO QCD LONGITUDINAL STRUCTURE FUNCTION.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTEG/ NTOT,NPASS
      DIMENSION XPQ(-6:6)
      DATA PI/3.14159/
      NTOT=NTOT+1
      CALL LNSTRF(Z,Q2,XPQ)
      FLQINT=0.
      DO 10  I=-LST(12),LST(12)
        IF(I.EQ.0) GOTO 10
        FLQINT=FLQINT+QC(IABS(I))**2*XPQ(I)
   10 CONTINUE
      FLQINT=4./3.*PARL(25)/PI*(X/Z)**2*FLQINT/Z
      NPASS=NPASS+1

      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE FLTABL

C...INTEGRATES THE LONGITUDINAL STRUCTURE FUNCTION, STORE ON GRID
C...IN X, Q**2.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LINTEG/ NTOT,NPASS
      COMMON /FLGRID/ NFX,NFQ,XR(2),QR(2),FLQT(41,16),FLGT(41,16),
     +FLMT(41,16)
      EXTERNAL FLQINT,FLGINT,FLTINT

      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      LHT=LST(11)/100
      IF(LST(3).GE.3) WRITE(6,10000) LST(11),LQCD,LTM,LHT
      IF(LQCD.LT.1.AND.LTM.LT.1) RETURN
      CALL LTIMEX(T1)
      DO 10 IX=1,NFX
        DO 10 IQ=1,NFQ
          FLQT(IX,IQ)=0.
          FLGT(IX,IQ)=0.
   10 FLMT(IX,IQ)=0.
      QR(1)=Q2MIN
      XR(1)=XMIN
      XR(2)=XMAX
      DO 60  IX=1,NFX
        X=10**(ALOG10(XR(1))+(ALOG10(XR(2))-ALOG10(XR(1)))*(IX-1)/(NFX-
     +  1))
        QR(2)=X*PARL(21)
        IF(QR(1).GT.QR(2)) GOTO 60
        LQ=0
        DO 50  IQ=1,NFQ
          Q2=10**(ALOG10(QR(1))+(ALOG10(QR(2))-ALOG10(QR(1)))* (IQ-1)/
     +    (NFQ-1))
CTEST IF(LQ.GT.0) GOTO 500
          IF(Q2.GT.PARL(21)) LQ=LQ+1
          Y=Q2/(PARL(21)*X)
          IF(Y.LT.0.0.OR.Y.GT.1.0) LQ=LQ+1
          PARL(25)=ULALPS(Q2)
          IF(LQCD.EQ.1) THEN
C...QUARK PART.
            ACCUR=PARL(11)
            IT=0
   20       IT=IT+1
            NTOT=0
            NPASS=0
            EPS=ACCUR
            CALL GADAP(X,1.,FLQINT,EPS,FLQ)
            IF(FLQ.LT.1) THEN
              ACCUR=FLQ*PARL(11)
              IF(IT.LT.2) GOTO 20
            ENDIF
            FLQT(IX,IQ)=FLQ
C...GLUON PART.
            ACCUR=PARL(11)
            IT=0
   30       IT=IT+1
            NTOT=0
            NPASS=0
            EPS=ACCUR
            CALL GADAP(X,1.,FLGINT,EPS,FLG)
            IF(FLG.LT.1.) THEN
              ACCUR=FLG*PARL(11)
              IF(IT.LT.2) GOTO 30
            ENDIF
            FLGT(IX,IQ)=FLG
          ENDIF
          IF(LTM.EQ.1) THEN
C...TARGET MASS  PART.
            ACCUR=PARL(11)
            IT=0
   40       IT=IT+1
            NTOT=0
            NPASS=0
            EPS=ACCUR
            CALL GADAP(X,1.,FLTINT,EPS,FLM)
            IF(FLM.LT.1) THEN
              ACCUR=FLM*PARL(11)
              IF(IT.LT.2) GOTO 40
            ENDIF
            FLMT(IX,IQ)=FLM
          ENDIF
   50   CONTINUE
   60 CONTINUE
   70 CONTINUE
      CALL LTIMEX(T2)
      IF(LST(3).GE.3) WRITE(6,10100) T2-T1
      RETURN

10000 FORMAT(' INITIALISATION FOR FL; QCD, TARGET MASS, HIGHER TWIST: ',
     +/,' LST(11) =',I5,' --> LQCD, LTM, LHT =',3I3)
10100 FORMAT(' FL INTEGRATIONS PERFORMED IF LQCD=1 AND/OR LTM=1, ',
     +'RESULTS ON GRID.'/,' TIME FOR FL INTEGRATIONS IS ',F7.1,' SEC.')
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      FUNCTION FLTINT(Z)

C...INTEGRAND FOR TARGET MASS CORRECTION CONTRIBUTION TO
C...QUARK LONGITUDINAL STRUCTURE FUNCTION

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTEG/ NTOT,NPASS
      DIMENSION XPQ(-6:6)
      DATA PM2/0.8804/
      NTOT=NTOT+1
      CALL LNSTRF(Z,Q2,XPQ)
      FLTINT=0.
      DO 10  I=-LST(12),LST(12)
        IF(I.EQ.0) GOTO 10
        FLTINT=FLTINT+QC(IABS(I))**2*XPQ(I)
   10 CONTINUE
      FLTINT=4.*PM2/Q2*(X/Z)**2*X*FLTINT
      NPASS=NPASS+1

      RETURN
      END
*CMZ :  1.02/01 12/01/97  16.42.39  by  J. Brunner
*CMZ :  1.01/50 23/05/96  12.34.50  by  Piero Zucchelli
*-- Author :
      SUBROUTINE FORCED_DECAY(NUFORCE,ISTATUS)

      COMMON/NTUPL10/      NUTYPE,IPARENT,EPARENT,XDECAY,YDECAY,ZDECAY,
     +                     PXPAR,PYPAR,PZPAR,XDET,YDET,XL,PXNU,PYNU,PZNU,
     +                     NPROT

      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000), KFDP(2000,5)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      REAL BRATKP(10)
      REAL BRATK0(10)
      DATA ZDETC/82342.0/
      DATA DETX/250./
      DATA DETY/250./

      ISTATUS=3
      IF (NUFORCE.LE.50) THEN

        R1=0
        R2=0

* R1= other K CHARGED decays Br's
        DO IO=361,366
          BRATKP(IO-360)=BRAT(IO)
        ENDDO
        DO IO=361,364
          MDME(IO,1)=0
          R1=R1+BRAT(IO)
        ENDDO
        MDME(366,1)=0
        R1=R1+BRAT(366)

* R2= other K NEUTRAL decays Br's
        DO IO=958,965
          BRATK0(IO-957)=BRAT(IO)
        ENDDO
        DO IO=958,959
          MDME(IO,1)=0
          R2=R2+BRAT(IO)
        ENDDO
        DO IO=962,965
          MDME(IO,1)=0
          R2=R2+BRAT(IO)
        ENDDO

* k+ nu e decays branching ratio
        EDECKP=1-R1
        EDECK0=1-R2
        EMA=MAX(EDECKP,EDECK0)
        E1=EDECKP/EMA
        E2=EDECK0/EMA
* now restore BRAT for charged kaon
        BRAT(361)=1-E1
        MDME(361,1)=1
        DO IO=365,365
          BRAT(IO)=BRAT(IO)*E1/EDECKP
        ENDDO
        BRAT(958)=1-E2
        MDME(958,1)=1
        DO IO=960,961
          BRAT(IO)=BRAT(IO)*E2/EDECK0
        ENDDO
        XFACT=1/EMA


      ENDIF

      IF (ISFIRST.EQ.0) THEN
        ISFIRST=1
        WRITE(*,*)' +++IMPORTANT!! GBEAM ENHANCED MODE:'
        WRITE(*,*)' +++STATISTICAL AMPLIFICATION OF'
        WRITE(*,*)' +++NEUTRINOS OF TYPE ',NUFORCE
        WRITE(*,*)' +++BY FACTOR ',XFACT
        CALL LULIST(12)
      ENDIF


      IF (IPARENT.EQ.9) THEN
        LPARENT=-211
      ELSEIF(IPARENT.EQ.8) THEN
        LPARENT=211
      ELSEIF(IPARENT.EQ.11) THEN
        LPARENT=321
      ELSEIF(IPARENT.EQ.12) THEN
        LPARENT=-321
      ELSEIF(IPARENT.EQ.5) THEN
        LPARENT=-13
      ELSEIF(IPARENT.EQ.6) THEN
        LPARENT=13
      ELSEIF(IPARENT.EQ.10) THEN
        LPARENT=130
      ELSE
        WRITE(*,*)'WARNING: UNKNOWN GEANT PARENT ID=',IPARENT
        ISTATUS=1
      ENDIF

      IF (NUFORCE.EQ.51) THEN
        LUFORCE=14
      ELSEIF(NUFORCE.EQ.52) THEN
        LUFORCE=-14
      ELSEIF(NUFORCE.EQ.49) THEN
        LUFORCE=12
      ELSEIF(NUFORCE.EQ.50) THEN
        LUFORCE=-12
      ENDIF


   10 CONTINUE
      N=1
      K(1,1)=5
      K(1,2)=LPARENT

      P(1,1)=PXPAR
      P(1,2)=PYPAR
      P(1,3)=PZPAR
      P(1,5)=ULMASS(LPARENT)
      P(1,4)=SQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)

*      WRITE(*,*)'BEFORE:'
*      CALL LULIST(3)
      CALL LUDECY(1)
*      WRITE(*,*)'AFTER:'
*      CALL LULIST(3)

      DO I=2,N
        IF (K(I,2).EQ.LUFORCE) THEN
          TMPXL=ZDETC-ZDECAY
          TMPXDET=TMPXL*P(I,1)/P(I,3)+XDECAY
          TMPYDET=TMPXL*P(I,2)/P(I,3)+YDECAY
* following line is wrong!!! for debugging
          IF (ABS(TMPXDET).GT.DETX.OR.ABS(TMPYDET).GT.DETY) GOTO 10
          XL=TMPXL
          XDET=TMPXDET
          YDET=TMPYDET
          PXNU=P(I,1)
          PYNU=P(I,2)
          PZNU=P(I,3)
          NUTYPE=NUFORCE
          ISTATUS=0
        ENDIF
      ENDDO

      DO IO=958,965
        MDME(IO,1)=1
        BRAT(IO)=BRATK0(IO-957)
      ENDDO
      DO IO=361,366
        MDME(IO,1)=1
        BRAT(IO)=BRATKP(IO-360)
      ENDDO

      RETURN
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION FORM1(MNUM,QQ,S1,SDWA)
C     ==================================================================
C     FORMFACTORFOR F1 FOR 3 SCALAR FINAL STATE
C     R. FISHER, J. WESS AND F. WAGNER Z. PHYS C3 (1980) 313
C     H. GEORGI, WEAK INTERACTIONS AND MODERN PARTICLE THEORY,
C     THE BENJAMIN/CUMMINGS PUB. CO., INC. 1984.
C     R. DECKER, E. MIRKES, R. SAUER, Z. WAS KARLSRUHE PREPRINT TTP92-25
C     AND ERRATUM !!!!!!
C     ==================================================================
C
      COMPLEX FORM1,WIGNER,WIGFOR,FPIKM,BWIGM
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      WIGNER(A,B,C)= CMPLX(1.0,0.0)/CMPLX(A-B**2,B*C)
      IF     (MNUM.EQ.0) THEN
C ------------  3 PI HADRONIC STATE (A1)
        GAMAX=GAMA1*GFUN(QQ)/GFUN(AMA1**2)
        FORM1=AMA1**2*WIGNER(QQ,AMA1,GAMAX)*FPIKM(SQRT(S1),AMPI,AMPI)
      ELSEIF (MNUM.EQ.1) THEN
C ------------ K- PI- K+
        FORM1=BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
        GAMAX=GAMA1*GFUN(QQ)/GFUN(AMA1**2)
        FORM1=AMA1**2*WIGNER(QQ,AMA1,GAMAX)*FORM1
      ELSEIF (MNUM.EQ.2) THEN
C ------------ K0 PI- K0B
        FORM1=BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
        GAMAX=GAMA1*GFUN(QQ)/GFUN(AMA1**2)
        FORM1=AMA1**2*WIGNER(QQ,AMA1,GAMAX)*FORM1
      ELSEIF (MNUM.EQ.3) THEN
C ------------ K- K0 PI0
        FORM1=0.0
        GAMAX=GAMA1*GFUN(QQ)/GFUN(AMA1**2)
        FORM1=AMA1**2*WIGNER(QQ,AMA1,GAMAX)*FORM1
      ELSEIF (MNUM.EQ.4) THEN
C ------------ PI0 PI0 K-
        XM2=1.402
        GAM2=0.174
        FORM1=BWIGM(S1,AMKST,GAMKST,AMK,AMPI)
        FORM1=WIGFOR(QQ,XM2,GAM2)*FORM1
      ELSEIF (MNUM.EQ.5) THEN
C ------------ K- PI- PI+
        XM2=1.402
        GAM2=0.174
        FORM1=WIGFOR(QQ,XM2,GAM2)*FPIKM(SQRT(S1),AMPI,AMPI)
      ELSEIF (MNUM.EQ.6) THEN
        FORM1=0.0
      ELSEIF (MNUM.EQ.7) THEN
C -------------- ETA PI- PI0 FINAL STATE
        FORM1=0.0
      ENDIF
C
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION FORM2(MNUM,QQ,S1,SDWA)
C     ==================================================================
C     FORMFACTORFOR F2 FOR 3 SCALAR FINAL STATE
C     R. FISHER, J. WESS AND F. WAGNER Z. PHYS C3 (1980) 313
C     H. GEORGI, WEAK INTERACTIONS AND MODERN PARTICLE THEORY,
C     THE BENJAMIN/CUMMINGS PUB. CO., INC. 1984.
C     R. DECKER, E. MIRKES, R. SAUER, Z. WAS KARLSRUHE PREPRINT TTP92-25
C     AND ERRATUM !!!!!!
C     ==================================================================
C
      COMPLEX FORM2,WIGNER,WIGFOR,FPIKM,BWIGM
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      WIGNER(A,B,C)= CMPLX(1.0,0.0)/CMPLX(A-B**2,B*C)
      IF     (MNUM.EQ.0) THEN
C ------------  3 PI HADRONIC STATE (A1)
        GAMAX=GAMA1*GFUN(QQ)/GFUN(AMA1**2)
        FORM2=AMA1**2*WIGNER(QQ,AMA1,GAMAX)*FPIKM(SQRT(S1),AMPI,AMPI)
      ELSEIF (MNUM.EQ.1) THEN
C ------------ K- PI- K+
        GAMAX=GAMA1*GFUN(QQ)/GFUN(AMA1**2)
        FORM2=AMA1**2*WIGNER(QQ,AMA1,GAMAX)*FPIKM(SQRT(S1),AMPI,AMPI)
      ELSEIF (MNUM.EQ.2) THEN
C ------------ K0 PI- K0B
        GAMAX=GAMA1*GFUN(QQ)/GFUN(AMA1**2)
        FORM2=AMA1**2*WIGNER(QQ,AMA1,GAMAX)*FPIKM(SQRT(S1),AMPI,AMPI)
      ELSEIF (MNUM.EQ.3) THEN
C ------------ K- K0 PI0
        GAMAX=GAMA1*GFUN(QQ)/GFUN(AMA1**2)
        FORM2=AMA1**2*WIGNER(QQ,AMA1,GAMAX)*FPIKM(SQRT(S1),AMPI,AMPI)
      ELSEIF (MNUM.EQ.4) THEN
C ------------ PI0 PI0 K-
        XM2=1.402
        GAM2=0.174
        FORM2=BWIGM(S1,AMKST,GAMKST,AMK,AMPI)
        FORM2=WIGFOR(QQ,XM2,GAM2)*FORM2
      ELSEIF (MNUM.EQ.5) THEN
C ------------ K- PI- PI+
        XM2=1.402
        GAM2=0.174
        FORM2=BWIGM(S1,AMKST,GAMKST,AMK,AMPI)
        FORM2=WIGFOR(QQ,XM2,GAM2)*FORM2
C
      ELSEIF (MNUM.EQ.6) THEN
        XM2=1.402
        GAM2=0.174
        FORM2=WIGFOR(QQ,XM2,GAM2)*FPIKM(SQRT(S1),AMPI,AMPI)
C
      ELSEIF (MNUM.EQ.7) THEN
C -------------- ETA PI- PI0 FINAL STATE
        FORM2=0.0
      ENDIF
C
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION FORM3(MNUM,QQ,S1,SDWA)
C     ==================================================================
C     FORMFACTORFOR F3 FOR 3 SCALAR FINAL STATE
C     R. FISHER, J. WESS AND F. WAGNER Z. PHYS C3 (1980) 313
C     H. GEORGI, WEAK INTERACTIONS AND MODERN PARTICLE THEORY,
C     THE BENJAMIN/CUMMINGS PUB. CO., INC. 1984.
C     R. DECKER, E. MIRKES, R. SAUER, Z. WAS KARLSRUHE PREPRINT TTP92-25
C     AND ERRATUM !!!!!!
C     ==================================================================
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMPLEX FORM3
      IF (MNUM.EQ.6) THEN
        FORM3=CMPLX(0.0)
      ELSE
        FORM3=CMPLX(0.0)
      ENDIF
      FORM3=0
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION FORM4(MNUM,QQ,S1,S2,S3)
C     ==================================================================
C     FORMFACTORFOR F4 FOR 3 SCALAR FINAL STATE
C     R. DECKER, IN PREPARATION
C     R. DECKER, E. MIRKES, R. SAUER, Z. WAS KARLSRUHE PREPRINT TTP92-25
C     AND ERRATUM !!!!!!
C     ==================================================================
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMPLEX FORM4,WIGNER,FPIKM
      REAL*4 M
      WIGNER(A,B,C)=CMPLX(1.0,0.0) /CMPLX(A-B**2,B*C)
      IF (MNUM.EQ.0) THEN
C ------------  3 PI HADRONIC STATE (A1)
        G1=5.8
        G2=6.08
        FPIP=0.02
        AMPIP=1.3
        GAMPIP=0.3
        S=QQ
        G=GAMPIP
        XM1=AMPIZ
        XM2=AMRO
        M  =AMPIP
        IF (S.GT.(XM1+XM2)**2) THEN
          QS=SQRT(ABS((S -(XM1+XM2)**2)*(S -(XM1-XM2)**2)))/SQRT(S)
          QM=SQRT(ABS((M**2-(XM1+XM2)**2)*(M**2-(XM1-XM2)**2)))/M
          W=SQRT(S)
          GS=G*(M/W)**2*(QS/QM)**5
        ELSE
          GS=0.0
        ENDIF
        GAMX=GS*W/M
        FORM4=G1*G2*FPIP/AMRO**4/AMPIP**2
     +       *AMPIP**2*WIGNER(QQ,AMPIP,GAMX)
     +       *( S1*(S2-S3)*FPIKM(SQRT(S1),AMPIZ,AMPIZ)
     +         +S2*(S1-S3)*FPIKM(SQRT(S2),AMPIZ,AMPIZ) )
      ELSEIF (MNUM.EQ.1) THEN
C ------------  3 PI HADRONIC STATE (A1)
        G1=5.8
        G2=6.08
        FPIP=0.02
        AMPIP=1.3
        GAMPIP=0.3
        S=QQ
        G=GAMPIP
        XM1=AMPIZ
        XM2=AMRO
        M  =AMPIP
        IF (S.GT.(XM1+XM2)**2) THEN
          QS=SQRT(ABS((S -(XM1+XM2)**2)*(S -(XM1-XM2)**2)))/SQRT(S)
          QM=SQRT(ABS((M**2-(XM1+XM2)**2)*(M**2-(XM1-XM2)**2)))/M
          W=SQRT(S)
          GS=G*(M/W)**2*(QS/QM)**5
        ELSE
          GS=0.0
        ENDIF
        GAMX=GS*W/M
        FORM4=G1*G2*FPIP/AMRO**4/AMPIP**2
     +       *AMPIP**2*WIGNER(QQ,AMPIP,GAMX)
     +       *( S1*(S2-S3)*FPIKM(SQRT(S1),AMPIZ,AMPIZ)
     +         +S2*(S1-S3)*FPIKM(SQRT(S2),AMPIZ,AMPIZ) )
      ELSE
        FORM4=CMPLX(0.0,0.0)
      ENDIF
C ---- THIS FORMFACTOR IS SWITCHED OFF .. .
      FORM4=CMPLX(0.0,0.0)
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION FORM5(MNUM,QQ,S1,S2)
C     ==================================================================
C     FORMFACTORFOR F5 FOR 3 SCALAR FINAL STATE
C     G. KRAMER, W. PALMER, S. PINSKY, PHYS. REV. D30 (1984) 89.
C     G. KRAMER, W. PALMER             Z. PHYS. C25 (1984) 195.
C     R. DECKER, E. MIRKES, R. SAUER, Z. WAS KARLSRUHE PREPRINT TTP92-25
C     AND ERRATUM !!!!!!
C     ==================================================================
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMPLEX FORM5,WIGNER,FPIKM,FPIKMD,BWIGM
      WIGNER(A,B,C)=CMPLX(1.0,0.0)/CMPLX(A-B**2,B*C)
      IF     (MNUM.EQ.0) THEN
C ------------  3 PI HADRONIC STATE (A1)
        FORM5=0.0
      ELSEIF (MNUM.EQ.1) THEN
C ------------ K- PI- K+
        ELPHA=-0.2
        FORM5=FPIKMD(SQRT(QQ),AMPI,AMPI)/(1+ELPHA) *( FPIKM(SQRT(S2),
     +  AMPI,AMPI) +ELPHA*BWIGM(S1,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.2) THEN
C ------------ K0 PI- K0B
        ELPHA=-0.2
        FORM5=FPIKMD(SQRT(QQ),AMPI,AMPI)/(1+ELPHA) *( FPIKM(SQRT(S2),
     +  AMPI,AMPI) +ELPHA*BWIGM(S1,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.3) THEN
C ------------ K- K0 PI0
        FORM5=0.0
      ELSEIF (MNUM.EQ.4) THEN
C ------------ PI0 PI0 K-
        FORM5=0.0
      ELSEIF (MNUM.EQ.5) THEN
C ------------ K- PI- PI+
        ELPHA=-0.2
        FORM5=BWIGM(QQ,AMKST,GAMKST,AMPI,AMK)/(1+ELPHA)
     +       *(       FPIKM(SQRT(S1),AMPI,AMPI)
     +         +ELPHA*BWIGM(S2,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.6) THEN
C ------------ PI- K0B PI0
        ELPHA=-0.2
        FORM5=BWIGM(QQ,AMKST,GAMKST,AMPI,AMKZ)/(1+ELPHA)
     +       *(       FPIKM(SQRT(S2),AMPI,AMPI)
     +         +ELPHA*BWIGM(S1,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.7) THEN
C -------------- ETA PI- PI0 FINAL STATE
        FORM5=FPIKMD(SQRT(QQ),AMPI,AMPI)*FPIKM(SQRT(S1),AMPI,AMPI)
      ENDIF
C
      END
*CMZ :  1.00/00 10/08/94  16.29.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION FORMOM(XMAA,XMOM)
C     ==================================================================
C     FORMFACTORFOR PI-PI0 GAMMA FINAL STATE
C      R. DECKER, Z. PHYS C36 (1987) 487.
C     ==================================================================
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON /TESTA1/ KEYA1
      COMPLEX BWIGN,FORMOM
      DATA ICONT /1/
* THIS INLINE FUNCT. CALCULATES THE SCALAR PART OF THE PROPAGATOR
      BWIGN(XM,AM,GAMMA)=1./CMPLX(XM**2-AM**2,GAMMA*AM)
* HADRON CURRENT
      FRO  =0.266*AMRO**2
      ELPHA=- 0.1
      AMROP = 1.7
      GAMROP= 0.26
      AMOM  =0.782
      GAMOM =0.0085
      AROMEG= 1.0
      GCOUP=12.924
      GCOUP=GCOUP*AROMEG
      FQED  =SQRT(4.0*3.1415926535/137.03604)
      FORMOM=FQED*FRO**2/SQRT(2.0)*GCOUP**2*BWIGN(XMOM,AMOM,GAMOM)
     $     *(BWIGN(XMAA,AMRO,GAMRO)+ELPHA*BWIGN(XMAA,AMROP,GAMROP))
     $     *(BWIGN( 0.0,AMRO,GAMRO)+ELPHA*BWIGN( 0.0,AMROP,GAMROP))
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      COMPLEX FUNCTION FPIK(W)
C **********************************************************
C     PION FORM FACTOR
C **********************************************************
      COMPLEX BWIG
      REAL ROM,ROG,ROM1,ROG1,BETA1,PI,PIM,S,W
      EXTERNAL BWIG
      DATA  INIT /0/
C
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0 ) THEN
        INIT=1
        PI=3.141592654
        PIM=.140
        ROM=0.773
        ROG=0.145
        ROM1=1.370
        ROG1=0.510
        BETA1=-0.145
      ENDIF
C -----------------------------------------------
      S=W**2
      FPIK= (BWIG(S,ROM,ROG)+BETA1*BWIG(S,ROM1,ROG1))
     + /(1+BETA1)
      RETURN
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      COMPLEX FUNCTION FPIKM(W,XM1,XM2)
C **********************************************************
C     PION FORM FACTOR
C **********************************************************
      COMPLEX BWIGM
      REAL ROM,ROG,ROM1,ROG1,BETA1,PI,PIM,S,W
      EXTERNAL BWIG
      DATA  INIT /0/
C
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0 ) THEN
        INIT=1
        PI=3.141592654
        PIM=.140
        ROM=0.773
        ROG=0.145
        ROM1=1.370
        ROG1=0.510
        BETA1=-0.145
      ENDIF
C -----------------------------------------------
      S=W**2
      FPIKM=(BWIGM(S,ROM,ROG,XM1,XM2)+BETA1*BWIGM(S,ROM1,ROG1,XM1,XM2))
     + /(1+BETA1)
      RETURN
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      COMPLEX FUNCTION FPIKMD(W,XM1,XM2)
C **********************************************************
C     PION FORM FACTOR
C **********************************************************
      COMPLEX BWIGM
      REAL ROM,ROG,ROM1,ROG1,PI,PIM,S,W
      EXTERNAL BWIG
      DATA  INIT /0/
C
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0 ) THEN
        INIT=1
        PI=3.141592654
        PIM=.140
        ROM=0.773
        ROG=0.145
        ROM1=1.500
        ROG1=0.220
        ROM2=1.750
        ROG2=0.120
        BETA=6.5
        DELTA=-26.0
      ENDIF
C -----------------------------------------------
      S=W**2
      FPIKMD=(DELTA*BWIGM(S,ROM,ROG,XM1,XM2)
     +      +BETA*BWIGM(S,ROM1,ROG1,XM1,XM2)
     +      +     BWIGM(S,ROM2,ROG2,XM1,XM2))
     + /(1+BETA+DELTA)
      RETURN
      END
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      COMPLEX FUNCTION FPIKMK(W,XM1,XM2)
C **********************************************************
C     KAON FORM FACTOR
C **********************************************************
      COMPLEX BWIGM
      REAL ROM,ROG,ROM1,ROG1,BETA1,PI,PIM,S,W
      EXTERNAL BWIG
      DATA  INIT /0/
C
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0 ) THEN
        INIT=1
        PI=3.141592654
        PIM=.140
        ROM=0.773
        ROG=0.145
        ROM1=1.570
        ROG1=0.510
C     BETA1=-0.111
        BETA1=-0.221
      ENDIF
C -----------------------------------------------
      S=W**2
      FPIKMK=(BWIGM(S,ROM,ROG,XM1,XM2)+BETA1*BWIGM(S,ROM1,ROG1,XM1,XM2))
     + /(1+BETA1)
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION FPIRHO(W)
C **********************************************************
C     SQUARE OF PION FORM FACTOR
C **********************************************************
      COMPLEX FPIK
      FPIRHO=CABS(FPIK(W))**2
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION FPIRK(W)
C ----------------------------------------------------------
C     SQUARE OF PION FORM FACTOR
C ----------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C     COMPLEX FPIKMK
      COMPLEX FPIKM
      FPIRK=CABS(FPIKM(W,AMK,AMKZ))**2
C     FPIRK=CABS(FPIKMK(W,AMK,AMKZ))**2
      END
*CMZ :  1.01/44 05/01/96  16.53.44  by  Piero Zucchelli
*CMZ :  1.01/40 09/11/95  14.58.06  by  Piero Zucchelli
*-- Author :
      SUBROUTINE FZCLOS
*-----------------------------------------------------*
*                                                     *
*                   CLOSE FZ-FILE                     *
*                                                     *
*-----------------------------------------------------*
*KEEP,zebra.

      PARAMETER (NNQ=1000000)
*
      DIMENSION LQ(NNQ),IQ(NNQ),Q(NNQ)
      EQUIVALENCE (Q(1),IQ(1),LQ(9),JSTRUC(8))
      COMMON /QUEST/IQUEST(100)
      COMMON /XQSTOR/IXEVT,IFENCE(16),JGEEV,JSTRUC(99),JREFER(100),
     +DIV12(NNQ)
      COMMON /FZLUN/LUNFZ
      COMMON/MZIOALL/IOGENF

*KEND.
*
      CALL FZENDO(LUNFZ,'TX')
*
      END
*CMZ :  1.01/44 05/01/96  16.06.50  by  Piero Zucchelli
*CMZ :  1.01/40 10/11/95  19.03.13  by  Piero Zucchelli
*-- Author :
      SUBROUTINE FZINI
*-----------------------------------------------------*
*                                                     *
*               INITIALIZE ZEBRA                      *
*                                                     *
*-----------------------------------------------------*
*KEEP,zebra.

      PARAMETER (NNQ=1000000)
*
      DIMENSION LQ(NNQ),IQ(NNQ),Q(NNQ)
      EQUIVALENCE (Q(1),IQ(1),LQ(9),JSTRUC(8))
      COMMON /QUEST/IQUEST(100)
      COMMON /XQSTOR/IXEVT,IFENCE(16),JGEEV,JSTRUC(99),JREFER(100),
     +DIV12(NNQ)
      COMMON /FZLUN/LUNFZ
      COMMON/MZIOALL/IOGENF

*KEEP,info.
         COMMON/INFONEW/IRDATE,IRTIME
*KEND.
*
*---   INITIALISATION OF ZEBRA
*
      CALL DATIME(IRDATE,IRTIME)
C      CALL MZEBRA(-3)
      CALL MZEBRA( -1 )
*
*---   INITIALISATION OF DYNAMIC STORE
*
      NLIM=NNQ/2
      IXSTOR=10
      CALL MZSTOR(IXSTOR,'/XQSTOR/','.',IFENCE,JGEEV,JREFER(1),
     +            DIV12(1),DIV12(NLIM),DIV12(NNQ))

      NDIV=NNQ/10
      NDIVM=NDIV*5
      CALL MZDIV(IXSTOR,IXEVT,'EVT_DIV',NDIV,NDIVM,'.')
      CALL DZVERI('After init.',IXEVT,'CLSU')

      CALL FZOPN('jetta.rfz')
      CALL FZRUN(LUNFZ,99999,0,0)
      END
*CMZ :  1.01/45 08/01/96  11.11.53  by  Piero Zucchelli
*CMZ :  1.01/40 10/11/95  16.07.03  by  Piero Zucchelli
*-- Author :
      SUBROUTINE FZOPN(CHNAME)
*-----------------------------------------------------*
*                                                     *
*          OPEN FZ-FILE WITH NAME CHNAME              *
*                                                     *
*-----------------------------------------------------*
*KEEP,zebra.

      PARAMETER (NNQ=1000000)
*
      DIMENSION LQ(NNQ),IQ(NNQ),Q(NNQ)
      EQUIVALENCE (Q(1),IQ(1),LQ(9),JSTRUC(8))
      COMMON /QUEST/IQUEST(100)
      COMMON /XQSTOR/IXEVT,IFENCE(16),JGEEV,JSTRUC(99),JREFER(100),
     +DIV12(NNQ)
      COMMON /FZLUN/LUNFZ
      COMMON/MZIOALL/IOGENF

*KEND.
*
      CHARACTER CHNAME*(*)
      CHARACTER OPT*4
      LUNFZ = 17
      OPT = 'XLO'
      MED = 0
*
      IF (CHNAME(1:3).EQ.'exa') THEN
        CHNAME = '/dev/rmt0'
        OPT = OPT(1:3)//'T'
        MED = 1
      ENDIF

      CALL CFOPEN(LUNPTR,MED,8100,'w',0,CHNAME,ISTAT)
      IF (ISTAT.NE.0) THEN
        WRITE(LPUNIT,10020) ISTAT,LUNFZ
        STOP
      ENDIF
      IQUEST(1) = LUNPTR
      CALL FZFILE(LUNFZ,8100,OPT)
      IF (IQUEST(1).NE.0)  THEN
        WRITE (LPUNIT,10010) IQUEST(1),LUNFZ
        STOP
      ENDIF
10010 FORMAT (//' +++FILOPN - fatal error no. =',I5,' returned from',
     +' FZFILE, LUN = ',I5,' ++++++++')
10020 FORMAT (//' +++FILOPN - fatal error no. =',I5,' returned from',
     +' CFOPEN, LUN = ',I5,' ++++++++')
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C#######################################################################
C
C   ONE- AND TWO-DIMENSIONAL ADAPTIVE GAUSSIAN INTEGRATION ROUTINES.
C
C **********************************************************************

      SUBROUTINE GADAP(A0,B0,F,EPS,SUM)
C
C   PURPOSE           - INTEGRATE A FUNCTION F(X)
C   METHOD            - ADAPTIVE GAUSSIAN
C   USAGE             - CALL GADAP(A0,B0,F,EPS,SUM)
C   PARAMETERS  A0    - LOWER LIMIT (INPUT,REAL)
C               B0    - UPPER LIMIT (INPUT,REAL)
C               F     - FUNCTION F(X) TO BE INTEGRATED. MUST BE
C                       SUPPLIED BY THE USER. (INPUT,REAL FUNCTION)
C               EPS   - DESIRED RELATIVE ACCURACY. IF SUM IS SMALL EPS
C                       WILL BE ABSOLUTE ACCURACY INSTEAD. (INPUT,REAL)
C               SUM   - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT,REAL)
C   PRECISION         - SINGLE
C   REQ'D PROG'S      - F
C   AUTHOR            - T. JOHANSSON, LUND UNIV. COMPUTER CENTER, 1973
C   REFERENCE(S)      - THE AUSTRALIAN COMPUTER JOURNAL,3 P.126 AUG. -71
C
      COMMON/GADAP1/ NUM,IFU
      EXTERNAL F
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)
10000 FORMAT(16H GADAP:I TOO BIG)
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)
      IF(EPS.LT.1.0E-8) EPS=1.0E-8
      RED=1.3
      L=1
      I=1
      SUM=0.
      C=SQRT(15.)/5.
      A(1)=A0
      B(1)=B0
      F1(1)=F(0.5*(1+C)*A0+0.5*(1-C)*B0)
      F2(1)=F(0.5*(A0+B0))
      F3(1)=F(0.5*(1-C)*A0+0.5*(1+C)*B0)
      IFU=3
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)
   10 CONTINUE
      L=L+1
      N(L)=3
      EPS=EPS*RED
      A(I+1)=A(I)+C*(B(I)-A(I))
      B(I+1)=B(I)
      A(I+2)=A(I)+B(I)-A(I+1)
      B(I+2)=A(I+1)
      A(I+3)=A(I)
      B(I+3)=A(I+2)
      W1=A(I)+(B(I)-A(I))/5.
      U2=2.*W1-(A(I)+A(I+2))/2.
      F1(I+1)=F(A(I)+B(I)-W1)
      F2(I+1)=F3(I)
      F3(I+1)=F(B(I)-A(I+2)+W1)
      F1(I+2)=F(U2)
      F2(I+2)=F2(I)
      F3(I+2)=F(B(I+2)+A(I+2)-U2)
      F1(I+3)=F(A(I)+A(I+2)-W1)
      F2(I+3)=F1(I)
      F3(I+3)=F(W1)
      IFU=IFU+6
      IF(IFU.GT.5000) GOTO 40
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))
      SS=S(I+1)+S(I+2)+S(I+3)
      I=I+3
      IF(I.GT.300)GOTO 30
      SOLD=S(I-3)
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 10
      SUM=SUM+SS
      I=I-4
      N(L)=0
      L=L-1
   20 CONTINUE
      IF(L.EQ.1) GOTO 40
      N(L)=N(L)-1
      EPS=EPS/RED
      IF(N(L).NE.0) GOTO 10
      I=I-1
      L=L-1
      GOTO 20
   30 WRITE(6,10000)
   40 RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE GADAP2(A0,B0,FL,FU,F,EPS,SUM)
C
C   PURPOSE           - INTEGRATE A FUNCTION F(X,Y) OF TWO VARIABLES
C   METHOD            - ADAPTIVE GAUSSIAN IN BOTH DIRECTIONS
C   USAGE             - CALL GADAP2(A0,B0,FL,FU,F,EPS,SUM)
C   PARAMETERS  A0    - LOWER X-LIMIT (INPUT,REAL)
C               B0    - UPPER X-LIMIT (INPUT,REAL)
C               FL    - USER SUPPLIED FUNCTION FL(X) GIVING THE LOWER
C                       Y-LIMIT FOR A GIVEN X-VALUE
C                       (INPUT,REAL FUNCTION)
C               FU    - USER SUPPLIED FUNCTION FU(X) GIVING THE UPPER
C                       Y-LIMIT FOR A GIVEN X-VALUE
C                       (INPUT,REAL FUNCTION)
C               F     - USER SUPPLIED FUNCTION F(X,Y) TO BE INTEGRATED
C                       (INPUT,REAL FUNCTION)
C               EPS   - DESIRED ACCURACY (INPUT,REAL)
C               SUM   - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT,REAL)
C   PRECISION         - SINGLE
C   REQ'D PROG'S      - FL,FU,F,GADAPF
C   AUTHOR            - THOMAS JOHANSSON, LDC,1973
C
      COMMON/GADAP1/ NUM,IFU
      EXTERNAL F,FL,FU
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)
10000 FORMAT(16H GADAP:I TOO BIG)
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)
      IF(EPS.LT.1.0E-8) EPS=1.0E-8
      RED=1.4
      L=1
      I=1
      SUM=0.
      C=SQRT(15.)/5.
      A(1)=A0
      B(1)=B0
      X=0.5*(1+C)*A0+0.5*(1-C)*B0
      AY=FL(X)
      BY=FU(X)
      F1(1)=GADAPF(X,AY,BY,F,EPS)
      X=0.5*(A0+B0)
      AY=FL(X)
      BY=FU(X)
      F2(1)=GADAPF(X,AY,BY,F,EPS)
      X=0.5*(1-C)*A0+0.5*(1+C)*B0
      AY=FL(X)
      BY=FU(X)
      F3(1)=GADAPF(X,AY,BY,F,EPS)
      IFU=3
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)
   10 CONTINUE
      L=L+1
      N(L)=3
      EPS=EPS*RED
      A(I+1)=A(I)+C*(B(I)-A(I))
      B(I+1)=B(I)
      A(I+2)=A(I)+B(I)-A(I+1)
      B(I+2)=A(I+1)
      A(I+3)=A(I)
      B(I+3)=A(I+2)
      W1=A(I)+(B(I)-A(I))/5.
      U2=2.*W1-(A(I)+A(I+2))/2.
      X=A(I)+B(I)-W1
      AY=FL(X)
      BY=FU(X)
      F1(I+1)=GADAPF(X,AY,BY,F,EPS)
      F2(I+1)=F3(I)
      X=B(I)-A(I+2)+W1
      AY=FL(X)
      BY=FU(X)
      F3(I+1)=GADAPF(X,AY,BY,F,EPS)
      X=U2
      AY=FL(X)
      BY=FU(X)
      F1(I+2)=GADAPF(X,AY,BY,F,EPS)
      F2(I+2)=F2(I)
      X=B(I+2)+A(I+2)-U2
      AY=FL(X)
      BY=FU(X)
      F3(I+2)=GADAPF(X,AY,BY,F,EPS)
      X=A(I)+A(I+2)-W1
      AY=FL(X)
      BY=FU(X)
      F1(I+3)=GADAPF(X,AY,BY,F,EPS)
      F2(I+3)=F1(I)
      X=W1
      AY=FL(X)
      BY=FU(X)
      F3(I+3)=GADAPF(X,AY,BY,F,EPS)
      IFU=IFU+6
      IF(IFU.GT.5000) GOTO 40
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))
      SS=S(I+1)+S(I+2)+S(I+3)
      I=I+3
      IF(I.GT.300)GOTO 30
      SOLD=S(I-3)
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 10
      SUM=SUM+SS
      I=I-4
      N(L)=0
      L=L-1
   20 CONTINUE
      IF(L.EQ.1) GOTO 40
      N(L)=N(L)-1
      EPS=EPS/RED
      IF(N(L).NE.0) GOTO 10
      I=I-1
      L=L-1
      GOTO 20
   30 WRITE(6,10000)
   40 RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************


      FUNCTION GADAPF(X,A0,B0,F,EPS)
      COMMON/GADAP1/ NUM,IFU
      EXTERNAL F
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)
10000 FORMAT(16H GADAP:I TOO BIG)
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)
      IF(EPS.LT.1.0E-8) EPS=1.0E-8
      RED=1.4
      L=1
      I=1
      SUM=0.
      C=SQRT(15.)/5.
      A(1)=A0
      B(1)=B0
      F1(1)=F(X,0.5*(1+C)*A0+0.5*(1-C)*B0)
      F2(1)=F(X,0.5*(A0+B0))
      F3(1)=F(X,0.5*(1-C)*A0+0.5*(1+C)*B0)
      IFU=3
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)
   10 CONTINUE
      L=L+1
      N(L)=3
      EPS=EPS*RED
      A(I+1)=A(I)+C*(B(I)-A(I))
      B(I+1)=B(I)
      A(I+2)=A(I)+B(I)-A(I+1)
      B(I+2)=A(I+1)
      A(I+3)=A(I)
      B(I+3)=A(I+2)
      W1=A(I)+(B(I)-A(I))/5.
      U2=2.*W1-(A(I)+A(I+2))/2.
      F1(I+1)=F(X,A(I)+B(I)-W1)
      F2(I+1)=F3(I)
      F3(I+1)=F(X,B(I)-A(I+2)+W1)
      F1(I+2)=F(X,U2)
      F2(I+2)=F2(I)
      F3(I+2)=F(X,B(I+2)+A(I+2)-U2)
      F1(I+3)=F(X,A(I)+A(I+2)-W1)
      F2(I+3)=F1(I)
      F3(I+3)=F(X,W1)
      IFU=IFU+6
      IF(IFU.GT.5000) GOTO 40
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))
      SS=S(I+1)+S(I+2)+S(I+3)
      I=I+3
      IF(I.GT.300)GOTO 30
      SOLD=S(I-3)
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 10
      SUM=SUM+SS
      I=I-4
      N(L)=0
      L=L-1
   20 CONTINUE
      IF(L.EQ.1) GOTO 40
      N(L)=N(L)-1
      EPS=EPS/RED
      IF(N(L).NE.0) GOTO 10
      I=I-1
      L=L-1
      GOTO 20
   30 WRITE(6,10000)
   40 GADAPF=SUM
      EPS=EPS/RED
      RETURN
      END
*CMZ :  1.01/15 14/05/95  11.30.16  BY  PIERO ZUCCHELLI
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.35.13  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 08/07/94  16.07.34  BY  PIERO ZUCCHELLI
*CMZ :  1.00/01 17/03/92  18.04.44  BY  UNKNOWN
*-- AUTHOR :
      SUBROUTINE GBINIT
*KEEP,CDEBEAM.
C--
      COMMON /CONTRO/  BINIT,LUNB,NPNEUT,NPANTI,CPNORM,XPSOUR,SIGDIV
      COMMON /FLUXES/  FLUXD(80000),WEIGHT(8),SPECD(800),SPECN(800)
      COMMON /INPUT/   IALL(80000),NCOUNT(8)
      LOGICAL          BINIT


C-
*KEND.
      DIMENSION IBUF(16)
      LUNB   = 10
      NPNEUT = 800000
      NPANTI = 800000
      CPNORM = 1.E9
      XPSOUR = 67860.
      SIGDIV = 0.0004
      BINIT  = .FALSE.
C
C- READ MONTE CARLO BEAM SPECTRA FROM FILE
C
      DO I=1,5000
        READ(UNIT=LUNB,FMT=101,END=10) (IBUF(J),J=1,16)
        DO J=1,16
          II = (I-1)*16+J
          IALL(II) = IBUF(J)
        ENDDO
      ENDDO
   10 CONTINUE
C
C- INTEGRATE SPECTRA
C
      DO I=1,8
        NCOUNT(I) = 0
        DO J=1,10000
          II = (I-1)*10000+J
          NCOUNT(I) = NCOUNT(I)+IALL(II)
        ENDDO
      ENDDO
C
C- NORMALISE SPECTRA
C
      CNEUT = CPNORM/1.E3/FLOAT(NPNEUT)
      CANTI = CPNORM/1.E3/FLOAT(NPANTI)
      DO I=1,40000
        FLUXD(I) = CNEUT*FLOAT(IALL(I))
        FLUXD(40000+I) = CANTI*FLOAT(IALL(40000+I))
      ENDDO
C
C- NORMALISE INTEGRALS, PREPARE CUMULATIVE DISTRIBUTION
C
      DO I=1,4
        WEIGHT(I) = CNEUT*NCOUNT(I)
        WEIGHT(4+I) = CANTI*NCOUNT(4+I)
      ENDDO
      DO I=2,4
        WEIGHT(I) = WEIGHT(I)+WEIGHT(I-1)
        WEIGHT(4+I) = WEIGHT(4+I)+WEIGHT(3+I)
      ENDDO
C
C- INTEGRATE OVER RADIUS, PREPARE CUMULATIVE SPECTRA
C
      CALL VZERO(SPECD,800)
      CALL VZERO(SPECN,800)
      DO I=1,8
        DO J=1,100
          II = (I-1)*100+J
          DO K=1,100
            JJ = (I-1)*10000+J+(K-1)*100
            SPECD(II) = SPECD(II)+FLUXD(JJ)
            SPECN(II) = SPECD(II)
          ENDDO
        ENDDO
      ENDDO
      DO I=1,8
        DO J=2,100
          II = (I-1)*100+J
          SPECN(II) = SPECN(II)+SPECN(II-1)
        ENDDO
      ENDDO
C
      BINIT = .TRUE.
      RETURN
101   FORMAT(16I5)
      END
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.35.13  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/07/94  14.12.00  BY  PIERO ZUCCHELLI
*CMZ :  1.00/01 17/03/92  18.04.44  BY  UNKNOWN
*-- AUTHOR :
      SUBROUTINE GBSPEC(BEAM,IFLAV,RADIUS,SPEC)
*KEEP,CDEBEAM.
C--
      COMMON /CONTRO/  BINIT,LUNB,NPNEUT,NPANTI,CPNORM,XPSOUR,SIGDIV
      COMMON /FLUXES/  FLUXD(80000),WEIGHT(8),SPECD(800),SPECN(800)
      COMMON /INPUT/   IALL(80000),NCOUNT(8)
      LOGICAL          BINIT


C-
*KEND.
      CHARACTER*4 BEAM
      DIMENSION SPEC(100)
      COMMON/MAXSPEC/RMAXSPEC,RINTSPEC
C
      IF(.NOT.BINIT) CALL GBINIT
C
      IF(BEAM.EQ.'NEUT') THEN
        IBEAM=0
      ELSEIF(BEAM.EQ.'ANTI') THEN
        IBEAM=4
      ELSE
        GOTO 10
      ENDIF
C
      IF(IFLAV.LT.1.OR.IFLAV.GT.4)               GOTO 10
      IF(RADIUS.LE.0.0.OR.RADIUS.GE.300.0)       GOTO 10
      IFLAV = IFLAV+IBEAM
C
      CALL VZERO(SPEC,100)
      RAD2   = RADIUS**2/900.+1.
      IRADIU = RAD2
      FRAC   = RAD2-FLOAT(IRADIU)*900.
      DO I=1,100
        DO J=1,IRADIU-1
          IP = (IFLAV-1)*10000+I+(J-1)*100
          SPEC(I) = SPEC(I)+FLUXD(IP)
        ENDDO
        II = (IFLAV-1)*10000+I+(IRADIU-1)*100
        IF(FRAC.GT.0.0) SPEC(I) = SPEC(I)+FRAC*FLUXD(II)
      ENDDO
C
      RINTSPEC=0
      DO I=1,100
        RMAXSPEC=MAX(RMAXSPEC,SPEC(I))
        RINTSPEC=RINTSPEC+SPEC(I)
      END DO
*     WRITE(*,*)'RMAXSPEC=',RMAXSPEC
      IFLAV = IFLAV-IBEAM
      RETURN
C
   10 WRITE(6,10000)
10000 FORMAT(1X,' GBSPEC: ERROR IN INPUT VARIABLES!')
      RETURN
      END
*CMZ :  1.02/05 13/01/97  15.02.17  by  P. Zucchelli
*CMZ :  1.02/04 13/01/97  14.41.19  by  P. Zucchelli
*CMZ :  1.02/00 12/01/97  16.15.37  by  J. Brunner
*CMZ :  1.01/50 17/04/96  21.39.11  by  Piero Zucchelli
*CMZ :  1.01/41 12/12/95  17.03.06  by  Piero Zucchelli
*CMZ :  1.01/39 20/10/95  14.40.20  by  Piero Zucchelli
*CMZ :  1.01/38 18/10/95  18.27.55  by  Piero Zucchelli
*CMZ :  1.01/37 18/10/95  18.21.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/34 25/07/95  12.04.18  BY  PIERO ZUCCHELLI
*CMZ :  1.01/31 02/06/95  20.18.22  BY  PIERO ZUCCHELLI
*CMZ :  1.01/30 02/06/95  20.08.19  BY  PIERO ZUCCHELLI
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.35.13  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 23/09/94  12.01.45  BY  PIERO ZUCCHELLI
*-- AUTHOR :    PIERO ZUCCHELLI   11/09/94
      SUBROUTINE GENTABLE(LFILE,LEPIN,ENERGY_FIX,PPZ,INTERACTION)
*KEEP,JETTA.
C--
         PARAMETER (ICENTO=100)
         PARAMETER (ISIZ=93)
         PARAMETER (IOF1=32)
         PARAMETER (IOF2=83)
         PARAMETER (LUX_LEVEL=4)
         INTEGER*4 JTAU(100),JPRI(100),JSTRO(100)
         REAL*4 FTUPLE(ISIZ)
         COMMON/JETTAGL/JTAU,JPRI,JSTRO
         COMMON/NTUPLA/FTUPLE,ISFIRST
         COMMON/BEAM/SPEC(ICENTO)
         COMMON /MAXSPEC/RMAXSPEC,RINTSPEC
         COMMON/SAV/XMINSAV(ICENTO),XMAXSAV(ICENTO),YMINSAV(ICENTO),
     &   YMAXSAV(ICENTO),Q2MINSAV(ICENTO),Q2MAXSAV(ICENTO),
     &   W2MINSAV(ICENTO),W2MAXSAV(ICENTO),PARIMAX(ICENTO),
     &   PPSAVE(ICENTO,3,4,5),PARICOR(ICENTO),INDEX,SIGMASAV(ICENTO),
     &   XMSIGMA,XSECT

*KEEP,KEYS.
        COMMON/CFREAD/SPACE(5000)
 	COMMON/MYKEYS/IEVAR,IF5CC,INEUT,IINTE,IFERM,IFLAT,ICOUN,
     &  REFIX,IMUDO,IKAT1,IKAT2,IKAT3,IKAT4,IKAT5,IKAT6,
     &  INEVT,IJAK1,IJAK2,IITDK,RPTAU,RXK0D,NTGR,IDIMUON,IGLU,
     &  IQDEN,LOME(2),IFILES,ICCHA,ISEED,IDSUBS,IONEM,EHAC,IPSEL,
     &  IHIST


*KEND.
      REAL VECT(3),GKIN(3),G4MES(4)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/RUNCOM/IMODE

      IMODE=0
      CALL GETNEU(IPNUM,NEUTYPE,VECT,GKIN,
     + MESTYPE,G4MES,NEUFORCE,IMODE)
      IF (LST(17).EQ.0) THEN
        CALL LINIT(LFILE,LEPIN,ENERGY_FIX,PPZ,INTERACTION)
        WRITE(55,*)ENERGY_FIX,PARL(23),XMIN,XMAX,YMIN,YMAX,Q2MIN,Q2MAX,1.
        CLOSE(55)
      ELSE
* 5 PERCENT ON DSIGMA/DE LOOKS REASONABLE
*         PARL(15)=0.05
        RMAXSPEC=0.
        DO II=1,100
          I=II
          ENE=(I-1)*3. + 1.5

          CALL CATS

          CALL LINIT(LFILE,LEPIN,ENE,PPZ,INTERACTION)

          PARICOR(I)=PARI(32)
          PARIMAX(I)=PARI(LST(23))*1.5
          XMINSAV(I)=XMIN
          XMAXSAV(I)=XMAX
          YMINSAV(I)=YMIN
          YMAXSAV(I)=YMAX
          Q2MINSAV(I)=Q2MIN
          Q2MAXSAV(I)=Q2MAX
          W2MINSAV(I)=W2MIN
          W2MAXSAV(I)=W2MAX
          SIGMASAV(I)=PARL(23)
          XMSIGMA=MAX(PARL(23),XMSIGMA)
          DO 10 IA=1,2
            DO 10 JA=1,5
   10     PPSAVE(I,3,IA,JA)=PSAVE(3,IA,JA)

          IF (LST(37).EQ.0) RMAXSPEC=MAX(RMAXSPEC,SPEC(I))
          WRITE(*,*)'X-SECTION AT',ENE,' GEV=',PARL(23),
     +    'PB  BEAM WEIGHTED:',SPEC(I)
          WRITE(55,*)ENE,PARL(23),XMIN,XMAX,YMIN,
     +    YMAX,Q2MIN,Q2MAX,SPEC(I),PARI(LST(23))
        END DO
        CLOSE(55)
        CALL LINIT(LFILE,LEPIN,100.,PPZ,INTERACTION)
      ENDIF
      RETURN
      END
*CMZ :  1.02/09 14/01/97  15.14.45  by  P. Zucchelli
*CMZ :  1.02/08 14/01/97  11.55.37  by  P. Zucchelli
*CMZ :  1.02/07 14/01/97  11.11.46  by  P. Zucchelli
*CMZ :  1.02/06 13/01/97  18.51.57  by  P. Zucchelli
*CMZ :  1.02/04 13/01/97  14.54.54  by  P. Zucchelli
*CMZ :  1.02/00 12/01/97  16.23.16  by  J. Brunner
*CMZ :  1.01/51 04/07/96  10.20.58  by  Piero Zucchelli
*CMZ :  1.01/50 23/05/96  12.34.50  by  Piero Zucchelli
*-- Author :
C======================================================================
C======================================================================
C
      SUBROUTINE GETHNEU(IPNUM,NEUTYPE,VECT,GKIN,
     + MESTYPE,G4MES,NEUFORCE,IMODE)

C  IPNUM=protons number (back)
C  NEUTYPE=neutrino type (back)
C  VECT=position of creation of neutrino (back)
C  GKIN=neutrino 3-momentum
C  MESTYPE=parent meson type (back)
C  G4MES=meson 4 momentum
C  NEUFORCE=to force particular decays (to be implemented)
C  IMODE=0 to start (input), 2 to get neutrinos (input),
C  it is set to 4 by the routine when the input files are completed;
C  to continue, set it to 2 again.
C
C
      PARAMETER(NBIN=100)
      CHARACTER*80 FILNAME
      CHARACTER CHVAR(16)*8
C
      DIMENSION VECT(3),GKIN(3),YARR(100)
      DIMENSION G4MES(4),IHIS(4)
C
      COMMON/NTUPL10/      NUTYPE,IPARENT,EPARENT,XDECAY,YDECAY,ZDECAY,
     +                     PXPAR,PYPAR,PZPAR,XDET,YDET,XL,PXNU,PYNU,PZNU,
     +                     NPROT
      COMMON/RUNCOM/IMODEOLD,IFILES,IRUN
C
      DATA IHIS/3001,4001,1001,2001/
      DATA CHVAR/'NUTYPE',   'IPARENT',   'EPARENT',   'XDECAY ',
     +           'YDECAY',   'ZDECAY',    'PXPAR',     'PYPAR',
     +           'PZPAR',    'XDET',      'YDET',      'XL',
     +           'PXNU',     'PYNU',      'PZNU',      'NPROT'/
C
C     -----------------------------------------------------------------
      IF (IMODE.EQ.0) THEN
C
C       Global initialization (IMODE=0)
C       Start with Neutrino Run Number 1
C
        IRUN =1
        ICOUNT =0
        NUMBER =0
        IDIFF =0
        IFILES =0
        NBASE=0
C
C       Zero Neutrino Counter and Return
C
        RETURN
      ENDIF
C     -----------------------------------------------------------------
C     -----------------------------------------------------------------
      IF (IMODE.EQ.2) THEN
C
C       Event Processing (IMODE=2)
C       ---------------------------------------------------------------
C
        IF(ICOUNT.EQ.0)THEN
          WRITE(FILNAME,1)'../beam/histos.rz'
    1          FORMAT(A)
C
C

          LREC = 0
          CALL HROPEN(2,'BDIR',FILNAME,' ',LREC,ISTAT)
* PZ LOOP
          INDEX=NEUFORCE-48
          IHI=IHIS(INDEX)
*e mo`?
* something like:

          CALL HCDIR('//BDIR',' ')
          CALL HRIN(IHI,9999,0)
          CALL HUNPAK (IHI,YARR,' ',0)

* put IHI histogram into YARR
          CALL HISPRE(YARR,NBIN)
          IF (ISTAT.NE.0) THEN
            STOP 'HROPEN ERROR'
          ENDIF


*up to here
          IF (IMODE.NE.4) THEN
            WRITE(6,55) FILNAME
          ENDIF
   55    FORMAT(1X,'OPENING HISTOS FILE ',A)
C
        ENDIF
C
* here play with dices...

        NUTYPE=NEUFORCE
* here extract gkin(3) according to the histogram distribution....

        CALL HISRAN(YARR,NBIN,0.,3.,XRAN)
        GKIN(3)=XRAN
        IF (NUTYPE.EQ.53) NUTYPE=49
        IF (NUTYPE.EQ.54) NUTYPE=50
c
        IF(IERR.NE.0)THEN
          PRINT *, 'ERROR', ICOUNT
        ELSE
          ICOUNT=ICOUNT+1
        ENDIF

        IF (NEUFORCE.NE.0) THEN
*          CALL FORCED_DECAY(NEUFORCE,ISTATUS)
        ENDIF


C
        G4MES(1)=0.
        G4MES(2)=0.
        G4MES(3)=0.
        G4MES(4)=0.
        MESTYPE=0.
        VECT(1) =0.
        VECT(2) =0.
        VECT(3) =0.
        GKIN(1) =0.
        GKIN(2) =0.
        GKIN(3) =XRAN
        IPNUM = NBASE+1
        NEUTYPE = NUTYPE


C
C
C
C        IF(ICOUNT.LE.1.AND.IMODE.EQ.2)THEN
C          WRITE(6,500)NPROT,NUTYPE, VECT(1),VECT(2),VECT(3), GKIN(1),
C     +    GKIN(2),GKIN(3)
C        END IF
C
  500          FORMAT(1x,'POT=',I10,1X,
     +           'NEUT=',I3 ,1X,' x =',E15.9
     +                              ,1X,' y =',E15.9
     +                              ,1X,' z =',E15.9
     +                              ,1X,' px=',E15.9
     +                              ,1X,' py=',E15.9
     +                              ,1X,' pz=',E15.9)
C
      ENDIF
      IF (IMODE.EQ.4) THEN
        WRITE(*,*)' END OF NEUTRINO STATISTICS ...REWINDING...'
        WRITE(*,*)' TO CONTINUE, you have to SET IMODE=2'
      ENDIF
C
      RETURN
      END
*CMZ :  1.02/09 14/01/97  15.14.45  by  P. Zucchelli
*CMZ :  1.02/06 13/01/97  17.18.47  by  P. Zucchelli
*CMZ :  1.02/00 12/01/97  16.23.16  by  J. Brunner
*CMZ :  1.01/51 04/07/96  10.20.58  by  Piero Zucchelli
*CMZ :  1.01/50 23/05/96  12.34.50  by  Piero Zucchelli
*-- Author :
C======================================================================
C======================================================================
C
      SUBROUTINE GETNEU(IPNUM,NEUTYPE,VECT,GKIN,
     + MESTYPE,G4MES,NEUFORCE,IMODE)

C  IPNUM=protons number (back)
C  NEUTYPE=neutrino type (back)
C  VECT=position of creation of neutrino (back)
C  GKIN=neutrino 3-momentum
C  MESTYPE=parent meson type (back)
C  G4MES=meson 4 momentum
C  NEUFORCE=to force particular decays (to be implemented)
C  IMODE=0 to start (input), 2 to get neutrinos (input),
C  it is set to 4 by the routine when the input files are completed;
C  to continue, set it to 2 again.
C
C
      CHARACTER*80 FILNAME
      CHARACTER CHVAR(16)*8
C
      DIMENSION VECT(3),GKIN(3)
      DIMENSION G4MES(4)
C
      COMMON/NTUPL10/      NUTYPE,IPARENT,EPARENT,XDECAY,YDECAY,ZDECAY,
     +                     PXPAR,PYPAR,PZPAR,XDET,YDET,XL,PXNU,PYNU,PZNU,
     +                     NPROT
      COMMON/RUNCOM/IMODEOLD,IFILES,IRUN
C
      DATA CHVAR/'NUTYPE',   'IPARENT',   'EPARENT',   'XDECAY ',
     +           'YDECAY',   'ZDECAY',    'PXPAR',     'PYPAR',
     +           'PZPAR',    'XDET',      'YDET',      'XL',
     +           'PXNU',     'PYNU',      'PZNU',      'NPROT'/
C
C     -----------------------------------------------------------------
      IF (IMODE.EQ.0) THEN
C
C       Global initialization (IMODE=0)
C       Start with Neutrino Run Number 1
C
        IRUN =1
        ICOUNT =0
        NUMBER =0
        IDIFF =0
        IFILES =0
        NBASE=0
C
C       Zero Neutrino Counter and Return
C
        RETURN
      ENDIF
C     -----------------------------------------------------------------
C     -----------------------------------------------------------------
      IF (IMODE.EQ.2) THEN
C
C       Event Processing (IMODE=2)
C       ---------------------------------------------------------------
C
        IF(ICOUNT.EQ.0)THEN
C
C if end-of file
C
C         OPEN next RZ  file of neutrino
C
          IF(IRUN.LT.10)THEN
            WRITE(FILNAME,1)IRUN
    1          FORMAT('../beam/neutrino',I1,'.rz')
          ELSE IF (irun.LT.100)THEN
            WRITE(FILNAME,2)IRUN
    2          FORMAT('../beam/neutrino',I2,'.rz')
          ELSE IF (irun.LT.1000)THEN
            WRITE(FILNAME,3)IRUN
    3          FORMAT('../beam/neutrino',I3,'.rz')
          END IF
C
C
          LREC = 0
          CALL HROPEN(2,'BDIR',FILNAME,'X',LREC,ISTAT)
* PZ LOOP
          IF (ISTAT .NE. 0) THEN
            IMODE=4
            IRUN=1
            WRITE(FILNAME,1)IRUN
            LREC = 0
            CALL HROPEN(2,'BDIR',FILNAME,'X',LREC,ISTAT)
            IF (ISTAT.NE.0) THEN
              STOP 'HROPEN ERROR'
            ENDIF
          ENDIF


          CALL HRIN(0, 999, 0)
          CALL HNOENT(1,IMAX)
          IF (IMODE.NE.4) THEN
            WRITE(6,55) IMAX,FILNAME
          ENDIF
   55       FORMAT(1X,'OPENING ',I8,' EVENTS FROM FILE ',A)
C
C       ---------------------------------------------------------------
C      COMMON/NTUPL10/      NUTYPE,IPARENT,EPARENT,XDECAY,YDECAY,ZDECAY,
C     +                     PXPAR,PYPAR,PZPAR,XDET,YDET,XL,PXNU,PYNU,PZNU,
C     +                     NPROT
          CALL HBNAME(1,' ',0,'$CLEAR')
          CALL HBNAME(1, 'XNUMU', NUTYPE,'$SET:NUTYPE')
          CALL HBNAME(1, 'XNUMU', IPARENT,'$SET:IPARENT')
          CALL HBNAME(1, 'XNUMU', EPARENT,'$SET:EPARENT')
          CALL HBNAME(1, 'XNUMU', XDECAY,'$SET:XDECAY')
          CALL HBNAME(1, 'XNUMU', YDECAY,'$SET:YDECAY')
          CALL HBNAME(1, 'XNUMU', ZDECAY,'$SET:ZDECAY')
          CALL HBNAME(1, 'XNUMU', PXPAR,'$SET:PXPAR')
          CALL HBNAME(1, 'XNUMU', PYPAR,'$SET:PYPAR')
          CALL HBNAME(1, 'XNUMU', PZPAR,'$SET:PZPAR')
          CALL HBNAME(1, 'XNUMU', XDET, '$SET:XDET')
          CALL HBNAME(1, 'XNUMU', YDET, '$SET:YDET')
          CALL HBNAME(1, 'XNUMU', XL, '$SET:XL')
          CALL HBNAME(1, 'XNUMU', PXNU, '$SET:PXNU')
          CALL HBNAME(1, 'XNUMU', PYNU, '$SET:PYNU')
          CALL HBNAME(1, 'XNUMU', PZNU, '$SET:PZNU')
          CALL HBNAME(1, 'XNUMU', NPROT, '$SET:NPROT')
        ENDIF
C
C       READ next neutrinos parameters
C

C
        CALL HGNTV(1, CHVAR, 16, ICOUNT+1, IERR)
* gbeam patch!!!
        IF (NUTYPE.EQ.53) NUTYPE=49
        IF (NUTYPE.EQ.54) NUTYPE=50
c
        IF(IERR.NE.0)THEN
          PRINT *, 'ERROR', ICOUNT
        ELSE
          ICOUNT=ICOUNT+1
        ENDIF

        IF (NEUFORCE.NE.0) THEN
*          CALL FORCED_DECAY(NEUFORCE,ISTATUS)
        ENDIF


C
        G4MES(1)=PXPAR
        G4MES(2)=PYPAR
        G4MES(3)=PZPAR
        G4MES(4)=EPAR
        MESTYPE=IPARENT
        VECT(1) =XDECAY
        VECT(2) =YDECAY
        VECT(3) =ZDECAY
        GKIN(1) = PXNU
        GKIN(2) = PYNU
        GKIN(3) = PZNU
        IPNUM = NBASE+NPROT
        NEUTYPE = NUTYPE


C
C
C
C        IF(ICOUNT.LE.1.AND.IMODE.EQ.2)THEN
C          WRITE(6,500)NPROT,NUTYPE, VECT(1),VECT(2),VECT(3), GKIN(1),
C     +    GKIN(2),GKIN(3)
C        END IF
C
  500          FORMAT(1x,'POT=',I10,1X,
     +           'NEUT=',I3 ,1X,' x =',E15.9
     +                              ,1X,' y =',E15.9
     +                              ,1X,' z =',E15.9
     +                              ,1X,' px=',E15.9
     +                              ,1X,' py=',E15.9
     +                              ,1X,' pz=',E15.9)
C
C

        IF(ICOUNT.EQ.IMAX)THEN
          ICOUNT=0
          NBASE=NBASE+NPROT
          IRUN=IRUN+1
          CALL HREND('BDIR')
          WRITE(6,56)FILNAME
   56         FORMAT(1X,'CLOSING FILE',1X,A)
C
        ENDIF
C
      ENDIF
      IF (IMODE.EQ.4) THEN
        WRITE(*,*)' END OF NEUTRINO STATISTICS ...REWINDING...'
        WRITE(*,*)' TO CONTINUE, you have to SET IMODE=2'
      ENDIF
C
      RETURN
      END
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION GFUN(QKWA)
C ****************************************************************
C     G-FUNCTION USED TO INRODUCE ENERGY DEPENDENCE IN A1 WIDTH
C ****************************************************************
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      IF (QKWA.LT.(AMRO+AMPI)**2) THEN
        GFUN=4.1*(QKWA-9*AMPIZ**2)**3 *(1.-3.3*(QKWA-9*AMPIZ**2)+5.8*
     +  (QKWA-9*AMPIZ**2)**2)
      ELSE
        GFUN=QKWA*(1.623+10.38/QKWA-9.32/QKWA**2+0.65/QKWA**3)
      ENDIF
      END
*CMZ :  1.01/17 14/05/95  11.47.38  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 14/08/94  03.47.49  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE INIMAS
C ----------------------------------------------------------------------
C     INITIALISATION OF MASSES
C
C     CALLED BY : KORALZ
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
C IN-COMING / OUT-GOING  FERMION MASSES
      AMTAU  = 1.7771
C      AMNUTA = 0.010
      AMNUTA = 0.0
      AMEL   = 0.0005111
      AMNUE  = 0.0
      AMMU   = 0.105659
      AMNUMU = 0.0
C
C MASSES USED IN TAU DECAYS
      AMPIZ  = 0.134964
      AMPI   = 0.139568
      AMRO   = 0.773
      GAMRO  = 0.145
CC    GAMRO  = 0.666
      AMA1   = 1.251
      GAMA1  = 0.599
      AMK    = 0.493667
      AMKZ   = 0.49772
      AMKST  = 0.8921
      GAMKST = 0.0513
C
      RETURN
      END
*CMZ :  1.00/00 09/08/94  17.43.59  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE INIPHY(XK00)
C ----------------------------------------------------------------------
C     INITIALISATION OF PARAMETERS
C     USED IN QED AND/OR GSW ROUTINES
C ----------------------------------------------------------------------
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
      REAL*8 PI8,XK00
C
      PI8    = 4.D0*DATAN(1.D0)
      ALFINV = 137.03604D0
      ALFPI  = 1D0/(ALFINV*PI8)
      XK0=XK00
      END
*CMZ :  1.01/50 22/05/96  18.06.09  by  Piero Zucchelli
*CMZ :  1.01/26 29/05/95  19.08.20  BY  PIERO ZUCCHELLI
*CMZ :  1.01/25 29/05/95  16.05.52  BY  PIERO ZUCCHELLI
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/11 13/05/95  19.10.40  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 09/08/94  17.43.59  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE INITDK
C ----------------------------------------------------------------------
C     INITIALISATION OF TAU DECAY PARAMETERS  AND ROUTINES
C
C     CALLED BY : KORALZ
C ----------------------------------------------------------------------
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / TAUBRA / GAMPRT(30),JLIST(30),NCHAN
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / MDECOMP /IDFFIN(9,NMODE),MULPIK(NMODE)
     +                ,NAMES
      CHARACTER NAMES(NMODE)*31
      REAL*4 PI
C
C LIST OF BRANCHING RATIOS
CAM NORMALISED TO E NU NUTAU CHANNEL
CAM                  ENU   MUNU   PINU  RHONU   A1NU   KNU    K*NU   PI'
CAM   DATA JLIST  /    1,     2,     3,     4,     5,     6,     7,
CAM   DATA GAMPRT /1.000,0.9730,0.6054,1.2432,0.8432,0.0432,O.O811,0.616
CAM
CAM  MULTIPION DECAYS
C
C    CONVENTIONS OF PARTICLES NAMES
C                 K-,P-,K+,  K0,P-,KB,  K-,P0,K0
C                  3, 1,-3  , 4, 1,-4  , 3, 2, 4  ,
C                 P0,P0,K-,  K-,P-,P+,  P-,KB,P0
C                  2, 2, 3  , 3, 1,-1  , 1,-4, 2  ,
C                 ET,P-,P0   P-,P0,GM
C                  9, 1, 2  , 1, 2, 8
C
      DIMENSION NOPIK(6,NMODE),NPIK(NMODE)
CAM   OUTGOING MULTIPLICITY AND FLAVORS OF MULTI-PION /MULTI-K MODES
      DATA   NPIK  /                4,                    4,
     +                              5,                    5,
     +                              6,                    6,
     +                              3,                    3,
     +                              3,                    3,
     +                              3,                    3,
     +                              3,                    3,
     +                              2                         /
      DATA  NOPIK / -1,-1, 1, 2, 0, 0,     2, 2, 2,-1, 0, 0,
     +              -1,-1, 1, 2, 2, 0,    -1,-1,-1, 1, 1, 0,
     +              -1,-1,-1, 1, 1, 2,    -1,-1, 1, 2, 2, 2,
     +              -3,-1, 3, 0, 0, 0,    -4,-1, 4, 0, 0, 0,
     +              -3, 2,-4, 0, 0, 0,     2, 2,-3, 0, 0, 0,
     +              -3,-1, 1, 0, 0, 0,    -1, 4, 2, 0, 0, 0,
     +               9,-1, 2, 0, 0, 0,    -1, 2, 8, 0, 0, 0,
     +              -3, 4, 0, 0, 0, 0                         /
C LIST OF BRANCHING RATIOS
      NCHAN = NMODE + 7
      DO 10 I = 1,30
        IF (I.LE.NCHAN) THEN
          JLIST(I) = I
          IF(I.EQ. 1) GAMPRT(I) = 1.00000
          IF(I.EQ. 2) GAMPRT(I) = 0.97980
          IF(I.EQ. 3) GAMPRT(I) = 0.64960
* EX 1.3405,
          IF(I.EQ. 4) GAMPRT(I) = 1.3405
          IF(I.EQ. 5) GAMPRT(I) = 1.2
          IF(I.EQ. 6) GAMPRT(I) = 0.0397
          IF(I.EQ. 7) GAMPRT(I) = 0.0696
* CHANGED FOR INCREASING 3 PROG DECAYS
          IF(I.EQ. 8) GAMPRT(I) = 0.0835
          IF(I.EQ. 9) GAMPRT(I) = 0.0170
          IF(I.EQ.10) GAMPRT(I) = 0.0641
          IF(I.EQ.11) GAMPRT(I) = 0.00286
          IF(I.EQ.12) GAMPRT(I) = 0.0043
          IF(I.EQ.13) GAMPRT(I) = 0.0042
          IF(I.EQ.14) GAMPRT(I) = 0.0061
          IF(I.EQ.15) GAMPRT(I) = 0.0056
          IF(I.EQ.16) GAMPRT(I) = 0.0005
          IF(I.EQ.17) GAMPRT(I) = 0.0059
          IF(I.EQ.18) GAMPRT(I) = 0.0321
          IF(I.EQ.19) GAMPRT(I) = 0.0320
          IF(I.EQ.20) GAMPRT(I) = 0.0110
          IF(I.EQ.21) GAMPRT(I) = 0.0031
          IF(I.EQ.22) GAMPRT(I) = 0.0181
          IF(I.EQ. 8) NAMES(I-7)='  TAU-  --> 2PI-,  PI0,  PI+   '
          IF(I.EQ. 9) NAMES(I-7)='  TAU-  --> 3PI0,        PI-   '
          IF(I.EQ.10) NAMES(I-7)='  TAU-  --> 2PI-,  PI+, 2PI0   '
          IF(I.EQ.11) NAMES(I-7)='  TAU-  --> 3PI-, 2PI+,        '
          IF(I.EQ.12) NAMES(I-7)='  TAU-  --> 3PI-, 2PI+,  PI0   '
          IF(I.EQ.13) NAMES(I-7)='  TAU-  --> 2PI-,  PI+, 3PI0   '
          IF(I.EQ.14) NAMES(I-7)='  TAU-  -->  K-, PI-,  K+      '
          IF(I.EQ.15) NAMES(I-7)='  TAU-  -->  K0, PI-, K0B      '
          IF(I.EQ.16) NAMES(I-7)='  TAU-  -->  K-,  K0, PI0      '
          IF(I.EQ.17) NAMES(I-7)='  TAU-  --> PI0, PI0,  K-      '
          IF(I.EQ.18) NAMES(I-7)='  TAU-  -->  K-, PI-, PI+      '
          IF(I.EQ.19) NAMES(I-7)='  TAU-  --> PI-, K0B, PI0      '
          IF(I.EQ.20) NAMES(I-7)='  TAU-  --> ETA, PI-, PI0      '
          IF(I.EQ.21) NAMES(I-7)='  TAU-  --> PI-, PI0, GAM      '
          IF(I.EQ.22) NAMES(I-7)='  TAU-  -->  K-,  K0           '
        ELSE
          JLIST(I) = 0
          GAMPRT(I) = 0.
        ENDIF
   10 CONTINUE
      DO I=1,NMODE
        MULPIK(I)=NPIK(I)
        DO J=1,MULPIK(I)
          IDFFIN(J,I)=NOPIK(J,I)
        ENDDO
      ENDDO
C
C
C --- COEFFICIENTS TO FIX RATIO OF:
C --- A1 3CHARGED/ A1 1CHARGED 2 NEUTRALS MATRIX ELEMENTS (MASLESS LIM.)
C --- PROBABILITY OF K0 TO BE KS
C --- PROBABILITY OF K0B TO BE KS
C --- RATIO OF COEFFICIENTS FOR K*--> K0 PI-
C --- ALL COEFFICENTS SHOULD BE IN THE RANGE (0.0,1.0)
C --- THEY MEANING IS PROBABILITY OF THE FIRST CHOICE ONLY IF ONE
C --- NEGLECTS MASS-PHASE SPACE EFFECTS
      BRA1=0.5
      BRK0=0.5
      BRK0B=0.5
      BRKS=0.6667
C
C --- REMAINING CONSTANTS
      PI =4.*ATAN(1.)
      GFERMI = 1.16637E-5
      CCABIB = 0.975
      GV     = 1.0
      GA     =-1.0
C ZW 13.04.89 HERE WAS AN ERROR
      SCABIB = SQRT(1.-CCABIB**2)
      GAMEL  = GFERMI**2*AMTAU**5/(192*PI**3)
C
C      CALL DEXAY(-1)
C
      RETURN
      END

*CMZ :  1.01/50 22/05/96  18.06.09  by  Piero Zucchelli
*CMZ :  1.01/22 29/05/95  15.21.25  BY  PIERO ZUCCHELLI
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/11 13/05/95  19.10.40  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 09/08/94  17.43.59  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE INITDK_NEW
C ----------------------------------------------------------------------
C     INITIALISATION OF TAU DECAY PARAMETERS  AND ROUTINES
C
C     CALLED BY : KORALZ
C ----------------------------------------------------------------------
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / TAUBRA / GAMPRT(30),JLIST(30),NCHAN
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / MDECOMP /IDFFIN(9,NMODE),MULPIK(NMODE)
     +                ,NAMES
      CHARACTER NAMES(NMODE)*31
      REAL*4 PI
C
C LIST OF BRANCHING RATIOS
CAM NORMALISED TO E NU NUTAU CHANNEL
CAM                  ENU   MUNU   PINU  RHONU   A1NU   KNU    K*NU   PI'
CAM   DATA JLIST  /    1,     2,     3,     4,     5,     6,     7,
CAM   DATA GAMPRT /1.000,0.9730,0.6054,1.2432,0.8432,0.0432,O.O811,0.616
CAM
CAM  MULTIPION DECAYS
C
C    CONVENTIONS OF PARTICLES NAMES
C                 K-,P-,K+,  K0,P-,KB,  K-,P0,K0
C                  3, 1,-3  , 4, 1,-4  , 3, 2, 4  ,
C                 P0,P0,K-,  K-,P-,P+,  P-,KB,P0
C                  2, 2, 3  , 3, 1,-1  , 1,-4, 2  ,
C                 ET,P-,P0   P-,P0,GM
C                  9, 1, 2  , 1, 2, 8
C
      DIMENSION NOPIK(6,NMODE),NPIK(NMODE)
CAM   OUTGOING MULTIPLICITY AND FLAVORS OF MULTI-PION /MULTI-K MODES
      DATA   NPIK  /                4,                    4,
     +                              5,                    5,
     +                              6,                    6,
     +                              3,                    3,
     +                              3,                    3,
     +                              3,                    3,
     +                              3,                    3,
     +                              2                         /
      DATA  NOPIK / -1,-1, 1, 2, 0, 0,     2, 2, 2,-1, 0, 0,
     +              -1,-1, 1, 2, 2, 0,    -1,-1,-1, 1, 1, 0,
     +              -1,-1,-1, 1, 1, 2,    -1,-1, 1, 2, 2, 2,
     +              -3,-1, 3, 0, 0, 0,    -4,-1, 4, 0, 0, 0,
     +              -3, 2,-4, 0, 0, 0,     2, 2,-3, 0, 0, 0,
     +              -3,-1, 1, 0, 0, 0,    -1, 4, 2, 0, 0, 0,
     +               9,-1, 2, 0, 0, 0,    -1, 2, 8, 0, 0, 0,
     +              -3, 4, 0, 0, 0, 0                         /
C LIST OF BRANCHING RATIOS
      NCHAN = NMODE + 7
      DO 10 I = 1,30
        IF (I.LE.NCHAN) THEN
          JLIST(I) = I
          IF(I.EQ. 1) GAMPRT(I) = 1.00000
          IF(I.EQ. 2) GAMPRT(I) = 0.98001
          IF(I.EQ. 3) GAMPRT(I) = 0.64964
          IF(I.EQ. 4) GAMPRT(I) = 1.39922
          IF(I.EQ. 5) GAMPRT(I) = 0.8432
          IF(I.EQ. 6) GAMPRT(I) = 0.03720
          IF(I.EQ. 7) GAMPRT(I) = 0.08051
          IF(I.EQ. 8) GAMPRT(I) = 0.0835
          IF(I.EQ. 9) GAMPRT(I) = 0.0170
          IF(I.EQ.10) GAMPRT(I) = 0.0641
          IF(I.EQ.11) GAMPRT(I) = 0.0286
          IF(I.EQ.12) GAMPRT(I) = 0.0043
          IF(I.EQ.13) GAMPRT(I) = 0.0042
          IF(I.EQ.14) GAMPRT(I) = 0.01222
          IF(I.EQ.15) GAMPRT(I) = 0.0056
          IF(I.EQ.16) GAMPRT(I) = 0.0005
          IF(I.EQ.17) GAMPRT(I) = 0.0059
          IF(I.EQ.18) GAMPRT(I) = 0.0321
          IF(I.EQ.19) GAMPRT(I) = 0.0320
          IF(I.EQ.20) GAMPRT(I) = 0.0110
          IF(I.EQ.21) GAMPRT(I) = 0.0031
          IF(I.EQ.22) GAMPRT(I) = 0.0181
          IF(I.EQ. 8) NAMES(I-7)='  TAU-  --> 2PI-,  PI0,  PI+   '
          IF(I.EQ. 9) NAMES(I-7)='  TAU-  --> 3PI0,        PI-   '
          IF(I.EQ.10) NAMES(I-7)='  TAU-  --> 2PI-,  PI+, 2PI0   '
          IF(I.EQ.11) NAMES(I-7)='  TAU-  --> 3PI-, 2PI+,        '
          IF(I.EQ.12) NAMES(I-7)='  TAU-  --> 3PI-, 2PI+,  PI0   '
          IF(I.EQ.13) NAMES(I-7)='  TAU-  --> 2PI-,  PI+, 3PI0   '
          IF(I.EQ.14) NAMES(I-7)='  TAU-  -->  K-, PI-,  K+      '
          IF(I.EQ.15) NAMES(I-7)='  TAU-  -->  K0, PI-, K0B      '
          IF(I.EQ.16) NAMES(I-7)='  TAU-  -->  K-,  K0, PI0      '
          IF(I.EQ.17) NAMES(I-7)='  TAU-  --> PI0, PI0,  K-      '
          IF(I.EQ.18) NAMES(I-7)='  TAU-  -->  K-, PI-, PI+      '
          IF(I.EQ.19) NAMES(I-7)='  TAU-  --> PI-, K0B, PI0      '
          IF(I.EQ.20) NAMES(I-7)='  TAU-  --> ETA, PI-, PI0      '
          IF(I.EQ.21) NAMES(I-7)='  TAU-  --> PI-, PI0, GAM      '
          IF(I.EQ.22) NAMES(I-7)='  TAU-  -->  K-,  K0           '
        ELSE
          JLIST(I) = 0
          GAMPRT(I) = 0.
        ENDIF
   10 CONTINUE
      DO I=1,NMODE
        MULPIK(I)=NPIK(I)
        DO J=1,MULPIK(I)
          IDFFIN(J,I)=NOPIK(J,I)
        ENDDO
      ENDDO
C
C
C --- COEFFICIENTS TO FIX RATIO OF:
C --- A1 3CHARGED/ A1 1CHARGED 2 NEUTRALS MATRIX ELEMENTS (MASLESS LIM.)
C --- PROBABILITY OF K0 TO BE KS
C --- PROBABILITY OF K0B TO BE KS
C --- RATIO OF COEFFICIENTS FOR K*--> K0 PI-
C --- ALL COEFFICENTS SHOULD BE IN THE RANGE (0.0,1.0)
C --- THEY MEANING IS PROBABILITY OF THE FIRST CHOICE ONLY IF ONE
C --- NEGLECTS MASS-PHASE SPACE EFFECTS
      BRA1=0.5
      BRK0=0.5
      BRK0B=0.5
      BRKS=0.6667
C
C --- REMAINING CONSTANTS
      PI =4.*ATAN(1.)
      GFERMI = 1.16637E-5
      CCABIB = 0.975
      GV     = 1.0
      GA     =-1.0
C ZW 13.04.89 HERE WAS AN ERROR
      SCABIB = SQRT(1.-CCABIB**2)
      GAMEL  = GFERMI**2*AMTAU**5/(192*PI**3)
C
C      CALL DEXAY(-1)
C
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/11 13/05/95  18.27.58  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.39  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE JAKER(JAK)
C     *********************
C
C **********************************************************************
C                                                                      *
C           *********TAUOLA LIBRARY: VERSION 2.5 ********              *
C           **************JUNE     1994******************              *
C           **      AUTHORS: S.JADACH, Z.WAS        *****              *
C           **  R. DECKER, M. JEZABEK, J.H.KUEHN,   *****              *
C           ********AVAILABLE FROM: WASM AT CERNVM ******              *
C           *******PUBLISHED IN COMP. PHYS. COMM.********              *
C           *** PREPRINT CERN-TH-5856 SEPTEMBER 1990 ****              *
C           *** PREPRINT CERN-TH-6195 OCTOBER   1991 ****              *
C           *** PREPRINT CERN-TH-6793 NOVEMBER  1992 ****              *
C **********************************************************************
C
C ----------------------------------------------------------------------
C SUBROUTINE JAKER,
C CHOOSES DECAY MODE ACCORDING TO LIST OF BRANCHING RATIOS
C JAK=1 ELECTRON MODE
C JAK=2 MUON MODE
C JAK=3 PION MODE
C JAK=4 RHO  MODE
C JAK=5 A1   MODE
C JAK=6 K    MODE
C JAK=7 K*   MODE
C JAK=8 NPI  MODE
C
C     CALLED BY : DEXAY
C ----------------------------------------------------------------------
      COMMON / TAUBRA / GAMPRT(30),JLIST(30),NCHAN
      COMMON/BERI/JALLY,JEIN
C      REAL   CUMUL(20)
      REAL   CUMUL(30)
      INTEGER JALLY(30)
C
      IF(NCHAN.LE.0.OR.NCHAN.GT.30) GOTO 30
      CALL RANMAR(RRR,1)
      SUM=0
      DO 10 I=1,NCHAN
        SUM=SUM+GAMPRT(I)
   10 CUMUL(I)=SUM
      DO 20 I=NCHAN,1,-1
        IF(RRR.LT.CUMUL(I)/CUMUL(NCHAN)) JI=I
   20 CONTINUE
      JAK=JLIST(JI)
      JEIN=JAK
      JALLY(JAK)=JALLY(JAK)+1
      RETURN
   30 PRINT 10000
10000 FORMAT(' ----- JAKER: WRONG NCHAN')
      STOP
      END
*CMZ :          13/03/97  16.01.10  by  Unknown
*CMZ :  1.02/09 14/01/97  15.08.05  by  P. Zucchelli
*CMZ :  1.01/37 01/08/95  15.05.38  BY  PIERO ZUCCHELLI
*CMZ :  1.01/36 31/07/95  17.54.39  BY  PIERO ZUCCHELLI
*CMZ :  1.01/34 12/07/95  11.32.39  BY  PIERO ZUCCHELLI
*CMZ :  1.01/33 12/07/95  10.47.16  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.35.13  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 08/09/94  12.44.55  BY  PIERO ZUCCHELLI
*CMZ :  1.01/00 08/09/94  09.48.55  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 15/08/94  07.06.13  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 20/07/94  12.08.20  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE JETMC
*KEEP,CDEBEAM.
C--
      COMMON /CONTRO/  BINIT,LUNB,NPNEUT,NPANTI,CPNORM,XPSOUR,SIGDIV
      COMMON /FLUXES/  FLUXD(80000),WEIGHT(8),SPECD(800),SPECN(800)
      COMMON /INPUT/   IALL(80000),NCOUNT(8)
      LOGICAL          BINIT


C-
*KEEP,POLAR.
C--
      COMMON /POLARIZ/POL(4000,3)
      REAL POLARX(4)
*KEEP,LUDAT1.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUJETS.
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      SAVE /LUJETS/
*KEEP,FOREFI.
C--
	INTEGER*4 IEVT
         COMMON/FOREFICASS/IEVT


*KEEP,JETTA.
C--
         PARAMETER (ICENTO=100)
         PARAMETER (ISIZ=93)
         PARAMETER (IOF1=32)
         PARAMETER (IOF2=83)
         PARAMETER (LUX_LEVEL=4)
         INTEGER*4 JTAU(100),JPRI(100),JSTRO(100)
         REAL*4 FTUPLE(ISIZ)
         COMMON/JETTAGL/JTAU,JPRI,JSTRO
         COMMON/NTUPLA/FTUPLE,ISFIRST
         COMMON/BEAM/SPEC(ICENTO)
         COMMON /MAXSPEC/RMAXSPEC,RINTSPEC
         COMMON/SAV/XMINSAV(ICENTO),XMAXSAV(ICENTO),YMINSAV(ICENTO),
     &   YMAXSAV(ICENTO),Q2MINSAV(ICENTO),Q2MAXSAV(ICENTO),
     &   W2MINSAV(ICENTO),W2MAXSAV(ICENTO),PARIMAX(ICENTO),
     &   PPSAVE(ICENTO,3,4,5),PARICOR(ICENTO),INDEX,SIGMASAV(ICENTO),
     &   XMSIGMA,XSECT

*KEND.

      CALL VZERO(JTAU,100)
      CALL VZERO(JPRI,100)
      CALL VZERO(JSTRO,100)

      I=4
      DO I=4,N
        IPA2=I
   11   CONTINUE
        IF(K(IPA2,2).EQ.92) JSTRO(I)=1
        IGPA=K(IPA2,3)
        IF (IGPA.NE.0) THEN
          IPA2=IGPA
          GOTO 11
        ENDIF

* SEARCH FOR UNDECAYED PARTICLES
        IF(K(I,4).EQ.0.AND.K(I,1).LE.10.AND.K(I,3).NE.0.AND. ABS(K(I,2)
     +  ).NE.12.AND.ABS(K(I,2)).NE.14.AND. ABS(K(I,2)).NE.16) THEN

*            WRITE(88,*)'PT:',IEVT,I
* LET'S FILL IN THE INTERESTING CALORIMETER INFORMATIONS:
          ICHARGE=LUCHGE(K(I,2))
          ISMISSED=0
          FTUPLE(65)=FTUPLE(65)+P(I,4)
          FTUPLE(68)=FTUPLE(68)+1
          IF (ICHARGE.EQ.0) FTUPLE(71)=FTUPLE(71)+1
          IF (ABS(K(I,2)).EQ.11.OR.ABS(K(I,2)).EQ.22) THEN
* PHOTON,E ARE CONSIDERED EM ENERGY DEPOSIT
            FTUPLE(66)=FTUPLE(66)+P(I,4)
            FTUPLE(69)=FTUPLE(69)+1
            IF (ICHARGE.EQ.0) FTUPLE(72)=FTUPLE(72)+1
            ISMISSED=1
          ENDIF

          IF (ABS(K(I,2)).GE.100.and.k(i,1).le.10.and.
     +   k(i,1).ge.1) THEN
* STABLE PARTICLES WITH QUARKS RELEASE HADRONIC ENERGY, SO...
            FTUPLE(67)=FTUPLE(67)+P(I,4)
            FTUPLE(70)=FTUPLE(70)+1
            IF (ICHARGE.EQ.0) FTUPLE(73)=FTUPLE(73)+1
            ISMISSED=1
          ENDIF

          IF (ISMISSED.EQ.0.AND.ABS(K(I,2)).NE.13) THEN
            WRITE(*,*)' MISSED EM/HAD/MU PARTICLE:',K(I,2)
          ENDIF

* NOW SEARCH FOR TAU ANCESTOR
          IPA=I
   10     CONTINUE
          IGRANDPA=K(IPA,3)
          IF (IGRANDPA.NE.0) THEN
            IPA=IGRANDPA
            GOTO 10
          ENDIF
          ETOT=P(I,4)
*       WRITE(*,*)'I,IPA,CHG E=',I,IPA,LUCHGE(K(I,2))/3,ETOT

          IF (IPA.EQ.1) THEN
            JTAU(I)=1
          ELSE
            JPRI(I)=1
          ENDIF

        ENDIF

      END DO

      RETURN
      END
*CMZ :  1.02/09 14/01/97  15.58.27  by  P. Zucchelli
*CMZ :  1.01/52 19/11/96  17.37.58  by  Piero Zucchelli
*CMZ :  1.01/51 05/08/96  09.01.24  by  Piero Zucchelli
*CMZ :  1.01/50 23/05/96  12.34.50  by  Piero Zucchelli
*CMZ :  1.01/49 29/01/96  16.02.04  by  Piero Zucchelli
*CMZ :  1.01/48 15/01/96  15.30.50  by  Piero Zucchelli
*CMZ :  1.01/47 11/01/96  10.11.37  by  Piero Zucchelli
*CMZ :  1.01/46 09/01/96  11.47.50  by  Piero Zucchelli
*CMZ :  1.01/45 08/01/96  14.22.55  by  Piero Zucchelli
*CMZ :  1.01/43 15/12/95  18.03.32  by  Piero Zucchelli
*CMZ :  1.01/40 11/12/95  12.39.22  by  Piero Zucchelli
*CMZ :  1.01/39 06/11/95  14.45.27  by  Piero Zucchelli
*CMZ :  1.01/37 21/09/95  12.47.07  BY  PIERO ZUCCHELLI
*CMZ :  1.01/36 31/07/95  18.03.46  BY  PIERO ZUCCHELLI
*CMZ :  1.01/35 26/07/95  15.13.48  BY  PIERO ZUCCHELLI
*CMZ :  1.01/34 25/07/95  11.31.42  BY  PIERO ZUCCHELLI
*CMZ :  1.01/33 12/07/95  11.18.32  BY  PIERO ZUCCHELLI
*CMZ :  1.01/21 27/05/95  15.58.38  BY  PIERO ZUCCHELLI
*CMZ :  1.01/20 21/05/95  14.50.37  BY  PIERO ZUCCHELLI
*CMZ :  1.01/19 16/05/95  08.23.30  BY  UNKNOWN
*CMZ :  1.01/18 14/05/95  16.16.54  BY  PIERO ZUCCHELLI
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/13 14/05/95  11.22.22  BY  PIERO ZUCCHELLI
*CMZ :  1.01/11 14/05/95  11.17.15  BY  PIERO ZUCCHELLI
*CMZ :  1.01/10 13/05/95  10.26.44  BY  PIERO ZUCCHELLI
*CMZ :  1.01/09 20/04/95  12.52.49  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/07 05/03/95  11.34.14  BY  PIERO ZUCCHELLI
*CMZ :  1.01/06 05/03/95  11.20.39  BY  PIERO ZUCCHELLI
*CMZ :  1.01/03 04/03/95  23.35.00  BY  PIERO ZUCCHELLI
*CMZ :  1.01/02 04/03/95  21.38.55  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 04/03/95  21.18.20  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 21/08/94  11.05.12  BY  PIERO ZUCCHELLI
*-- AUTHOR :    PIERO ZUCCHELLI   29/07/94

      SUBROUTINE JETTOUT
*KEEP,FOREFI.
C--
	INTEGER*4 IEVT
         COMMON/FOREFICASS/IEVT


*KEEP,JETTA.
C--
         PARAMETER (ICENTO=100)
         PARAMETER (ISIZ=93)
         PARAMETER (IOF1=32)
         PARAMETER (IOF2=83)
         PARAMETER (LUX_LEVEL=4)
         INTEGER*4 JTAU(100),JPRI(100),JSTRO(100)
         REAL*4 FTUPLE(ISIZ)
         COMMON/JETTAGL/JTAU,JPRI,JSTRO
         COMMON/NTUPLA/FTUPLE,ISFIRST
         COMMON/BEAM/SPEC(ICENTO)
         COMMON /MAXSPEC/RMAXSPEC,RINTSPEC
         COMMON/SAV/XMINSAV(ICENTO),XMAXSAV(ICENTO),YMINSAV(ICENTO),
     &   YMAXSAV(ICENTO),Q2MINSAV(ICENTO),Q2MAXSAV(ICENTO),
     &   W2MINSAV(ICENTO),W2MAXSAV(ICENTO),PARIMAX(ICENTO),
     &   PPSAVE(ICENTO,3,4,5),PARICOR(ICENTO),INDEX,SIGMASAV(ICENTO),
     &   XMSIGMA,XSECT

*KEEP,LUDAT1.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUJETS.
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      SAVE /LUJETS/
*KEEP,KEYS.
        COMMON/CFREAD/SPACE(5000)
 	COMMON/MYKEYS/IEVAR,IF5CC,INEUT,IINTE,IFERM,IFLAT,ICOUN,
     &  REFIX,IMUDO,IKAT1,IKAT2,IKAT3,IKAT4,IKAT5,IKAT6,
     &  INEVT,IJAK1,IJAK2,IITDK,RPTAU,RXK0D,NTGR,IDIMUON,IGLU,
     &  IQDEN,LOME(2),IFILES,ICCHA,ISEED,IDSUBS,IONEM,EHAC,IPSEL,
     &  IHIST


*KEEP,zebra.

      PARAMETER (NNQ=1000000)
*
      DIMENSION LQ(NNQ),IQ(NNQ),Q(NNQ)
      EQUIVALENCE (Q(1),IQ(1),LQ(9),JSTRUC(8))
      COMMON /QUEST/IQUEST(100)
      COMMON /XQSTOR/IXEVT,IFENCE(16),JGEEV,JSTRUC(99),JREFER(100),
     +DIV12(NNQ)
      COMMON /FZLUN/LUNFZ
      COMMON/MZIOALL/IOGENF

*KEND.
      PARAMETER(MAXV=20)
      PARAMETER(MAXT=50)
      REAL*4 VERT(MAXV,3)
      REAL*4 WA59(MAXT,9),VEC(8)
      REAL*4 TMPAR(3),VSTR(3),PM(3),PH(3)
      INTEGER NTRV(MAXV)
      INTEGER LBEA(MAXV)
      INTEGER NTBEAM(MAXV)
      INTEGER IPART(MAXV,MAXT)
      INTEGER LPART(MAXV,MAXT)
      INTEGER LSTRFROM(MAXV,MAXT)
      INTEGER LSTRDA(MAXV,MAXT)
      INTEGER KTRK(MAXV,MAXT)
      REAL*4 PLAB(MAXV,MAXT,3)
      REAL*4 UBUFT(MAXV,MAXT,7)
      INTEGER UNOUNO
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON/SBEAM/ PNUMBER,NEUTYPE,VECT(3),GKIN(3),MESTYPE,G4MES(4)
      REAL*4 BEAMFT(7)
      INTEGER*4 DALUAEF(200)
      COMMON/DALUA/DALUAEF

      CALL VZERO(LSTRFROM,MAXV*MAXT)
      CALL VZERO(LSTRDA,MAXV*MAXT)


      LOUT=12
      IF (UNOUNO.EQ.0) THEN
        UNOUNO=1
*       OPEN(UNIT=LOUT,ERR=9191,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      ENDIF

      EMIN=.000000
* 5 MEV CUTOFF,SORRY ABOUT THIS....


      IDEVT=IEVT
      LNU =LUTOGE(K(1,2))
      LLEP=LUTOGE(K(4,2))
      LCHA=0
      NTRK=0
      NVTX=0
      ENU=P(1,4)
      ILASTD=0

      CURVZ=-1.

      LCHA=0.
      CALL VZERO(DALUAEF,200)

      DO I=1,N

        IF (K(I,2).EQ.92) THEN
          STR=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
          DO IB=1,3
            VSTR(IB)=P(I,IB)/STR
          ENDDO
        ENDIF


        IF ((K(I,1).LE.10.OR.(K(I,1).GT.10.AND.
     +  K(I,4).GT.0.AND.V(I,5).NE.0).AND.P(I,4)
     +  .GT.EMIN).OR.I.EQ.4) THEN

          LTMP =LUTOGE(K(I,2))

*DETECT CHARM DECAY, KEEP INDEX IN ILASTD
          IF((ABS(K(I,2)).GT.400.AND.ABS(K(I,2)).LT.500)
     + .OR.ABS(K(I,2)).EQ.4122) THEN
            LCHA=LUTOGE(K(I,2))
            FTUPLE(IOF1+7)=FTUPLE(IOF1+7)+1
            ILASTD=I
          ENDIF

* ASSUME 1ST ORDER  THAT LUND DECAY VERTEX ARE ORDERED BY Z DIRECTION
          IF (V(I,3).NE.CURVZ) THEN
* HAVE WE FOUND A NEW DECAY VERTEX ?
            ISFOUND=1
* TO BE SURE, CHECK ALL THE VERTEX FOUND
            DO JJ=1,NVTX
              IF (V(I,3).EQ.VERT(JJ,3)) THEN
* SORRY , I WAS WRONG
                ISFOUND=0
              ENDIF
            END DO
            IF (ISFOUND.EQ.1) THEN
* THIS IS A NEW VERTEX....
              NVTX=NVTX+1
* UPDATE CURVZ
              CURVZ=V(I,3)
* FILL VERT INFORMATIONS
              DO J=1,3
                VERT(NVTX,J)=V(I,J)
              END DO
* AND TAKE A "CANDIDATE" PARENT
              IPARENT=K(I,3)

   10         CONTINUE

              IF (IPARENT.NE.0) THEN
* IF IT'S NOT DIRECTLY COMING FROM THE BEAM,
                IDT=K(IPARENT,2)
* NOW WE HAVE THE ID OF THE PARENT CANDIDATE
* AND WE ASK IT TO BE A TRACK (I.E. CTAU>0) WHICH HAS ACTUALLY DECAYED
                IF (K(IPARENT,4).GT.0..AND.
     +         V(IPARENT,5).NE.0.AND.K(IPARENT,1).GT.10) THEN
                  LBEA(NVTX)=LUTOGE(IDT)
                ELSE
* GO UP AND LOOK FOR IT
                  IPARENT=K(IPARENT,3)
                  GOTO 10
                ENDIF


              ELSE
                LBEA(NVTX)=4
              ENDIF
* DEBUG
              IF (IEVT.GE.LOME(1).AND.IEVT.LE.LOME(2)) THEN
                WRITE(*,*)' VERTEX PARENT FOR PARTICLE NUMBER',
     +          I,' = ', RENT,IPA
              ENDIF
              NTBEAM(NVTX)=IPARENT
            ENDIF
          ENDIF

        ENDIF
      END DO

* COUNT TRACKS IN EACH VERTEX
* LOOP ON VERTEXES
      DO J=1,NVTX
        NTRV(J)=0
        DO I=1,N
          IF (K(I,2).EQ.92) ISTRINGA=I
          IF ((K(I,1).LE.10.OR.(K(I,1).GT.10.AND.
     +    K(I,4).GT.0.AND.V(I,5).NE.0)).AND.P(I,
     +    4).GT.EMIN.OR.I.EQ.4) THEN
            IF (V(I,3).EQ.VERT(J,3)) THEN
* COUNT TOTAL TRACKS AND TRACKS FROM THAT VERTEX
              NTRV(J)=NTRV(J)+1
              NTRK=NTRK+1
* PUT IN DALUAEF ARRAY THE EFICASS POSITION FOR THAT PARTICLE
              DALUAEF(I)=NTRK
              KTRK(J,NTRV(J))=NTRK
* TAKE GEANT PARTICLE CODE
              IPART(J,NTRV(J))=LUTOGE(K(I,2))
              LPART(J,NTRV(J))=K(I,2)
              LSTRFROM(J,NTRV(J))=JSTRO(I)
              IF (K(I,3).EQ.ISTRINGA) THEN
                LSTRDA(J,NTRV(J))=1
              ENDIF

* TAKE PARTICLE MOMENTUM
              DO KK=1,3
                PLAB(J,NTRV(J),KK)=P(I,KK)
              END DO
* FILL PT, P, E
              UBUFT(J,NTRV(J),1)=SQRT(P(I,1)**2+P(I,2)**2)
              UBUFT(J,NTRV(J),2)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
              UBUFT(J,NTRV(J),3)=P(I,4)
              PPP=UBUFT(J,NTRV(J),2)
* BETAI=PI/P
              IF (PPP.GT.0) THEN
                DO KK=1,3
                  UBUFT(J,NTRV(J),3+KK)=P(I,KK)/PPP
                END DO
              ELSE
                DO KK=1,3
                  UBUFT(J,NTRV(J),3+KK)=0.
                END DO
              ENDIF
* FILL 0 THE TRACK LENGTH
              UBUFT(J,NTRV(J),7)=0.
* CALCULATE TRACK LENGTH: GAMMA BETA C TAU
              IF (K(I,4).NE.0.AND.V(I,5).NE.0) THEN
                IF (P(I,5).GT.0) THEN
                  UBUFT(J,NTRV(J),7)=V(I,5)*PPP/P(I,5)/10.
                ELSE
                  UBUFT(J,NTRV(J),7)=0.
                ENDIF
              ENDIF
* VERY IMPORTANT FOR NEW JETSET CUT IN DECAYS
              IF (K(I,4).EQ.0)  UBUFT(J,NTRV(J),7)=0.
              IF (IEVT.GE.LOME(1).AND.IEVT.LE.LOME(2)) THEN
                WRITE(*,*)'TRACK ',I,'FROM VTX ',J,' PARENTID=',LBEA(J)
              ENDIF
            ENDIF
          ENDIF
        END DO
      END DO


* SIMULATE OUTPUT

      TEMULX=RNDMM(ISEED)*140.-70.
      TEMULY=RNDMM(ISEED)*140.-70.
      FTUPLE(IOF2+1)=TEMULX
      FTUPLE(IOF2+2)=TEMULY
      XBARI=0
      YBARI=0
      XMBARI=0
      YMBARI=0
      XBARIA=0
      YBARIA=0
      PLAA=0

      NMU=0
      NHAD=0

      CALL MZBOOK(IXEVT,JGEGE,JGEEV,-3,'GEGE',NVTX,NVTX,7,3,0)
      CALL MZBOOK(IXEVT,JGEBE,JGEEV,-4,'GEBE',0,0,14,3,0)
      JGEGE=LQ(JGEEV-3)
      Q(JGEGE+1)=IDEVT
      Q(JGEGE+2)=LNU
      Q(JGEGE+3)=LLEP
      Q(JGEGE+4)=LCHA
      Q(JGEGE+5)=NVTX
      Q(JGEGE+6)=NTRK
      Q(JGEGE+7)=ENU

      WRITE(LOUT)IDEVT,LNU,LLEP,LCHA,NVTX+1,NTRK,ENU
      IF (IEVT.GE.LOME(1).AND.IEVT.LE.LOME(2)) THEN
        WRITE(*,*)'OUTPUT-EVT',IDEVT,LNU,LLEP,LCHA,NVTX+1,NTRK,ENU
      ENDIF

      JGEBE=LQ(JGEEV-4)
      DO IJU=1,3
        Q(JGEBE+IJU)=VECT(IJU)
        Q(JGEBE+IJU+3)=GKIN(IJU)
      ENDDO
      Q(JGEBE+7)=NEUTYPE
      Q(JGEBE+8)=XSECT
      Q(JGEBE+9)=PNUMBER

      BEAMFT(1)=MESTYPE
      BEAMFT(2)=G4MES(1)
      BEAMFT(3)=G4MES(2)
      BEAMFT(4)=G4MES(3)
      BEAMFT(5)=G4MES(4)
      BEAMFT(6)=XSECT
      BEAMFT(7)=PNUMBER

      DO IRI=1,7
        Q(JGEBE+6+IRI)=BEAMFT(IRI)
      ENDDO

      WRITE(LOUT)1,1,0,0,(VECT(KK),KK=1,3)
      WRITE(LOUT)1,1,1,NEUTYPE,(GKIN(KK),KK=1,3), (BEAMFT(KK),KK=1,7)
      IF (IEVT.GE.LOME(1).AND.IEVT.LE.LOME(2)) THEN
        WRITE(*,*)'OUTPUT-VERT-NEUTRINO ',1,1,0,0,(VECT(KK),KK=1,3)
        WRITE(*,*)'OUTPUT-PART-NEUTRINO ' ,1,1,1,NEUTYPE,(GKIN(KK),KK=
     +  1,3), (BEAMFT(KK),KK=1,7)
      ENDIF

      DO J=1,NVTX
        ITRUENTBEAM=DALUAEF(NTBEAM(J))
        IF (NTRV(J).EQ.0) THEN
          WRITE(*,*)'DANGER! VERTEX WITH 0 TRACKS:NVTX=',NVTX
          WRITE(*,*)'VERT,IDEVT=',VERT(J,3),IDEVT
          CALL LULIST(3)
          WRITE(*,*)'END DANGER! '
        ENDIF
        DO JJ=1,3
          TMPAR(JJ)=VERT(J,JJ)/10.
        END DO
        JGEGE=LQ(JGEEV-3)
        CALL MZBOOK(IXEVT,JGEVT,JGEGE,-J,'GEVT',NTRV(J),NTRV(J),7,3,0)
        Q(JGEVT+1)=J
        Q(JGEVT+2)=NTRV(J)
        Q(JGEVT+3)=ITRUENTBEAM
        Q(JGEVT+4)=LBEA(J)
        DO JJ=1,3
          Q(JGEVT+4+JJ)=TMPAR(JJ)
        END DO
        WRITE(LOUT)J+1,NTRV(J),ITRUENTBEAM,LBEA(J),
     +  (TMPAR(JJ),JJ=1,3)
        IF  (IEVT.GE.LOME(1).AND.IEVT.LE.LOME(2)) THEN
          WRITE(*,*)'OUTPUT-VERT=',J,NTRV(J),ITRUENTBEAM,LBEA(J),
     +    (TMPAR(JJ),JJ=1,3)
        ENDIF

* NOW LOOP ON SINGLE VERTEX HERE
        DO I=1,NTRV(J)
          JGEGE=LQ(JGEEV-3)
          JGEVT=LQ(JGEGE-J)
          CALL MZBOOK(IXEVT,JGETR,JGEVT,-I,'GETR',0,0,14,3,0)
          Q(JGETR+1)=KTRK(J,I)
          Q(JGETR+2)=I
          Q(JGETR+3)=J
          Q(JGETR+4)=IPART(J,I)
          DO KK=1,3
            Q(JGETR+4+KK)=PLAB(J,I,KK)
          ENDDO
          DO KK=1,7
            Q(JGETR+7+KK)=UBUFT(J,I,KK)
          ENDDO
          WRITE(LOUT)KTRK(J,I)+1,I,J+1,IPART(J,I),(PLAB(J,I,KK),KK=1,3),
     +    (UBUFT(J,I,KK),KK=1,7)
          IF (IEVT.GE.LOME(1).AND.IEVT.LE.LOME(2)) THEN
            WRITE(*,*)'OUTPUT-TRACK=',KTRK(J,I),I,J,IPART(J,I), (PLAB(J
     +      ,I,KK),KK=1,3), (UBUFT(J,I,KK),KK=1,7)
          ENDIF
          IF (IPART(J,I).EQ.34.OR.
     +    IPART(J,I).EQ.33) THEN
            FTUPLE(80)=UBUFT(J,I,7)
          ENDIF

          IF (IPART(J,I).EQ.35.OR.IPART(J,I).EQ.36.OR.
     +    IPART(J,I).EQ.37.OR.IPART(J,I).EQ.38.OR.
     +    IPART(J,I).EQ.39.OR.IPART(J,I).EQ.40.OR.
     +    IPART(J,I).EQ.41.OR.IPART(J,I).EQ.42.OR.
     +    IPART(J,I).EQ.43.OR.IPART(J,I).EQ.44.OR.
     +    IPART(J,I).EQ.45.OR.IPART(J,I).EQ.46
     +    ) THEN
            FTUPLE(IOF1+10)=UBUFT(J,I,7)
            FTUPLE(IOF1+11)=IPART(J,I)
            FTUPLE(IOF1+12)=UBUFT(J,I,2)
            FTUPLE(IOF1+13)=UBUFT(J,I,3)
            FTUPLE(IOF1+15)=PLAB(J,I,1)
            FTUPLE(IOF1+16)=PLAB(J,I,2)
            FTUPLE(IOF1+17)=PLAB(J,I,3)
          ENDIF

          IF(UBUFT(J,I,7).EQ.0) THEN

            PLA=SQRT(PLAB(J,I,1)**2+PLAB(J,I,2)**2+PLAB(J,I,3)**2)

* CHECK PT IN FRAGMENTATION


* CONSTRUCT Z VARIABLE

            ETHIS=UBUFT(J,I,3)

            IF  (LSTRFROM(J,I).GT.0) THEN
              PTSR=PLAB(J,I,1)*VSTR(1)+ PLAB(J,I,2)*VSTR(2)+PLAB(J,I,3)
     +        *VSTR(3)
              PTSTRFIN=SQRT(PLA**2-PTSR**2)
              THCAS=RNDMM(ISEED)*3.141*2
              CALL HFILL(1011,PTSTRFIN*COS(THCAS),0.,1.)
              CALL HFILL(1011,PTSTRFIN*SIN(THCAS),0.,1.)
            ENDIF
            IF  (LSTRDA(J,I).GT.0) THEN
              PTSR=PLAB(J,I,1)*VSTR(1)+ PLAB(J,I,2)*VSTR(2)+PLAB(J,I,3)
     +        *VSTR(3)
              PTSTRFIN=SQRT(PLA**2-PTSR**2)
              THCAS=RNDMM(ISEED)*3.141*2

* CHECK ON D(Z) AS IN ALLASIA ET AL., Z PHYS C 24, 119-131 (1984)
* CHECK ON D(X_F) AS IN ALLEN ET AL., NUCL PHYS. B214 (1983) 369-391
              ZZZ=ETHIS/U
              XF=2*PTSR/SQRT(W2)
              ICHG=LUCHGE(LPART(J,I))
              IF (IPART(J,I).EQ.8) THEN
                IF (W2.GT.3) THEN
                  CALL HFILL(1015,XF,0.,1.)
                ENDIF
              ENDIF
              IF (IPART(J,I).EQ.9) THEN
                IF (W2.GT.3) THEN
                  CALL HFILL(1016,XF,0.,1.)
                ENDIF
              ENDIF
              IF (W2.GT.5.AND.Q2.GT.1) THEN
                IF (ICHG.GT.0) THEN
                  CALL HFILL(1013,ZZZ,0.,1.)
                ELSE
                  CALL HFILL(1014,ZZZ,0.,1.)
                ENDIF
              ENDIF
              CALL HFILL(1012,PTSTRFIN*COS(THCAS),0.,1.)
              CALL HFILL(1012,PTSTRFIN*SIN(THCAS),0.,1.)
            ENDIF
          ENDIF
        END DO
      END DO


      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LAZIMU(XP,ZP)

C...CHOOSE AZIMUTHAL ANGLE (PHI) ACCORDING TO QCD MATRIX ELEMENTS.
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)

      J=LST(24)-1
      SGN=SIGN(1.,2.5-LST(24))
      IFL=LST(25)
      I=IABS(IFL)
      IH=1
      IF(LST(30).EQ.1) IH=2

      IF(LST(23).EQ.2) THEN
        A=PARI(24)*DQCD(0,J,1,XP,ZP,Y)+PARI(25)*DQCD(0,J,2,XP,ZP,Y)
     &  -LST(30)*ISIGN(1,IFL)*PARI(26)*DQCD(0,J,3,XP,ZP,Y)
        B=DQCD(1,J,1,XP,ZP,Y)
     &  +SGN*LST(30)*ISIGN(1,IFL)*DQCD(1,J,3,XP,ZP,Y)
        C=DQCD(2,J,1,XP,ZP,Y)
      ELSE
        A=(EWQC(1,IH,I)+EWQC(2,IH,I))*(PARI(24)*DQCD(0,J,1,XP,ZP,Y)+
     &    PARI(25)*DQCD(0,J,2,XP,ZP,Y))
     &    -LST(30)*ISIGN(1,IFL)*(EWQC(1,IH,I)-EWQC(2,IH,I))
     &    *PARI(26)*DQCD(0,J,3,XP,ZP,Y)
        B=(EWQC(1,IH,I)+EWQC(2,IH,I))*DQCD(1,J,1,XP,ZP,Y)
     &    +SGN*LST(30)*ISIGN(1,IFL)*(EWQC(1,IH,I)-EWQC(2,IH,I))
     &    *DQCD(1,J,3,XP,ZP,Y)
        C=(EWQC(1,IH,I)+EWQC(2,IH,I))*DQCD(2,J,1,XP,ZP,Y)
      ENDIF
      PHIMAX=ABS(A)+ABS(B)+ABS(C)
   10 PHI=6.2832*RLU(0)
      IF(A+B*COS(PHI)+C*COS(2.*PHI).LT.RLU(0)*PHIMAX) GOTO 10
      CALL LUROBO(0.,PHI,0.,0.,0.)

      RETURN
      END
*CMZ :  1.01/33 11/07/95  18.31.20  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.35.13  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 29/07/94  17.22.12  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 20/07/94  12.26.37  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE JETTA
*KEEP,CDEBEAM.
C--
      COMMON /CONTRO/  BINIT,LUNB,NPNEUT,NPANTI,CPNORM,XPSOUR,SIGDIV
      COMMON /FLUXES/  FLUXD(80000),WEIGHT(8),SPECD(800),SPECN(800)
      COMMON /INPUT/   IALL(80000),NCOUNT(8)
      LOGICAL          BINIT


C-
*KEEP,POLAR.
C--
      COMMON /POLARIZ/POL(4000,3)
      REAL POLARX(4)
*KEEP,LUDAT1.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUJETS.
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      SAVE /LUJETS/
*KEEP,JETTA.
C--
         PARAMETER (ICENTO=100)
         PARAMETER (ISIZ=93)
         PARAMETER (IOF1=32)
         PARAMETER (IOF2=83)
         PARAMETER (LUX_LEVEL=4)
         INTEGER*4 JTAU(100),JPRI(100),JSTRO(100)
         REAL*4 FTUPLE(ISIZ)
         COMMON/JETTAGL/JTAU,JPRI,JSTRO
         COMMON/NTUPLA/FTUPLE,ISFIRST
         COMMON/BEAM/SPEC(ICENTO)
         COMMON /MAXSPEC/RMAXSPEC,RINTSPEC
         COMMON/SAV/XMINSAV(ICENTO),XMAXSAV(ICENTO),YMINSAV(ICENTO),
     &   YMAXSAV(ICENTO),Q2MINSAV(ICENTO),Q2MAXSAV(ICENTO),
     &   W2MINSAV(ICENTO),W2MAXSAV(ICENTO),PARIMAX(ICENTO),
     &   PPSAVE(ICENTO,3,4,5),PARICOR(ICENTO),INDEX,SIGMASAV(ICENTO),
     &   XMSIGMA,XSECT

*KEEP,FOREFI.
C--
	INTEGER*4 IEVT
         COMMON/FOREFICASS/IEVT


*KEND.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U

* FIRST, CALCULATE MEAN CHARGED MULTIPLICITIES

      ICHGPRI=0
      INEUPRI=0
      ICHGTAU=0
      INEUTAU=0
      ICHGCHA=0

      ENETAUCHG=0
      ENETAUNEU=0
      ENEPRICHG=0
      ENEPRINEU=0


      DO I=1,N

        IF (K(I,2).EQ.92) THEN
          FTUPLE(31)=P(I,4)
          FTUPLE(32)=P(I,5)
*        WRITE(*,*)'W2 AND WSTR2:',W2,P(I,5)**2
        ENDIF

        IF(JTAU(I).EQ.1) THEN
          IF(LUCHGE(K(I,2)).NE.0) THEN
            ICHGTAU=ICHGTAU+1
            ENETAUCHG=ENETAUCHG+P(I,4)
          ELSE
            INEUTAU=INEUTAU+1
            ENETAUNEU=ENETAUNEU+P(I,4)
          ENDIF
        ENDIF

        IF(JPRI(I).EQ.1) THEN
          IF(LUCHGE(K(I,2)).NE.0) THEN
            ICHGPRI=ICHGPRI+1
            ENEPRICHG=ENEPRICHG+P(I,4)
          ELSE
            INEUPRI=INEUPRI+1
            ENEPRINEU=ENEPRINEU+P(I,4)
          ENDIF
        ENDIF

      END DO

*     WRITE(*,*)'ICHGPRI,ENEPRICHG=',ICHGPRI,ENEPRICHG

      FTUPLE(20)=ICHGTAU
      FTUPLE(21)=INEUTAU
      FTUPLE(22)=ICHGPRI
      FTUPLE(23)=INEUPRI
      FTUPLE(24)=ENETAUCHG
      FTUPLE(25)=ENETAUNEU
      FTUPLE(26)=ENEPRICHG
      FTUPLE(27)=ENEPRINEU

      RETURN
      END
*CMZ :          04/03/97  12.54.46  by  Unknown
*CMZ :  1.01/22 27/05/95  16.18.36  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.26  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      BLOCK DATA LEPTOD

C...GIVE SENSIBLE DEFAULT VALUES TO SWITCHES AND PARAMETERS.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LFLMIX/ CABIBO(4,4)
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON /LGRID/ NXX,NWW,XX(20),WW(15),PQG(20,15,3),PQQB(20,15,2),
     &QGMAX(20,15,3),QQBMAX(20,15,2),YCUT(20,15),XTOT(20,15),NP
      COMMON /FLGRID/ NFX,NFQ,XR(2),QR(2),FLQT(41,16),FLGT(41,16),
     &FLMT(41,16)
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON /LMINUI/ XKIN(4),UKIN(4),WKIN(4),AIN(4),BIN(4),
     &MAXFIN,RELUP,RELERR,RELER2,FCNMAX
      COMMON /LMINUC/ NAMKIN(4),NAM(30)
      CHARACTER*10 NAMKIN,NAM

C...LEPTOU: CUTS, BASIC SWITCHES AND PARAMETERS.
      DATA CUT/1.E-04,1.,0.,1.,4.,1.E+08,5.,1.E+08,.1,1.E+08,.1,1.E+08,
     &0.,3.1416/
      DATA LST/0,1,5,1,3,1,1,12,5,1,0,4,5,1,1,1,0,2,3,21*0/
      DATA PARL/1.,1.,0.44,0.75,0.226,0.,0.,0.015,2.,0.,0.01,4.,
     &0.001,0.44,0.01,7.29735E-03,1.16637E-05,0.044,0.03,1.,10*0./
C...INTERNALLY USED VARIABLES.
      DATA PARI/40*0./
      DATA QC/-.33333,.66667,-.33333,.66667,-.33333,.66667,
     &        -.33333,.66667/
      DATA CABIBO/.95,.05,2*0.,.05,.948,.002,2*0.,.002,.998,4*0.,1./
      DATA OPTX/1.,3*0./,OPTY/1.,3*0./,OPTQ2/1.,3*0./,OPTW2/1.,3*0./
      DATA NXX,NWW/20,15/
      DATA PQG,PQQB,QGMAX,QQBMAX/3000*0./,YCUT/300*0./,XTOT/300*0./
      DATA NFX,NFQ/41,16/,FLQT,FLGT,FLMT/1968*0./
      DATA XKIN/1.,2.,3.,4./,UKIN,WKIN,AIN,BIN/16*0./,MAXFIN/2000/
      DATA RELUP,RELERR,RELER2/0.1,0.05,0.05/
      DATA NAMKIN/'         X','          ','          ','          '/
      DATA IPY/
     1 0,     0,     2,     2,     6,     1,     1,     6,     3,     1,
     2 3,     1,     1,     2,     1,     1,     4,     1,     1,     1,
     3 0,     1,     1,     1,     1,     1,     1,     0,     0,     0,
     4 1,     2,     1,     1,    30,    33,     1,     1,     7,     0,
     5 0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     6 0,     0,     0,     1,   100,     0,     0,     0,     0,     0,
     7 0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     8 0,     0,     0,     0,     0,     0,     0,     0,     0,     0/
      DATA (PYPAR(I),I=1,40)/
     1   7.299E-03,   2.290E-01,   2.000E-01,   2.500E-01,   4.000E+00,
     1   1.000E+00,   4.400E-01,   4.400E-01,   7.500E-02,   0.000E+00,
     2   2.000E+00,   2.000E+00,   1.000E+00,   0.000E+00,   3.000E+00,
     2   1.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   1.000E+00,
     3   2.500E-01,   1.000E+00,   2.000E+00,   1.000E-03,   1.000E+00,
     3   1.000E+00,   1.000E+00,  -2.000E-02,  -1.000E-02,   0.000E+00,
     4   0.000E+00,   1.600E+00,   0.500E+00,   0.200E+00,   3.894E-01,
     4   1.000E+00,   3.300E-01,   6.600E-01,   0.000E+00,   1.000E+00/
      DATA (PYPAR(I),I=41,80)/
     5   2.260E+00,   1.000E+04,   1.000E-04,   0.000E+00,   0.000E+00,
     5   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     6   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     6   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     7   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     7   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     8   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,
     8   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00/
      DATA PYVAR/80*0./
      END
*CMZ :          06/03/97  15.28.28  by  Unknown
*CMZ :  1.01/50 20/03/96  12.38.57  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 23/09/94  12.10.13  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 19/08/94  11.05.13  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 20/07/94  12.12.49  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LEPTOX

C...SELECT PROCESS AND CHOOSE KINEMATICAL VARIABLES (X,Y; X,Q2; X,W2)
C...ACCORDING TO THE DIFFERENTIAL CROSS SECTION.

*KEEP,JETTA.
C--
         PARAMETER (ICENTO=100)
         PARAMETER (ISIZ=93)
         PARAMETER (IOF1=32)
         PARAMETER (IOF2=83)
         PARAMETER (LUX_LEVEL=4)
         INTEGER*4 JTAU(100),JPRI(100),JSTRO(100)
         REAL*4 FTUPLE(ISIZ)
         COMMON/JETTAGL/JTAU,JPRI,JSTRO
         COMMON/NTUPLA/FTUPLE,ISFIRST
         COMMON/BEAM/SPEC(ICENTO)
         COMMON /MAXSPEC/RMAXSPEC,RINTSPEC
         COMMON/SAV/XMINSAV(ICENTO),XMAXSAV(ICENTO),YMINSAV(ICENTO),
     &   YMAXSAV(ICENTO),Q2MINSAV(ICENTO),Q2MAXSAV(ICENTO),
     &   W2MINSAV(ICENTO),W2MAXSAV(ICENTO),PARIMAX(ICENTO),
     &   PPSAVE(ICENTO,3,4,5),PARICOR(ICENTO),INDEX,SIGMASAV(ICENTO),
     &   XMSIGMA,XSECT

*KEND.

      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON /FLINFO/ RFLQ,RFLG,RFLM,RFLT
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/WLIST/WW1,WW2,WW3,WW5
      DIMENSION PQH(17,2),PNT(2,2),XPQ(-6:6)
      DOUBLE PRECISION DARI27,DARI28
      DATA DARI27,DARI28/2*0.D0/
      DATA W2LOW,W2UPP,YLOW,YUPP,Q2LOW,Q2UPP/6*0./

      W2LOW=0
      W2UPP=0
      YLOW=0
      YUPP=0
      Q2LOW=0
      Q2UPP=0

      DO 30 IH=1,2
        DO 10 I=1,2
   10   PNT(I,IH)=0.
        DO 20 I=1,8
          EWQC(1,IH,I)=0.
   20   EWQC(2,IH,I)=0.
        DO 30 I=1,17
   30 PQH(I,IH)=0.
      DO 40 I=1,17
   40 PQ(I)=0.

      LST(21)=0
      NCUT=0
      RML=P(4,5)
      EE=P(1,4)
      RMM=P(2,5)
      HEL=LST(30)

      S=PARL(21)+PSAVE(3,1,5)**2+PSAVE(3,2,5)**2
      PM2=PSAVE(3,2,5)**2
      IF(LST(2).NE.1) THEN

* NOT SO USEFUL AS I THOUGHT
        AAA=0.5
        BBB=0.5
        IF (PARL(21).GT.0) THEN
          AAA= 0.5*(1.-RML**2/(PARL(21)*X)-2.* (RML*RMM/PARL(21))**2)
     +    /(1.+X*RMM**2/PARL(21))
          BBB=0.5*SQRT((1.-RML**2/(PARL(21)*X))**2-(2.*RML*RMM/PARL(21)
     +    )**2) /(1.+X*RMM**2/PARL(21))
        ELSE
          WRITE(*,*)'WARNING: 2PK==0'
        ENDIF

* EE-RMM OR EE+RMM ? THIS IS THE QUESTION.....

        TXLOW=RML**2/2./RMM/(EE+RMM)

        Q2LOW=MAX(Q2MIN,X*YMIN*S,(W2MIN-PM2)*X/MAX(1.-X,1.E-22))
        Q2UPP=MIN(Q2MAX,X*YMAX*S,(W2MAX-PM2)*X/MAX(1.-X,1.E-22))
        YLOWCMP=MAX(YMIN,Q2MIN/MAX(S*X,1.E-22), (W2MIN-PM2)/MAX(S*(1.-
     +  X),1.E-22))
        YUPPCMP=MIN(YMAX,Q2MAX/MAX(S*X,1.E-22), (W2MAX-PM2)/MAX(S*(1.-
     +  X),1.E-22))
        YLOW=MAX(YMIN,Q2MIN/MAX(S*X,1.E-22), (W2MIN-PM2)/MAX(S*(1.-X),
     +  1.E-22))
        YUPP=MIN(YMAX,Q2MAX/MAX(S*X,1.E-22), (W2MAX-PM2)/MAX(S*(1.-X),
     +  1.E-22))



        W2LOW=MAX(W2MIN,(1.-X)*YMIN*S+PM2,Q2MIN*(1.-X)/MAX(X,1.E-22)+
     +  PM2)
        W2UPP=MIN(W2MAX,(1.-X)*YMAX*S+PM2,Q2MAX*(1.-X)/MAX(X,1.E-22)+
     +  PM2)
        GOTO 70
      ENDIF

      IF(PARI(28).LT.0.5) THEN
C...FOR FIRST CALL, RESET DOUBLE PRECISION COUNTERS.
        DARI27=0.D0
        DARI28=0.D0
      ENDIF
   50 DARI28=DARI28+1.D0
      PARI(28)=DARI28
   60 CONTINUE
*      WRITE(*,*)'OPTX=',OPTX
*      WRITE(*,*)'OPTY=',OPTY
C...CHOOSE X ACCORDING TO THE DISTRIBUTION
C...HX(X) =  A + B/X + C/X**2 + D/X**3. IN DETAIL
C...HQ=OPTX(1)/(XMAX-XMIN) + 1/LN(XMAX/XMIN)*OPTX(2)/X
C...   +XMIN*XMAX/(XMAX-XMIN)*OPTX(3)/X**2
C...   +2*(XMIN*XMAX)**2/(XMAX**2-XMIN**2)*OPTX(4)/X**3
      WHICH=(OPTX(1)+OPTX(2)+OPTX(3)+OPTX(4))*RLU(0)
      IF(WHICH.LE.OPTX(1)) THEN
        X=XMIN+RLU(0)*(XMAX-XMIN)
      ELSEIF(WHICH.LE.(OPTX(1)+OPTX(2))) THEN
        X=XMIN*(XMAX/XMIN)**RLU(0)
      ELSEIF(WHICH.LE.(OPTX(1)+OPTX(2)+OPTX(3))) THEN
        X=XMIN*XMAX/(XMAX+RLU(0)*(XMIN-XMAX))
      ELSE
        X=SQRT((XMIN*XMAX)**2/(XMAX**2+RLU(0)*(XMIN**2-XMAX**2)))
      ENDIF
      IF(LST(31).EQ.1) THEN
C...CHOOSE Q**2 ACCORDING TO THE DISTRIBUTION
C...HQ(Q2) =  A + B/(Q2) + C/(Q2)**2 + D/(Q2)**3. IN DETAIL
C...HQ=OPTQ2(1)/(Q2MAX-Q2MIN) + 1/LN(Q2MAX/Q2MIN)*OPTQ2(2)/Q2
C...   +Q2MIN*Q2MAX/(Q2MAX-Q2MIN)*OPTQ2(3)/Q2**2
C...   +2*(Q2MIN*Q2MAX)**2/(Q2MAX**2-Q2MIN**2)*OPTQ2(4)/Q2**3
        Q2LOW=MAX(Q2MIN,X*YMIN*S,(W2MIN-PM2)*X/(1.-X))
        Q2UPP=MIN(Q2MAX,X*YMAX*S,(W2MAX-PM2)*X/(1.-X))
        IF(Q2UPP.LT.Q2LOW) GOTO 60
        WHICH=(OPTQ2(1)+OPTQ2(2)+OPTQ2(3)+OPTQ2(4))*RLU(0)
        IF(WHICH.LE.OPTQ2(1)) THEN
          Q2=Q2LOW+RLU(0)*(Q2UPP-Q2LOW)
        ELSEIF(WHICH.LE.(OPTQ2(1)+OPTQ2(2))) THEN
          Q2=Q2LOW*(Q2UPP/Q2LOW)**RLU(0)
        ELSEIF(WHICH.LE.(OPTQ2(1)+OPTQ2(2)+OPTQ2(3))) THEN
          Q2=Q2LOW*Q2UPP/(Q2UPP+RLU(0)*(Q2LOW-Q2UPP))
        ELSE
          Q2=SQRT((Q2LOW*Q2UPP)**2/(Q2UPP**2+RLU(0)*(Q2LOW**2-Q2UPP**2)
     +    ))
        ENDIF
        Y=Q2/(PARL(21)*X)
        IF(Y.LT.YMIN.OR.Y.GT.YMAX) GOTO 50
      ELSEIF(LST(31).EQ.2) THEN
C...CHOOSE Y ACCORDING TO THE DISTRIBUTION
C...HY(Y) =  A + B/Y + C/Y**2 + D/Y**3. IN DETAIL
C...HY=OPTY(1)/(YMAX-YMIN) + 1/LN(YMAX/YMIN)*OPTY(2)/Y
C...   +YMIN*YMAX/(YMAX-YMIN)*OPTY(3)/Y**2
C...   +2*(YMIN*YMAX)**2/(YMAX**2-YMIN**2)*OPTY(4)/Y**3
        YLOW=MAX(YMIN,Q2MIN/(S*X),(W2MIN-PM2)/(S*(1.-X)))
        YUPP=MIN(YMAX,Q2MAX/(S*X),(W2MAX-PM2)/(S*(1.-X)))
        IF(YUPP.LT.YLOW) GOTO 60
        WHICH=(OPTY(1)+OPTY(2)+OPTY(3)+OPTY(4))*RLU(0)
        IF(WHICH.LE.OPTY(1)) THEN
          Y=YLOW+RLU(0)*(YUPP-YLOW)
        ELSEIF(WHICH.LE.(OPTY(1)+OPTY(2))) THEN
          Y=YLOW*(YUPP/YLOW)**RLU(0)
        ELSEIF(WHICH.LE.(OPTY(1)+OPTY(2)+OPTY(3))) THEN
          Y=YLOW*YUPP/(YUPP+RLU(0)*(YUPP-YLOW))
        ELSE
          Y=SQRT((YLOW*YUPP)**2/(YUPP**2+RLU(0)*(YLOW**2-YUPP**2)))
        ENDIF
        Q2=X*Y*PARL(21)
        IF(Q2.LT.Q2MIN.OR.Q2.GT.Q2MAX) GOTO 50
      ELSEIF(LST(31).EQ.3) THEN
C...CHOOSE W**2 ACCORDING TO THE DISTRIBUTION
C...HW(W2) =  A + B/(W2) + C/(W2)**2 + D/(W2)**3. IN DETAIL
C...HW=OPTW2(1)/(W2MAX-W2MIN) + 1/LN(W2MAX/W2MIN)*OPTW2(2)/W2
C...   +W2MIN*W2MAX/(W2MAX-W2MIN)*OPTW2(3)/W2**2
C...   +2*(W2MIN*W2MAX)**2/(W2MAX**2-W2MIN**2)*OPTW2(4)/W2**3
        W2LOW=MAX(W2MIN,(1.-X)*YMIN*S+PM2,Q2MIN*(1.-X)/X+PM2)
        W2UPP=MIN(W2MAX,(1.-X)*YMAX*S+PM2,Q2MAX*(1.-X)/X+PM2)
        IF(W2UPP.LT.W2LOW) GOTO 60
        WHICH=(OPTW2(1)+OPTW2(2)+OPTW2(3)+OPTW2(4))*RLU(0)
        IF(WHICH.LE.OPTW2(1)) THEN
          W2=W2LOW+RLU(0)*(W2UPP-W2LOW)
        ELSEIF(WHICH.LE.(OPTW2(1)+OPTW2(2))) THEN
          W2=W2LOW*(W2UPP/W2LOW)**RLU(0)
        ELSEIF(WHICH.LE.(OPTW2(1)+OPTW2(2)+OPTW2(3))) THEN
          W2=W2LOW*W2UPP/(W2UPP+RLU(0)*(W2LOW-W2UPP))
        ELSE
          W2=SQRT((W2LOW*W2UPP)**2/(W2UPP**2+RLU(0)*(W2LOW**2-W2UPP**2)
     +    ))
        ENDIF
        Y=(W2-P(2,5)**2)/((1.-X)*PARL(21))
        Q2=X*Y*PARL(21)
        IF(Y.LT.YMIN.OR.Y.GT.YMAX) GOTO 50
        IF(Q2.LT.Q2MIN.OR.Q2.GT.Q2MAX) GOTO 50
      ENDIF
c          write(*,*) x,y

             aa=lkinem(lst(2))
c          write(*,*) aa,lst(2)
   70 IF(LKINEM(LST(2)).NE.0) THEN
        NCUT=NCUT+1
        IF(LST(2).EQ.1) THEN
          IF(NCUT.LE.9999) GOTO 50
          IF(LST(3).GE.1) then
          do kk=1,14
c           write(*,*) kk,cut(kk)
          end do
c          write(*,*) lst(11), lst(12),lst(13),lst(14)
c          print*,'yes here'
          WRITE(6,10000)
          end if
        ENDIF
        LST(21)=1
        RETURN
      ENDIF

c      print*,' survived in leptox'
      PARI(24)=(1.+(1.-Y)**2)/2.
      PARI(25)=1.-Y
      PARI(26)=(1.-(1.-Y)**2)/2.
      CALL LNSTRF(X,Q2,XPQ)
C...LEPTON HELICITY STATE, ONLY ONE CONTRIBUTES IN SOME CASES.
      IH=1
      IF(PARL(6).GT.+0.99) IH=2
   80 LST(30)=SIGN(1.,IH-1.5)
      PQH(17,IH)=0.
      PNT(1,IH)=0.
      PNT(2,IH)=0.
      IF(LST(23).EQ.2) THEN
C...CHARGED CURRENT: ZERO CROSS-SECTION FOR ONE HELICITY STATE.
        IF(KSAVE(1).LT.0.AND.IH.EQ.1
     +  .OR.KSAVE(1).GT.0.AND.IH.EQ.2) GOTO 110
*LST(30)=- O +1 A SECONDA DI HELICITA'  LEFT/RIGHT
        YQ=PARI(24)-LST(30)*PARI(26)
        YQB=PARI(24)+LST(30)*PARI(26)
        IF(PARI(11).GT.1.E-06) THEN
          IF(K(3,2).LT.0) THEN
            PNT(1,IH)=(1.-PARI(11))*PARI(13)*YQ
            PNT(2,IH)=PARI(11)*PARI(12)*YQ
          ELSE
            PNT(1,IH)=(1.-PARI(11))*PARI(12)*YQ
            PNT(2,IH)=PARI(11)*PARI(13)*YQ
          ENDIF
        ENDIF
        DO 90  I=1,LST(12)
          IF(K(3,2)*QC(I).LT.0) THEN
            PQH(I,IH)=XPQ(I)*YQ
          ELSE
            PQH(I+LST(12),IH)=XPQ(-I)*YQB
          ENDIF
   90   CONTINUE
      ELSE
C...NEUTRAL CURRENT: ELECTROMAGNETIC OR WEAK OR BOTH WITH INTERFERENCE.
        GFQ2=Q2/(PMAS(23,1)**2+Q2)*SQRT(2.)*PARL(17)*PMAS(23,1)**2/
     +  (3.1415927*PARL(16))
C...CORRECTION TO OBTAIN Q**2 DEPENDENT ALPHA-EM, IF DESIRED.
        AEMCOR=1.
        IF(LST(18).GE.2) AEMCOR=ULALEM(Q2)/PARL(16)
        II=3-IH
        ZLEP=ZL(IH,ILEP+2*INU)
        DO 100 I=1,MAX(LST(12),LST(13))
          A=(-IG*QC(I)*AEMCOR+IZ*GFQ2*ZLEP*ZQ(IH,I))**2
          B=(-IG*QC(I)*AEMCOR+IZ*GFQ2*ZLEP*ZQ(II,I))**2
C...SAVE HELICITY-DEPENDENT ELECTROWEAK QUARK COUPLINGS FOR LATER USE.
          EWQC(1,IH,I)=A
          EWQC(2,IH,I)=B
          IF(I.GT.LST(12)) GOTO 100
          FYQ=(A+B)*PARI(24)+(A-B)*PARI(26)
          PQH(I,IH)=XPQ(I)*FYQ
          IF(I.LE.2.AND.PARI(11).GT.1.E-06) THEN
            PNT(1,IH)=PNT(1,IH)+(1.-PARI(11))*PARI(11+I)*FYQ
            PNT(2,IH)=PNT(2,IH)+PARI(11)*PARI(14-I)*FYQ
          ENDIF
          PQH(I+LST(12),IH)=XPQ(-I)*((A+B)*PARI(24)-(A-B)*PARI(26))
  100   CONTINUE
      ENDIF
  110 CONTINUE
      DO 120 I=1,LST(12)
  120 PQH(17,IH)=PQH(17,IH)+PQH(I,IH)+PQH(I+LST(12),IH)

      IF(ABS(PARL(6)).LT.0.99.AND.IH.EQ.1) THEN
        IH=2
        GOTO 80
      ENDIF

      IF (LST(32).NE.0.AND.LST(23).EQ.2) THEN

        F1=0.
        F5=0.
        F2CC=0.
        F3CC=0.
*              D      U~      S      C~
        F2CC=XPQ(1)+XPQ(-2)+XPQ(3)+XPQ(-4)
        F3CC=(XPQ(1)-XPQ(-2)+XPQ(3)-XPQ(-4))/X



        IF(X.NE.0) THEN
          F5=F2CC/X
          F1=F2CC*0.5/X
        ELSE
          WRITE(*,*)'WARNING:X==0'
        ENDIF

* NOW CALCULATE FULL CROS SECTION, EXCEPT
* FOR COMFAC AND PARL(19) FACTORS
* KEEPING IN ACCOUNT BOTH NUCLEON AND LEPTON MASS

        RML=P(4,5)
        EE=P(1,4)
        RMM=P(2,5)
        HEL=LST(30)


        A1=( X*Y + RML**2 /PARL(21)  ) *Y
        A2=(1.-Y)- ( RMM**2*X*Y/PARL(21) + (RML*RMM/PARL(21))**2 )
        A3=( X*Y*(1.-Y/2.) - RML**2/(2.*PARL(21))*Y )
        A5=-RML**2/PARL(21)


        RML=0.
        RMM=0.0000001

        C1=( X*Y + RML**2 /PARL(21)  ) *Y
        C2=(1.-Y)- ( RMM**2*X*Y/PARL(21) + (RML*RMM/PARL(21))**2 )
        C3=( X*Y*(1.-Y/2.) - RML**2/(2.*PARL(21))*Y )
        C5=-RML**2/PARL(21)





*        WRITE(*,*)'A1,A2,A3,A5=',A1,A2,A3,A5

        PP=A1*F1+A2*F2CC+A3*F3CC+A5*F5
        PPPP=C1*F1+C2*F2CC+C3*F3CC+C5*F5
        PPP=A1*F1+A2*F2CC+A3*F3CC
        P17=(1.-PARL(6))/2.*PQH(17,1)+(1.+PARL(6))/2.*PQH(17,2)


        P17OK=(1.-PARL(6))/2.*PPPP+(1.+PARL(6))/2.*PQH(17,2)
        IF( ABS(P17OK-P17).GT.0.00005) THEN
          WRITE(*,*)'TAU LEPTON X-SECTION WRONG',P17OK,P17
        ENDIF

        IF(PP.LT.0.) PP=0.
* BY HAND, BUT A GOOD THING IS TO CHECKIT...

        IF(X.LT.TXLOW) PP=0.


* NOW W CALCULATIONS FOR TAU POLARIZATION

        RMM=P(2,5)
        WW1=F2CC*(1./U+0.5/X/RMM)
        WW2=F2CC/U
        WW3=-F2CC/U/X*SQRT(1+ (Q2/U)**2 )
        WW5=F2CC/X/U
*        WRITE(*,*)'W1,W2,W3,W5=',WW1,WW2,WW3,WW5

*     WRITE(*,*)'3F=',PPP,' 5F=',PP,' STD=',P17OK,' P17=',P17
*     WRITE(*,*)'Y=',Y,' 5F=',PP,' X=',X,' P17=',P17
* CRUCIAL!!!
* HERE STAYS THE FINAL JUMP INTO THE DARKNESS.....

        PQH(17,1)=PP


      ENDIF


      FLQ=0.
      FLG=0.
      FLM=0.
      FLT=0.
      IF(LST(23).EQ.1.AND.LST(11).NE.0.AND.LST(2).NE.-3) THEN
C-CHECK: IF(LST(23).EQ.1.AND.LST(11).NE.0) THEN
        LQCD=MOD(LST(11),10)
        LTM=MOD(LST(11)/10,10)
        LHT=LST(11)/100
C...INCLUDE QCD, TARGET MASS AND/OR HIGHER TWIST CONTR. TO LONG. STR FCN
C...FL FROM INTERPOLATION.
        IF(LQCD.EQ.1.OR.LTM.EQ.1) CALL FLIPOL(FLQ,FLG,FLM)
C...EVENT SIMULATION: IF REQUESTED, GET FL BY EVENT-BY-EVENT INTEGRATION
        IF(LST(2).GT.0.AND.
     +  (LQCD.EQ.2.OR.LTM.EQ.2)) CALL FLINTG(FLQ,FLG,FLM)
        IF(LTM.GE.1.OR.LHT.GE.1) THEN
          F2EM=0.
          DO 130 I=1,LST(12)
  130     F2EM=F2EM+QC(I)**2*(XPQ(I)+XPQ(-I))
          IF(LTM.GE.1) FLM=FLM-2.*X**2*PSAVE(3,2,5)**2/Q2*F2EM
          IF(LHT.GE.1) FLT=8.*PARL(19)/Q2*F2EM
        ENDIF
        DO 140 IH=1,2
          PQH17=PQH(17,IH)
C...NOTE FACTOR 2 AT THE END, SINCE PQH(IH,17) CONTAINS OVERALL FACTOR 2
          PQH(17,IH)=PQH(17,IH)-Y**2*(FLQ+FLG+FLM+FLT)
          DO 140 I=1,16
  140   PQH(I,IH)=PQH(I,IH)*PQH(17,IH)/PQH17
      ENDIF

      DO 150 I=1,17
  150 PQ(I)=(1.-PARL(6))/2.*PQH(I,1)+(1.+PARL(6))/2.*PQH(I,2)

C...RELATIVE CONTRIBUTION FROM LONGITUDINAL STR. FCN. AND HIGHER TWIST.
      RFLQ=-Y**2*FLQ/PQ(17)
      RFLG=-Y**2*FLG/PQ(17)
      RFLM=-Y**2*FLM/PQ(17)
      RFLT=-Y**2*FLT/PQ(17)

C...COMMON FACTOR FOR MATRIX ELEMENTS.
      IF(LST(31).EQ.1) THEN
        IF(LST(23).EQ.2) THEN
          COMFAC=1./X/(1.+Q2/PMAS(24,1)**2)**2
        ELSE
          COMFAC=1./X/Q2**2
        ENDIF
      ELSEIF(LST(31).EQ.2) THEN
        IF(LST(23).EQ.2) THEN
          COMFAC=1./(1.+Q2/PMAS(24,1)**2)**2*PARL(21)
        ELSE
          COMFAC=1./Q2**2*PARL(21)
        ENDIF
      ELSEIF(LST(31).EQ.3) THEN
        IF(LST(23).EQ.2) THEN
          COMFAC=1./X/(1.+Q2/PMAS(24,1)**2)**2  * X/(1.-X)
        ELSE
          COMFAC=1./X/Q2**2 * X/(1.-X)
        ENDIF
      ENDIF
C-CHECK: MOVE CHANGE OF COMFAC TO BELOW??
C...PREPARE FOR Q2 WEIGHTING.
C     WEIGHT=1/Q2**2
      WEIGHT=1.D0
      COMFAC=COMFAC/WEIGHT

c      print*,'before cec lst(2)'
      IF(LST(2).LE.-2) RETURN
c          print*,' after'
      HX=OPTX(1)/(XMAX-XMIN) + 1./ALOG(XMAX/XMIN)*OPTX(2)/X
     ++XMIN*XMAX/(XMAX-XMIN)*OPTX(3)/X**2
     ++2*(XMIN*XMAX)**2/(XMAX**2-XMIN**2)*OPTX(4)/X**3
      XFACT=OPTX(1)+OPTX(2)+OPTX(3)+OPTX(4)
      IF(LST(31).EQ.1) THEN
        HQ2=OPTQ2(1)/(Q2UPP-Q2LOW)
     +  +1./ALOG(Q2UPP/Q2LOW)*OPTQ2(2)/Q2
     +  +Q2LOW*Q2UPP/(Q2UPP-Q2LOW)*OPTQ2(3)/Q2**2
     +  +2*(Q2LOW*Q2UPP)**2/(Q2UPP**2-Q2LOW**2)*OPTQ2(4)/Q2**3
        Q2FACT=OPTQ2(1)+OPTQ2(2)+OPTQ2(3)+OPTQ2(4)
        COMFAC=COMFAC*XFACT*Q2FACT/HX/HQ2
      ELSEIF(LST(31).EQ.2) THEN
        HY=OPTY(1)/(YUPP-YLOW)+1./ALOG(YUPP/YLOW)*OPTY(2)/Y
     +  +YLOW*YUPP/(YUPP-YLOW)*OPTY(3)/Y**2
     +  +2*(YLOW*YUPP)**2/(YUPP**2-YLOW**2)*OPTY(4)/Y**3
        YFACT=OPTY(1)+OPTY(2)+OPTY(3)+OPTY(4)
        COMFAC=COMFAC*XFACT*YFACT/HX/HY
      ELSEIF(LST(31).EQ.3) THEN
        HW2=OPTW2(1)/(W2UPP-W2LOW)
     +  +1./ALOG(W2UPP/W2LOW)*OPTW2(2)/W2
     +  +W2LOW*W2UPP/(W2UPP-W2LOW)*OPTW2(3)/W2**2
     +  +2*(W2LOW*W2UPP)**2/(W2UPP**2-W2LOW**2)*OPTW2(4)/W2**3
        W2FACT=OPTW2(1)+OPTW2(2)+OPTW2(3)+OPTW2(4)
        COMFAC=COMFAC*XFACT*W2FACT/HX/HW2
      ENDIF
      IF(LST(2).LE.0) RETURN

C-CHECK: MOVE CHANGE OF COMFAC TO HERE?
      SIGL=(1.-PARL(6))/2.*PQH(17,1)
      SIGR=(1.+PARL(6))/2.*PQH(17,2)
      SIGMA=SIGL+SIGR
      IF(LST(2).EQ.1) THEN
C...WHEN CHOSING (X,Y), REJECT ACCORDING TO MAXIMUM OF "CROSS-SECTION",
C...UPDATE CROSS-SECTION ESTIMATE.
        DARI27=DARI27+DBLE(SIGMA)*DBLE(COMFAC)*WEIGHT
        PARI(27)=DARI27
        VIOL=SIGMA*COMFAC/PARI(LST(23))
        IF(VIOL.GT.PARI(32)) THEN
          PARI(32)=VIOL
          IF(PARI(32).GT.1.) THEN
            PARI(LST(23))=PARI(LST(23))*PARI(32)
            IF(LST(3).GE.1) WRITE(6,10100) PARI(32),INT(PARI(30)+1),
     +      PARI(LST(23)),X,Y,Q2,W2
            PARI(32)=1.
          ENDIF
        ENDIF
        IF(VIOL.LT.RLU(0)) GOTO 50
        PARL(24)=PARI(31)*DARI27/DARI28
      ENDIF

      IF(ABS(PARL(6)).LT.0.99) THEN
C...CHOOSE HELICITY OF INCOMING LEPTON.
        IH=1
        IF(RLU(0)*SIGMA.GT.SIGL) IH=2
      ENDIF
      LST(30)=SIGN(1.,IH-1.5)

C...CHOOSE TARGET NUCLEON, PROTON OR NEUTRON.
      LST(22)=1
      K(2,2)=2212
      IF(PARI(11).GT.1.E-06) THEN
        IF(RLU(0).LT.(PARI(11)*(PQH(17,IH)-PNT(1,IH)-PNT(2,IH))+
     +  PNT(2,IH))/PQH(17,IH)) THEN
          LST(22)=2
          K(2,2)=2112
        ENDIF
      ENDIF
      RCROSS=PARI(31)*PQ(17)*COMFAC
      FTUPLE(1)=RCROSS
      FTUPLE(2)=X
      FTUPLE(3)=Y

c        print*,' end of leptox'

      RETURN
10000 FORMAT(' WARNING: LEPTOX IS LOOPING, CANNOT FIND ALLOWED ',
     +'PHASE SPACE POINT DUE TO CUTS,',/,
     +10X,'CHECK, IN PARTICULAR, CUT(11) TO CUT(14)')
10100 FORMAT(' WARNING: MAXIMUM VIOLATED BY A FACTOR ',F7.3,
     +' IN EVENT ',I7,/,' MAXIMUM INCREASED BY THIS FACTOR TO ',E12.3,
     +/,' POINT OF VIOLATION: X, Y, Q**2, W**2 = ',4G10.3)
      END
*CMZ :  1.01/50 29/02/96  12.09.02  by  Piero Zucchelli
*CMZ :  1.01/45 08/01/96  11.11.52  by  Piero Zucchelli
*CMZ :  1.01/43 15/12/95  18.01.11  by  Piero Zucchelli
*CMZ :  1.01/41 12/12/95  16.06.05  by  Piero Zucchelli
*CMZ :  1.01/40 08/12/95  16.56.00  by  Piero Zucchelli
*CMZ :  1.01/36 26/07/95  17.35.22  BY  PIERO ZUCCHELLI
*CMZ :  1.01/30 02/06/95  19.54.53  BY  PIERO ZUCCHELLI
*CMZ :  1.01/29 02/06/95  19.51.39  BY  PIERO ZUCCHELLI
*CMZ :  1.01/28 02/06/95  18.17.14  BY  PIERO ZUCCHELLI
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LFLAV(IFL,IFLR)
*KEEP,KEYS.
        COMMON/CFREAD/SPACE(5000)
 	COMMON/MYKEYS/IEVAR,IF5CC,INEUT,IINTE,IFERM,IFLAT,ICOUN,
     &  REFIX,IMUDO,IKAT1,IKAT2,IKAT3,IKAT4,IKAT5,IKAT6,
     &  INEVT,IJAK1,IJAK2,IITDK,RPTAU,RXK0D,NTGR,IDIMUON,IGLU,
     &  IQDEN,LOME(2),IFILES,ICCHA,ISEED,IDSUBS,IONEM,EHAC,IPSEL,
     &  IHIST


*KEEP,PERROR.
*-- AUTHOR :    PIERO ZUCCHELLI   01/09/94
        PARAMETER (CHARMSENS=10000)
 	COMMON/MYERR/ICRACK


*KEND.
C...CHOOSE FLAVOUR OF STRUCK QUARK AND THE
C...CORRESPONDING FLAVOUR OF THE TARGET REMNANT JET.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON /LFLMIX/ CABIBO(4,4)

      LST(21)=0
      IF(LST(24).EQ.3) THEN
        NFL=LST(13)
      ELSE
        NFL=LST(12)
      ENDIF

      IMAXDIMU=CHARMSENS

   10 R=RLU(0)*PQ(17)
      PSUB=0.
      DO 20 I=1,2*NFL
        IFL=I
        PSUB=PSUB+PQ(I)
        IF(R.LE.PSUB) GOTO 30
   20 CONTINUE
   30 CONTINUE
      IF(IFL.GT.NFL) IFL=NFL-IFL
      LST(25)=IFL
      IFLR=-IFL

      IF(LST(23).EQ.2) THEN
C...WEAK CHARGED CURRENT, CHANGE THE FLAVOUR OF THE STRUCK
C...QUARK USING GENERALIZED CABIBBO MIXING MATRIX.
        IFLA=IABS(IFL)
        J1=(IFLA+1)/2
        M1=MOD(IFLA,2)
        M2=MOD(IFLA+1,2)
        R=RLU(0)
        PSUB=0.
        DO 40  J=1,4
          J2=J
          PSUB=PSUB+CABIBO(M1*J2+M2*J1,M2*J2+M1*J1)
          IF(R.LT.PSUB) GOTO 50
   40   CONTINUE
   50   IFL=2*J2-M2
        IF(LST(25).LT.0) IFL=-IFL
      ENDIF

      IFLA=IABS(IFL)
      IFLRA=IABS(IFLR)

* PIEROZ PATCH FOR DIMUONS
      IF (IDIMUON.GE.1.AND.IMAXDIMU.EQ.0) THEN
        WRITE(*,*)'SKIPPING CHARM PRODUCTION: CRACK ALARM'
        ICRACK=1
        GOTO 3434
* THIS IS PROGRAMMING!!!
      ENDIF
      IF (IDIMUON.GE.1.AND.IFLA.NE.4) THEN
        RSMALL=RNDMM(ISEED)
        IF (RSMALL.GT.0.05) GOTO 10
      ENDIF
      IF(IFLA.GE.4.OR.IFLRA.GE.4) THEN
C...THRESHOLD FUNCTION FOR HEAVY QUARKS OF FLAVOUR IFLA AND IFLRA.
        IF(1.-(.938+PMAS(LUCOMP(IFLA),1)+PMAS(LUCOMP(IFLRA),1)
     +  +2.*PMAS(1,1))**2/W2.LT.RLU(0)) THEN
          IMAXDIMU=IMAXDIMU-1
          GOTO(10,60 ,60 ) LST(24)
        ENDIF
      ENDIF

 3434 CONTINUE
C...REMNANT FLAVOUR TAKEN CARE OF LATER FOR QQBAR EVENT AND ME+PS CASE
      IF(LST(24).EQ.3) RETURN
      IF(LST(8).GT.10.AND.LST(8).NE.19) RETURN

C...WITH LST(14)=0/1(DEFAULT) BARYON PRODUCTION FROM THE TARGET REMNANT
C...IS EXCLUDED/INCLUDED.
      IF(LST(14).EQ.0) RETURN
      IF(IFLR.EQ.-2) THEN
        IF(LST(22).EQ.1) THEN
          IFLR=2101
          IF(RLU(0).GT.PARL(4)) IFLR=2103
        ELSE
          IFLR=1103
        ENDIF
      ELSEIF(IFLR.EQ.-1) THEN
        IF(LST(22).EQ.1) THEN
          IFLR=2203
        ELSE
          IFLR=2101
          IF(RLU(0).GT.PARL(4)) IFLR=2103
        ENDIF
      ENDIF
      RETURN

   60 LST(21)=1
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 24/07/94  16.02.59  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 19/07/94  17.08.36  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LFRAME(IFR,IPH)

C...MAKE TRANSFORMATION FROM HADRONIC CM FRAME TO LAB FRAME.

      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
*KEEP,POLAR.
C--
      COMMON /POLARIZ/POL(4000,3)
      REAL POLARX(4)
*KEND.
      DOUBLE PRECISION DTHETA,DPHI,DBETA,DBETADBG(3)
      INTEGER IFR,IPH,IFRAME,IPHI

*      WRITE(*,*)'ENTER IFR,PHI,LST28,LST29',IFR,IPH,LST(28),LST(29)

      IFRAME=IFR
      IPHI=IPH
      IF(IFRAME.LT.1.OR.IFRAME.GT.4.OR.IPHI.LT.0.OR.IPHI.GT.1)
     +GOTO 100
      IF(IFRAME.EQ.1) IPHI=0
      N=N+1
      DO 10 J=1,5
   10 P(N,J)=PB(J)

   20 CONTINUE
      IF(IPHI.NE.LST(29)) THEN
        IFRAME=2
      ELSE
        IFRAME=IFR
      ENDIF
      IF((IFRAME.EQ.LST(28)).AND.(IPHI.EQ.LST(29))) THEN
        DO 30 J=1,5
   30   PB(J)=P(N,J)
        N=N-1
        RETURN
      ENDIF

*      WRITE(*,*)'IFR,PHI,LST28,LST29,LST6',
*     &IFRAME,IPHI,LST(28),LST(29),LST(6)
      GOTO(40 ,50 ,70 ,90 ), LST(28)
      GOTO 100

   40 IF(IFRAME.GE.2) THEN
        CALL LUDBRB(0,0,STHETA(2),SPHI(2),0.D0,0.D0,0.D0)

        CALL LUDBRB(0,0,0.,0.,DBETA(2,1),DBETA(2,2),DBETA(2,3))
        LST(28)=2
      ELSE
        GOTO 100
      ENDIF
      GOTO 20

   50 IF(LST(6).NE.0.AND.IPHI.NE.LST(29)) THEN
        CALL LUDBRB(0,0,0.,SIGN(PHIR,FLOAT(IPHI-LST(29))),0.D0,0.D0,
     +  0.D0)
        LST(29)=IPHI
      ENDIF

      IF(IFRAME.EQ.1) THEN
        CALL LUDBRB(0,0,0.,0.,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3))
        CALL LUDBRB(0,0,-STHETA(2),0.,0.D0,0.D0,0.D0)
        LST(28)=1
      ELSEIF(IFRAME.GE.3) THEN
        IF(LST(17).EQ.0) THEN
          CALL LUDBRB(0,0,0.,0.,0.D0,0.D0,DBETA(1,3))
          IF(PSAVE(3,1,3).LT.0.) THEN
            DO 60  I=1,N
   60       P(I,3)=-P(I,3)
          ENDIF
        ELSE
          IF(DBETA(1,3).EQ.0) THEN
            DO J=1,3
              DBETADBG(J)= (DBLE(PSAVE(3,1,J))+DBLE(PSAVE(3,2,J)))/
     +        (DBLE(PSAVE(3,1,4))+DBLE(PSAVE(3,2,4)))
            END DO
            CALL LUDBRB(0,0,STHETA(1),SPHI(1),0.D0,0.D0,0.D0)
            CALL LUDBRB(0,0,0.,0.,DBETADBG(1),DBETADBG(2),DBETADBG(3))
          ELSE
            CALL LUDBRB(0,0,STHETA(1),SPHI(1),0.D0,0.D0,0.D0)
            CALL LUDBRB(0,0,0.,0.,DBETA(1,1),DBETA(1,2),DBETA(1,3))
          ENDIF
        ENDIF
        LST(28)=3
      ENDIF
      GOTO 20

   70 IF(IFRAME.LE.2) THEN
        IF(LST(17).EQ.0) THEN
          IF(PSAVE(3,1,3).LT.0.) THEN
            DO 80  I=1,N
   80       P(I,3)=-P(I,3)
          ENDIF
          CALL LUDBRB(0,0,0.,0.,0.D0,0.D0,-DBETA(1,3))
        ELSE
          CALL LUDBRB(0,0,0.,0.,-DBETA(1,1),-DBETA(1,2),-DBETA(1,3))
          CALL LUDBRB(0,0,0.,-SPHI(1),0.D0,0.D0,0.D0)
          CALL LUDBRB(0,0,-STHETA(1),0.,0.D0,0.D0,0.D0)
        ENDIF
        LST(28)=2
      ELSEIF(IFRAME.EQ.4) THEN
        THEBOS=PLU(N,13)
        PHIBOS=PLU(N,15)
        CALL LUDBRB(0,0,0.,-PHIBOS,0.D0,0.D0,0.D0)
        CALL LUDBRB(0,0,-THEBOS,0.,0.D0,0.D0,0.D0)
        LST(28)=4
      ENDIF
      GOTO 20

   90 IF(IFRAME.LE.3) THEN
        CALL LUDBRB(0,0,THEBOS,PHIBOS,0.D0,0.D0,0.D0)
        LST(28)=3
      ENDIF
      GOTO 20

  100 WRITE(*,10000) IFRAME,IPHI,LST(28),LST(29)
10000 FORMAT(' BAD VARIABLES IN SUBROUTINE LFRAME: IFRAME,IPHI,',
     +'LST(28),LST(29) =',4I5)
      RETURN
      END
*CMZ :  1.01/40 11/12/95  19.19.33  by  Piero Zucchelli
*CMZ :  1.01/39 19/10/95  11.02.39  by  Piero Zucchelli
*CMZ :  1.01/16 14/05/95  11.46.23  BY  PIERO ZUCCHELLI
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 04/03/95  17.16.12  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 19/08/94  11.04.46  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 17/07/94  18.30.49  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LINIT(LFILE,LEPIN,PLZ,PPZ,INTER)

C...INITIALIZE FOR AN INCOMING LEPTON (TYPE LEPIN, MOMENTUM PZ=PLZ)
C...AND TARGET NUCLEON (MOMENTUM PZ=PPZ) TO INTERACT VIA INTER.
C...FIND MAXIMUM OF DIFFERENTIAL CROSS SECTION, CALCULATE QCD EVENT
C...PROBABILITIES OR READ THEM FROM LOGICAL FILE LFILE (IF >0).
C...NUMERICAL INTEGRATION TO OBTAIN TOTAL CROSS-SECTION.

      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LGRID/ NXX,NWW,XX(20),WW(15),PQG(20,15,3),PQQB(20,15,2),
     +QGMAX(20,15,3),QQBMAX(20,15,2),YCUT(20,15),XTOT(20,15),NP
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
      COMMON /LMINUI/ XKIN(4),UKIN(4),WKIN(4),AIN(4),BIN(4),
     +MAXFIN,RELUP,RELERR,RELER2,FCNMAX
      COMMON /LMINUC/ NAMKIN(4),NAM(30)
      COMMON /LPFLAG/ LST3
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      CHARACTER*10 NAMKIN,NAM
      DIMENSION LSTW(40),PARLW(30)
      DOUBLE PRECISION DTHETA,DPHI,DBETA
      DATA PI/3.1415927/,NCALL/0/

      NCALL=NCALL+1
      LST3=LST(3)

* S SUPPRESSION ACCORDING TO JONES ET AL,
* Z PHYS C 27, 43-52 (1985)

      PARJ(2)=0.203


* WARNING !!! UNDOCUMENTED VARIATION IN Q0 MASS SCALE
      PYPAR(12)=0.2*2

      PARJ(32)=0.5
* (D=1. GeV) is, with quark masses added, used to define
* the minimum allowable energy of a
* colour-singlet jet system


* TAU MASS CORRECTION ACCORDING TO PDG

      IF (NCALL.EQ.1) THEN
        CALL LUGIVE('PMAS(15,1)=1.777')
        WRITE(*,*)' ACTUAL CTAU LIFETIME IN MM:',PMAS(15,4)
        CALL LUGIVE('PMAS(15,4)=0.0886')
      ENDIF



* NOW SOME PERSONAL PATCHES TO LUSHOW:
*       Q2=M^2/4
*     MSTJ(44)=1
* LAMBDA=0.25 ACCORDING TO THE STRUCTURE FUNCTIONS VALUE
      PARJ(81)=0.25
* INVARIANT MASS CUTOFF FOR PARTON SHOWERS
*        PARJ(82)=0.5
* INVARIANT MASS CUTOFF FOR PHOTON EMISSION
*        PARJ(83)=0.5




      IF(LST(32).NE.0) THEN
        IF (LEPIN.NE.16) THEN
          WRITE(*,*)'***WARNING: LST(32) OPTION TESTED ONLY WITH TAU'
        ENDIF
      ENDIF


      IF(LST(8).LT.2) THEN
C...DEFAULT FRAGMENTATION PARAMETERS SUITABLE FOR PARTON SHOWER CASE,
C...RESET WHEN USING MATRIX ELEMENTS OR NO QCD.
        PARJ(21)=0.4
        PARJ(33)=1.1
        PARJ(41)=1.
        PARJ(42)=0.7
        PARJ(43)=1.
        PARJ(44)=0.7
      ELSE
C...RESET PYTHIA PARAMETERS FROM LEPTO PARAMETERS.
        IF(MOD(LST(8),10).EQ.3.OR.MOD(LST(8),10).EQ.5) IPY(13)=0
        IF(MOD(LST(8),10).EQ.4.OR.MOD(LST(8),10).EQ.5) IPY(14)=0
        IPY(8)=LST(12)
      ENDIF

      IF(LST(18).GE.1) THEN
C...W, Z MASSES FROM THETA-WEINBERG, FERMI CONSTANT GF AND RAD. CORR.
        PMAS(24,1)=SQRT(PI*PARL(16)/(SQRT(2.)*PARL(17)*PARL(5)*
     +  (1.-PARL(18))))
        PMAS(23,1)=PMAS(24,1)/SQRT(1.-PARL(5))
      ENDIF
C...COUPLINGS BETWEEN Z0 AND LEFT/RIGHT-HANDED LEPTONS AND QUARKS.
      ZL(1,1)=-.5+PARL(5)
      ZL(1,2)=PARL(5)
      ZL(2,1)=ZL(1,2)
      ZL(2,2)=ZL(1,1)
      ZL(1,3)=0.5
      ZL(2,3)=0.
      ZL(1,4)=0.
      ZL(2,4)=0.5
      DO 10 IFL=1,8
        ZQ(1,IFL)=SIGN(0.5,QC(IFL))-QC(IFL)*PARL(5)
   10 ZQ(2,IFL)=-QC(IFL)*PARL(5)

C...SET INITIAL STATE.
      LST(23)=INTER
      KSAVE(1)=LEPIN
      KSAVE(2)=2212
      K(1,1)=21
      K(1,2)=KSAVE(1)
      K(1,3)=0
      K(1,4)=0
      K(1,5)=0
      K(2,1)=21
      K(2,2)=KSAVE(2)
      K(2,3)=0
      K(2,4)=0
      K(2,5)=0
      P(1,1)=0.
      P(1,2)=0.
      P(1,3)=PLZ
      P(1,5)=ULMASS(KSAVE(1))
      P(1,4)=SQRT(P(1,3)**2+P(1,5)**2)
      P(2,1)=0.
      P(2,2)=0.
      P(2,3)=PPZ
      P(2,5)=ULMASS(KSAVE(2))
      P(2,4)=SQRT(P(2,3)**2+P(2,5)**2)
      N=2
      LST(28)=3
C...SAVE MOMENTUM VECTORS OF INCOMING PARTICLES
      DO 20 I=1,2
        DO 20 J=1,5
   20 PSAVE(3,I,J)=P(I,J)
C...DOT-PRODUCT OF INITIAL PARTICLES, CMS ENERGY
      PARL(21)=2.*(DBLE(P(1,4))*DBLE(P(2,4))-DBLE(P(1,3))*DBLE(P(2,3)))
      ROOTS=SQRT((DBLE(P(1,4))+DBLE(P(2,4)))**2
     +          -(DBLE(P(1,3))+DBLE(P(2,3)))**2)
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10000)
     +LEPIN,(P(1,J),J=1,3),PARL(1),PARL(2),(P(2,J),J=1,3),INTER,ROOTS
      IF(PLZ*PPZ.GT.0.1) THEN
        WRITE(6,10100)
        STOP
      ENDIF

      IF(PSAVE(3,1,3).LT.0.) THEN
C...FLIP EVENT TO HAVE INITIAL LEPTON ALONG +Z AXIS
        P(1,3)=-P(1,3)
        P(2,3)=-P(2,3)
      ENDIF
C...BOOST PARAMETERS TO CMS OF INCOMING PARTICLES
      DBETA(1,1)=0.D0
      DBETA(1,2)=0.D0
      DBETA(1,3)=(DBLE(P(1,3))+DBLE(P(2,3)))/(DBLE(P(1,4))+DBLE(P(2,4)))
      SPHI(1)=0.D0
      STHETA(1)=0.D0
      IF(LST(17).NE.0) THEN
C...FOR VARYING BEAM ENERGIES, TRANSFORM TO CMS, LEPTON ALONG +Z AXIS.
        CALL LUDBRB(0,0,0.,0.,0.D0,0.D0,-DBETA(1,3))
        SPHI(1)=ULANGL(P(1,1),P(1,2))
        CALL LUDBRB(0,0,0.,-SPHI(1),0.D0,0.D0,0.D0)
        STHETA(1)=ULANGL(P(1,3),P(1,1))
        CALL LUDBRB(0,0,-STHETA(1),0.,0.D0,0.D0,0.D0)
        LST(28)=2
      ENDIF

C...EFFECTIVE LIMITS ON KINEMATIC VARIABLES X, Y, Q**2, W**2
      PM2=P(2,5)**2
      S=ROOTS**2
      CUT(1)=MAX(CUT(1),0.)
      CUT(2)=MIN(CUT(2),1.)
      CUT(3)=MAX(CUT(3),0.)
      CUT(4)=MIN(CUT(4),1.)
      CUT(5)=MAX(CUT(5),0.)
      CUT(6)=MIN(CUT(6),S)
      CUT(7)=MAX(CUT(7),0.)
      CUT(8)=MIN(CUT(8),S)
      CUT(9)=MAX(CUT(9),0.)
      CUT(10)=MIN(CUT(10),S/(2.*P(2,5)))
      XMIN =CUT(1)
      XMAX =CUT(2)
      YMIN =CUT(3)
      YMAX =CUT(4)
      Q2MIN=CUT(5)
      Q2MAX=CUT(6)
      W2MIN=CUT(7)
      W2MAX=CUT(8)
      UMIN =CUT(9)
      UMAX =CUT(10)
      DO 30 I=1,2
        IF(LST(32).NE.0) THEN
          XMIN=MAX(XMIN,Q2MIN/(S*YMAX),Q2MIN/(2.*P(2,5)*CUT(10)),
     +    1.-(W2MAX-PM2)/MAX(S*YMIN,1.E-22), 1.-(W2MAX-PM2)/MAX(2.*P(2,
     +    5)*UMIN,1.E-22))
          XMIN=MAX(XMIN,XMIN, ULMASS(LEPIN-1)**2/(2.*ULMASS(2212)*(PLZ+
     +    ULMASS(2212))) )
          XMINCMP=MAX(XMIN,Q2MIN/(S*YMAX),Q2MIN/(2.*P(2,5)*CUT(10)),
     +    1.-(W2MAX-PM2)/MAX(S*YMIN,1.E-22), 1.-(W2MAX-PM2)/MAX(2.*P(2,
     +    5)*UMIN,1.E-22))
*      WRITE(*,*)'XMIN, XMINCMP=',XMIN,XMINCMP
        ELSE
          XMIN=MAX(XMIN,Q2MIN/(S*YMAX),Q2MIN/(2.*P(2,5)*CUT(10)),
     +    1.-(W2MAX-PM2)/MAX(S*YMIN,1.E-22), 1.-(W2MAX-PM2)/MAX(2.*P(2,
     +    5)*UMIN,1.E-22))
        ENDIF
        XMAX=MIN(XMAX,Q2MAX/MAX(S*YMIN,1.E-22), Q2MAX/MAX(2.*P(2,5)*
     +  UMIN,1.E-22), 1.-(W2MIN-PM2)/(S*YMAX),1.-(W2MIN-PM2)/(2.*P(2,5)
     +  *UMAX))
        YMIN=MAX(YMIN,Q2MIN/(S*XMAX),(W2MIN-PM2)/(S*(1.-XMIN)), (W2MIN-
     +  PM2+Q2MIN)/S,2.*P(2,5)*UMIN/S)
        YMAX=MIN(YMAX,Q2MAX/MAX(S*XMIN,1.E-22), (W2MAX-PM2)/MAX(S*(1.-
     +  XMAX),1.E-22), (W2MAX-PM2+Q2MAX)/S,2.*P(2,5)*UMAX/S)
        Q2MIN=MAX(Q2MIN,S*XMIN*YMIN,S*YMIN-W2MAX+PM2, 2.*P(2,5)*UMIN*
     +  XMIN,(W2MIN-PM2)*XMIN/(1.-XMIN))
        Q2MAX=MIN(Q2MAX,S*XMAX*YMAX,S*YMAX-W2MIN+PM2, 2.*P(2,5)*UMAX*
     +  XMAX,(W2MAX-PM2)*XMAX/MAX(1.-XMAX,1.E-22))
        W2MIN=MAX(W2MIN,S*(1.-XMAX)*YMIN+PM2,Q2MIN*(1.-XMAX)/XMAX+PM2,
     +  S*YMIN-Q2MAX+PM2,2.*P(2,5)*UMIN*(1.-XMAX)+PM2)
        W2MAX=MIN(W2MAX,S*(1.-XMIN)*YMAX+PM2, Q2MAX*(1.-XMIN)/MAX(XMIN,
     +  1.E-22)+PM2, S*YMAX-Q2MIN+PM2,2.*P(2,5)*UMAX*(1.-XMIN)+PM2)
C     UMIN=MAX(UMIN,....)
C     UMAX=MIN(UMAX,....)
   30 CONTINUE




      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10200)
     +CUT,XMIN,XMAX,YMIN,YMAX,Q2MIN,Q2MAX,W2MIN,W2MAX,UMIN,UMAX
      IF(XMAX.LT.XMIN.OR.YMAX.LT.YMIN.OR.Q2MAX.LT.Q2MIN.OR.
     +W2MAX.LT.W2MIN) THEN
        IF(LST(3).GE.1) WRITE(6,10300)
        IF(LST(3).GE.2) THEN
*          WRITE(6,11600)
*          STOP
        ENDIF
      ENDIF


      PARI(11)=(PARL(1)-PARL(2))/PARL(1)
      KSAVE(4)=LEPIN
      ILEP=1
      IF(LEPIN.LT.0) ILEP=2
      INU=0
      IF(IABS(LEPIN).EQ.12.OR.IABS(LEPIN).EQ.14.OR.
     +IABS(LEPIN).EQ.16) INU=1
      IF(INU.EQ.1) THEN
C...SET FULL POLARISATION FOR INCOMING NEUTRINO.
        PARL(6)=-1.
        IF(LEPIN.LT.0) PARL(6)=1.
      ENDIF
      IF(LST(23).EQ.1.AND.INU.EQ.0) THEN
C...ELECTROMAGNETIC INTERACTION.
        KSAVE(3)=22
        IG=1
        IZ=0
      ELSEIF(LST(23).EQ.2) THEN
C...WEAK CHARGED CURRENT, ONLY ONE HELICITY STATE CONTRIBUTES.
        IF(KSAVE(1).LT.0.AND.PARL(6).LT.-0.99
     +  .OR.KSAVE(1).GT.0.AND.PARL(6).GT.0.99) THEN
          IF(LST(3).GE.1) WRITE(6,10400) LEPIN,PARL(6)
          IF(LST(3).GE.2) THEN
            WRITE(6,11600)
            STOP
          ENDIF
        ENDIF
        IF(MOD(IABS(LEPIN),2).EQ.0) THEN
          KSAVE(3)=ISIGN(24,LEPIN)
          KSAVE(4)=ISIGN(IABS(LEPIN)-1,LEPIN)
        ELSE
          KSAVE(3)=ISIGN(24,-LEPIN)
          KSAVE(4)=ISIGN(IABS(LEPIN)+1,LEPIN)
        ENDIF
      ELSEIF(LST(23).EQ.3.OR.(LST(23).EQ.4.AND.INU.EQ.1)) THEN
C...WEAK NEUTRAL CURRENT.
        KSAVE(3)=23
        IG=0
        IZ=1
      ELSEIF(LST(23).EQ.4.AND.INU.EQ.0) THEN
C...NEUTRAL CURRENT, ELECTROMAGNETIC AND WEAK WITH INTERFERENCE.
        KSAVE(3)=23
        IG=1
        IZ=1
      ELSE
        IF(LST(3).GE.1) WRITE(6,10500) INTER,LEPIN
        IF(LST(3).GE.2) THEN
          WRITE(6,11600)
          STOP
        ENDIF
      ENDIF

C...CHOICE OF INDEPENDENT VARIABLES.
      IF(LST(1).EQ.0) THEN
        LST(31)=1
        IF(INTER.EQ.2.OR.INTER.EQ.3) LST(31)=2
      ELSE
        LST(31)=IABS(LST(1))
      ENDIF
      IF(LST(31).LT.1.OR.LST(31).GT.3) THEN
        IF(LST(3).GE.1) WRITE(6,10600) LST(1),LST(31)
        IF(LST(3).GE.2) THEN
          WRITE(6,11600)
          STOP
        ENDIF
      ENDIF
      IF(LST(1).LT.0) THEN
C...USER-DEFINED OPTIMIZATION PARAMETERS.
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10700)
     +  OPTX,OPTY,OPTQ2,OPTW2
      ELSE
C...SET OPTIMIZATION PARAMETERS.
        DO 40 I=1,4
          OPTX(I)=0.
          OPTY(I)=0.
          OPTQ2(I)=0.
   40   OPTW2(I)=0.
        IF(INTER.EQ.1) THEN
          OPTX(2)=1.
          OPTY(1)=1.
          OPTQ2(3)=1.
          OPTW2(3)=1.
        ELSEIF(INTER.EQ.4) THEN
          OPTX(1)=0.1
          OPTX(2)=1.
          OPTY(1)=1.
          OPTQ2(1)=0.5
          OPTQ2(2)=0.5
          OPTQ2(3)=1.
          OPTW2(1)=0.5
          OPTW2(2)=0.5
          OPTW2(3)=1.
        ELSE
          OPTX(1)=1.
          OPTY(1)=1.
          OPTQ2(1)=1.
          OPTW2(1)=1.
        ENDIF
      ENDIF

C...INITIALIZE MONTE CARLO ESTIMATE OF CROSS SECTION.
      PARL(24)=0.
      PARI(27)=0.
      PARI(28)=0.
      PARI(29)=0.
      PARI(30)=0.
      PARI(32)=0.
      IF(LST(23).EQ.2) THEN
C...CONSTANT FACTOR GF**2/PI FOR CC, TRANSFORMATION TO PICOBARN.
        PARI(31)=PARL(17)**2/PI*0.39E+09
      ELSE
C...CONSTANT FACTOR 2*PI*ALPHA**2 FOR NC, TRANSFORMATION TO PICOBARN.
        PARI(31)=2.*PI*PARL(16)**2*0.39E+09
      ENDIF
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10800)
     +(I,LST(I),LST(I+10),PARL(I),PARL(I+10),I=1,10)

C...SET UP GRID WITH LONGITUDINAL STRUCTURE FUNCTION, QCD OR TARGET MASS
      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      IF(INTER.EQ.1.AND.LST(11).NE.0) CALL FLTABL

C...GET INTEGRATED CROSS-SECTION.
      PARL(23)=0.
      IF(LST(10).GT.0) CALL LXSECT
      IF(LQCD.EQ.2.OR.LTM.EQ.2) THEN
        WRITE(6,10900)
        IF(LQCD.EQ.2) WRITE(6,11000)
        IF(LTM .EQ.2) WRITE(6,11100)
        WRITE(6,11200)
      ENDIF

      IF(LST(2).EQ.1) THEN
C...FIND MAX VALUE OF DIFFERENTIAL CROSS SECTION FOR REJECTION.
        UKIN(1)=(XMAX+XMIN)/2.
        WKIN(1)=0.8*(XMAX-XMIN)/2.
        AIN(1)=XMIN
        BIN(1)=XMAX
        IF(LST(31).EQ.1) THEN
          UKIN(2)=(Q2MAX+Q2MIN)/2.
          WKIN(2)=0.8*(Q2MAX-Q2MIN)/2.
          AIN(2)=Q2MIN
          BIN(2)=Q2MAX
          NAMKIN(2)='      Q**2'
        ELSEIF(LST(31).EQ.2) THEN
          UKIN(2)=(YMAX+YMIN)/2.
          WKIN(2)=0.8*(YMAX-YMIN)/2.
          AIN(2)=YMIN
          BIN(2)=YMAX
          NAMKIN(2)='         Y'
        ELSEIF(LST(31).EQ.3) THEN
          UKIN(2)=(W2MAX+W2MIN)/2.
          WKIN(2)=0.8*(W2MAX-W2MIN)/2.
          AIN(2)=W2MIN
          BIN(2)=W2MAX
          NAMKIN(2)='      W**2'
        ENDIF
C...MAXIMUM OBTAINED BY MINIMIZING -(DIFF. X-SECTION).
        CALL LTIMEX(TI1)
        CALL LMINEW
        CALL LTIMEX(TI2)
*        WRITE(*,*)'MIN. TIME:',TI2-TI1,' S'
        PARI(LST(23))=FCNMAX*1.5
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,11300)
     +  PARI(LST(23)),TI2-TI1
      ENDIF

      IF(LFILE.GT.0) THEN
C...READ QCD WEIGHTS FROM FILE.
        READ(LFILE) LSTW,PARLW,NXX,NWW,NP,XX,WW
        IPMAX=2
        IF(LSTW(17).NE.0) IPMAX=3
        READ(LFILE) (((PQG(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,NP),
     +  (((PQQB(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,NP),
     +  (((QGMAX(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,IPMAX),
     +  (((QQBMAX(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,MIN(2,IPMAX)),
     +  YCUT
        IF(NP.NE.1) READ(LFILE) XTOT
        CLOSE(LFILE)
C...RESET PARAMETERS FOR MATRIX ELEMENT INTEGRATION.
        PARL(8)=PARLW(8)
        PARL(9)=PARLW(9)
        PARL(11)=PARLW(11)
        PARL(12)=PARLW(12)
        PARL(13)=PARLW(13)
C...CHECK CURRENT PARAMETER VALUES AGAINST THOSE USED WHEN
C...CALCULATING WEIGHTS.
        IF(LST(12).NE.LSTW(12).OR.LST(13).NE.LSTW(13)
     +  .OR.LST(15).NE.LSTW(15).OR.LST(16).NE.LSTW(16)
     +  .OR.LST(17).NE.LSTW(17).OR.LST(23).NE.LSTW(23)
     +  .OR.ABS(PARL(1)-PARLW(1)).GT.0.1.OR.ABS(PARL(2)-PARLW(2)).GT.0.1
     +  .OR.ABS(PARL(5)-PARLW(5)).GT.0.01
     +  .OR.ABS(PARL(6)-PARLW(6)).GT.0.1) THEN
          IF(LST(3).GE.1) WRITE(6,11400) LST(12),LSTW(12),LST(13),
     +    LSTW(13),LST(15), LSTW(15),LST(16),LSTW(16),LST(17),LSTW(17),
     +    LST(23),LSTW(23), PARL(1),PARLW(1),PARL(2),PARLW(2),PARL(5),
     +    PARLW(5),PARL(6), PARLW(6)
          IF(LST(3).GE.2) THEN
            WRITE(6,11600)
            STOP
          ENDIF
        ENDIF
      ELSEIF(LST(8).EQ.1.OR.LST(8)/10.EQ.1.OR.MOD(LST(8),10).EQ.9) THEN
C...CALCULATE WEIGHTS IF 1ST ORDER QCD REQUESTED.
        CALL LTIMEX(TI1)
        CALL LWEITS(LFILE)
        CALL LTIMEX(TI2)
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,11500)
     +  TI2-TI1
      ENDIF

C...RESET COUNTERS TO ZERO FOR MONTE CARLO ESTIMATE OF CROSS SECTION.
      PARI(27)=0.
      PARI(28)=0.
      PARI(29)=0.
      PARI(30)=0.
      LST(20)=0
      RETURN

10000 FORMAT('1',//,5X,'THE LUND MONTE CARLO FOR DEEP INELASTIC LEPTON-'
     +,'NUCLEON SCATTERING',/,5X,65('='),//,
     +25X,'LEPTO VERSION 6.1, MAY 4, 1992',//,
     +' LEPTON: TYPE =',I3,5X,'MOMENTUM (PX,PY,PZ) =',3F8.1,
     +' GEV',//,' TARGET: A, Z =',2F3.0,2X,
     +'MOMENTUM (PX,PY,PZ) =',3F8.1,' GEV',//,
     +' INTERACTION :',I3,14X,' CMS ENERGY =',1PG12.4,' GEV',/)
10100 FORMAT(' WARNING: LEPTON AND NUCLEON MOMENTA IN SAME DIRECTION',
     +' NOT ALLOWED.',/,10X,'EXECUTION STOPPED.')
10200 FORMAT(/,' USER APPLIED CUTS (+ PHASE SPACE) : ',1P,
     +      G12.4,' <   X   < ',G12.4,
     +/,37X,G12.4,' <   Y   < ',G12.4,
     +/,37X,G12.4,' < Q**2  < ',G12.4,
     +/,37X,G12.4,' < W**2  < ',G12.4,
     +/,37X,G12.4,' <  NU   < ',G12.4,
     +/,37X,G12.4,' <  E''   < ',G12.4,
     +/,37X,G12.4,' < THETA < ',G12.4,/,
     +/,       ' EFFECTIVE RANGES (FROM ABOVE CUTS): ',
     +      G12.4,' <   X   < ',G12.4,
     +/,37X,G12.4,' <   Y   < ',G12.4,
     +/,37X,G12.4,' < Q**2  < ',G12.4,
     +/,37X,G12.4,' < W**2  < ',G12.4,
     +/,37X,G12.4,' <  NU   < ',G12.4)
10300 FORMAT(' WARNING: EFFECTIVE UPPER LIMIT OF KINEMATICAL ',
     +'VARIABLE(S) SMALLER THAN CORRESPONDING LOWER LIMIT.')
10400 FORMAT(' WARNING: WEAK CHARGED CURRENT CROSS SECTION ZERO FOR ',
     +'SPECIFIED LEPTON HELICITY; LEPIN, PARL(6) =',I3,F5.2)
10500 FORMAT(' WARNING: UNRECOGNIZED INTERACTION IN LINIT CALL: ',
     +'INTER = ',I5,'  FOR LEPTON LEPIN =',I5)
10600 FORMAT(' WARNING: UNALLOWED VALUE OF LST(1) =',I3,
     +' AND/OR LST(31) =',I3)
10700 FORMAT(/,' USER-DEFINED OPTIMIZATION PARAMETERS:',
     +/,5X,'OPTX(1...4)  =',4G11.3,/,5X,'OPTY(1...4)  =',4G11.3,
     +/,5X,'OPYQ2(1...4) =',4G11.3,/,5X,'OPTW2(1...4) =',4G11.3,/)
10800 FORMAT(/,' PARAMETER VALUES:',//,9X,'I',4X,'LST(I)',1X,
     +'LST(I+10)',8X,'PARL(I)',5X,'PARL(I+10)',1P,
     +/,5X,55('-'),10(/,3I10,2G15.4),/)
10900 FORMAT(' WARNING: CROSS SECTION, PARL(23), EXCLUDES FL (SEE ',
     +'LST(11)) FROM:')
11000 FORMAT(10X,'QCD, SINCE EVALUATED EVENT BY EVENT FOR LQCD=2')
11100 FORMAT(10X,'TM , SINCE EVALUATED EVENT BY EVENT FOR LTM =2')
11200 FORMAT(' CROSS SECTION IN PARL(24) INCLUDES THESE CONTRIBUTIONS.')
11300 FORMAT(' MAX OF DIFFERENTIAL CROSS SECTION (FOR WEIGHTING) =',
     +E12.4,/,' OBTAINED IN ',F7.2,' SECONDS.',/)
11400 FORMAT(//,' WARNING: CURRENT PARAMETER VALUES DO NOT MATCH ',
     +'WITH THOSE USED WHEN CALCULATING QCD WEIGHTS.',//,15X,
     +'CURRENT VALUE     VALUE FOR WEIGHTS',/,
     +/,'     LST(12)   ',I12,10X,I12,
     +/,'     LST(13)   ',I12,10X,I12,
     +/,'     LST(15)   ',I12,10X,I12,
     +/,'     LST(16)   ',I12,10X,I12,
     +/,'     LST(17)   ',I12,10X,I12,
     +/,'     LST(23)   ',I12,10X,I12,
     +/,'     PARL(1)   ',E12.4,10X,E12.4,
     +/,'     PARL(2)   ',E12.4,10X,E12.4,
     +/,'     PARL(5)   ',E12.4,10X,E12.4,
     +/,'     PARL(6)   ',E12.4,10X,E12.4)
11500 FORMAT(/,' TIME FOR CALCULATING QCD WEIGHTS =',F5.1,' SECONDS',/)
11600 FORMAT(' EXECUTION STOPPED ',/)
      END
*CMZ :          03/03/97  19.17.26  by  Unknown
*CMZ :  1.01/22 27/05/95  16.02.43  BY  PIERO ZUCCHELLI
*CMZ :  1.01/21 27/05/95  15.58.38  BY  PIERO ZUCCHELLI
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 22/07/94  17.46.01  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      FUNCTION LKINEM(L)

C...CALCULATE KINEMATICAL VARIABLES AND REJECT (OPTIONALLY) IF OUTSIDE
C...REQUIRED LIMITS.

      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
      DOUBLE PRECISION DTHETA,DPHI,DBETA,DE,DPZ,DPT,DETOT
C     print*,' in lkinem'
      LKINEM=1
      IF(L.EQ.-3) THEN
C...X,W KNOWN FROM LWEITS, NO CUTS APPLIED.
        U=(W2-P(2,5)**2)/(2.*P(2,5)*(1.-X))
        Q2=2.*P(2,5)*U*X
        Y=Q2/(PARL(21)*X)
        GOTO 20
      ENDIF
C...X,Y GIVEN.
      PARL(22)=Y*PARL(21)
      Q2=X*PARL(22)
      U=PARL(22)/(2.*P(2,5))
      W2=PARL(22)*(1.-X)+P(2,5)**2
      P(4,5)=ULMASS(K(4,2))
      IF(P(4,5)/SQRT(PARL(21)).LT.0.001) THEN
C...SIMPLER FORMULAE FOR EFFECTIVELY MASSLESS SCATTERED LEPTON.
        DE=DBLE(P(1,4))*(1.-DBLE(Y))+DBLE(X)*DBLE(Y)*DBLE(ABS(P(2,3)))
        DPZ=DE-DBLE(X)*DBLE(Y)*(DBLE(P(2,4))+DBLE(ABS(P(2,3))))
      ELSE
C...FORMULAE FOR MASSIVE SCATTERED LEPTON.
        DE=DBLE(P(1,4))+(DBLE(ABS(P(2,3)))*(DBLE(Q2)+DBLE(P(4,5))**2)/
     +  (2.D0*DBLE(P(1,4)))-DBLE(PARL(22))/2.D0)/
     +  (DBLE(P(2,4))+DBLE(ABS(P(2,3))))
        DPZ=DBLE(P(1,4))-(DBLE(P(2,4))*(DBLE(Q2)+DBLE(P(4,5))**2)/
     +  (2.D0*DBLE(P(1,4)))+DBLE(PARL(22))/2.D0)/
     +  (DBLE(P(2,4))+DBLE(ABS(P(2,3))))
      ENDIF
      DPT=DE**2-DPZ**2-DBLE(P(4,5))**2
      IF(DPT.LT.0.D0) RETURN
      DPT=SQRT(DPT)
      P(4,1)=DPT
      P(4,2)=0.
      P(4,3)=DPZ
      P(4,4)=DE
      P(3,1)=-DPT
      P(3,2)=0.
      P(3,3)=DBLE(P(1,3))-DPZ
      P(3,4)=DBLE(P(1,4))-DE
      P(3,5)=-SQRT(Q2)
      K(3,3)=1
      K(4,3)=1
      N=4
      IF(L.EQ.3) GOTO 20
C      print*,' cut cut cut'
      IF(X.LT.CUT(1).OR.X.GT.CUT(2))  RETURN
      IF(Y.LT.CUT(3).OR.Y.GT.CUT(4))  RETURN
      IF(Q2.LT.CUT(5).OR.Q2.GT.CUT(6))  RETURN
      IF(W2.LT.CUT(7).OR.W2.GT.CUT(8))  RETURN
      IF(U.LT.CUT(9).OR.U.GT.CUT(10))  RETURN
      IF(LST(17).EQ.0) THEN
C            print*,'survived'
        IF(P(4,4).LT.CUT(11).OR.P(4,4).GT.CUT(12)) THEN
          WRITE(*,*)'CUTTING TOO LOW LEPTON ENERGY',P(4,4),CUT(11)
          RETURN
        ENDIF
        THETAL=PLU(4,13)
C       THETAL=ACOS((P(1,1)*P(4,1)+P(1,2)*P(4,2)+P(1,3)*P(4,3))
C    &  /SQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2)/
C    &  SQRT(P(4,1)**2+P(4,2)**2+P(4,3)**2))
      ELSE
C...NO CUTS ON ENERGY, ANGLE FOR INITIALISATION OF VARYING ENERGY MODE
        IF(LST(20).NE.0) GOTO 20
C...TRANSFORM SCATTERED LEPTON BACK TO LAB SYSTEM TO MAKE CUT
C...IN ENERGY AND ANGLE (DEFINED AS SPACE ANGLE TO INCOMING LEPTON).
        DO 10  J=1,5
          K(6,J)=K(4,J)
   10   P(6,J)=P(4,J)
*        WRITE(*,*)'LKINEM BOOST'
        CALL LUDBRB(6,6,STHETA(1),SPHI(1),0.D0,0.D0,0.D0)
        CALL LUDBRB(6,6,0.,0.,DBETA(1,1),DBETA(1,2),DBETA(1,3))
        IF(P(6,4).LT.CUT(11).OR.P(6,4).GT.CUT(12))  RETURN
        THETAL=ACOS((PSAVE(3,1,1)*P(6,1)+PSAVE(3,1,2)*P(6,2)+
     +  PSAVE(3,1,3)*P(6,3))
     +  /SQRT(PSAVE(3,1,1)**2+PSAVE(3,1,2)**2+PSAVE(3,1,3)**2)/
     +  SQRT(P(6,1)**2+P(6,2)**2+P(6,3)**2))
      ENDIF
      IF(THETAL.LT.CUT(13).OR.THETAL.GT.CUT(14))  then
C            print*,' problem with thetal'
         RETURN
      end if
 20   LKINEM=0
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C#######################################################################
C
C   THE FOLLOWING ROUTINES ARE SLIGHTLY MODIFIED MINIMIZATION ROUTINES
C   FROM THE MINUIT PROGRAM PACKAGE.
C
C **********************************************************************

      SUBROUTINE LMCMND

C...THIS IS THE MINUIT ROUTINE COMAND.
CC        GETS IFORMATION FROM /LMINUI/  AND TAKES APPROPRIATE ACTION,
CC        EITHER DIRECTLY BY SKIPPING TO THE CORRESPONDING CODE IN
CC        LMCMND, OR BY SETTING UP A CALL TO A SUBROUTINE
CC
      COMMON /LMINUI/ XKIN(4),UKIN(4),WKIN(4),AIN(4),BIN(4),
     &MAXFIN,RELUP,RELERR,RELER2,FCNMAX
      COMMON /LPFLAG/ LST3

      COMMON
     1/LMMINE/ ERP(30)  ,ERN(30)
     2/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     3/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     4/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     5/LMVARI/ V(15,15)
     7/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     7/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     C/LMCASC/ Y(16)    ,JH       ,JL
     F/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     G/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     J/LMVART/ VT(15,15)
      COMMON
     6/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     8/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     9/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     A/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     B/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC
      FVAL3 = 2.0*AMIN+1.0
C                                        . . . . . . . . . . ERROR DEF
      WORD7(1)=RELUP*ABS(AMIN)
      UP = WORD7(1)
      IF (UP .LE. 0.)  UP = 1.0
      IF (ISW(2) .GE. 1)  CALL LMPRIN(1,AMIN)
      WORD7(1)=MAXFIN
      WORD7(2)=RELERR*UP
      NFCNMX = WORD7(1) + 0.5
      IF (NFCNMX .LE. 0)  NFCNMX = 1000
      EPSI = WORD7(2)
      IF (EPSI .LE. 0.)  EPSI = 0.1 * UP
      NEWMIN = 0
      ITAUR = 0
      ISW(1) = 0
      CALL LMSIMP
      IF(ABS(DIRIN(1)).LE.ABS(EPSMAC*X(1)).AND.
     +   ABS(DIRIN(2)).LE.ABS(EPSMAC*X(2))) THEN
        IF(LST3.GE.1) WRITE(6,10000)
        GOTO 10
      ENDIF
      WORD7(1)=MAXFIN
      RELERR=RELER2*RELERR
      WORD7(2)=RELERR*UP
      NFCNMX = WORD7(1) + 0.5
      IF (NFCNMX .LE. 0)  NFCNMX = 1000
      EPSI = WORD7(2)
      IF (EPSI .LE. 0.)  EPSI = 0.1 * UP
      CALL LMSIMP
   10 FCNMAX=ABS(AMIN)
      IF(ISW(1).GE.1) THEN
        IF(LST3.GE.1) WRITE(6,10100)
        FCNMAX=FCNMAX*1.25
      ENDIF
      FMAX=ABS(AMIN)
C                                        . . . . . . . . . . END, EXIT
      WORD7(1)=0.
   20 IT = WORD7(1) + 0.5
      IF (FVAL3 .EQ. AMIN .OR. IT .GT. 0) RETURN
      IFLAG = 3
      CALL LSIGMX(NPAR,GIN,F,U,IFLAG)
      NFCN = NFCN + 1
      IF(LST3.GE.1.AND.ABS(F).GT.FMAX) WRITE(6,10200) F
      RETURN

10000 FORMAT(' WARNING: STEPSIZES ARE LESS THAN MACHINE ACCURACY ',
     &'TIMES VARIABLE VALUES. NO FURTHER MINIMIZATION ATTEMPTED.')
10100 FORMAT(' WARNING: SIMPLEX MINIMIZATION HAS NOT CONVERGED ',
     &'PROPERLY.',/,10X,'RETURNED MAXIMUM INCREASED BY A FACTOR 1.25.')
10200 FORMAT(' WARNING FROM LMCMND: FUNCTION AT MINIMUM, ',E12.4,
     &'  IS SMALLER THAN STORED MINIMUM.')

      END
*CMZ :  1.01/50 22/05/96  12.22.19  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LMEPS

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON /PYPROC/ ISUB,KFL(3,2),XPY(2),SH,TH,UH,Q2PY,XSEC(0:40)
      COMMON /PYINT1/ XQPY(2,-6:6)
      DOUBLE PRECISION DTHETA,DPHI,DBETA
      DOUBLE PRECISION DPQ2,DPB(3),DPA(3),DCTHET,DROBO(5)
      DIMENSION KS(9,5),PS(9,5),ROBO(5),XPQ(-6:6)
      DATA MG/2/
      SAVE KS,PS

C     CALL GULIST(100,2)
C...SAVE EVENT RECORD IN HADRONIC CMS
      DO 10 I=1,7
        DO 10 J=1,5
          KS(I,J)=K(I,J)
   10 PS(I,J)=P(I,J)
C...REARRANGE EVENT RECORD TO PYSSPB STANDARD
      IP2=6
      IF(LST(24).EQ.3) IP2=7
      DO 20 J=1,5
        K(3,J)=0.
        P(3,J)=0.
        K(4,J)=0
        P(4,J)=0.
        K(5,J)=KS(3,J)
        P(5,J)=PS(3,J)
        K(7,J)=KS(4,J)
        P(7,J)=PS(4,J)
        K(8,J)=KS(5,J)
        P(8,J)=PS(5,J)
        K(9,J)=KS(4,J)
        P(9,J)=PS(4,J)
        K(10,J)=KS(IP2,J)
   20 P(10,J)=PS(IP2,J)
      K(5,3)=3
      K(6,3)=4
      K(7,3)=5
      K(8,3)=6
      K(9,3)=5
      K(10,3)=6
      DO 30 I=5,10
   30 K(I,1)=21
      K(9,1)=0
C...INCOMING PARTON = OUTGOING 2 PARTON - BOSON FOURVECTORS
      DO 40 J=1,4
   40 P(6,J)=P(8,J)+P(10,J)-P(5,J)
      P(6,5)=0.
      K(6,2)=LST(25)
      IF(LST(24).EQ.3) K(6,2)=21
      N=10
C     CALL GULIST(101,2)

      XR=X
      DPQ2=DBLE(Q2)
C...PARTONS WITH COLOUR INFORMATION IN HADRONIC CMS FRAME.
      DO 50  I=11,27
        DO 50  J=1,5
          K(I,J)=0
          P(I,J)=0.
   50 V(I,J)=0.
      NS=20
      DO 60  J=1,5
        K(NS+1,J)=K(5,J)
        P(NS+1,J)=P(5,J)
        K(NS+3,J)=K(6,J)
        P(NS+3,J)=P(6,J)
        K(NS+5,J)=K(8,J)
        P(NS+5,J)=P(8,J)
        K(NS+6,J)=K(10,J)
   60 P(NS+6,J)=P(10,J)
C...OLD STANDARD CONTINUATION LINES
      K(NS+2,1)=-1
      K(NS+2,3)=NS+1
      K(NS+4,1)=-1
      K(NS+4,3)=NS+3
      P(NS+4,3)=27
      P(NS+4,4)=27
C...ORIGIN AND COLOUR INFO FOR INCOMING PARTON
      K(NS+3,1)=13
      K(NS+3,3)=2
      K(NS+3,4)=27
      K(NS+3,5)=27
C...COLOUR INFO FOR TWO OUTGOING PARTONS
      K(NS+5,1)=3
      K(NS+6,1)=3
      IF(K(NS+6,2).EQ.21) THEN
C...QG-EVENT
        IF(K(NS+5,2).GT.0) THEN
          K(NS+5,4)=(NS+6)*MSTU(5)
          K(NS+5,5)=(NS+7)*MSTU(5)
          K(NS+6,4)=(NS+7)*MSTU(5)
          K(NS+6,5)=(NS+5)*MSTU(5)
        ELSE
          K(NS+5,4)=(NS+7)*MSTU(5)
          K(NS+5,5)=(NS+6)*MSTU(5)
          K(NS+6,4)=(NS+5)*MSTU(5)
          K(NS+6,5)=(NS+7)*MSTU(5)
        ENDIF
      ELSE
C...QQBAR-EVENT
        K(NS+5,4)=(NS+7)*MSTU(5)
        K(NS+5,5)=(NS+7)*MSTU(5)
        K(NS+6,4)=(NS+7)*MSTU(5)
        K(NS+6,5)=(NS+7)*MSTU(5)
      ENDIF
C...EFFECTIVE OUTGOING PARTON = SUM OF BOTH OUTGOING PARTONS
      K(NS+7,1)=14
      K(NS+7,3)=3
      IF(LST(24).EQ.2) THEN
        K(NS+7,2)=K(NS+5,2)
        IF(K(NS+7,2).EQ.21) WRITE(6,*) ' WARNING: K(NS+7,2)=',K(NS+7,2)
        IF(K(NS+7,2).GT.0) THEN
          K(NS+7,4)=(NS+3)*MSTU(5)+26
          K(NS+7,5)=(NS+3)*MSTU(5)+25
        ELSE
          K(NS+7,4)=(NS+3)*MSTU(5)+25
          K(NS+7,5)=(NS+3)*MSTU(5)+26
        ENDIF
      ELSE
        K(NS+7,2)=21
        IF(K(NS+5,2).GT.0) THEN
          K(NS+7,4)=(NS+3)*MSTU(5)+25
          K(NS+7,5)=(NS+3)*MSTU(5)+26
        ELSE
          K(NS+7,4)=(NS+3)*MSTU(5)+26
          K(NS+7,5)=(NS+3)*MSTU(5)+25
        ENDIF
      ENDIF
      DO 70  J=1,4
   70 P(NS+7,J)=P(8,J)+P(10,J)
      P(NS+7,5)=SQRT(P(NS+7,4)**2-P(NS+7,1)**2-P(NS+7,2)**2-
     +P(NS+7,3)**2)
      N=NS+7
C     CALL GULIST(103,2)

C...SCALE FOR BREMSSTRAHLUNG ETC.
      Q2PY=Q2
      IPY(40)=10
      IPY(47)=N
C...SAVE QUANTITIES FOR LATER USE.
      XPY(1)=1.
      XPY(2)=XR
      CALL PYSTFU(K(2,2),XR,Q2,XPQ)
      DO 80  IFL=-6,6
   80 XQPY(2,IFL)=XPQ(IFL)
      IF(LST(23).EQ.1) THEN
        ISUB=39
        IPY(11)=1
      ELSEIF(LST(23).EQ.3) THEN
        ISUB=39
        IPY(11)=2
      ELSEIF(LST(23).EQ.4) THEN
        ISUB=39
        IPY(11)=3
      ELSEIF(LST(23).EQ.2) THEN
        ISUB=40
      ENDIF
      IF(ISUB.EQ.39.AND.IPY(11).EQ.1) THEN
        KFL(2,1)=22
      ELSEIF(ISUB.EQ.39.AND.IPY(11).EQ.2) THEN
        KFL(2,1)=23
      ELSEIF(ISUB.EQ.39.AND.IPY(11).EQ.3) THEN
        KFL(2,1)=23
      ELSEIF(ISUB.EQ.40) THEN
        KFL(2,1)=-24
      ENDIF
      KFL(2,2)=K(6,2)
      KFL(1,1)=KFL(2,1)
      KFL(1,2)=KFL(2,2)
      IF(ISUB.EQ.39) KFL(3,1)=K(1,2)
      IF(ISUB.EQ.40) KFL(3,1)=K(1,2)+ISIGN(1,K(1,2))
      KFL(3,2)=K(27,2)
      PYVAR(2)=PARL(21)
      PYVAR(1)=SQRT(PYVAR(2))
      PYVAR(3)=P(1,5)
      PYVAR(4)=P(2,5)
      PYVAR(5)=PYVAR(1)/2.
      IPY(41)=K(1,2)
      IPY(42)=K(2,2)
      IPY(48)=0

C...GENERATE TIMELIKE PARTON SHOWER (IF REQUIRED)
      IF(IPY(13).EQ.1) THEN
        CALL LSCALE(1,QMAX)
        CALL LUSHOW(25,26,QMAX)
      ENDIF
      IT=25
      IF(N.GE.27) IT=27
      NS=N
C     CALL GULIST(104,2)

C...GENERATE SPACELIKE PARTON SHOWER (IF REQUIRED)
      IPU1=0
      IPU2=23
      IF(XPY(2)*(1.+(P(IT,5)**2+PYPAR(22))/P(21,5)**2).GT.0.999) THEN
        LST(21)=47
        RETURN
      ENDIF
      IF(IPY(14).GE.1) THEN
        CALL PYSSPB(IPU1,IPU2)
      ELSE
        DO 90  I=NS+1,NS+4
          DO 90  J=1,5
            K(I,J)=0
            P(I,J)=0.
   90   V(I,J)=0.
        K(NS+1,1)=11
        K(NS+1,2)=KFL(2,1)
        K(NS+1,3)=21
        DO 100 J=1,5
  100   P(NS+1,J)=P(21,J)
        K(NS+2,1)=-1
        K(NS+2,3)=NS+1
        K(NS+3,1)=13
        K(NS+3,2)=KFL(2,2)
        K(NS+3,3)=23
        K(NS+3,4)=23
        K(NS+3,5)=23
        P(NS+3,3)=(P(IT,5)**2+Q2)*(P(21,4)-P(21,3))/(2.*Q2)
        P(NS+3,4)=-P(NS+3,3)
        K(NS+4,1)=-1
        K(NS+4,3)=NS+3
        P(NS+4,3)=23
        P(NS+4,4)=23
        P(24,1)=NS+3
        P(24,2)=NS+3
        K(23,4)=K(23,4)+(NS+3)*MSTU(5)
        K(23,5)=K(23,5)+(NS+3)*MSTU(5)
        IPU1=0
        IPU2=NS+3
        N=N+4
      ENDIF
C     CALL GULIST(105,2)

C...ROTATE AND BOOST OUTGOING PARTON SHOWER
      IF(N.GT.31) THEN
        K(N+1,1)=0
        DO 110 J=1,4
  110   P(N+1,J)=P(NS+1,J)+P(NS+3,J)
        IF(P(N+1,4).LE.1.01*P(IT,5)) THEN
          LST(21)=50
          RETURN
        ENDIF
        ROBO(1)=ULANGL(P(IT,3),SQRT(P(IT,1)**2+P(IT,2)**2))
        ROBO(2)=ULANGL(P(IT,1),P(IT,2))
        CALL LUDBRB(25,NS,0.,-ROBO(2),0.D0,0.D0,0.D0)
        CALL LUDBRB(25,NS,-ROBO(1),0.,0.D0,0.D0,0.D0)
        DROBO(5)=-(P(IT,3)*P(IT,4)-P(N+1,4)*SQRT(P(N+1,4)**2-
     +  P(IT,4)**2+P(IT,3)**2))/(P(IT,3)**2+P(N+1,4)**2)
        CALL LUDBRB(25,NS,0.,0.,0.D0,0.D0,DROBO(5))
        ROBO(1)=ULANGL(P(N+1,3),SQRT(P(N+1,1)**2+P(N+1,2)**2))
        ROBO(2)=ULANGL(P(N+1,1),P(N+1,2))
        CALL LUDBRB(25,NS,ROBO(1),ROBO(2),0.D0,0.D0,0.D0)
      ENDIF
C     CALL GULIST(106,2)

      Q2PY=Q2
C...HADRON REMNANT AND PRIMORDIAL KT
      IPY(47)=N
      CALL PYREMM(IPU1,IPU2)
      IF(IPY(48).EQ.1) THEN
        LST(21)=48
        RETURN
      ENDIF
C     CALL GULIST(107,2)

C...REARRANGE PARTONS ALONG STRINGS
      MSTU(24)=0
      CALL LUPREP(0)
      IF(MSTU(24).NE.0) THEN
C       CALL GULIST(188,2)
        IF(LST(3).GE.1) WRITE(6,*) ' LUPREP ERROR MSTU(24)= ',MSTU(24)
        LST(21)=49
        RETURN
      ENDIF
C     CALL GULIST(109,2)

C...CLEAN UP EVENT RECORD -> ORDER:
C...1=INC. LEPTON; 2=INC. NUCLEON; 3=EXCH BOSON; 4=SCAT. LEPTON;
C...5=INC. PARTON BEFORE INITIAL SHOWER; 6=INC. PARTON AT HARD SCATTERING
C...AFTER SHOWER; 7,8=FIRST,SECOND PARTON FROM HARD SCATTERING
C...BEFORE FINAL SHOWER
      LST(26)=7
      DO 120 J=1,5
        K(N+1,J)=K(4,J)
  120 P(N+1,J)=P(4,J)
      DO 130 J=1,5
        K(3,J)=K(5,J)
        P(3,J)=P(5,J)
        K(4,J)=K(9,J)
        P(4,J)=P(9,J)
        K(5,J)=K(N+1,J)
        P(5,J)=P(N+1,J)
        K(6,J)=K(NS+3,J)
        P(6,J)=P(NS+3,J)
C     K(7,J)=K(IT,J)
C     P(7,J)=P(IT,J)
        K(7,J)=K(25,J)
        P(7,J)=P(25,J)
        K(8,J)=K(26,J)
        P(8,J)=P(26,J)
  130 CONTINUE
      K(3,3)=1
      K(4,3)=1
      K(6,1)=21
      K(6,3)=5
      K(6,4)=0
      K(6,5)=0
      K(7,1)=21
      K(7,3)=6
      K(7,4)=0
      K(7,5)=0
      K(8,1)=21
      K(8,3)=6
      K(8,4)=0
      K(8,5)=0
C...ACTIVATE LINE WITH SCATTERED LEPTON.
      K(4,1)=1
C...DEACTIVATE OBSOLETE LINES 9, 10, 21, NS+1 (EXTRA LINES WITH BOSON)
      K(9,1)=0
      K(10,1)=0
      K(21,1)=0
      IF(K(NS+1,2).EQ.K(3,2)) K(NS+1,1)=0
C...ZERO IRRELEVANT LINES WITH K(I,1)<0
      DO 150 I=1,N
        IF(K(I,1).LT.0) THEN
          DO 140 J=1,5
            K(I,J)=0
  140     P(I,J)=0.
        ENDIF
  150 CONTINUE
C     CALL GULIST(110,2)
C...DELETE INTERNAL PARTON LINES, I.E. WITH K(I,1)=13,14
      IF(MOD(LST(4)/10,10).EQ.0) THEN
        CALL LTIMEX(T1)
        CALL LUEDIT(14)
        CALL LTIMEX(T2)
C       CALL GULIST(111,2)
      ENDIF
C...DELETE EMPTY LINES
      CALL LTIMEX(T1)
      CALL LUEDIT(12)
      CALL LTIMEX(T2)
C     CALL GULIST(112,2)

      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 20/09/94  14.38.50  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LMIDAT

C...THIS IS THE MINUIT ROUTINE MIDATA.
CC        GETS PARAMETERS FROM /LMINUI/  AND /LMINUC/
CC        AND SETS UP THE STARTING PARAMETER LISTS.
CC        CONTROL THEN PASSES TO LMCMND FOR READING THE COMMAND "CARDS".
CC

      COMMON
     +/LMMINE/ ERP(30)  ,ERN(30)
     +/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     +/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     +/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     +/LMVARI/ V(15,15)
     +/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     +/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     +/LMCASC/ Y(16)    ,JH       ,JL
     +/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     +/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     +/LMVART/ VT(15,15)
      COMMON
     +/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     +/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     +/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     +/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     +/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC
      COMMON /LMINUI/ XKIN(4),UKIN(4),WKIN(4),AIN(4),BIN(4),
     +MAXFIN,RELUP,RELERR,RELER2,FCNMAX
      COMMON /LMINUC/ NAMKIN(4),NAM(30)
      COMMON /LPFLAG/ LST3
      CHARACTER*10 NAMKIN,NAM,NAMK,BLANK
      CHARACTER XTITLE*60
      REAL LMPINT
      DATA BLANK/'          '/
      DATA XTITLE/' FIND MINIMUM OF -(DIFFERENTIAL CROSS SECTION)'/
      DATA MNINIT/0/,IFATAL,NINT/0,0/
C                                        . INITIALIZE NEW DATA BLOCK . .

      MNINIT=0
      IFATAL=0
      NINT=0

      IF (MNINIT .EQ. 0)  NBLOCK=0
      MNINIT = 1
      NBLOCK = NBLOCK + 1
      VERSN = 11.79
      IF(LST3.GE.5) THEN
        WRITE (ISYSWR,10200) MAXINT,MAXEXT,VERSN,NBLOCK
        WRITE (ISYSWR,10300)
      ENDIF
      DO 10 I= 1, 7
   10 ISW(I) = 0
      SIGMA = 0.
      CALL LTIMEX(TIME)
      IF(LST3.GE.5) THEN
        WRITE (ISYSWR,11200) XTITLE,TIME,EPSMAC
        WRITE (ISYSWR,10300)
      ENDIF
      NPFIX = 0
      NINT = 0
      NU = 0
      NPAR = 0
      IFATAL = 0
      IF(LST3.GE.5) WRITE (ISYSWR,10300)
      DO 20  I= 1, MAXEXT
        U(I) = 0.0
        NAM(I) = BLANK
        ERP(I) = 0.0
        ERN(I) = 0.0
        LCODE(I) = 0
   20 LCORSP (I) = 0
      UP = 1.0
      ISW(5) = 1
      IUNIT = ISYSRD
C                                        . . . READ PARAMETER CARDS . .
      ENTRY LMIDA2
      DO 150 I= 1, 200
        IF(I.GE.5) GOTO 160
        XK=XKIN(I)
        NAMK=NAMKIN(I)
        UK=UKIN(I)
        WK=WKIN(I)
        A=AIN(I)
        B=BIN(I)
        K = XK + 0.1
        NU = MAX0(NU,K)
        IF (K .LE. 0) GO TO 160
        IF (K .LE. MAXEXT) GO TO 30
        IFATAL = IFATAL + 1
        IF(LST3.GE.1) THEN
          WRITE (ISYSWR,10700) K,MAXEXT
          WRITE (ISYSWR,10000) K,NAMK,UK,WK,A,B
        ENDIF
        GO TO 150
   30   CONTINUE
        IF(NAM(K).EQ.BLANK) GO TO 40
C         PREVIOUSLY DEFINED PARAMETER IS BEING REDEFINED
        IF(LST3.GE.1) WRITE(ISYSWR,10500)
        IF(WERR(K).GT..0) NINT=NINT-1
   40   CONTINUE
        NAM(K) = NAMK
        U(K) = UK
        WERR(K) = WK
        IF (WK .GT. 0.0) GO TO 50
C                                        . . . FIXED PARAMETER . . . .
        IF(LST3.GE.5) WRITE (ISYSWR, 10000) K,NAMK,UK
        LCODE(K) = 0
        GO TO 140
C                                        . . . VARIABLE PARAMETER . . .
   50   IF(LST3.GE.5) WRITE (ISYSWR, 10000) K,NAMK,UK,WK,A,B
        NINT = NINT + 1
        ISW(2) = 0
        IF (A) 80 ,60 ,80
   60   IF (B) 80 ,70 ,80
   70   LCODE(K) = 1
        GO TO 140
   80   IF (B-A) 100,90 ,110
   90   IFATAL = IFATAL + 1
        IF(LST3.GE.1) WRITE (ISYSWR,10800)
        GO TO 110
  100   SAV = B
        B = A
        A = SAV
        IF(LST3.GE.1) WRITE (ISYSWR,10100)
  110   ALIM(K) = A
        BLIM(K) = B
        LCODE(K) = 4
        IF ((B-U(K))*(U(K)-A)) 120,130,140
  120   IFATAL = IFATAL + 1
        IF(LST3.GE.1) WRITE (ISYSWR,10900)
        GO TO 140
  130   IF(LST3.GE.1) WRITE (ISYSWR,10400)
  140   CONTINUE
  150 CONTINUE
      IFATAL = IFATAL + 1
      IF(LST3.GE.1) WRITE (ISYSWR,11000)
C                                       . . . END PARAMETER CARDS
C                                       . .   . STOP IF FATAL ERROR
  160 IF(LST3.GE.5) WRITE (ISYSWR,10300)
      IF (NINT .LE. MAXINT)  GO TO 170
      IF(LST3.GE.1) WRITE (ISYSWR,10600) NINT,MAXINT
      IFATAL = IFATAL + 1
  170 IF (IFATAL .LE. 0)  GO TO 180
      IF(LST3.GE.1) WRITE (ISYSWR,11100) IFATAL
      IF(LST3.GE.2) STOP
C                                       CALCULATE STEP SIZES DIRIN
  180 NPAR = 0
      DO 190 K= 1, NU
        IF (LCODE(K) .LE. 0) GO TO 190
        NPAR = NPAR + 1
        LCORSP(K) = NPAR
        SAV = U(K)
        X(NPAR) = LMPINT(SAV,K)
        XT(NPAR) = X(NPAR)
        SAV2 = SAV + WERR(K)
        VPLU = LMPINT(SAV2,K) - X(NPAR)
        SAV2 = SAV - WERR(K)
        VMINU = LMPINT(SAV2,K) - X(NPAR)
        DIRIN(NPAR) = 0.5 * (ABS(VPLU) +ABS(VMINU))
        G2(NPAR) = 2.0 / DIRIN(NPAR)**2
        GSTEP(NPAR) = DIRIN(NPAR)
        IF (LCODE(K) .GT. 1) GSTEP(NPAR) = -GSTEP(NPAR)
  190 CONTINUE
      SIGMA = 1.0E10
      IUNIT = ISYSRD
      RETURN
C...     THE FORMAT BELOW IS MACHINE-DEPENDENT. (A10) , (A6,4X) , ETC.
10000 FORMAT (I10,2X,A10,2X,2G12.6,2X,2G12.6)
10100 FORMAT(' WARNING           - ABOVE LIMITS HAVE BEEN REVERSED.')
10200 FORMAT (1H1/42X,21(1H*)/42X,21H*   D506   MINUIT   */42X,
     +12H* DIMENSIONS, I3, 1H/, I3, 2H */   42X,
     +'*  MODIFICATION OF  *',/,42X,
     +11H*   VERSION ,F6.2,4H   */42X,16H* DATA BLOCK NO. ,I3,2H *)
10300 FORMAT (4X,96(1H*))
10400 FORMAT(' WARNING           - ABOVE PARAMETER IS AT LIMIT ')
10500 FORMAT(' WARNING  *******  - PARAMETER REQUESTED ON FOLLOWING',
     +' CARD HAS ALREADY APPEARED.  PREVIOUS VALUES IGNORED.')
10600 FORMAT('0   TOO MANY VARIABLE PARAMETERS.  YOU REQUEST',I5/,
     +'   THIS VERSION OF MINUIT IS ONLY DIMENSIONED FOR',I4//)
10700 FORMAT('0FATAL ERROR. PARAMETER NUMBER',I11,' GREATER THAN ',
     +'ALLOWED MAXIMUM',I4)
10800 FORMAT(' FATAL ERROR. UPPER AND LOWER LIMITS ARE EQUAL.')
10900 FORMAT(' FATAL ERROR. PARAMETER OUTSIDE LIMITS',/)
11000 FORMAT('0FATAL ERROR. MORE THAN 200 PARAMETER CARDS',/)
11100 FORMAT(/I5,' FATAL ERRORS ON PARAMETER CARDS.  ABORT.',//)
11200 FORMAT(5X,A60,5X,'TIME',F8.3,' SECONDS',/,70X,'MACH. PREC.=',
     +E10.2)
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LMINEW

C...THIS IS THE MINUIT ROUTINE MINNEW.
CC        THIS IS THE MAIN PROGRAM, DISGUISED AS A SUBROUTINE FOR
CC        REASONS OF COMPATIBILITY BETWEEN SYSTEMS.    IT INITIALIZES
CC        SOME CONSTANTS IN COMMON (INCLUDING THE LOGICAL I/O UNIT NOS.)
CC        THEN VERIFIES THAT FCN GIVES THE SAME VALUE WHEN CALLED
CC        TWICE WITH THE SAME ARGUMENTS, AND PASSES CONTROL TO LMCMND.
CC

      COMMON /LPFLAG/ LST3
      COMMON
     +/LMMINE/ ERP(30)  ,ERN(30)
     +/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     +/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     +/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     +/LMVARI/ V(15,15)
     +/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     +/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     +/LMCASC/ Y(16)    ,JH       ,JL
     +/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     +/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     +/LMVART/ VT(15,15)
      COMMON
     +/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     +/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     +/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     +/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     +/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC

C         UNIT NUMBERS FOR CARD READER, PRINTER, PUNCH
C
      ISYSRD = 5
      ISYSWR = 6
      ISYSPU = 7
      MAXINT=15
      MAXEXT=30
C                   DETERMINE MACHINE ACCURACY EPSMAC
      EPSMAC = 0.5
      DO 10 I= 1, 100
        EPSMAC = EPSMAC * 0.5
        IF ((1.0+EPSMAC) .EQ. 1.0) GO TO 20
   10 CONTINUE
      EPSMAC = 1.0E-6
   20 EPSMAC = 2.0 * EPSMAC
C                             . . . . . . . . .
   30 CONTINUE
      NFCN = 1
      CALL LMIDAT
      CALL LMINTO(X)
      IF(LST3.GE.5) WRITE (ISYSWR,10000)
10000 FORMAT (/,'0FIRST ENTRY TO FCN ')
      CALL LSIGMX(NPAR,GIN,AMIN,U,1)
      CALL LSIGMX(NPAR,GIN,AMIN,U,4)
      CALL LMPRIN(1,AMIN)
      CALL LSIGMX(NPAR,GIN,F   ,U,4)
      IF  (F .NE. AMIN)  GO TO 40
      NFCN = 3
      CALL LMCMND
      RETURN
   40 CONTINUE
      IF(LST3.GE.1) WRITE (ISYSWR,10100) AMIN, F
      IF(LST3.GE.2) STOP
10100 FORMAT('0FOR THE ABOVE VALUES OF THE PARAMETERS, FCN IS TIME-',
     +'DEPENDENT',/,'0F = ',E22.14,' FOR FIRST CALL',/,' F =',E22.14,
     +' FOR SECOND')
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LMINTO(PINT)

C...THIS IS THE MINUIT ROUTINE INTOEX.
CC        TRANSFORMS FROM INTERNAL COORDINATES (PINT) TO EXTERNAL
CC        PARAMETERS (U).   THE MINIMIZING ROUTINES WHICH WORK IN
CC        INTERNAL COORDINATES CALL THIS ROUTINE BEFORE CALLING FCN.

      COMMON
     +/LMMINE/ ERP(30)  ,ERN(30)
     +/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     +/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     +/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     +/LMVARI/ V(15,15)
     +/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     +/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     +/LMCASC/ Y(16)    ,JH       ,JL
     +/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     +/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     +/LMVART/ VT(15,15)
      COMMON
     +/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     +/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     +/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     +/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     +/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC

      DIMENSION PINT(2)
      DO 30  I= 1, NU
        J = LCORSP(I)
        IF ( J ) 30 ,30 ,10
   10   CONTINUE
        IF (LCODE(I) .EQ. 1) GO TO 20
        AL = ALIM(I)
        U(I) = AL + 0.5 *(SIN(PINT(J)) +1.0) * (BLIM(I) -AL)
        GO TO 30
   20   U(I) = PINT(J)
   30 CONTINUE
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      REAL FUNCTION LMPINT(PEXTI,I)

C...THIS IS THE MINUIT ROUTINE PINTF.
CC        CALCULATES THE INTERNAL PARAMETER VALUE LMPINT CORRESPONDING
CC        TO THE EXTERNAL VALUE PEXTI FOR PARAMETER I.
CC
      COMMON
     1/LMMINE/ ERP(30)  ,ERN(30)
     2/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     3/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     4/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     5/LMVARI/ V(15,15)
     7/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     7/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     C/LMCASC/ Y(16)    ,JH       ,JL
     F/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     G/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     J/LMVART/ VT(15,15)
      COMMON
     6/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     8/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     9/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     A/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     B/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC
      COMMON /LPFLAG/ LST3
      DATA BIG, SMALL / 1.570796326795 , -1.570796326795 /
      IGO = LCODE(I)
      GO TO (10 ,20 ,30 ,40 ),IGO
C--       IGO = 1  MEANS NO LIMITS
   10 LMPINT = PEXTI
      GO TO 120
   20 CONTINUE
   30 CONTINUE
C--       IGO = 4  MEANS THERE ARE TWO LIMITS
   40 ALIMI = ALIM(I)
      BLIMI = BLIM(I)
      IF (PEXTI-ALIMI)  50 ,100,70
   50 A = SMALL
   60 LMPINT = A
      PEXTI = ALIMI + 0.5* (BLIMI-ALIMI) *(SIN(A) +1.0)
      LIMSET=1
      IF(LST3.GE.1) WRITE (ISYSWR,10000) I
      GO TO 120
   70 IF (BLIMI-PEXTI)  80 ,110,90
   80 A = BIG
      GO TO 60
   90 YY=2.0*(PEXTI-ALIMI)/(BLIMI-ALIMI) - 1.0
      LMPINT = ATAN(YY/SQRT(1.0- YY**2) )
      GO TO 120
  100 LMPINT = SMALL
      GO TO 120
  110 LMPINT = BIG
  120 RETURN
10000 FORMAT(' WARNING - VARIABLE',I3,' HAS BEEN BROUGHT BACK IN',
     +'SIDE LIMITS BY LMPINT.')
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LMPRIN  (IKODE,FVAL)

C...THIS IS THE MINUIT ROUTINE MPRINT.
CC        PRINTS THE VALUES OF THE PARAMETERS AT THE TIME OF THE CALL.
CC        ALSO PRINTS OTHER RELEVANT INFORMATION SUCH AS FUNCTION VALUE,
CC        ESTIMATED DISTANCE TO MINIMUM, PARAMETER ERRORS, STEP SIZES.
CC        ACCORDING TO THE VALUE OF IKODE,THE PRINTOUT IS LONG FORMAT,
CC        SHORT FORMAT, OR MINOS FORMAT (0,1,2)
CC

      COMMON
     +/LMMINE/ ERP(30)  ,ERN(30)
     +/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     +/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     +/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     +/LMVARI/ V(15,15)
     +/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     +/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     +/LMCASC/ Y(16)    ,JH       ,JL
     +/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     +/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     +/LMVART/ VT(15,15)
      COMMON
     +/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     +/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     +/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     +/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     +/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC
      COMMON /LMINUC/ NAMKIN(4),NAM(30)
      COMMON /LPFLAG/ LST3
      CHARACTER*10 NAMKIN,NAM
C                                        . GET TIME AND PRINT HEADINGS .
      CALL LTIMEX(TI)
      IF(LST3.GE.5) WRITE (ISYSWR,10000)
      E = SIGMA
      KOUNT = 0
C                                        . . . LOOP OVER PARAMETERS . .
      DO 110 I= 1, NU
        IF(NAM(I).EQ.'          ') GOTO 110
   10   L = LCORSP(I)
        IF (L .EQ. 0) GO TO 80
C              VARIABLE PARAMETER.  CALCULATE EXTERNAL ERROR IF V EXISTS
        IF (ISW(2) .LT. 1) GO TO 30
        DX = SQRT(ABS(V(L,L)*UP))
        IF (LCODE(I) .LE. 1) GO TO 20
        AL = ALIM(I)
        BA = BLIM(I) - AL
        DU1 = AL + 0.5 *(SIN(X(L)+DX) +1.0) * BA - U(I)
        DU2 = AL + 0.5 *(SIN(X(L)-DX) +1.0) * BA - U(I)
        IF (DX .GT. 1.0) DU1 = BA
        DX = 0.5 * (ABS(DU1) + ABS(DU2))
   20   WERR(I) = DX
   30   X1 = X(L)
        X2 = DIRIN(L)
        IF (IKODE .LT. 2) GO TO 40
        X1 = ERP(I)
        X2 = ERN(I)
   40   IF (KOUNT) 50,50,60
   50   KOUNT = 1
        IF(LST3.GE.5) WRITE (ISYSWR,10100) FVAL,NFCN,TI,E, L,I,NAM(I),
     +  U(I),WERR(I),X1,X2
        GO TO 70
   60   IF(LST3.GE.5) WRITE (ISYSWR,10200) L,I,NAM(I),U(I),WERR(I),X1,
     +  X2
   70   IF (LCODE(I) .LE. 1) GO TO 110
        IF(LST3.GE.1.AND. ABS(COS(X(L))) .LT. 0.001) WRITE (ISYSWR,
     +  10400)
        GO TO 110
C                           FIXED PARAMETER.  PRINT ONLY IF IKODE .GT.0
   80   IF (IKODE .EQ. 0) GO TO 110
        IF (KOUNT) 90,90,100
   90   KOUNT = 1
        IF(LST3.GE.5) WRITE (ISYSWR,10100) FVAL,NFCN,TI,E, L,I,NAM(I),
     +  U(I)
        GO TO 110
  100   IF(LST3.GE.5) WRITE (ISYSWR,10300) I,NAM(I),U(I)
  110 CONTINUE
      IF(LST3.GE.5.AND. IKODE.GE.1 .AND.ISW(2).GE.1) WRITE (ISYSWR,
     +10500) UP
      RETURN
10000 FORMAT(/ 4X,'FCN VALUE',5X,'CALLS',4X,'TIME',4X,' EDM  ',4X ,
     +'INT.EXT. PARAMETER     VALUE         ERROR      INTERN.VALUE  ',
     +'INT.STEP SIZE')
10100 FORMAT(E15.7,I7,F9.2,E11.2,I6,I4,1X,A10,4E14.5)
10200 FORMAT(1H ,41X,I6,I4,1X,A10,4E14.5)
10300 FORMAT(1H  ,47X  ,I4,1X,A10,4E14.5)
10400 FORMAT(1H ,52X  ,'WARNING -   - ABOVE PARAMETER IS AT LIMIT.')
10500 FORMAT(/45X,'ERRORS CORRESPOND TO FUNCTION CHANGE OF ',E12.4)
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 12/09/94  16.18.19  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LMRAZZ(YNEW,PNEW)

C...THIS IS THE MINUIT ROUTINE RAZZIA.
CC        CALLED ONLY BY SIMPLEX (AND IMPROV) TO ADD A NEW POINT
CC        AND REMOVE AN OLD ONE FROM THE CURRENT SIMPLEX, AND GET THE
CC        ESTIMATED DISTANCE TO MINIMUM.
CC
      COMMON
     +/LMMINE/ ERP(30)  ,ERN(30)
     +/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     +/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     +/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     +/LMVARI/ V(15,15)
     +/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     +/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     +/LMCASC/ Y(16)    ,JH       ,JL
     +/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     +/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     +/LMVART/ VT(15,15)
      COMMON
     +/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     +/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     +/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     +/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     +/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC
      COMMON /LPFLAG/ LST3
      DIMENSION PNEW(15)
      DO 10 I=1,NPAR
   10 P(I,JH)=PNEW(I)
      Y(JH)=YNEW
      IF(YNEW.GE.AMIN) GO TO 30
      DO 20 I=1,NPAR
   20 X(I)=PNEW(I)
      CALL LMINTO(X)
      AMIN=YNEW
      JL=JH
   30 CONTINUE
      JH=1
      NPARP1=NPAR+1
   40 DO 50 J=2,NPARP1
        IF (Y(J) .GT. Y(JH)) JH = J
   50 CONTINUE
      SIGMA = Y(JH) - Y(JL)
      IF (SIGMA .LE. 0.)  GO TO 90
      US = 1.0/SIGMA
      DO 70 I= 1, NPAR
        PBIG = P(I,1)
        PLIT = PBIG
        DO 60 J= 2, NPARP1
          IF (P(I,J) .GT. PBIG) PBIG = P(I,J)
          IF (P(I,J) .LT. PLIT) PLIT = P(I,J)
   60   CONTINUE
        DIRIN(I) = PBIG - PLIT
        IF (ITAUR .LT. 1 ) V(I,I) = 0.5*(V(I,I) +US*DIRIN(I)**2)
   70 CONTINUE
   80 RETURN
   90 IF(LST3.GE.1.AND.MOD(ITOO,10).EQ.0) THEN
        WRITE (ISYSWR, 10000) NPAR
        ITOO=ITOO+1
      ENDIF
      GO TO 80
10000 FORMAT('0***** FUNCTION VALUE DOES NOT SEEM TO DEPEND ON ANY ',
     +'OF THE',I3,' VARIABLE PARAMETERS',/15X ,'VERIFY THAT STEP SIZES',
     +' ARE BIG ENOUGH AND CHECK FCN LOGIC.',/1X,81(1H*)/1X,81(1H*)//)
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 20/09/94  17.12.04  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LMSIMP

C...THIS IS THE MINUIT ROUTINE SIMPLEX.
CC        PERFORMS A MINIMIZATION USING THE SIMPLEX METHOD OF NELDER
CC        AND MEAD (REF. -- COMP. J. 7,308 (1965)).
      COMMON /LMINUI/ XKIN(4),UKIN(4),WKIN(4),AIN(4),BIN(4),
     +MAXFIN,RELUP,RELERR,RELER2,FCNMAX
      COMMON /LPFLAG/ LST3
      COMMON
     +/LMMINE/ ERP(30)  ,ERN(30)
     +/LMPARI/ X(15)    ,XT(15)   ,DIRIN(15) ,MAXINT     ,NPAR
     +/LMPARE/ U(30)              ,WERR(30)  ,MAXEXT     ,NU
     +/LMLIMI/ ALIM(30) ,BLIM(30) ,LCODE(30) ,LCORSP(30) ,LIMSET
     +/LMVARI/ V(15,15)
     +/LMFIX / IPFIX(15),XS(15)   ,XTS(15)   ,DIRINS(15) ,NPFIX
     +/LMFIX2/ GRDS(15) ,G2S(15)  ,GSTEPS(15),ABERFS(15)
     +/LMCASC/ Y(16)    ,JH       ,JL
     +/LMDERI/ GIN(30)  ,GRD(15)  ,G2(15)    ,GSTEP(15)  ,ABERF(15)
     +/LMSIMV/ P(15,16) ,PSTAR(15),PSTST(15) ,PBAR(15)   ,PRHO(15)
     +/LMVART/ VT(15,15)
      COMMON
     +/LMUNIT/ ISYSRD   ,ISYSWR   ,ISYSPU
     +/LMTITL/ TITLE(13),DATE(2)  ,ISW(7)    ,NBLOCK
     +/LMCONV/ EPSI     ,APSI     ,VTEST     ,NSTEPQ     ,NFCN ,NFCNMX
     +/LMCARD/ CWORD    ,CWORD2   ,CWORD3    ,WORD7(7)
     +/LMMINI/ AMIN     ,UP       ,NEWMIN    ,ITAUR      ,SIGMA,EPSMAC

      DATA ALPHA,BETA,GAMMA,RHOMIN,RHOMAX / 1.0, 0.5, 2.0, 4.0, 8.0/
      ALPHA=1.0
      BETA=0.5
      GAMMA=2.0
      RHOMIN=4.0
      RHOMAX=8.0

      IF (NPAR .LE. 0)  RETURN
      NPFN=NFCN
      NPARP1=NPAR+1
      RHO1 = 1.0 + ALPHA
      RHO2 = RHO1 + ALPHA*GAMMA
      WG = 1.0/FLOAT(NPAR)
      IFLAG=4
      IF(LST3.GE.5) WRITE(ISYSWR,10000) EPSI
      DO 10 I= 1, NPAR
        IF (ISW(2) .GE. 1) DIRIN(I) = SQRT(V(I,I)*UP)
        IF (ABS(DIRIN(I)) .LT. 1.0E-10*ABS(X(I))) DIRIN(I)=1.0E-8*X(I)
        IF(ITAUR.LT. 1) V(I,I) = DIRIN(I)**2/UP
   10 CONTINUE
      IF (ITAUR .LT. 1)  ISW(2) = 1
C**       CHOOSE THE INITIAL SIMPLEX USING SINGLE-PARAMETER SEARCHES
   20 CONTINUE
      YNPP1 = AMIN
      JL = NPARP1
      Y(NPARP1) = AMIN
      ABSMIN = AMIN
      DO 70 I= 1, NPAR
        AMING = AMIN
        PBAR(I) = X(I)
        BESTX = X(I)
        KG = 0
        NS = 0
        NF = 0
   30   X(I) = BESTX + DIRIN(I)
        CALL LMINTO(X)
        CALL LSIGMX(NPAR,GIN, F, U, 4)
        NFCN = NFCN + 1
        IF (F .LE. AMING) GO TO 40
C         FAILURE
        IF (KG .EQ. 1) GO TO 50
        KG = -1
        NF = NF + 1
        DIRIN(I) = DIRIN(I) * (-0.4)
        IF (NF .LT. 3) GO TO 30
        NS = 6
C         SUCCESS
   40   BESTX = X(I)
        DIRIN(I) = DIRIN(I) * 3.0
        AMING = F
        KG = 1
        NS = NS + 1
        IF (NS .LT. 6) GO TO 30
C         LOCAL MINIMUM FOUND IN ITH DIRECTION
   50   Y(I) = AMING
        IF (AMING .LT. ABSMIN) JL = I
        IF (AMING .LT. ABSMIN) ABSMIN = AMING
        X(I) = BESTX
        DO 60 K= 1, NPAR
   60   P(K,I) = X(K)
   70 CONTINUE
      JH = NPARP1
      AMIN=Y(JL)
      CALL LMRAZZ(YNPP1,PBAR)
      DO 80 I= 1, NPAR
   80 X(I) = P(I,JL)
      CALL LMINTO(X)
      DO 90 I=1,NPAR
   90 IF(ABS(DIRIN(I)).LE.ABS(EPSMAC*X(I))) DIRIN(I)=4.*EPSMAC*X(I)
      IF (ISW(5) .GE. 1)  CALL LMPRIN(0,AMIN)
      SIGMA = SIGMA * 10.
      SIG2 = SIGMA
      NCYCL=0
C                                        . . . . .  START MAIN LOOP
  100 CONTINUE
C...CHANGE IN SIMPLX; ERROR REDEFINED FOR SECOND CALL TO LMSIMP.
      UP=RELUP*ABS(AMIN)
      EPSI=RELERR*UP
      IF (SIG2 .LT. EPSI .AND. SIGMA.LT.EPSI) GO TO 220
      SIG2 = SIGMA
      IF ((NFCN-NPFN) .GT. NFCNMX) GO TO 230
C         CALCULATE NEW POINT * BY REFLECTION
      DO 120 I= 1, NPAR
        PB = 0.
        DO 110 J= 1, NPARP1
  110   PB = PB + WG * P(I,J)
        PBAR(I) = PB - WG * P(I,JH)
  120 PSTAR(I)=(1.+ALPHA)*PBAR(I)-ALPHA*P(I,JH)
      CALL LMINTO(PSTAR)
      CALL LSIGMX(NPAR,GIN,YSTAR,U,4)
      NFCN=NFCN+1
      IF(YSTAR.GE.AMIN) GO TO 190
C         POINT * BETTER THAN JL, CALCULATE NEW POINT **
      DO 130 I=1,NPAR
  130 PSTST(I)=GAMMA*PSTAR(I)+(1.-GAMMA)*PBAR(I)
      CALL LMINTO(PSTST)
      CALL LSIGMX(NPAR,GIN,YSTST,U,4)
      NFCN=NFCN+1
C         TRY A PARABOLA THROUGH PH, PSTAR, PSTST.  MIN = PRHO
      Y1 = (YSTAR-Y(JH)) * RHO2
      Y2 = (YSTST-Y(JH)) * RHO1
      RHO = 0.5 * (RHO2*Y1 -RHO1*Y2) / (Y1 -Y2)
      IF (RHO .LT. RHOMIN) GO TO 160
      IF (RHO .GT. RHOMAX)  RHO = RHOMAX
      DO 140 I= 1, NPAR
  140 PRHO(I) = RHO*PBAR(I) + (1.0-RHO)*P(I,JH)
      CALL LMINTO(PRHO)
      CALL LSIGMX(NPAR,GIN,YRHO, U,4)
      NFCN = NFCN + 1
      IF (YRHO .LT. Y(JL) .AND. YRHO .LT. YSTST) GO TO 150
      IF (YSTST .LT. Y(JL)) GO TO 170
      IF (YRHO .GT. Y(JL)) GO TO 160
C         ACCEPT MINIMUM POINT OF PARABOLA, PRHO
  150 CALL LMRAZZ (YRHO,PRHO)
      GO TO 180
  160 IF (YSTST .LT. Y(JL)) GO TO 170
      CALL LMRAZZ(YSTAR,PSTAR)
      GO TO 180
  170 CALL LMRAZZ(YSTST,PSTST)
  180 NCYCL=NCYCL+1
      IF (ISW(5) .LT. 2) GO TO 100
      IF (ISW(5) .GE. 3 .OR. MOD(NCYCL, 10) .EQ. 0) CALL LMPRIN(0,AMIN)
      GO TO 100
C         POINT * IS NOT AS GOOD AS JL
  190 IF (YSTAR .GE. Y(JH)) GO TO 200
      JHOLD = JH
      CALL LMRAZZ(YSTAR,PSTAR)
      IF (JHOLD .NE. JH) GO TO 100
C         CALCULATE NEW POINT **
  200 DO 210 I=1,NPAR
  210 PSTST(I)=BETA*P(I,JH)+(1.-BETA)*PBAR(I)
      CALL LMINTO (PSTST)
      CALL LSIGMX(NPAR,GIN,YSTST,U,4)
      NFCN=NFCN+1
      IF(YSTST.GT.Y(JH)) GO TO 20
C     POINT ** IS BETTER THAN JH
      IF (YSTST .LT. AMIN) GO TO 170
      CALL LMRAZZ(YSTST,PSTST)
      GO TO 100
C                                        . . . . . .  END MAIN LOOP
  220 IF(LST3.GE.5) WRITE(ISYSWR,10100)
      GO TO 240
  230 IF(LST3.GE.5) WRITE(ISYSWR,10200)
      ISW(1) = 1
  240 DO 260 I=1,NPAR
        PB = 0.
        DO 250 J=1,NPARP1
  250   PB = PB + WG * P(I,J)
  260 PBAR(I) = PB - WG * P(I,JH)
      CALL LMINTO(PBAR)
      CALL LSIGMX(NPAR,GIN,YPBAR,U,IFLAG)
      NFCN=NFCN+1
      IF (YPBAR .LT. AMIN)  CALL LMRAZZ(YPBAR,PBAR)
      CALL LMINTO(X)
      IF (NFCNMX+NPFN-NFCN .LT. 3*NPAR) GO TO 270
      IF (SIGMA .GT. 2.0*EPSI) GO TO 20
  270 CALL LMPRIN(1-ITAUR, AMIN)
      RETURN
10000 FORMAT(' START SIMPLEX MINIMIZATION          ',8X   ,'CON',
     +'VERGENCE CRITERION -- ESTIMATED DISTANCE TO MINIMUM (EDM) .LT.',
     +E10.2 )
10100 FORMAT(1H ,'SIMPLEX MINIMIZATION HAS CONVERGED')
10200 FORMAT(1H ,'SIMPLEX TERMINATES WITHOUT CONVERGENCE')
      END
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C######################################################################
C
C   VARIOUS ROUTINES TO GIVE STRUCTURE FUNCTION PARAMETRIZATIONS.
C   ALL BUT LNSTRF CAN BE USED SEPARATELY WITHOUT INITIALIZATION.
C
C ********************************************************************

      SUBROUTINE LNSTRF(X,Q2,XPQ)

C...STRUCTURE FUNCTION PER NUCLEON FOR A PROTON/NEUTRON MIXTURE
C...ACCORDING TO DEFINED NUCLEUS.

      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      DIMENSION XPQ(-6:6)

      CALL PYSTFU(2212,X,Q2,XPQ)

      IF(PARI(11).LE.1.E-06) RETURN
      XDV=XPQ(1)-XPQ(-1)
      XUV=XPQ(2)-XPQ(-2)
C...FOR NUCLEAR TARGET, MIX U- AND D-VALENCE DISTRIBUTIONS.
      XPQ(1)=(1.-PARI(11))*XDV+PARI(11)*XUV + XPQ(-1)
      XPQ(2)=(1.-PARI(11))*XUV+PARI(11)*XDV + XPQ(-2)
C...SAVE D AND U VALENCE IN PROTON
      PARI(12)=XDV
      PARI(13)=XUV

      RETURN
      END
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LPRIKT(S,PT,PHI)

C...SIZE (PT) AND AZIMUTHAL ANGLE (PHI) OF PRIMORDIAL KT ACCORDING
C...TO A GAUSSIAN DISTRIBUTION.

      PT=S*SQRT(-ALOG(RLU(0)))
      PHI=6.2832*RLU(0)
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LPRWTS(NSTEP)

C...PRINTS PROBABILITIES FOR Q-, QG- AND QQBAR-EVENTS USING THE PRESENT
C...QCD WEIGHTS STORED IN COMMON BLOCK LGRID.

      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LGRID/ NXX,NWW,XX(20),WW(15),PQG(20,15,3),PQQB(20,15,2),
     +QGMAX(20,15,3),QQBMAX(20,15,2),YCUT(20,15),XTOT(20,15),NP
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ

      WMAX=SQRT(PARL(21)+PSAVE(3,1,5)**2+PSAVE(3,2,5)**2)
      WRITE(6,10000) PARL(11),LST(13),MSTU(112),PARU(112), PARL(8),
     +PARL(9),PARL(12),PARL(13)
      IF(NP.EQ.1) THEN
        WRITE(6,10100)
      ELSE
        WRITE(6,10200)
      ENDIF
      WRITE(6,10300) LST(19),NWW,NXX,WW,XX
      IF(WMAX.GT.WW(NWW)) WRITE(6,10400) WMAX,WW(NWW)
      WRITE(6,10500)

      LW=0
      DO 30  IW=1,NWW,MAX(1,NSTEP)
        W=WW(IW)
        IF(LW.GT.0) GOTO 40
        IF(W.GT.WMAX) LW=LW+1
        W2=W**2
        LX=0
        DO 20  IX=1,NXX,MAX(1,NSTEP)
          X=XX(IX)
          IF(LX.GT.0) GOTO 30
          U=(W2-PSAVE(3,2,5)**2)/(2.*PSAVE(3,2,5)*(1.-X))
          Q2=2.*PSAVE(3,2,5)*U*X
          Y=Q2/(PARL(21)*X)
          PARI(24)=(1.+(1.-Y)**2)/2.
          PARI(25)=1.-Y
          PARI(26)=(1.-(1.-Y)**2)/2.
          PARL(25)=ULALPS(Q2)
          IF(Y.GT.1.) LX=LX+1
          RQG=0.
          RQQB=0.
          DO 10  IP=1,NP
            IF(NP.EQ.1) THEN
              RQG=PQG(IX,IW,IP)
              RQQB=PQQB(IX,IW,IP)
            ELSE
              RQG=RQG+PQG(IX,IW,IP)*PARI(23+IP)/XTOT(IX,IW)
              IF(IP.LT.3) RQQB=RQQB+PQQB(IX,IW,IP)*PARI(23+IP)/XTOT(IX,
     +        IW)
            ENDIF
   10     CONTINUE
C...INCLUDE ALPHA-STRONG IN WEIGHT.
          RQG=RQG*PARL(25)
          RQQB=RQQB*PARL(25)
          IF(LST(39).EQ.-91) THEN
C...INCLUDE 3-JET CROSS SECTION IN DENOMINATOR
            QTOT=1.+RQG+RQQB
            RQG =RQG/QTOT
            RQQB=RQQB/QTOT
          ENDIF
          RQ=1.-RQG-RQQB
          WRITE(6,10600) W,X,Y,Q2,PARL(25),YCUT(IX,IW),RQ,RQG,RQQB
   20   CONTINUE
   30 CONTINUE
   40 CONTINUE
      RETURN

10000 FORMAT('1',/,5X,'SUMMARY OF QCD MATRIX ELEMENT INTEGRATION',
     +           /,5X,'-----------------------------------------',//,
     +/,' FOR GLUON RADIATION (QG-EVENT) AND BOSON-GLUON FUSION ',
     +'(QQ-EVENT) PROBABILITY.',
     +//,' REQUIRED PRECISION IN INTEGRATION, PARL(11) =',F8.4,
     +//,' HEAVIEST FLAVOUR PRODUCED IN BOSON-GLUON FUSION, LST(13) =',
     +I5,//,' ALPHA-STRONG PARAMETERS: # FLAVOURS, MSTU(112) =',I3,
     +'  QCD LAMBDA, PARU(112) =',F6.3,' GEV',
     +//,' CUTS ON MATRIX ELEMENTS:',
     +/,' PARL(8), PARL(9), PARL(12), PARL(13) =',4F8.4,/)
10100 FORMAT(' LEPTON ENERGY NOT ALLOWED TO VARY IN SIMULATION.',/)
10200 FORMAT(' LEPTON ENERGY ALLOWED TO VARY IN SIMULATION, ',/,
     +' Y IN TABLE BELOW CALCULATED ASSUMING MAX ENERGY.',/)
10300 FORMAT(' GRID CHOICE, LST(19) =',I3,5X,'# GRID POINTS IN W, X =',
     +2I5,/,' W-VALUES IN ARRAY WW:',/,10F8.1,/,5F8.1,
     +/,' X-VALUES IN ARRAY XX:',/,10F8.4,/,10F8.4,/)
10400 FORMAT(' MAX W OUTSIDE GRID, EXECUTION STOPPED ] WMAX, GRID-MAX ='
     +,2F12.1)
10500 FORMAT(//,6X,'W',7X,'X',7X,'Y',8X,'Q**2',3X,'ALPHA',
     +5X,'CUT',2X,'Q-EVENT',1X,'QG-EVENT',1X,'QQ-EVENT',
     +/,1X,77(1H-),/)
10600 FORMAT(F7.1,2F8.4,1PG12.3,0PF8.2,F8.4,3F9.4)
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LQCDPR(QG,QQB)

C...PROBABILITIES FOR HARD QCD EVENTS, FOR GIVEN X AND W, ARE OBTAINED
C...BY MAKING A LINEAR INTERPOLATION OF THE VALUES ON THE X-W GRID.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LGRID/ NXX,NWW,XX(20),WW(15),PQG(20,15,3),PQQB(20,15,2),
     +QGMAX(20,15,3),QQBMAX(20,15,2),YCUT(20,15),XTOT(20,15),NP
      DATA NOUT,NABOVE/2*0/,NWARN/10/

      QG=0.
      QQB=0.
      W=SQRT(W2)
C...IF W IS VERY SMALL OR X CLOSE TO 1, SET QCD WEIGHTS TO ZERO.
      IF(WW(1).LT.6..AND.W.LT.WW(1)) RETURN
      IF(X.GT.XX(NXX)) RETURN

      XP=X
      IF(X.LT.XX(1).OR.X.GT.XX(NXX).OR.
     +W.LT.WW(1).OR.W.GT.WW(NWW)) THEN
C...X AND/OR W OUTSIDE LIMITS OF GRID, WRITE WARNING FOR
C...FIRST NWARN CASES.
        NOUT=NOUT+1
        IF(LST(3).GE.1.AND.NOUT.LE.NWARN) WRITE(6,10000) X,W,
     +  INT(PARI(29)),NWARN
        IF(X.LT.XX(1)) XP=XX(1)
        IF(X.GT.XX(NXX)) XP=XX(NXX)
        IF(W.LT.WW(1)) W=WW(1)
        IF(W.GT.WW(NWW)) W=WW(NWW)
      ENDIF

      IH=1
      IF(LST(30).EQ.1) IH=2
      IX=0
   10 IX=IX+1
      IF(XP.GT.XX(IX+1)) GOTO 10
      IW=0
   20 IW=IW+1
      IF(W.GT.WW(IW+1)) GOTO 20
      WD=(W-WW(IW))/(WW(IW+1)-WW(IW))
      XD=(XP-XX(IX))/(XX(IX+1)-XX(IX))

      DO 30  IP=1,NP
        X1P=(PQG(IX+1,IW,IP)-PQG(IX,IW,IP))*XD+PQG(IX,IW,IP)
        X2P=(PQG(IX+1,IW+1,IP)-PQG(IX,IW+1,IP))*XD+PQG(IX,IW+1,IP)
        QGIP=(X2P-X1P)*WD+X1P
        IF(NP.EQ.1) THEN
          QG=QGIP
          PARI(15)=MAX(QGMAX(IX,IW,IH),QGMAX(IX+1,IW+1,IH), QGMAX(IX+1,
     +    IW,IH),QGMAX(IX,IW+1,IH))
        ELSE
          QG=QG+PARI(23+IP)*QGIP
          PARI(14+IP)=MAX(QGMAX(IX,IW,IP),QGMAX(IX+1,IW+1,IP), QGMAX(IX
     +    +1,IW,IP),QGMAX(IX,IW+1,IP))
        ENDIF
        IF(IP.EQ.3) GOTO 30
        X1P=(PQQB(IX+1,IW,IP)-PQQB(IX,IW,IP))*XD+PQQB(IX,IW,IP)
        X2P=(PQQB(IX+1,IW+1,IP)-PQQB(IX,IW+1,IP))*XD+PQQB(IX,IW+1,IP)
        QQBIP=(X2P-X1P)*WD+X1P
        IF(NP.EQ.1) THEN
          QQB=QQBIP
          PARI(18)=MAX(QQBMAX(IX,IW,IH),QQBMAX(IX+1,IW+1,IH), QQBMAX(IX
     +    +1,IW,IH),QQBMAX(IX,IW+1,IH))
        ELSE
          QQB=QQB+PARI(23+IP)*QQBIP
          PARI(17+IP)=MAX(QQBMAX(IX,IW,IP),QQBMAX(IX+1,IW+1,IP),
     +    QQBMAX(IX+1,IW,IP),QQBMAX(IX,IW+1,IP))
        ENDIF
   30 CONTINUE

      IF(NP.NE.1) THEN
C...GET TOTAL X-SECTION FROM INTERPOLATION TO BE USED FOR NORMALIZATION.
        X1P=(XTOT(IX+1,IW)-XTOT(IX,IW))*XD+XTOT(IX,IW)
        X2P=(XTOT(IX+1,IW+1)-XTOT(IX,IW+1))*XD+XTOT(IX,IW+1)
        PQ17=(X2P-X1P)*WD+X1P
        QG=QG/PQ17
        QQB=QQB/PQ17
      ENDIF

C...GET VALUE OF Y-CUT, IE MINIMUM SCALED INVARIANT MASS SQUARED.
      PARL(27)=MAX(YCUT(IX,IW),YCUT(IX+1,IW+1),
     +YCUT(IX+1,IW),YCUT(IX,IW+1))

C...INCLUDE ALPHA-STRONG IN WEIGHT.
      QG=QG*PARL(25)
      QQB=QQB*PARL(25)
      IF(LST(39).EQ.-91) THEN
C...INCLUDE 3-JET CROSS SECTION IN DENOMINATOR
        QTOT=1.+QG+QQB
        QG =QG/QTOT
        QQB=QQB/QTOT
      ENDIF
      IF(QG+QQB.GT.1) THEN
C...SUM OF QCD EVENT PROBABILITIES LARGER THAN UNITY, RESCALE TO UNITY
C...AND PRINT WARNING FOR FIRST NWARN CASES.
        NABOVE=NABOVE+1
        IF(LST(3).GE.1.AND.NABOVE.LE.NWARN) WRITE(6,10100) QG,QQB,X,W,
     +  INT(PARI(29)),NWARN
        QGQQB=QG+QQB
        QG=QG/QGQQB
        QQB=QQB/QGQQB
      ENDIF

10000 FORMAT(' WARNING: X=',F7.4,' OR W=',F6.1,' OUTSIDE QCD GRID',
     +' IN EVENT NO.',I8,/,10X,
     +'WEIGHT ON LIMIT OF GRID USED. ONLY FIRST',I5,' WARNINGS PRINTED')
10100 FORMAT(' WARNING: SUM OF QCD PROBABILITIES LARGER THAN UNITY ',
     +' QG, QQB =',2F8.4,/10X,'OCCURS AT X, W =',F8.4,F6.1,
     +' IN EVENT NO.',I8,/,10X,
     +'WEIGHTS RESCALED TO UNIT SUM. ONLY FIRST',I5,' WARNINGS PRINTED')
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LQGEV

C...GENERATE QUARK-GLUON JET EVENT, CHOOSE XP AND ZP ACCORDING TO QCD
C...MATRIX ELEMENTS AND APPLY CUTS FOR SOFT AND COLLINEAR GLUONS.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)

      LST(24)=2
      W=SQRT(W2)
      J1=MSTU(1)
      J2=MSTU(1)+1
      J3=MSTU(1)+2
      J4=MSTU(1)+3

      CALL LXP(XP,IFAIL)
      IF(IFAIL.NE.0) GOTO 30

C...CHOOSE FLAVOUR OF SCATTERED QUARK AND TARGET REMNANT.
   10 CALL LFLAV(IFL,IFLR)
      IF(LST(21).NE.0) RETURN
      CALL LZP(XP,ZP,IFAIL)
      IF(IFAIL.NE.0) GOTO 30
      AMIFL=ULMASS(IFL)
      AMIFLR=ULMASS(IFLR)

      IF(LST(14).EQ.0.OR.IFLR.GT.10
     +.OR.(LST(8).GE.2.AND.MOD(LST(8),10).NE.9)) THEN
        IF(W.LT.AMIFL+AMIFLR+PARJ(32)) GOTO 30
        IF(LQMCUT(XP,ZP,AMIFL,0.,AMIFLR).NE.0) GOTO 30
        CALL LU3ENT(J1,IFL,21,IFLR,W,PARI(21),PARI(23))
        K(MSTU(1)+2,3)=2
        CALL LUROBO(ACOS(-P(J3,3)/SQRT(P(J3,3)**2+P(J3,1)**2)),
     +  0.,0.,0.,0.)
      ELSE
C...TARGET REMNANT IS NOT A SIMPLE DIQUARK, SPECIAL TREATMENT NEEDED.
        IF(W.LT.AMIFL+AMIFLR+1.+PARJ(32)) GOTO 30
        IF(LQMCUT(XP,ZP,AMIFL,0.,1.).NE.0) GOTO 30
        IFLRO=IFLR
        NREMH=0
   20   NREMH=NREMH+1
        IF(NREMH.GT.100) GOTO 30
        CALL LREMH(IFLRO,IFLR,K2,XT)
        AMK2=ULMASS(K2)
        AMIFLR=ULMASS(IFLR)
        P(J1,5)=AMIFL
        P(J2,5)=0.
        CALL LPRIKT(PARL(14),PT,PHI)
        PT2=PT**2
        TM2K2=AMK2**2+PT2
        TMIFLR=AMIFLR**2+PT2
        P(J3,5)=SQRT(TM2K2/XT+TMIFLR/(1.-XT))
        IF(LQMCUT(XP,ZP,AMIFL,0.,P(J3,5)).NE.0) GOTO 20
        MSTU(10)=1
        CALL LU3ENT(J1,IFL,21,IFLR,W,PARI(21),PARI(23))
        K(MSTU(1)+2,3)=2
        MSTU(10)=2
        CALL LUROBO(ACOS(-P(J3,3)/SQRT(P(J3,3)**2+P(J3,1)**2)),
     +  0.,0.,0.,0.)
        EPZ=P(J3,4)-P(J3,3)
        P(J3,1)=PT*COS(PHI)
        P(J3,2)=PT*SIN(PHI)
        P(J3,3)=-0.5*((1.-XT)*EPZ-TMIFLR/(1.-XT)/EPZ)
        P(J3,4)= 0.5*((1.-XT)*EPZ+TMIFLR/(1.-XT)/EPZ)
        P(J3,5)=AMIFLR
        P(J4,1)=-P(J3,1)
        P(J4,2)=-P(J3,2)
        P(J4,3)=-0.5*(XT*EPZ-TM2K2/XT/EPZ)
        P(J4,4)= 0.5*(XT*EPZ+TM2K2/XT/EPZ)
        P(J4,5)=AMK2
        K(J4,1)=1
        K(J4,2)=K2
        K(J4,3)=2
        K(J4,4)=0
        K(J4,5)=0
        N=J4
        IF((P(J3,4)+P(J2,4)/2.)**2-(P(J3,1)+P(J2,1)/2.)**2-P(J3,2)**2
     +  -(P(J3,3)+P(J2,3)/2.)**2.LT.(AMIFLR+2.5*PARJ(32))**2) GOTO 20
      ENDIF

      CALL LAZIMU(XP,ZP)
      LST(21)=0
      RETURN

   30 LST(21)=1
      RETURN
      END
*CMZ :  1.00/00 24/07/94  15.46.44  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LQEV

C...GENERATE AN ORDINARY 2-JET EVENT, Q-EVENT.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)

      LST(24)=1
      W=SQRT(W2)

C...CHOOSE FLAVOUR OF SCATTERED QUARK AND TARGET REMNANT.
   10 CALL LFLAV(IFL,IFLR)
      IF(LST(21).NE.0) GOTO 10

      GOTO 20
C...ENTRY USED FOR ARIADNE
      ENTRY LQEVAR(IFLAR,IFLRAR)
      IFL=IFLAR
      IFLR=IFLRAR
      LST(24)=1
      W=SQRT(W2)

   20 CONTINUE
      IF(LST(14).EQ.0.OR.IFLR.GT.10
     +.OR.(LST(8).GE.2.AND.MOD(LST(8),10).NE.9)) THEN
C...CHECK IF ENERGY IN JET SYSTEM IS ENOUGH FOR FRAGMENTATION.
C...PARJ(32) DEFAULTS TO 1 GEV
        IF(W.LT.ULMASS(IFL)+ULMASS(IFLR)+PARJ(32)) GOTO 10
        CALL LU2ENT(MSTU(1),IFL,IFLR,W)
        K(MSTU(1)+1,3)=2
*        WRITE(*,*)'POPPING LU2ENT'
      ELSE
*        WRITE(*,*)'REMNANTS TREATEMENT'
C...TARGET REMNANT IS NOT A SIMPLE DIQUARK, SPECIAL TREATMENT NEEDED.
        AMIFL=ULMASS(IFL)
        IF(W.LT.AMIFL+ULMASS(IFLR)+0.9+PARJ(32)) GOTO 10
        IFLRO=IFLR
        NREMH=0
   30   NREMH=NREMH+1
        IF(NREMH.GT.100) GOTO 40
        CALL LREMH(IFLRO,IFLR,K2,XT)
        AMK2=ULMASS(K2)
        AMIFLR=ULMASS(IFLR)
C...GIVE BALANCING PT TO IFLQ AND IFLQQ.
        CALL LPRIKT(PARL(14),PT,PHI)
        PT2=PT**2
        TM2K2=AMK2**2+PT2
        EK2=.5*(XT*W+TM2K2/XT/W)
        PZK2=-.5*(XT*W-TM2K2/XT/W)
        EPZ=W-TM2K2/XT/W
        WT=(1.-XT)*W*EPZ-PT2
C...CHECK IF ENERGY IN JET SYSTEM IS ENOUGH FOR FRAGMENTATION.
        IF(WT.LT.(AMIFL+AMIFLR+PARJ(32))**2) GOTO 30
        WT=SQRT(WT+PT2)
        TMIFLR=AMIFLR**2+PT2
        EIFL=.5*(WT+(AMIFL**2-TMIFLR)/WT)
        EIFLR=.5*(WT+(-AMIFL**2+TMIFLR)/WT)
        THER=ULANGL(-SQRT(EIFLR**2-TMIFLR),PT)
C...FORM JET SYSTEM.
        CALL LU1ENT(-MSTU(1),IFL,EIFL,0.,0.)
        CALL LU1ENT(MSTU(1)+1,IFLR,EIFLR,THER,PHI)
        CALL LUDBRB(MSTU(1),0,0.,0.,0.D0,0.D0,
     +  (DBLE(EPZ)-(1.D0-DBLE(XT))*DBLE(W))/
     +  (DBLE(EPZ)+(1.D0-DBLE(XT))*DBLE(W)))
        THEK2=ULANGL(PZK2,PT)
C...ADD FORMED "TARGET" PARTICLE.
        MSTU(10)=1
        P(MSTU(1)+2,5)=AMK2
        CALL LU1ENT(MSTU(1)+2,K2,EK2,THEK2,PHI+3.1415927)
        MSTU(10)=2
        K(MSTU(1)+1,3)=2
        K(MSTU(1)+2,3)=2
      ENDIF

      LST(21)=0
      RETURN

   40 LST(21)=1
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      FUNCTION LQMCUT(XP,ZP,AM1,AM2,AM3)

C...APPLY CUTS, IF NECESSARY, ON THE EVENT CONFIGURATION
C...OBTAINED FROM QCD MATRIX ELEMENTS.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      DATA S12,S23,S13/3*0./

      IF(LST(24).EQ.2) THEN
        S12=Q2*(1.-XP)/XP
        S23=Q2*(XP-X)*(1.-ZP)/X/XP+AM2**2+AM3**2
        S13=Q2*(XP-X)*ZP/X/XP+AM1**2+AM3**2
      ELSEIF(LST(24).EQ.3) THEN
        S13=Q2*(1.-XP)/XP
        S23=Q2*(XP-X)*(1.-ZP)/X/XP+AM2**2+AM3**2
        S12=Q2*(XP-X)*ZP/X/XP+AM1**2+AM2**2
        IF(S13.LT.(AM1+AM3)**2) GOTO 10
      ENDIF

      W=SQRT(W2)
      X1=1.-(S23-AM1**2)/W2
      X3=1.-(S12-AM3**2)/W2
      X2=2.-X1-X3
      PARI(21)=X1
      PARI(22)=X2
      PARI(23)=X3
      IF(X1.GT.1..OR.X2.GT.1..OR.X3.GT.1.) GOTO 10
      IF(X1*W/2..LT.AM1.OR.X2*W/2..LT.AM2.OR.X3*W/2..LT.AM3)  GOTO 10
      PA1=SQRT((0.5*X1*W)**2-AM1**2)
      PA2=SQRT((0.5*X2*W)**2-AM2**2)
      PA3=SQRT((0.5*X3*W)**2-AM3**2)
      IF(ABS((PA3**2-PA1**2-PA2**2)/(2.*PA1*PA2)).GE.1.) GOTO 10
      IF(ABS((PA2**2-PA1**2-PA3**2)/(2.*PA1*PA3)).GE.1.) GOTO 10
      LQMCUT=0
      RETURN

   10 LQMCUT=1
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LQQBEV

C...GENERATE BOSON-GLUON FUSION EVENT, CHOOSE XP AND ZP ACCORDING TO
C...QCD MATRIX ELEMENTS AND APPLY CUTS FOR SOFTNESS AND COLLINEARNESS.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)

      LST(24)=3
      W=SQRT(W2)
      J1=MSTU(1)
      J2=MSTU(1)+1
      J3=MSTU(1)+2
      J4=MSTU(1)+3

      CALL LXP(XP,IFAIL)
      IF(IFAIL.NE.0) GOTO 50

C...CHOOSE FLAVOUR OF PRODUCED QUARK-ANTIQUARK PAIR.
   10 CALL LFLAV(IFL1,IFL3)
      IF(LST(21).NE.0) RETURN
      CALL LZP(XP,ZP,IFAIL)
      IF(IFAIL.NE.0) GOTO 50
      IFL1A=IABS(IFL1)
      IFL3A=IABS(IFL3)
      AMIFL1=ULMASS(IFL1)
      AMIFL3=ULMASS(IFL3)

      IF(LST(14).EQ.0.OR.(LST(8).GE.2.AND.MOD(LST(8),10).NE.9)) THEN
C...IF BARYON PRODUCTION FROM TARGET REMNANT IS NEGLECTED THE
C...TARGET REMNANT IS APPROXIMATED BY A GLUON.
        IF(W.LT.AMIFL1+AMIFL3+PARJ(32)) GOTO 50
        IF(LQMCUT(XP,ZP,AMIFL1,0.,AMIFL3).NE.0) GOTO 50
        CALL LU3ENT(J1,IFL1,21,IFL3,W,PARI(21),PARI(23))
        K(MSTU(1)+1,3)=2
        CALL LUROBO(-ACOS(-P(J2,3)/SQRT(P(J2,3)**2+P(J2,1)**2)),
     +  0.,0.,0.,0.)
        GOTO 40
      ENDIF

      IF(W.LT.AMIFL1+AMIFL3+0.9+2.*PARJ(32)) GOTO 50
      IF(LQMCUT(XP,ZP,AMIFL1,1.,AMIFL3).NE.0) GOTO 50
      P(J1,5)=AMIFL1
      P(J3,5)=AMIFL3
C...CHOOSE TARGET VALENCE QUARK/DIQUARK TO FORM JET SYSTEM WITH
C...PRODUCED ANTIQUARK/QUARK.
      IFLR2=INT(1.+LST(22)/3.+RLU(0))
      IF(IFLR2.EQ.LST(22)) THEN
        IFLR1=2101
        IF(RLU(0).GT.PARL(4)) IFLR1=2103
      ELSE
        IFLR1=1000*IFLR2+100*IFLR2+3
      ENDIF
      IFLR2=3-IFLR2
      AMR1=ULMASS(IFLR1)
      AMR2=ULMASS(IFLR2)
      NREMH=0
   20 NREMH=NREMH+1
      IF(NREMH.GT.100) GOTO 50
      CALL LREMH(0,IFLR1,IFLR2,XT)
      CALL LPRIKT(PARL(14),PT,PHI)
      PT2=PT**2
      TM2R1=AMR1**2+PT2
      TM2R2=AMR2**2+PT2
      P(J2,5)=SQRT(TM2R1/(1.-XT)+TM2R2/XT)
      IF(LQMCUT(XP,ZP,AMIFL1,P(J2,5),AMIFL3).NE.0) GOTO 20
      MSTU(10)=1
      CALL LU3ENT(J1,IFL1,21,IFL3,W,PARI(21),PARI(23))
      MSTU(10)=2
      CALL LUROBO(-ACOS(-P(J2,3)/SQRT(P(J2,3)**2+P(J2,1)**2)),
     +0.,0.,0.,0.)
      EPZ=P(J2,4)-P(J2,3)
      IF(IFL1.GT.0) THEN
        IR1=J2
        IR2=J4
      ELSE
        IR1=J4
        IR2=J2
      ENDIF
      P(IR1,1)=PT*COS(PHI)
      P(IR1,2)=PT*SIN(PHI)
      P(IR1,3)=-0.5*((1.-XT)*EPZ-TM2R1/(1.-XT)/EPZ)
      P(IR1,4)= 0.5*((1.-XT)*EPZ+TM2R1/(1.-XT)/EPZ)
      P(IR1,5)=AMR1
      P(IR2,1)=-P(IR1,1)
      P(IR2,2)=-P(IR1,2)
      P(IR2,3)=-0.5*(XT*EPZ-TM2R2/XT/EPZ)
      P(IR2,4)= 0.5*(XT*EPZ+TM2R2/XT/EPZ)
      P(IR2,5)=AMR2
      K(IR1,1)=1
      K(IR1,2)=IFLR1
      K(IR2,1)=1
      K(IR2,2)=IFLR2
      K(J3,1)=2
      DO 30  I=J1,J4
        DO 30  J=3,5
   30 K(I,J)=0
      N=J4
      K(IR1,3)=2
      K(IR2,3)=2
      IF((P(J1,4)+P(J2,4))**2-(P(J1,1)+P(J2,1))**2-(P(J1,3)+P(J2,3))**2
     +  -P(J2,2)**2.LT.(P(J1,5)+P(J2,5)+PARJ(32))**2) GOTO 20
      IF((P(J3,4)+P(J4,4))**2-(P(J3,1)+P(J4,1))**2-(P(J3,3)+P(J4,3))**2
     +  -P(J4,2)**2.LT.(P(J3,5)+P(J4,5)+PARJ(32))**2) GOTO 20

   40 CONTINUE

      CALL LAZIMU(XP,ZP)
      LST(21)=0
      RETURN

   50 LST(21)=1
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LREMH(IFLRO,IFLR,K2,Z)

C...GIVES FLAVOUR AND ENERGY-MOMENTUM FRACTION Z FOR THE PARTICLE
C...TO BE PRODUCED OUT OF THE TARGET REMNANT WHEN THAT IS NOT A
C...SIMPLE DIQUARK.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)

C...FLAVOURS FIXED WHEN CALLING FROM PYREMM OR LQQBEV
      IF(IFLRO.EQ.0) GOTO 20

C...SPLIT TARGET REMNANT QQQQ -> QQQ + Q OR QQQQBAR -> QQBAR + QQ
C...Q (QBAR) IS THE PARTNER TO THE STRUCK SEA QUARK
C...QQQ ARE THE NUCLEON VALENCE QUARKS FROM WHICH A QUARK Q OR A
C...DIQUARK QQ IS CHOSEN AT RANDOM TO FORM A JET SYSTEM WITH THE
C...SCATTERED SEA ANTIQUARK OR QUARK, RESPECTIVELY, THE OTHER PARTON
C...FORMS A BARYON QQQ OR MESON QQBAR, RESPECTIVELY.
   10 IFLQ=INT(1.+LST(22)/3.+RLU(0))
      IF(IFLQ.EQ.LST(22)) THEN
        IFLQQ=2101
        IF(RLU(0).GT.PARL(4)) IFLQQ =2103
      ELSE
        IFLQQ=1000*IFLQ+100*IFLQ+3
      ENDIF
      IFLQ=3-IFLQ

C...CHOOSE FLAVOUR OF HADRON AND PARTON FOR JET SYSTEM
      IF(IFLRO.GT.0) THEN
        CALL LUKFDI(IFLQQ,IFLRO,IDUM,K2)
        IF(K2.EQ.0) GOTO 10
        IFLR=IFLQ
      ELSE
        CALL LUKFDI(IFLQ,IFLRO,IDUM,K2)
        IF(K2.EQ.0) GOTO 10
        IFLR=IFLQQ
      ENDIF

C...ENTRY FOR USE FROM PYSSPB, FLAVOURS GIVEN, CHOOSE E-P FRACTION
   20 KSP=IFLR
C...SPLIT ENERGY-MOMENTUM OF TARGET REMNANT ACCORDING TO FUNCTIONS P(Z)
C...Z=E-PZ FRACTION FOR QQ (Q) FORMING JET-SYSTEM WITH STRUCK Q (QBAR)
C...1-Z=E-PZ FRACTION FOR QQBAR (QQQ) HADRON
C...MQ=MASS OF (LIGHT) PARTON REMNANT Q (QQ) IN JET SYSTEM
C...MQ=MASS OF PRODUCED (HEAVY FLAVOUR) HADRON
      AMSP=ULMASS(KSP)
      AMK2=ULMASS(K2)
C...OLD LEPTO TREATMENT
C...P(Z)=(A+1)(1-Z)**A WITH <Z>=1/(A+2)=1/3 SINCE A=1 FIXED
      Z=1.-SQRT(RLU(0))
C...FLIP IF BARYON PRODUCED
      KC2=IABS(LUCOMP(K2))
      IF(KC2.GE.301.AND.KC2.LE.400) Z=1.-Z
      IF(LST(14).EQ.2) THEN
C...UPDATE OF LEPTO TREATMENT
C...P(Z)=(A+1)(1-Z)**A WITH <Z>=1/(A+2)=MQ/(MQ+MQ) --> A=A(MQ,MQ)
        A=(AMSP+AMK2)/AMSP - 2.
        Z=RLU(0)**(1./(A+1.))
      ELSEIF(LST(14).EQ.3) THEN
C...USING PETERSON FRAGMENTATION FUNCTION
C...P(Z)=N/(Z(1-1/Z-C/(1-Z))**2)  WHERE C=(MQ/MQ)**2  (FC=-C)
        FC=-(AMSP/AMK2)**2
   30   Z=RLU(0)
        IF(-4.*FC*Z*(1.-Z)**2.LT.RLU(0)*((1.-Z)**2-FC*Z)**2) GOTO 30
      ENDIF
      LST(27)=1
      K2A=IABS(K2)
      IF((K2A.GE.1.AND.K2A.LE.8).OR.K2A.EQ.21.OR.LUCOMP(K2A).EQ.90)
     +LST(27)=2

      RETURN
      END
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LSCALE(INFIN,QMAX)

C...GIVE MAXIMUM VIRTUALITY OF PARTONS IN PARTON SHOWERS.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)

      QMAX=0.1
      IF(LST(8).GE.2.AND.LST(8).LE.5) THEN
C...PARTON SHOWERS WITHOUT MATRIX ELEMENTS MATCHING
        IF(LST(9).EQ.1) THEN
          QMAX=Q2
        ELSEIF(LST(9).EQ.2) THEN
          QMAX=W2
        ELSEIF(LST(9).EQ.3) THEN
          QMAX=SQRT(W2*Q2)
        ELSEIF(LST(9).EQ.4) THEN
          QMAX=Q2*(1.-X)
        ELSEIF(LST(9).EQ.5) THEN
          QMAX=Q2*(1.-X)*MAX(1.,LOG(1./MAX(1.E-06,X)))
        ELSEIF(LST(9).EQ.9) THEN
          QMAX=W2**(2./3.)
        ELSE
          WRITE(6,*) ' WARNING, LSCALE: LST(9)=',LST(9),' NOT ALLOWED'
        ENDIF
      ELSEIF(LST(8).GT.10.AND.LST(24).EQ.1.AND.LST(8).NE.19) THEN
C...PARTON SHOWERS ADDED TO Q-EVENT FROM 1ST ORDER MATRIX ELEMENTS
C...SCALE GIVEN BY Y_CUT*W**2
        QMAX=PARL(27)*W2
      ELSEIF(LST(8).GT.10.AND.LST(8).NE.19) THEN
C...PARTON SHOWERS ADDED TO QG-/QQBAR-EVENT FROM 1ST ORDER MATRIX ELEMENTS
C...SCALE GIVEN BY INVARIANT MASS OF FINAL PARTON PAIR
        QMAX=P(27,5)**2
        IF(INFIN.LT.0) QMAX=MAX(ABS(-Q2-2.*FOUR(25,21)),
     &  ABS(-Q2-2.*FOUR(26,21)))
      ENDIF
      IF(INFIN.LT.0) THEN
        QMAX=SQRT(PYPAR(26)*QMAX)
      ELSE
        QMAX=SQRT(PYPAR(25)*QMAX)
      ENDIF

      RETURN
      END
*CMZ :  1.02/02 12/01/97  17.48.59  by  P. Zucchelli
*CMZ :  1.02/01 12/01/97  16.40.21  by  J. Brunner
*CMZ :  1.01/50 22/05/96  12.24.00  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  18.41.15  BY  PIERO ZUCCHELLI
*CMZ :  1.01/06 05/03/95  10.59.58  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 28/07/94  17.53.54  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C#######################################################################
C
C  THE FOLLOWING ROUTINES FOR PARTON CASCADES WERE MADE TOGETHER
C  WITH M. BENGTSSON AND T. SJOSTRAND (Z. PHYS. C37 (1988) 465,
C  NUCL. PHYS. B301 (1988) 554). CONTAIN MODIFICATIONS OF
C  ROUTINES IN PYTHIA 4.8 (SJOSTRAND, BENGTSSON, CPC 46 (1987) 43).
C
C **********************************************************************

      SUBROUTINE LSHOWR(ICALL)

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LBOOST/ DBETA(2,3),STHETA(2),SPHI(2),PB(5),PHIR
* ADDED BY ME
*      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON /PYPROC/ ISUB,KFL(3,2),XPY(2),SH,TH,UH,Q2PY,XSEC(0:40)
      COMMON /PYINT1/ XQPY(2,-6:6)
*      COMMON/PYINT1/MINT(400),VINT(400)
*KEEP,FOREFI.
C--
	INTEGER*4 IEVT
         COMMON/FOREFICASS/IEVT


*KEND.

      DOUBLE PRECISION DTHETA,DPHI,DBETA
      DOUBLE PRECISION DPQ2,DPB(3),DPA(3),DCTHET,DROBO(5)
      DIMENSION KS(9,5),PS(9,5),ROBO(5),XPQ(-6:6)
      SAVE KS,PS
      DATA MG/2/
      IF(ICALL.EQ.0) THEN
C...INITIALIZE CASCADE FOR EACH EVENT, SAVE EVENT RECORD IN OVERALL CMS.
        DO 10 I=1,9
          DO 10 J=1,5
            KS(I,J)=0
   10   PS(I,J)=0.
        DO 20 J=1,5
          KS(1,J)=K(1,J)
          PS(1,J)=P(1,J)
          KS(2,J)=K(2,J)
          PS(2,J)=P(2,J)
          KS(5,J)=K(3,J)
          PS(5,J)=P(3,J)
          KS(7,J)=K(4,J)
   20   PS(7,J)=P(4,J)
        KS(5,3)=3
        KS(7,1)=21
        KS(7,3)=5
*        WRITE(*,*)'LUJET SAVING CALLED'
*        CALL LULIST(1)
        RETURN
      ENDIF

C     CALL GULIST(1,2)
C...APPLY PARTON CASCADE ON QPM EVENT.
C...SAVE INCOMING AND OUTGOING QUARK AS WELL AS SCATTERED LEPTON.
      KS(6,1)=21
      KS(6,2)=LST(25)
      KS(6,3)=4
      KS(8,1)=21
      KS(8,2)=K(5,2)
      KS(8,3)=6
      KS(9,1)=0
      KS(9,2)=K(4,2)
      KS(9,3)=5
      DO 30  J=1,5
        PS(6,J)=0.
        PS(8,J)=P(5,J)
   30 PS(9,J)=P(4,J)
      XR=X
      DPQ2=DBLE(Q2)
      PMA1=0.
      PS(6,5)=PMA1
      PMA2=PS(8,5)
      DPB(1)=0.5D0*(DPQ2*(1D0/XR-1D0)+DBLE(PS(1,5))**2-
     +ULMASS(IABS(KS(7,2)))**2)/(PS(1,4)+PS(2,4))
      DPB(2)=DSQRT(DPB(1)**2+DPQ2)
      DCTHET=(DBLE(PS(2,4))*DPB(1)-DPQ2/(2D0*XR))/(DBLE(PS(2,3))*
     +DPB(2))
      DPA(1)=(DPB(2)*DCTHET)**2-DPB(1)**2
      DPA(2)=DPQ2-DBLE(PMA1)**2+DBLE(PMA2)**2
      PS(6,4)=-(DPA(2)*DPB(1)-DPB(2)*DCTHET*DSQRT(DPA(2)**2+4D0*
     +DBLE(PMA1)**2*DPA(1)))/(2D0*DPA(1))
      PS(6,3)=-SQRT((PS(6,4)+PMA1)*(PS(6,4)-PMA1))
C...PARTONS WITH COLOUR INFORMATION IN HADRONIC CMS FRAME.
      DO 40  I=10,26
        DO 40  J=1,5
          K(I,J)=0
          P(I,J)=0.
   40 V(I,J)=0.
      NS=20
      K(NS+1,1)=21
      K(NS+1,2)=K(3,2)
      K(NS+1,3)=3
      K(NS+2,1)=-1
      K(NS+2,3)=NS+1
      K(NS+3,2)=KS(6,2)
      DO 50  J=1,5
   50 P(NS+1,J)=P(3,J)
      K(NS+3,1)=13
      K(NS+3,3)=2
      IF(MG.EQ.1) THEN
        PW2=W2
        DPA(3)=DSQRT(DPA(2)**2+4D0*DPQ2*DBLE(PMA1)**2)
        DPB(1)=(1D0/DBLE(XR)-2D0)*DPQ2/(2D0*SQRT(PW2))
        DPB(2)=DSQRT(DPB(1)**2+DPQ2)
        P(NS+3,4)=(DPB(2)*DPA(3)-DPB(1)*DPA(2))/(2D0*DPQ2)
        P(NS+3,3)=-P(NS+3,4)
      ENDIF
      P(NS+3,5)=0.
      K(NS+4,1)=-1
      K(NS+4,3)=NS+3
      K(NS+3,4)=NS+5
      K(NS+3,5)=NS+5
      P(NS+4,3)=NS+5
      P(NS+4,4)=NS+5
      K(NS+5,1)=3
      K(NS+5,3)=8
      K(NS+5,2)=KS(8,2)
      K(NS+6,1)=-1
      K(NS+6,3)=NS+5
      DO 60  J=1,4
        IF(MG.EQ.1) THEN
          P(NS+5,J)=P(NS+1,J)+P(NS+3,J)
        ELSE
          P(NS+5,J)=P(5,J)
          P(NS+3,J)=P(NS+5,J)-P(NS+1,J)
        ENDIF
   60 CONTINUE
      P(NS+5,5)=PMA2
      P(NS+6,1)=NS+3
      P(NS+6,2)=NS+3
      K(NS+5,4)=(NS+3)*MSTU(5)
      K(NS+5,5)=(NS+3)*MSTU(5)
      N=NS+6
C     CALL GULIST(2,2)
C...COPY SAVED RECORD IN OVERALL CMS TO LINE 1 THROUGH 9.
C...LINES 1,2,5,6,7 IN EP CMS, 8,9 IN HADRONIC CMS
      DO 70  I=1,9
        DO 70  J=1,5
          K(I,J)=KS(I,J)
   70 P(I,J)=PS(I,J)
*      WRITE(*,*)'1,2,5,6,7 EP CMS, 8,9, HAD CMS'
*      CALL LULIST(1)

C     CALL GULIST(3,2)
C...SCALE FOR BREMSSTRAHLUNG ETC.
      Q2PY=Q2
      IPY(40)=8
      IPY(47)=N
C...SAVE QUANTITIES FOR LATER USE.
      XPY(1)=1.
      XPY(2)=XR
      CALL PYSTFU(K(2,2),XR,Q2,XPQ)
      DO 80  IFL=-6,6
   80 XQPY(2,IFL)=XPQ(IFL)
      IF(LST(23).EQ.1) THEN
        ISUB=39
        IPY(11)=1
      ELSEIF(LST(23).EQ.3) THEN
        ISUB=39
        IPY(11)=2
      ELSEIF(LST(23).EQ.4) THEN
        ISUB=39
        IPY(11)=3
      ELSEIF(LST(23).EQ.2) THEN
        ISUB=40
      ENDIF
      IF(ISUB.EQ.39.AND.IPY(11).EQ.1) THEN
        KFL(2,1)=22
      ELSEIF(ISUB.EQ.39.AND.IPY(11).EQ.2) THEN
        KFL(2,1)=23
      ELSEIF(ISUB.EQ.39.AND.IPY(11).EQ.3) THEN
        KFL(2,1)=23
      ELSEIF(ISUB.EQ.40) THEN
        KFL(2,1)=-24
      ENDIF
      IFL1=K(6,2)
      IFL2=K(8,2)
      KFL(2,2)=IFL1
      KFL(1,1)=KFL(2,1)
      KFL(1,2)=KFL(2,2)
      IF(ISUB.EQ.39) KFL(3,1)=K(1,2)
      IF(ISUB.EQ.40) KFL(3,1)=K(1,2)+ISIGN(1,K(1,2))
      KFL(3,2)=IFL2
      PYVAR(2)=(P(1,4)+P(2,4))**2
      PYVAR(1)=SQRT(PYVAR(2))
      PYVAR(3)=P(1,5)
      PYVAR(4)=P(2,5)
      PYVAR(5)=P(1,3)
      IPY(41)=K(1,2)
      IPY(42)=K(2,2)
      IPY(48)=0

C...GENERATE TIMELIKE PARTON SHOWER (IF REQUIRED)
      IF(IPY(13).EQ.1) THEN
        CALL LSCALE(1,QMAX)
        QMAX=MIN(QMAX,P(25,4))
*        WRITE(*,*)' GENERATE TIME-LIKE SHOWER WITH QMAX=',QMAX
*        MSTP(22)=3
*        WRITE(*,*)'PRE LUSHOW'
*        CALL LULIST(1)
        CALL LUSHOW(25,0,QMAX)

*        WRITE(*,*)' AFTER LUSHOW'
*        CALL LULIST(1)
      ENDIF
      IT=25
      IF(N.GE.27) IT=27
      NS=N

*      WRITE(*,*)'TEST BOZZO:NS,N=',NS,N
C     CALL GULIST(4,2)

C...GENERATE SPACELIKE PARTON SHOWER (IF REQUIRED)
      IPU1=0
      IPU2=23
      IF(XPY(2)*(1.+(P(IT,5)**2+PYPAR(22))/P(21,5)**2).GT.0.999) THEN
        WRITE(*,*)'21-47 ERROR'
        LST(21)=47
        RETURN
      ENDIF
      IF (NS.EQ.26) WRITE(*,*)'MAYBE BOZZO...',IEVT+1
      IF(IPY(14).GE.1) THEN
        WRITE(*,*)'SPACE-LIKE SHOWER?',IEVT+1
        CALL PYSSPB(IPU1,IPU2)
      ELSE
        DO 90  I=NS+1,NS+4
          DO 90  J=1,5
            K(I,J)=0
            P(I,J)=0.
   90   V(I,J)=0.
        K(NS+1,1)=11
        K(NS+1,2)=KFL(2,1)
        K(NS+1,3)=21
        DO 100 J=1,5
  100   P(NS+1,J)=P(21,J)
        K(NS+2,1)=-1
        K(NS+2,3)=NS+1
        K(NS+3,1)=13
        K(NS+3,2)=KFL(2,2)
        K(NS+3,3)=23
        K(NS+3,4)=23
        K(NS+3,5)=23
        P(NS+3,3)=(P(IT,5)**2+Q2)*(P(21,4)-P(21,3))/(2.*Q2)
        P(NS+3,4)=-P(NS+3,3)
        K(NS+4,1)=-1
        K(NS+4,3)=NS+3
        P(NS+4,3)=23
        P(NS+4,4)=23
        P(24,1)=NS+3
        P(24,2)=NS+3
        K(23,4)=K(23,4)+(NS+3)*MSTU(5)
        K(23,5)=K(23,5)+(NS+3)*MSTU(5)
        IPU1=0
        IPU2=NS+3
        N=N+4
      ENDIF

*       CALL LULIST(1)

C     CALL GULIST(5,2)

C...ROTATE AND BOOST OUTGOING PARTON SHOWER
      IF(N.GT.30) THEN
*        WRITE(*,*)'ROTATE AND BOOST?'
        K(N+1,1)=0
        DO 110 J=1,4
  110   P(N+1,J)=P(NS+1,J)+P(NS+3,J)
        IF(P(N+1,4).LE.1.01*P(IT,5)) THEN
          LST(21)=50
          RETURN
        ENDIF
        ROBO(1)=ULANGL(P(IT,3),SQRT(P(IT,1)**2+P(IT,2)**2))
        ROBO(2)=ULANGL(P(IT,1),P(IT,2))
        CALL LUDBRB(25,NS,0.,-ROBO(2),0.D0,0.D0,0.D0)
        CALL LUDBRB(25,NS,-ROBO(1),0.,0.D0,0.D0,0.D0)
        DROBO(5)=-(P(IT,3)*P(IT,4)-P(N+1,4)*SQRT(P(N+1,4)**2-
     +  P(IT,4)**2+P(IT,3)**2))/(P(IT,3)**2+P(N+1,4)**2)
        CALL LUDBRB(25,NS,0.,0.,0.D0,0.D0,DROBO(5))
        ROBO(1)=ULANGL(P(N+1,3),SQRT(P(N+1,1)**2+P(N+1,2)**2))
        ROBO(2)=ULANGL(P(N+1,1),P(N+1,2))
        CALL LUDBRB(25,NS,ROBO(1),ROBO(2),0.D0,0.D0,0.D0)
      ENDIF
C     CALL GULIST(6,2)

      Q2PY=Q2
C...HADRON REMNANT AND PRIMORDIAL KT
      IPY(47)=N
*      WRITE(*,*)' THEN CALL PYREMM IPU1,IPU2',IPU1,IPU2
      CALL PYREMM(IPU1,IPU2)
      IF(IPY(48).EQ.1) THEN
*        WRITE(*,*)'IPY(48) ERROR'
        LST(21)=48
        RETURN
      ENDIF
C     CALL GULIST(7,2)

C...TRANSFORM LINE 1,2 AND 5-7 TO HADRONIC CMS FRAME.
*      WRITE(*,*)'1-2,5-7 BOOST'
      CALL LUDBRB(1,2,0.,0.,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3))
      CALL LUDBRB(1,2,-STHETA(2),0.,0.D0,0.D0,0.D0)
      CALL LUDBRB(5,7,0.,0.,-DBETA(2,1),-DBETA(2,2),-DBETA(2,3))
      CALL LUDBRB(5,7,-STHETA(2),0.,0.D0,0.D0,0.D0)
C     CALL GULIST(8,2)
*      WRITE(*,*)'1-2,5-7 TO HADRON CMS'
*      CALL LULIST(1)

C...REARRANGE PARTONS ALONG STRINGS
      MSTU(24)=0
      CALL LUPREP(0)
*      WRITE(*,*)'CALL LUPREP'
*      CALL LULIST(1)

      IF(MSTU(24).NE.0) THEN
C       CALL GULIST(88,2)
        IF(LST(3).GE.1) WRITE(6,*) ' LUPREP ERROR MSTU(24)= ',MSTU(24)
        LST(21)=49
        RETURN
      ENDIF
C     CALL GULIST(9,2)

C...CLEAN UP EVENT RECORD -> ORDER:
C...1=INC. LEPTON; 2=INC. NUCLEON; 3=EXCH BOSON; 4=SCAT. LEPTON;
C...5=INC. PARTON BEFORE INITIAL SHOWER; 6=INC. QUARK AT BOSON VERTEX
C...AFTER SHOWER; 7=SCAT. QUARK AT BOSON VERTEX BEFORE FINAL SHOWER
      LST(26)=7
      DO 120 J=1,5
        K(N+1,J)=K(4,J)
  120 P(N+1,J)=P(4,J)
      DO 130 J=1,5
        K(3,J)=K(5,J)
        P(3,J)=P(5,J)
        K(4,J)=K(9,J)
        P(4,J)=P(9,J)
        K(5,J)=K(N+1,J)
        P(5,J)=P(N+1,J)
C     K(7,J)=K(8,J)
C     P(7,J)=P(8,J)
        K(6,J)=K(NS+3,J)
        P(6,J)=P(NS+3,J)
        K(7,J)=K(IT,J)
        P(7,J)=P(IT,J)
  130 CONTINUE
      K(3,3)=1
      K(4,3)=1
      K(6,1)=21
      K(6,3)=5
      K(6,4)=0
      K(6,5)=0
      K(7,1)=21
      K(7,3)=6
      K(7,4)=0
      K(7,5)=0
C...ACTIVATE LINE WITH SCATTERED LEPTON.
      K(4,1)=1
C...DEACTIVATE OBSOLETE LINES 8, 9 AND 21, NS+1 (EXTRA LINES WITH BOSON)
      K(8,1)=0
      K(9,1)=0
      K(21,1)=0
      IF(K(NS+1,2).EQ.K(3,2)) K(NS+1,1)=0
C...ZERO IRRELEVANT LINES WITH K(I,1)<0
      DO 150 I=1,N
        IF(K(I,1).LT.0) THEN
          DO 140 J=1,5
            K(I,J)=0
  140     P(I,J)=0.
        ENDIF
  150 CONTINUE

*      WRITE(*,*)'AFTER CLEANUP'
*      CALL LULIST(1)

C     CALL GULIST(10,2)
C...DELETE INTERNAL PARTON LINES, I.E. WITH K(I,1)=13,14
      IF(MOD(LST(4)/10,10).EQ.0) THEN
        CALL LTIMEX(T1)
        CALL LUEDIT(14)
*      WRITE(*,*)'AFTER LUEDIT=14'
*      CALL LULIST(1)
        CALL LTIMEX(T2)
C       CALL GULIST(11,2)
      ENDIF
C...DELETE EMPTY LINES
      CALL LTIMEX(T1)
      CALL LUEDIT(12)
*      WRITE(*,*)'AFTER LUEDIT=12'
*      CALL LULIST(1)
      CALL LTIMEX(T2)
C     CALL GULIST(12,2)

      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 20/09/94  14.43.49  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LSIGMX(NPAR,DERIV,DIFSIG,XF,IFLAG)

C...CALCULATES THE NEGATIVE OF THE DIFFERENTIAL CROSS-SECTION.
C...IN THE GENERATION PROCEDURE THE MAXIMUM OF THE DIFFERENTIAL CROSS-
C...SECTION IS NEEDED FOR WEIGHTING PURPOSES. THIS MAXIMUM IS FOUND BY
C...MINIMIZING THE NEGATIVE DIFFERENTIAL CROSS-SECTION USING THE MINUIT
C...ROUTINES WHICH ARE THEN CALLING THIS ROUTINE.
C...MORE PRECISLY, ONLY THE PART OF THE CROSS-SECTION FORMULA WHICH IS
C...NEEDED FOR THE WEIGHTING PROCEDURE IS INCLUDED HERE.
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      DIMENSION DERIV(30),XF(30)
      COMMON/LINPATCH/NCALLS,NCALL

      DUMMY=NPAR+DERIV(1)
      IF(IFLAG.EQ.1) NCALLS=0
      IF(IFLAG.EQ.2) WRITE(6,10000)

      DIFSIG=1.E+10
      NCALLS=NCALLS+1
      X=XF(1)
      IF(X.LT.XMIN.OR.X.GT.XMAX) RETURN
      S=PARL(21)
      PM2=PSAVE(3,2,5)**2
      IF(LST(31).EQ.1) THEN
        Q2=XF(2)
        Y=Q2/(PARL(21)*X)
        W2=(1.-X)*Y*PARL(21)+PSAVE(3,2,5)**2
      ELSEIF(LST(31).EQ.2) THEN
        Y=XF(2)
        Q2=Y*X*PARL(21)
        W2=(1.-X)*Y*PARL(21)+PSAVE(3,2,5)**2
      ELSEIF(LST(31).EQ.3) THEN
        W2=XF(2)
        Y=(W2-PSAVE(3,2,5)**2)/((1.-X)*PARL(21))
        Q2=X*Y*PARL(21)
      ENDIF
      Q2LOW=MAX(Q2MIN,X*YMIN*S,(W2MIN-PM2)*X/(1.-X))
      Q2UPP=MIN(Q2MAX,X*YMAX*S,(W2MAX-PM2)*X/(1.-X))
      YLOW=MAX(YMIN,Q2MIN/(S*X),(W2MIN-PM2)/(S*(1.-X)))
      YUPP=MIN(YMAX,Q2MAX/(S*X),(W2MAX-PM2)/(S*(1.-X)))
      W2LOW=MAX(W2MIN,(1.-X)*YMIN*S+PM2,Q2MIN*(1.-X)/X+PM2)
      W2UPP=MIN(W2MAX,(1.-X)*YMAX*S+PM2,Q2MAX*(1.-X)/X+PM2)
      IF(Q2.LT.Q2LOW.OR.Q2.GT.Q2UPP) RETURN
      IF(Y.LT.YLOW.OR.Y.GT.YUPP) RETURN
      IF(W2.LT.W2LOW.OR.W2.GT.W2UPP) RETURN
      LST2=LST(2)
      LST(2)=-1
      CALL LEPTO
      LST(2)=LST2
      IF(LST(21).NE.0) RETURN
      DIFSIG=-PQ(17)*COMFAC

      IF(LST(3).GE.4.AND.IFLAG.EQ.3) WRITE(6,10100) NCALLS,DIFSIG,X,Y,
     +Q2,W2
      RETURN

10000 FORMAT(' WARNING: IFLAG = 2 IN CALL TO LSIGMX, WHICH DOES NOT '
     &,'CALCULATE DERIVATIVES.')
10100 FORMAT(/,5X,'TERMINATING ENTRY IN LSIGMX AFTER ',I5,' CALLS.',/,
     &5X,'BEST ESTIMATE OF MINIMUM FOUND TO BE ',E12.4,/,
     &5X,'LOCATED AT X, Y, Q**2, W**2 = ',4G10.3,/)

      END
*CMZ :  1.00/00 04/07/94  15.02.26  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C   14/06/92 206150131  MEMBER NAME  LEPTO61  (LUND)     M  FVS
C######################################################################C
C                                                                      C
C  THE LUND MONTE CARLO FOR DEEP INELASTIC LEPTON-NUCLEON SCATTERING   C
C                                                                      C
C                          L E P T O                                   C
C                                                                      C
C                   VERSION 6.1, MAY 4, 1992                           C
C                                                                      C
C   AUTHOR:   GUNNAR INGELMAN          PART-TIME ALSO AT               C
C             DESY THEORY GROUP        DEPT OF RADIATION SCIENCES      C
C             (ROOM 202  BLDG 2A)      UPPSALA UNIVERSITY              C
C             NOTKESTRASSE 85          BOX 535                         C
C             D-2000 HAMBURG 52, FRG   S-751 21 UPPSALA, SWEDEN        C
C             PHONE: +49(40)8998-2795  +46(18)18-3884                  C
C             TELEFAX:          -2777           -3833                  C
C             E-MAIL:   USE  INGELMAN@DESYVAX  FORWARD SET TO          C
C                     T00ING@DHHDESY3  INGELMAN@TSL.UU.SE              C
C                                                                      C
C   CONTRIBUTIONS ON PARTON CASCADES: M. BENGTSSON, T. SJOSTRAND       C
C                                                                      C
C   AVAILABILITY: ON REQUEST OR FROM DESY IBM AND VAX/VMS SYSTEMS:     C
C   DESY IBM LIBRARY      VXDESY DIRECTORY         CONTENT             C
C   T00ING.LUND(MEMBER)   DISK$T__:[INGELMAN.LUND]                     C
C              LEPTOINF            LEPTO.INFO      INFO, NEWS, UPDATES C
C              LEPTOTEX            LEPTO.TEX       MANUAL IN LATEX     C
C              LEPTO61             LEPTO61.FOR     SOURCE CODE         C
C              LEPTODEM            LEPTODEM.FOR    DEMO PROGRAM        C
C                                  LEPTODEM.COM    DEMO COMMAND FILE   C
C   T00ING.OBJECT(LEPTO61)         LEPTO61.OBJ     OBJECT CODE         C
C                                                                      C
C   MANUAL: G. INGELMAN, UPPSALA PREPRINT TSL/ISV 92-0065 AND          C
C   IN PROC. `PHYSICS AT HERA', EDS. W. BUCHMUELLER, G. INGELMAN,      C
C   DESY HAMBURG 1992, VOL. 3, P. 1366                                 C
C                                                                      C
C   PLEASE REPORT ANY PROBLEMS OR SUGGESTIONS FOR IMPROMEVENTS.        C
C                                                                      C
C######################################################################C

      SUBROUTINE LTIMEX(TIME)
C...INTERFACE ROUTINE TO TRANSFER A CALL TO SOME MACHINE-DEPENDENT
C...ROUTINE TO GET THE EXECUTION TIME USED SINCE JOB STARTED.
C...NICE, BUT NOT NECESSARY INFORMATION. CAN ALSO BE CALLED BY USER.

      TIME=0.
C...USE OF CERN LIBRARY ROUTINE Z007, REPLACE/DELETE IF NOT AVAILABLE.
      CALL TIMEX(TIME)
      RETURN
      END
*CMZ :  1.00/00 09/08/94  17.43.59  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION LUNPIK(ID,ISGN)
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      IDENT=ID*ISGN
      IF      (IDENT.EQ. 1) THEN
        IPKDEF=-211
      ELSEIF  (IDENT.EQ.-1) THEN
        IPKDEF= 211
      ELSEIF  (IDENT.EQ. 2) THEN
        IPKDEF=111
      ELSEIF  (IDENT.EQ.-2) THEN
        IPKDEF=111
      ELSEIF  (IDENT.EQ. 3) THEN
        IPKDEF=-321
      ELSEIF  (IDENT.EQ.-3) THEN
        IPKDEF= 321
      ELSEIF  (IDENT.EQ. 4) THEN
C
C K0 --> K0_LONG (IS 130) / K0_SHORT (IS 310) = 1/1
        CALL RANMAR(XIO,1)
        IF (XIO.GT.BRK0) THEN
          IPKDEF= 130
        ELSE
          IPKDEF= 310
        ENDIF
      ELSEIF  (IDENT.EQ.-4) THEN
C
C K0B--> K0_LONG (IS 130) / K0_SHORT (IS 310) = 1/1
        CALL RANMAR(XIO,1)
        IF (XIO.GT.BRK0B) THEN
          IPKDEF= 130
        ELSE
          IPKDEF= 310
        ENDIF
      ELSEIF  (IDENT.EQ. 8) THEN
        IPKDEF= 22
      ELSEIF  (IDENT.EQ.-8) THEN
        IPKDEF= 22
      ELSEIF  (IDENT.EQ. 9) THEN
        IPKDEF= 221
      ELSEIF  (IDENT.EQ.-9) THEN
        IPKDEF= 221
      ELSE
        PRINT *, 'STOP IN IPKDEF, WRONG IDENT=',IDENT
        STOP
      ENDIF
      LUNPIK=IPKDEF
      END
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.35.13  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 19/07/94  17.12.54  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C*********************************************************************

      SUBROUTINE LUROBO(THE,PHI,BEX,BEY,BEZ)

C...PURPOSE: TO PERFORM ROTATIONS AND BOOSTS.
      IMPLICIT DOUBLE PRECISION(D)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
*KEEP,POLAR.
C--
      COMMON /POLARIZ/POL(4000,3)
      REAL POLARX(4)
*KEND.
      SAVE /LUJETS/,/LUDAT1/
      DIMENSION ROT(3,3),PR(3),VR(3),DP(4),DV(4),POR(3)

C...FIND RANGE OF ROTATION/BOOST. CONVERT BOOST TO DOUBLE PRECISION.
      IMIN=1
      IF(MSTU(1).GT.0) IMIN=MSTU(1)
      IMAX=N
      IF(MSTU(2).GT.0) IMAX=MSTU(2)
      DBX=BEX
      DBY=BEY
      DBZ=BEZ
      GOTO 30

C...ENTRY FOR SPECIFIC RANGE AND DOUBLE PRECISION BOOST.
      ENTRY LUDBRB(IMI,IMA,THE,PHI,DBEX,DBEY,DBEZ)
      IMIN=IMI
      IF(IMIN.LE.0) IMIN=1
      IMAX=IMA
      IF(IMAX.LE.0) IMAX=N
      DBX=DBEX
      DBY=DBEY
      DBZ=DBEZ

C...OPTIONAL RESETTING OF V (WHEN NOT SET BEFORE.)
      IF(MSTU(33).NE.0) THEN
        DO 20  I=MIN(IMIN,MSTU(4)),MIN(IMAX,MSTU(4))
          DO 10  J=1,5
            V(I,J)=0.
   10     CONTINUE
   20   CONTINUE
        MSTU(33)=0
      ENDIF

C...CHECK RANGE OF ROTATION/BOOST.
   30 IF(IMIN.GT.MSTU(4).OR.IMAX.GT.MSTU(4)) THEN
        CALL LUERRM(11,'(LUROBO:) RANGE OUTSIDE LUJETS MEMORY')
        RETURN
      ENDIF

C...ROTATE, TYPICALLY FROM Z AXIS TO DIRECTION (THETA,PHI).
      IF(THE**2+PHI**2.GT.1E-20) THEN
        ROT(1,1)=COS(THE)*COS(PHI)
        ROT(1,2)=-SIN(PHI)
        ROT(1,3)=SIN(THE)*COS(PHI)
        ROT(2,1)=COS(THE)*SIN(PHI)
        ROT(2,2)=COS(PHI)
        ROT(2,3)=SIN(THE)*SIN(PHI)
        ROT(3,1)=-SIN(THE)
        ROT(3,2)=0.
        ROT(3,3)=COS(THE)
        DO 60  I=IMIN,IMAX
          IF(K(I,1).LE.0) GOTO 60
          DO 40  J=1,3
            PR(J)=P(I,J)
            VR(J)=V(I,J)
            POR(J)=POL(I,J)
   40     CONTINUE
          DO 50  J=1,3
            POL(I,J)=ROT(J,1)*POR(1)+ROT(J,2)*POR(2)+ROT(J,3)*POR(3)
            P(I,J) =ROT(J,1)*PR(1)+ROT(J,2)*PR(2)+ROT(J,3)*PR(3)
            V(I,J) =ROT(J,1)*VR(1)+ROT(J,2)*VR(2)+ROT(J,3)*VR(3)
   50     CONTINUE
   60   CONTINUE




      ENDIF

C...BOOST, TYPICALLY FROM REST TO MOMENTUM/ENERGY=BETA.
      IF(DBX**2+DBY**2+DBZ**2.GT.1E-20) THEN
        DB=SQRT(DBX**2+DBY**2+DBZ**2)
        IF(DB.GT.0.99999999D0) THEN
C...RESCALE BOOST VECTOR IF TOO CLOSE TO UNITY.
          CALL LUERRM(3,'(LUROBO:) BOOST VECTOR TOO LARGE')
          DBX=DBX*(0.99999999D0/DB)
          DBY=DBY*(0.99999999D0/DB)
          DBZ=DBZ*(0.99999999D0/DB)
          DB=0.99999999D0
        ENDIF
        DGA=1D0/SQRT(1D0-DB**2)
        DO 80  I=IMIN,IMAX
          IF(K(I,1).LE.0) GOTO 80
          DO 70  J=1,4
            DP(J)=P(I,J)
            DV(J)=V(I,J)
   70     CONTINUE
          DBP=DBX*DP(1)+DBY*DP(2)+DBZ*DP(3)
          DGABP=DGA*(DGA*DBP/(1D0+DGA)+DP(4))
          P(I,1)=DP(1)+DGABP*DBX
          P(I,2)=DP(2)+DGABP*DBY
          P(I,3)=DP(3)+DGABP*DBZ
          P(I,4)=DGA*(DP(4)+DBP)
          DBV=DBX*DV(1)+DBY*DV(2)+DBZ*DV(3)
          DGABV=DGA*(DGA*DBV/(1D0+DGA)+DV(4))
          V(I,1)=DV(1)+DGABV*DBX
          V(I,2)=DV(2)+DGABV*DBY
          V(I,3)=DV(3)+DGABV*DBZ
          V(I,4)=DGA*(DV(4)+DBV)
   80   CONTINUE
      ENDIF

      RETURN
      END
*CMZ :  1.02/09 14/01/97  15.14.44  by  P. Zucchelli
*CMZ :  1.01/51 24/05/96  11.26.15  by  Piero Zucchelli
*-- Author :
C*********************************************************************

      SUBROUTINE LUSTRF(IP)
C...Purpose: to handle the fragmentation of an arbitrary colour singlet
C...jet system according to the Lund string fragmentation model.
      IMPLICIT DOUBLE PRECISION(D)
      PARAMETER (MAXPZTRY=1000, MAXPZTRYR=10)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/
      DIMENSION DPS(5),KFL(3),PMQ(3),PX(3),PY(3),GAM(3),IE(2),PR(2),
     +IN(9),DHM(4),DHG(4),DP(5,5),IRANK(2),MJU(4),IJU(3),PJU(5,5),
     +TJU(5),KFJH(2),NJS(2),KFJS(2),PJS(4,5),MSTU9T(8),PARU9T(8)

C...Function: four-product of two vectors.
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)
      DFOUR(I,J)=DP(I,4)*DP(J,4)-DP(I,1)*DP(J,1)-DP(I,2)*DP(J,2)-
     +DP(I,3)*DP(J,3)

C...Reset counters. Identify parton system.
      MSTJ(91)=0
      NSAV=N
      MSTU90=MSTU(90)
      NP=0
      KQSUM=0
      DO 100 J=1,5
        DPS(J)=0D0
  100 CONTINUE
      MJU(1)=0
      MJU(2)=0
      I=IP-1
  110 I=I+1
      IF(I.GT.MIN(N,MSTU(4)-MSTU(32))) THEN
        CALL LUERRM(12,'(LUSTRF:) failed to reconstruct jet system')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      IF(K(I,1).NE.1.AND.K(I,1).NE.2.AND.K(I,1).NE.41) GOTO 110
      KC=LUCOMP(K(I,2))
      IF(KC.EQ.0) GOTO 110
      KQ=KCHG(KC,2)*ISIGN(1,K(I,2))
      IF(KQ.EQ.0) GOTO 110
      IF(N+5*NP+11.GT.MSTU(4)-MSTU(32)-5) THEN
        CALL LUERRM(11,'(LUSTRF:) no more memory left in LUJETS')
        IF(MSTU(21).GE.1) RETURN
      ENDIF

C...Take copy of partons to be considered. Check flavour sum.
      NP=NP+1
      DO 120 J=1,5
        K(N+NP,J)=K(I,J)
        P(N+NP,J)=P(I,J)
        IF(J.NE.4) DPS(J)=DPS(J)+P(I,J)
  120 CONTINUE
      DPS(4)=DPS(4)+SQRT(DBLE(P(I,1))**2+DBLE(P(I,2))**2+
     +DBLE(P(I,3))**2+DBLE(P(I,5))**2)
      K(N+NP,3)=I
      IF(KQ.NE.2) KQSUM=KQSUM+KQ
      IF(K(I,1).EQ.41) THEN
        KQSUM=KQSUM+2*KQ
        IF(KQSUM.EQ.KQ) MJU(1)=N+NP
        IF(KQSUM.NE.KQ) MJU(2)=N+NP
      ENDIF
      IF(K(I,1).EQ.2.OR.K(I,1).EQ.41) GOTO 110
      IF(KQSUM.NE.0) THEN
        CALL LUERRM(12,'(LUSTRF:) unphysical flavour combination')
        IF(MSTU(21).GE.1) RETURN
      ENDIF

C...Boost copied system to CM frame (for better numerical precision).
      IF(ABS(DPS(3)).LT.0.99D0*DPS(4)) THEN
        MBST=0
        MSTU(33)=1
        CALL LUDBRB(N+1,N+NP,0.,0.,-DPS(1)/DPS(4),-DPS(2)/DPS(4),
     +  -DPS(3)/DPS(4))
      ELSE
        MBST=1
        HHBZ=SQRT(MAX(1D-6,DPS(4)+DPS(3))/MAX(1D-6,DPS(4)-DPS(3)))
        DO 130 I=N+1,N+NP
          HHPMT=P(I,1)**2+P(I,2)**2+P(I,5)**2
          IF(P(I,3).GT.0.) THEN
            HHPEZ=(P(I,4)+P(I,3))/HHBZ
            P(I,3)=0.5*(HHPEZ-HHPMT/HHPEZ)
            P(I,4)=0.5*(HHPEZ+HHPMT/HHPEZ)
          ELSE
            HHPEZ=(P(I,4)-P(I,3))*HHBZ
            P(I,3)=-0.5*(HHPEZ-HHPMT/HHPEZ)
            P(I,4)=0.5*(HHPEZ+HHPMT/HHPEZ)
          ENDIF
  130   CONTINUE
      ENDIF

C...Search for very nearby partons that may be recombined.
      NTRYR=0
      PARU12=PARU(12)
      PARU13=PARU(13)
      MJU(3)=MJU(1)
      MJU(4)=MJU(2)
      NR=NP
  140 IF(NR.GE.3) THEN
        PDRMIN=2.*PARU12
        DO 150 I=N+1,N+NR
          IF(I.EQ.N+NR.AND.IABS(K(N+1,2)).NE.21) GOTO 150
          I1=I+1
          IF(I.EQ.N+NR) I1=N+1
          IF(K(I,1).EQ.41.OR.K(I1,1).EQ.41) GOTO 150
          IF(MJU(1).NE.0.AND.I1.LT.MJU(1).AND.IABS(K(I1,2)).NE.21)
     +    GOTO 150
          IF(MJU(2).NE.0.AND.I.GT.MJU(2).AND.IABS(K(I,2)).NE.21) GOTO
     +    150
          PAP=SQRT((P(I,1)**2+P(I,2)**2+P(I,3)**2)*(P(I1,1)**2+ P(I1,2)
     +    **2+P(I1,3)**2))
          PVP=P(I,1)*P(I1,1)+P(I,2)*P(I1,2)+P(I,3)*P(I1,3)
          PDR=4.*(PAP-PVP)**2/MAX(1E-6,PARU13**2*PAP+2.*(PAP-PVP))
          IF(PDR.LT.PDRMIN) THEN
            IR=I
            PDRMIN=PDR
          ENDIF
  150   CONTINUE

C...Recombine very nearby partons to avoid machine precision problems.
        IF(PDRMIN.LT.PARU12.AND.IR.EQ.N+NR) THEN
          DO 160 J=1,4
            P(N+1,J)=P(N+1,J)+P(N+NR,J)
  160     CONTINUE
          P(N+1,5)=SQRT(MAX(0.,P(N+1,4)**2-P(N+1,1)**2-P(N+1,2)**2-
     +    P(N+1,3)**2))
          NR=NR-1
          GOTO 140
        ELSEIF(PDRMIN.LT.PARU12) THEN
          DO 170 J=1,4
            P(IR,J)=P(IR,J)+P(IR+1,J)
  170     CONTINUE
          P(IR,5)=SQRT(MAX(0.,P(IR,4)**2-P(IR,1)**2-P(IR,2)**2-
     +    P(IR,3)**2))
          DO 190 I=IR+1,N+NR-1
            K(I,2)=K(I+1,2)
            DO 180 J=1,5
              P(I,J)=P(I+1,J)
  180       CONTINUE
  190     CONTINUE
          IF(IR.EQ.N+NR-1) K(IR,2)=K(N+NR,2)
          NR=NR-1
          IF(MJU(1).GT.IR) MJU(1)=MJU(1)-1
          IF(MJU(2).GT.IR) MJU(2)=MJU(2)-1
          GOTO 140
        ENDIF
      ENDIF
      NTRYR=NTRYR+1

C...Reset particle counter. Skip ahead if no junctions are present;
C...this is usually the case!
      NRS=MAX(5*NR+11,NP)
      NTRY=0
  200 NTRY=NTRY+1
      IF(NTRY.GT.MAXPZTRY.AND.NTRYR.LE.MAXPZTRYR) THEN
        PARU12=4.*PARU12
        PARU13=2.*PARU13
        GOTO 140
      ELSEIF(NTRY.GT.MAXPZTRY) THEN
        CALL LUERRM(14,'(LUSTRF:) caught in infinite loop')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      I=N+NRS
      MSTU(90)=MSTU90
      IF(MJU(1).EQ.0.AND.MJU(2).EQ.0) GOTO 580
      DO 570 JT=1,2
        NJS(JT)=0
        IF(MJU(JT).EQ.0) GOTO 570
        JS=3-2*JT

C...Find and sum up momentum on three sides of junction. Check flavours.
        DO 220 IU=1,3
          IJU(IU)=0
          DO 210 J=1,5
            PJU(IU,J)=0.
  210     CONTINUE
  220   CONTINUE
        IU=0
        DO 240 I1=N+1+(JT-1)*(NR-1),N+NR+(JT-1)*(1-NR),JS
          IF(K(I1,2).NE.21.AND.IU.LE.2) THEN
            IU=IU+1
            IJU(IU)=I1
          ENDIF
          DO 230 J=1,4
            PJU(IU,J)=PJU(IU,J)+P(I1,J)
  230     CONTINUE
  240   CONTINUE
        DO 250 IU=1,3
          PJU(IU,5)=SQRT(PJU(IU,1)**2+PJU(IU,2)**2+PJU(IU,3)**2)
  250   CONTINUE
        IF(K(IJU(3),2)/100.NE.10*K(IJU(1),2)+K(IJU(2),2).AND. K(IJU(3),
     +  2)/100.NE.10*K(IJU(2),2)+K(IJU(1),2)) THEN
          CALL LUERRM(12,'(LUSTRF:) unphysical flavour combination')
          IF(MSTU(21).GE.1) RETURN
        ENDIF

C...Calculate (approximate) boost to rest frame of junction.
        T12=(PJU(1,1)*PJU(2,1)+PJU(1,2)*PJU(2,2)+PJU(1,3)*PJU(2,3))/
     +  (PJU(1,5)*PJU(2,5))
        T13=(PJU(1,1)*PJU(3,1)+PJU(1,2)*PJU(3,2)+PJU(1,3)*PJU(3,3))/
     +  (PJU(1,5)*PJU(3,5))
        T23=(PJU(2,1)*PJU(3,1)+PJU(2,2)*PJU(3,2)+PJU(2,3)*PJU(3,3))/
     +  (PJU(2,5)*PJU(3,5))
        T11=SQRT((2./3.)*(1.-T12)*(1.-T13)/(1.-T23))
        T22=SQRT((2./3.)*(1.-T12)*(1.-T23)/(1.-T13))
        TSQ=SQRT((2.*T11*T22+T12-1.)*(1.+T12))
        T1F=(TSQ-T22*(1.+T12))/(1.-T12**2)
        T2F=(TSQ-T11*(1.+T12))/(1.-T12**2)
        DO 260 J=1,3
          TJU(J)=-(T1F*PJU(1,J)/PJU(1,5)+T2F*PJU(2,J)/PJU(2,5))
  260   CONTINUE
        TJU(4)=SQRT(1.+TJU(1)**2+TJU(2)**2+TJU(3)**2)
        DO 270 IU=1,3
          PJU(IU,5)=TJU(4)*PJU(IU,4)-TJU(1)*PJU(IU,1)-TJU(2)*PJU(IU,2)-
     +    TJU(3)*PJU(IU,3)
  270   CONTINUE

C...Put junction at rest if motion could give inconsistencies.
        IF(PJU(1,5)+PJU(2,5).GT.PJU(1,4)+PJU(2,4)) THEN
          DO 280 J=1,3
            TJU(J)=0.
  280     CONTINUE
          TJU(4)=1.
          PJU(1,5)=PJU(1,4)
          PJU(2,5)=PJU(2,4)
          PJU(3,5)=PJU(3,4)
        ENDIF

C...Start preparing for fragmentation of two strings from junction.
        ISTA=I
        DO 550 IU=1,2
          NS=IJU(IU+1)-IJU(IU)

C...Junction strings: find longitudinal string directions.
          DO 310 IS=1,NS
            IS1=IJU(IU)+IS-1
            IS2=IJU(IU)+IS
            DO 290 J=1,5
              DP(1,J)=0.5*P(IS1,J)
              IF(IS.EQ.1) DP(1,J)=P(IS1,J)
              DP(2,J)=0.5*P(IS2,J)
              IF(IS.EQ.NS) DP(2,J)=-PJU(IU,J)
  290       CONTINUE
            IF(IS.EQ.NS) DP(2,4)=SQRT(PJU(IU,1)**2+PJU(IU,2)**2+PJU(IU,
     +      3)**2)
            IF(IS.EQ.NS) DP(2,5)=0.
            DP(3,5)=DFOUR(1,1)
            DP(4,5)=DFOUR(2,2)
            DHKC=DFOUR(1,2)
            IF(DP(3,5)+2.*DHKC+DP(4,5).LE.0.) THEN
              DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
              DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
              DP(3,5)=0D0
              DP(4,5)=0D0
              DHKC=DFOUR(1,2)
            ENDIF
            DHKS=SQRT(DHKC**2-DP(3,5)*DP(4,5))
            DHK1=0.5*((DP(4,5)+DHKC)/DHKS-1.)
            DHK2=0.5*((DP(3,5)+DHKC)/DHKS-1.)
            IN1=N+NR+4*IS-3
            P(IN1,5)=SQRT(DP(3,5)+2.*DHKC+DP(4,5))
            DO 300 J=1,4
              P(IN1,J)=(1.+DHK1)*DP(1,J)-DHK2*DP(2,J)
              P(IN1+1,J)=(1.+DHK2)*DP(2,J)-DHK1*DP(1,J)
  300       CONTINUE
  310     CONTINUE

C...Junction strings: initialize flavour, momentum and starting pos.
          ISAV=I
          MSTU91=MSTU(90)
  320     NTRY=NTRY+1
          IF(NTRY.GT.MAXPZTRY.AND.NTRYR.LE.MAXPZTRYR) THEN
            PARU12=4.*PARU12
            PARU13=2.*PARU13
            GOTO 140
          ELSEIF(NTRY.GT.MAXPZTRY) THEN
            CALL LUERRM(14,'(LUSTRF:) caught in infinite loop')
            IF(MSTU(21).GE.1) RETURN
          ENDIF
          I=ISAV
          MSTU(90)=MSTU91
          IRANKJ=0
          IE(1)=K(N+1+(JT/2)*(NP-1),3)
          IN(4)=N+NR+1
          IN(5)=IN(4)+1
          IN(6)=N+NR+4*NS+1
          DO 340 JQ=1,2
            DO 330 IN1=N+NR+2+JQ,N+NR+4*NS-2+JQ,4
              P(IN1,1)=2-JQ
              P(IN1,2)=JQ-1
              P(IN1,3)=1.
  330       CONTINUE
  340     CONTINUE
          KFL(1)=K(IJU(IU),2)
          PX(1)=0.
          PY(1)=0.
          GAM(1)=0.
          DO 350 J=1,5
            PJU(IU+3,J)=0.
  350     CONTINUE

C...Junction strings: find initial transverse directions.
          DO 360 J=1,4
            DP(1,J)=P(IN(4),J)
            DP(2,J)=P(IN(4)+1,J)
            DP(3,J)=0.
            DP(4,J)=0.
  360     CONTINUE
          DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
          DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
          DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4)
          DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4)
          DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4)
          IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1.
          IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1.
          IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1.
          IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1.
          DHC12=DFOUR(1,2)
          DHCX1=DFOUR(3,1)/DHC12
          DHCX2=DFOUR(3,2)/DHC12
          DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12)
          DHCY1=DFOUR(4,1)/DHC12
          DHCY2=DFOUR(4,2)/DHC12
          DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12
          DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2)
          DO 370 J=1,4
            DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
            P(IN(6),J)=DP(3,J)
            P(IN(6)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     +      DHCYX*DP(3,J))
  370     CONTINUE

C...Junction strings: produce new particle, origin.
  380     I=I+1
          IF(2*I-NSAV.GE.MSTU(4)-MSTU(32)-5) THEN
            CALL LUERRM(11,'(LUSTRF:) no more memory left in LUJETS')
            IF(MSTU(21).GE.1) RETURN
          ENDIF
          IRANKJ=IRANKJ+1
          K(I,1)=1
          K(I,3)=IE(1)
          K(I,4)=0
          K(I,5)=0

C...Junction strings: generate flavour, hadron, pT, z and Gamma.
  390     CALL LUKFDI(KFL(1),0,KFL(3),K(I,2))
          IF(K(I,2).EQ.0) GOTO 320
          IF(MSTJ(12).GE.3.AND.IRANKJ.EQ.1.AND.IABS(KFL(1)).LE.10.AND.
     +    IABS(KFL(3)).GT.10) THEN
            IF(RLU(0).GT.PARJ(19)) GOTO 390
          ENDIF
          P(I,5)=ULMASS(K(I,2))
          CALL LUPTDI(KFL(1),PX(3),PY(3))
          PR(1)=P(I,5)**2+(PX(1)+PX(3))**2+(PY(1)+PY(3))**2
          CALL LUZDIS(KFL(1),KFL(3),PR(1),Z)
          IF(IABS(KFL(1)).GE.4.AND.IABS(KFL(1)).LE.8.AND. MSTU(90)
     +    .LT.8) THEN
            MSTU(90)=MSTU(90)+1
            MSTU(90+MSTU(90))=I
            PARU(90+MSTU(90))=Z
          ENDIF
          GAM(3)=(1.-Z)*(GAM(1)+PR(1)/Z)
          DO 400 J=1,3
            IN(J)=IN(3+J)
  400     CONTINUE

C...Junction strings: stepping within or from 'low' string region easy.
          IF(IN(1)+1.EQ.IN(2).AND.Z*P(IN(1)+2,3)*P(IN(2)+2,3)* P(IN(1),
     +    5)**2.GE.PR(1)) THEN
            P(IN(1)+2,4)=Z*P(IN(1)+2,3)
            P(IN(2)+2,4)=PR(1)/(P(IN(1)+2,4)*P(IN(1),5)**2)
            DO 410 J=1,4
              P(I,J)=(PX(1)+PX(3))*P(IN(3),J)+(PY(1)+PY(3))*P(IN(3)+1,
     +        J)
  410       CONTINUE
            GOTO 500
          ELSEIF(IN(1)+1.EQ.IN(2)) THEN
            P(IN(2)+2,4)=P(IN(2)+2,3)
            P(IN(2)+2,1)=1.
            IN(2)=IN(2)+4
            IF(IN(2).GT.N+NR+4*NS) GOTO 320
            IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN
              P(IN(1)+2,4)=P(IN(1)+2,3)
              P(IN(1)+2,1)=0.
              IN(1)=IN(1)+4
            ENDIF
          ENDIF

C...Junction strings: find new transverse directions.
  420     IF(IN(1).GT.N+NR+4*NS.OR.IN(2).GT.N+NR+4*NS.OR. IN(1)
     +    .GT.IN(2)) GOTO 320
          IF(IN(1).NE.IN(4).OR.IN(2).NE.IN(5)) THEN
            DO 430 J=1,4
              DP(1,J)=P(IN(1),J)
              DP(2,J)=P(IN(2),J)
              DP(3,J)=0.
              DP(4,J)=0.
  430       CONTINUE
            DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
            DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
            DHC12=DFOUR(1,2)
            IF(DHC12.LE.1E-2) THEN
              P(IN(1)+2,4)=P(IN(1)+2,3)
              P(IN(1)+2,1)=0.
              IN(1)=IN(1)+4
              GOTO 420
            ENDIF
            IN(3)=N+NR+4*NS+5
            DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4)
            DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4)
            DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4)
            IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1.
            IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1.
            IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1.
            IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1.
            DHCX1=DFOUR(3,1)/DHC12
            DHCX2=DFOUR(3,2)/DHC12
            DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12)
            DHCY1=DFOUR(4,1)/DHC12
            DHCY2=DFOUR(4,2)/DHC12
            DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12
            DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2)
            DO 440 J=1,4
              DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
              P(IN(3),J)=DP(3,J)
              P(IN(3)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     +        DHCYX*DP(3,J))
  440       CONTINUE
C...Express pT with respect to new axes, if sensible.
            PXP=-(PX(3)*FOUR(IN(6),IN(3))+PY(3)*FOUR(IN(6)+1,IN(3)))
            PYP=-(PX(3)*FOUR(IN(6),IN(3)+1)+PY(3)*FOUR(IN(6)+1,IN(3)+1)
     +      )
            IF(ABS(PXP**2+PYP**2-PX(3)**2-PY(3)**2).LT.0.01) THEN
              PX(3)=PXP
              PY(3)=PYP
            ENDIF
          ENDIF

C...Junction strings: sum up known four-momentum, coefficients for m2.
          DO 470 J=1,4
            DHG(J)=0.
            P(I,J)=PX(1)*P(IN(6),J)+PY(1)*P(IN(6)+1,J)+PX(3)*P(IN(3),J)
     +      + PY(3)*P(IN(3)+1,J)
            DO 450 IN1=IN(4),IN(1)-4,4
              P(I,J)=P(I,J)+P(IN1+2,3)*P(IN1,J)
  450       CONTINUE
            DO 460 IN2=IN(5),IN(2)-4,4
              P(I,J)=P(I,J)+P(IN2+2,3)*P(IN2,J)
  460       CONTINUE
  470     CONTINUE
          DHM(1)=FOUR(I,I)
          DHM(2)=2.*FOUR(I,IN(1))
          DHM(3)=2.*FOUR(I,IN(2))
          DHM(4)=2.*FOUR(IN(1),IN(2))

C...Junction strings: find coefficients for Gamma expression.
          DO 490 IN2=IN(1)+1,IN(2),4
            DO 480 IN1=IN(1),IN2-1,4
              DHC=2.*FOUR(IN1,IN2)
              DHG(1)=DHG(1)+P(IN1+2,1)*P(IN2+2,1)*DHC
              IF(IN1.EQ.IN(1)) DHG(2)=DHG(2)-P(IN2+2,1)*DHC
              IF(IN2.EQ.IN(2)) DHG(3)=DHG(3)+P(IN1+2,1)*DHC
              IF(IN1.EQ.IN(1).AND.IN2.EQ.IN(2)) DHG(4)=DHG(4)-DHC
  480       CONTINUE
  490     CONTINUE

C...Junction strings: solve (m2, Gamma) equation system for energies.
          DHS1=DHM(3)*DHG(4)-DHM(4)*DHG(3)
          IF(ABS(DHS1).LT.1E-4) GOTO 320
          DHS2=DHM(4)*(GAM(3)-DHG(1))-DHM(2)*DHG(3)-DHG(4)* (P(I,5)**2-
     +    DHM(1))+DHG(2)*DHM(3)
          DHS3=DHM(2)*(GAM(3)-DHG(1))-DHG(2)*(P(I,5)**2-DHM(1))
          P(IN(2)+2,4)=0.5*(SQRT(MAX(0D0,DHS2**2-4.*DHS1*DHS3))/
     +    ABS(DHS1)- DHS2/DHS1)
          IF(DHM(2)+DHM(4)*P(IN(2)+2,4).LE.0.) GOTO 320
          P(IN(1)+2,4)=(P(I,5)**2-DHM(1)-DHM(3)*P(IN(2)+2,4))/ (DHM(2)+
     +    DHM(4)*P(IN(2)+2,4))

C...Junction strings: step to new region if necessary.
          IF(P(IN(2)+2,4).GT.P(IN(2)+2,3)) THEN
            P(IN(2)+2,4)=P(IN(2)+2,3)
            P(IN(2)+2,1)=1.
            IN(2)=IN(2)+4
            IF(IN(2).GT.N+NR+4*NS) GOTO 320
            IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN
              P(IN(1)+2,4)=P(IN(1)+2,3)
              P(IN(1)+2,1)=0.
              IN(1)=IN(1)+4
            ENDIF
            GOTO 420
          ELSEIF(P(IN(1)+2,4).GT.P(IN(1)+2,3)) THEN
            P(IN(1)+2,4)=P(IN(1)+2,3)
            P(IN(1)+2,1)=0.
            IN(1)=IN(1)+JS
            GOTO 820
          ENDIF

C...Junction strings: particle four-momentum, remainder, loop back.
  500     DO 510 J=1,4
            P(I,J)=P(I,J)+P(IN(1)+2,4)*P(IN(1),J)+P(IN(2)+2,4)*P(IN(2),
     +      J)
            PJU(IU+3,J)=PJU(IU+3,J)+P(I,J)
  510     CONTINUE
          IF(P(I,4).LT.P(I,5)) GOTO 320
          PJU(IU+3,5)=TJU(4)*PJU(IU+3,4)-TJU(1)*PJU(IU+3,1)- TJU(2)*
     +    PJU(IU+3,2)-TJU(3)*PJU(IU+3,3)
          IF(PJU(IU+3,5).LT.PJU(IU,5)) THEN
            KFL(1)=-KFL(3)
            PX(1)=-PX(3)
            PY(1)=-PY(3)
            GAM(1)=GAM(3)
            IF(IN(3).NE.IN(6)) THEN
              DO 520 J=1,4
                P(IN(6),J)=P(IN(3),J)
                P(IN(6)+1,J)=P(IN(3)+1,J)
  520         CONTINUE
            ENDIF
            DO 530 JQ=1,2
              IN(3+JQ)=IN(JQ)
              P(IN(JQ)+2,3)=P(IN(JQ)+2,3)-P(IN(JQ)+2,4)
              P(IN(JQ)+2,1)=P(IN(JQ)+2,1)-(3-2*JQ)*P(IN(JQ)+2,4)
  530       CONTINUE
            GOTO 380
          ENDIF

C...Junction strings: save quantities left after each string.
          IF(IABS(KFL(1)).GT.10) GOTO 320
          I=I-1
          KFJH(IU)=KFL(1)
          DO 540 J=1,4
            PJU(IU+3,J)=PJU(IU+3,J)-P(I+1,J)
  540     CONTINUE
  550   CONTINUE

C...Junction strings: put together to new effective string endpoint.
        NJS(JT)=I-ISTA
        KFJS(JT)=K(K(MJU(JT+2),3),2)
        KFLS=2*INT(RLU(0)+3.*PARJ(4)/(1.+3.*PARJ(4)))+1
        IF(KFJH(1).EQ.KFJH(2)) KFLS=3
        IF(ISTA.NE.I) KFJS(JT)=ISIGN(1000*MAX(IABS(KFJH(1)), IABS(KFJH(
     +  2)))+100*MIN(IABS(KFJH(1)),IABS(KFJH(2)))+ KFLS,KFJH(1))
        DO 560 J=1,4
          PJS(JT,J)=PJU(1,J)+PJU(2,J)+P(MJU(JT),J)
          PJS(JT+2,J)=PJU(4,J)+PJU(5,J)
  560   CONTINUE
        PJS(JT,5)=SQRT(MAX(0.,PJS(JT,4)**2-PJS(JT,1)**2-PJS(JT,2)**2-
     +  PJS(JT,3)**2))
  570 CONTINUE

C...Open versus closed strings. Choose breakup region for latter.
  580 IF(MJU(1).NE.0.AND.MJU(2).NE.0) THEN
        NS=MJU(2)-MJU(1)
        NB=MJU(1)-N
      ELSEIF(MJU(1).NE.0) THEN
        NS=N+NR-MJU(1)
        NB=MJU(1)-N
      ELSEIF(MJU(2).NE.0) THEN
        NS=MJU(2)-N
        NB=1
      ELSEIF(IABS(K(N+1,2)).NE.21) THEN
        NS=NR-1
        NB=1
      ELSE
        NS=NR+1
        W2SUM=0.
        DO 590 IS=1,NR
          P(N+NR+IS,1)=0.5*FOUR(N+IS,N+IS+1-NR*(IS/NR))
          W2SUM=W2SUM+P(N+NR+IS,1)
  590   CONTINUE
        W2RAN=RLU(0)*W2SUM
        NB=0
  600   NB=NB+1
        W2SUM=W2SUM-P(N+NR+NB,1)
        IF(W2SUM.GT.W2RAN.AND.NB.LT.NR) GOTO 600
      ENDIF

C...Find longitudinal string directions (i.e. lightlike four-vectors).
      DO 630 IS=1,NS
        IS1=N+IS+NB-1-NR*((IS+NB-2)/NR)
        IS2=N+IS+NB-NR*((IS+NB-1)/NR)
        DO 610 J=1,5
          DP(1,J)=P(IS1,J)
          IF(IABS(K(IS1,2)).EQ.21) DP(1,J)=0.5*DP(1,J)
          IF(IS1.EQ.MJU(1)) DP(1,J)=PJS(1,J)-PJS(3,J)
          DP(2,J)=P(IS2,J)
          IF(IABS(K(IS2,2)).EQ.21) DP(2,J)=0.5*DP(2,J)
          IF(IS2.EQ.MJU(2)) DP(2,J)=PJS(2,J)-PJS(4,J)
  610   CONTINUE
        DP(3,5)=DFOUR(1,1)
        DP(4,5)=DFOUR(2,2)
        DHKC=DFOUR(1,2)
        IF(DP(3,5)+2.*DHKC+DP(4,5).LE.0.) THEN
          DP(3,5)=DP(1,5)**2
          DP(4,5)=DP(2,5)**2
          DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2+DP(1,5)**2)
          DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2+DP(2,5)**2)
          DHKC=DFOUR(1,2)
        ENDIF
        DHKS=SQRT(DHKC**2-DP(3,5)*DP(4,5))
        DHK1=0.5*((DP(4,5)+DHKC)/DHKS-1.)
        DHK2=0.5*((DP(3,5)+DHKC)/DHKS-1.)
        IN1=N+NR+4*IS-3
        P(IN1,5)=SQRT(DP(3,5)+2.*DHKC+DP(4,5))
        DO 620 J=1,4
          P(IN1,J)=(1.+DHK1)*DP(1,J)-DHK2*DP(2,J)
          P(IN1+1,J)=(1.+DHK2)*DP(2,J)-DHK1*DP(1,J)
  620   CONTINUE
  630 CONTINUE

C...Begin initialization: sum up energy, set starting position.
      ISAV=I
      MSTU91=MSTU(90)
  640 NTRY=NTRY+1
      IF(NTRY.GT.MAXPZTRY.AND.NTRYR.LE.MAXPZTRYR) THEN
        PARU12=4.*PARU12
        PARU13=2.*PARU13
        GOTO 140
      ELSEIF(NTRY.GT.MAXPZTRY) THEN
        CALL LUERRM(14,'(LUSTRF:) caught in infinite loop')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      I=ISAV
      MSTU(90)=MSTU91
      DO 660 J=1,4
        P(N+NRS,J)=0.
        DO 650 IS=1,NR
          P(N+NRS,J)=P(N+NRS,J)+P(N+IS,J)
  650   CONTINUE
  660 CONTINUE
      DO 680 JT=1,2
        IRANK(JT)=0
        IF(MJU(JT).NE.0) IRANK(JT)=NJS(JT)
        IF(NS.GT.NR) IRANK(JT)=1
        IE(JT)=K(N+1+(JT/2)*(NP-1),3)
        IN(3*JT+1)=N+NR+1+4*(JT/2)*(NS-1)
        IN(3*JT+2)=IN(3*JT+1)+1
        IN(3*JT+3)=N+NR+4*NS+2*JT-1
        DO 670 IN1=N+NR+2+JT,N+NR+4*NS-2+JT,4
          P(IN1,1)=2-JT
          P(IN1,2)=JT-1
          P(IN1,3)=1.
  670   CONTINUE
  680 CONTINUE

C...Initialize flavour and pT variables for open string.
      IF(NS.LT.NR) THEN
        PX(1)=0.
        PY(1)=0.
        IF(NS.EQ.1.AND.MJU(1)+MJU(2).EQ.0) CALL LUPTDI(0,PX(1),PY(1))
        PX(2)=-PX(1)
        PY(2)=-PY(1)
        DO 690 JT=1,2
          KFL(JT)=K(IE(JT),2)
          IF(MJU(JT).NE.0) KFL(JT)=KFJS(JT)
          MSTJ(93)=1
          PMQ(JT)=ULMASS(KFL(JT))
          GAM(JT)=0.
  690   CONTINUE

C...Closed string: random initial breakup flavour, pT and vertex.
      ELSE
        KFL(3)=INT(1.+(2.+PARJ(2))*RLU(0))*(-1)**INT(RLU(0)+0.5)
        CALL LUKFDI(KFL(3),0,KFL(1),KDUMP)
        KFL(2)=-KFL(1)
        IF(IABS(KFL(1)).GT.10.AND.RLU(0).GT.0.5) THEN
          KFL(2)=-(KFL(1)+ISIGN(10000,KFL(1)))
        ELSEIF(IABS(KFL(1)).GT.10) THEN
          KFL(1)=-(KFL(2)+ISIGN(10000,KFL(2)))
        ENDIF
        CALL LUPTDI(KFL(1),PX(1),PY(1))
        PX(2)=-PX(1)
        PY(2)=-PY(1)
        PR3=MIN(25.,0.1*P(N+NR+1,5)**2)
  700   CALL LUZDIS(KFL(1),KFL(2),PR3,Z)
        ZR=PR3/(Z*P(N+NR+1,5)**2)
        IF(ZR.GE.1.) GOTO 700
        DO 710 JT=1,2
          MSTJ(93)=1
          PMQ(JT)=ULMASS(KFL(JT))
          GAM(JT)=PR3*(1.-Z)/Z
          IN1=N+NR+3+4*(JT/2)*(NS-1)
          P(IN1,JT)=1.-Z
          P(IN1,3-JT)=JT-1
          P(IN1,3)=(2-JT)*(1.-Z)+(JT-1)*Z
          P(IN1+1,JT)=ZR
          P(IN1+1,3-JT)=2-JT
          P(IN1+1,3)=(2-JT)*(1.-ZR)+(JT-1)*ZR
  710   CONTINUE
      ENDIF

C...Find initial transverse directions (i.e. spacelike four-vectors).
      DO 750 JT=1,2
        IF(JT.EQ.1.OR.NS.EQ.NR-1) THEN
          IN1=IN(3*JT+1)
          IN3=IN(3*JT+3)
          DO 720 J=1,4
            DP(1,J)=P(IN1,J)
            DP(2,J)=P(IN1+1,J)
            DP(3,J)=0.
            DP(4,J)=0.
  720     CONTINUE
          DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
          DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
          DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4)
          DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4)
          DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4)
          IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1.
          IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1.
          IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1.
          IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1.
          DHC12=DFOUR(1,2)
          DHCX1=DFOUR(3,1)/DHC12
          DHCX2=DFOUR(3,2)/DHC12
          DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12)
          DHCY1=DFOUR(4,1)/DHC12
          DHCY2=DFOUR(4,2)/DHC12
          DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12
          DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2)
          DO 730 J=1,4
            DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
            P(IN3,J)=DP(3,J)
            P(IN3+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     +      DHCYX*DP(3,J))
  730     CONTINUE
        ELSE
          DO 740 J=1,4
            P(IN3+2,J)=P(IN3,J)
            P(IN3+3,J)=P(IN3+1,J)
  740     CONTINUE
        ENDIF
  750 CONTINUE

C...Remove energy used up in junction string fragmentation.
      IF(MJU(1)+MJU(2).GT.0) THEN
        DO 770 JT=1,2
          IF(NJS(JT).EQ.0) GOTO 770
          DO 760 J=1,4
            P(N+NRS,J)=P(N+NRS,J)-PJS(JT+2,J)
  760     CONTINUE
  770   CONTINUE
      ENDIF

C...Produce new particle: side, origin.
  780 I=I+1
      IF(2*I-NSAV.GE.MSTU(4)-MSTU(32)-5) THEN
        CALL LUERRM(11,'(LUSTRF:) no more memory left in LUJETS')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      JT=1.5+RLU(0)
      IF(IABS(KFL(3-JT)).GT.10) JT=3-JT
      IF(IABS(KFL(3-JT)).GE.4.AND.IABS(KFL(3-JT)).LE.8) JT=3-JT
      JR=3-JT
      JS=3-2*JT
      IRANK(JT)=IRANK(JT)+1
      K(I,1)=1
      K(I,3)=IE(JT)
      K(I,4)=0
      K(I,5)=0

C...Generate flavour, hadron and pT.
  790 CALL LUKFDI(KFL(JT),0,KFL(3),K(I,2))
      IF(K(I,2).EQ.0) GOTO 640
      IF(MSTJ(12).GE.3.AND.IRANK(JT).EQ.1.AND.IABS(KFL(JT)).LE.10.AND.
     +IABS(KFL(3)).GT.10) THEN
        IF(RLU(0).GT.PARJ(19)) GOTO 790
      ENDIF
      P(I,5)=ULMASS(K(I,2))
      CALL LUPTDI(KFL(JT),PX(3),PY(3))
      PR(JT)=P(I,5)**2+(PX(JT)+PX(3))**2+(PY(JT)+PY(3))**2

C...Final hadrons for small invariant mass.
      MSTJ(93)=1
      PMQ(3)=ULMASS(KFL(3))
      PARJST=PARJ(33)
      IF(MSTJ(11).EQ.2) PARJST=PARJ(34)
      WMIN=PARJST+PMQ(1)+PMQ(2)+PARJ(36)*PMQ(3)
      IF(IABS(KFL(JT)).GT.10.AND.IABS(KFL(3)).GT.10) WMIN=
     +WMIN-0.5*PARJ(36)*PMQ(3)
      WREM2=FOUR(N+NRS,N+NRS)
      IF(WREM2.LT.0.10) GOTO 640
      IF(WREM2.LT.MAX(WMIN*(1.+(2.*RLU(0)-1.)*PARJ(37)),
     +PARJ(32)+PMQ(1)+PMQ(2))**2) GOTO 940

C...Choose z, which gives Gamma. Shift z for heavy flavours.
      CALL LUZDIS(KFL(JT),KFL(3),PR(JT),Z)
      IF(IABS(KFL(JT)).GE.4.AND.IABS(KFL(JT)).LE.8.AND.
     +MSTU(90).LT.8) THEN
        MSTU(90)=MSTU(90)+1
        MSTU(90+MSTU(90))=I
        PARU(90+MSTU(90))=Z
      ENDIF
      KFL1A=IABS(KFL(1))
      KFL2A=IABS(KFL(2))
      IF(MAX(MOD(KFL1A,10),MOD(KFL1A/1000,10),MOD(KFL2A,10),
     +MOD(KFL2A/1000,10)).GE.4) THEN
        PR(JR)=(PMQ(JR)+PMQ(3))**2+(PX(JR)-PX(3))**2+(PY(JR)-PY(3))**2
        PW12=SQRT(MAX(0.,(WREM2-PR(1)-PR(2))**2-4.*PR(1)*PR(2)))
        Z=(WREM2+PR(JT)-PR(JR)+PW12*(2.*Z-1.))/(2.*WREM2)
        PR(JR)=(PMQ(JR)+PARJST)**2+(PX(JR)-PX(3))**2+(PY(JR)-PY(3))**2
        IF((1.-Z)*(WREM2-PR(JT)/Z).LT.PR(JR)) GOTO 940
      ENDIF
      GAM(3)=(1.-Z)*(GAM(JT)+PR(JT)/Z)
      DO 800 J=1,3
        IN(J)=IN(3*JT+J)
  800 CONTINUE

C...Stepping within or from 'low' string region easy.
      IF(IN(1)+1.EQ.IN(2).AND.Z*P(IN(1)+2,3)*P(IN(2)+2,3)*
     +P(IN(1),5)**2.GE.PR(JT)) THEN
        P(IN(JT)+2,4)=Z*P(IN(JT)+2,3)
        P(IN(JR)+2,4)=PR(JT)/(P(IN(JT)+2,4)*P(IN(1),5)**2)
        DO 810 J=1,4
          P(I,J)=(PX(JT)+PX(3))*P(IN(3),J)+(PY(JT)+PY(3))*P(IN(3)+1,J)
  810   CONTINUE
        GOTO 900
      ELSEIF(IN(1)+1.EQ.IN(2)) THEN
        P(IN(JR)+2,4)=P(IN(JR)+2,3)
        P(IN(JR)+2,JT)=1.
        IN(JR)=IN(JR)+4*JS
        IF(JS*IN(JR).GT.JS*IN(4*JR)) GOTO 640
        IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN
          P(IN(JT)+2,4)=P(IN(JT)+2,3)
          P(IN(JT)+2,JT)=0.
          IN(JT)=IN(JT)+4*JS
        ENDIF
      ENDIF

C...Find new transverse directions (i.e. spacelike string vectors).
  820 IF(JS*IN(1).GT.JS*IN(3*JR+1).OR.JS*IN(2).GT.JS*IN(3*JR+2).OR.
     +IN(1).GT.IN(2)) GOTO 640
      IF(IN(1).NE.IN(3*JT+1).OR.IN(2).NE.IN(3*JT+2)) THEN
        DO 830 J=1,4
          DP(1,J)=P(IN(1),J)
          DP(2,J)=P(IN(2),J)
          DP(3,J)=0.
          DP(4,J)=0.
  830   CONTINUE
        DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
        DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
        DHC12=DFOUR(1,2)
        IF(DHC12.LE.1E-2) THEN
          P(IN(JT)+2,4)=P(IN(JT)+2,3)
          P(IN(JT)+2,JT)=0.
          IN(JT)=IN(JT)+4*JS
          GOTO 820
        ENDIF
        IN(3)=N+NR+4*NS+5
        DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4)
        DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4)
        DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4)
        IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1.
        IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1.
        IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1.
        IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1.
        DHCX1=DFOUR(3,1)/DHC12
        DHCX2=DFOUR(3,2)/DHC12
        DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12)
        DHCY1=DFOUR(4,1)/DHC12
        DHCY2=DFOUR(4,2)/DHC12
        DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12
        DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2)
        DO 840 J=1,4
          DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
          P(IN(3),J)=DP(3,J)
          P(IN(3)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     +    DHCYX*DP(3,J))
  840   CONTINUE
C...Express pT with respect to new axes, if sensible.
        PXP=-(PX(3)*FOUR(IN(3*JT+3),IN(3))+PY(3)*
     +  FOUR(IN(3*JT+3)+1,IN(3)))
        PYP=-(PX(3)*FOUR(IN(3*JT+3),IN(3)+1)+PY(3)*
     +  FOUR(IN(3*JT+3)+1,IN(3)+1))
        IF(ABS(PXP**2+PYP**2-PX(3)**2-PY(3)**2).LT.0.01) THEN
          PX(3)=PXP
          PY(3)=PYP
        ENDIF
      ENDIF

C...Sum up known four-momentum. Gives coefficients for m2 expression.
      DO 870 J=1,4
        DHG(J)=0.
        P(I,J)=PX(JT)*P(IN(3*JT+3),J)+PY(JT)*P(IN(3*JT+3)+1,J)+ PX(3)*
     +  P(IN(3),J)+PY(3)*P(IN(3)+1,J)
        DO 850 IN1=IN(3*JT+1),IN(1)-4*JS,4*JS
          P(I,J)=P(I,J)+P(IN1+2,3)*P(IN1,J)
  850   CONTINUE
        DO 860 IN2=IN(3*JT+2),IN(2)-4*JS,4*JS
          P(I,J)=P(I,J)+P(IN2+2,3)*P(IN2,J)
  860   CONTINUE
  870 CONTINUE
      DHM(1)=FOUR(I,I)
      DHM(2)=2.*FOUR(I,IN(1))
      DHM(3)=2.*FOUR(I,IN(2))
      DHM(4)=2.*FOUR(IN(1),IN(2))

C...Find coefficients for Gamma expression.
      DO 890 IN2=IN(1)+1,IN(2),4
        DO 880 IN1=IN(1),IN2-1,4
          DHC=2.*FOUR(IN1,IN2)
          DHG(1)=DHG(1)+P(IN1+2,JT)*P(IN2+2,JT)*DHC
          IF(IN1.EQ.IN(1)) DHG(2)=DHG(2)-JS*P(IN2+2,JT)*DHC
          IF(IN2.EQ.IN(2)) DHG(3)=DHG(3)+JS*P(IN1+2,JT)*DHC
          IF(IN1.EQ.IN(1).AND.IN2.EQ.IN(2)) DHG(4)=DHG(4)-DHC
  880   CONTINUE
  890 CONTINUE

C...Solve (m2, Gamma) equation system for energies taken.
      DHS1=DHM(JR+1)*DHG(4)-DHM(4)*DHG(JR+1)
      IF(ABS(DHS1).LT.1E-4) GOTO 640
      DHS2=DHM(4)*(GAM(3)-DHG(1))-DHM(JT+1)*DHG(JR+1)-DHG(4)*
     +(P(I,5)**2-DHM(1))+DHG(JT+1)*DHM(JR+1)
      DHS3=DHM(JT+1)*(GAM(3)-DHG(1))-DHG(JT+1)*(P(I,5)**2-DHM(1))
      P(IN(JR)+2,4)=0.5*(SQRT(MAX(0D0,DHS2**2-4.*DHS1*DHS3))/ABS(DHS1)-
     +DHS2/DHS1)
      IF(DHM(JT+1)+DHM(4)*P(IN(JR)+2,4).LE.0.) GOTO 640
      P(IN(JT)+2,4)=(P(I,5)**2-DHM(1)-DHM(JR+1)*P(IN(JR)+2,4))/
     +(DHM(JT+1)+DHM(4)*P(IN(JR)+2,4))

C...Step to new region if necessary.
      IF(P(IN(JR)+2,4).GT.P(IN(JR)+2,3)) THEN
        P(IN(JR)+2,4)=P(IN(JR)+2,3)
        P(IN(JR)+2,JT)=1.
        IN(JR)=IN(JR)+4*JS
        IF(JS*IN(JR).GT.JS*IN(4*JR)) GOTO 640
        IF(FOUR(IN(1),IN(2)).LE.1E-2) THEN
          P(IN(JT)+2,4)=P(IN(JT)+2,3)
          P(IN(JT)+2,JT)=0.
          IN(JT)=IN(JT)+4*JS
        ENDIF
        GOTO 820
      ELSEIF(P(IN(JT)+2,4).GT.P(IN(JT)+2,3)) THEN
        P(IN(JT)+2,4)=P(IN(JT)+2,3)
        P(IN(JT)+2,JT)=0.
        IN(JT)=IN(JT)+4*JS
        GOTO 820
      ENDIF

C...Four-momentum of particle. Remaining quantities. Loop back.
  900 DO 910 J=1,4
        P(I,J)=P(I,J)+P(IN(1)+2,4)*P(IN(1),J)+P(IN(2)+2,4)*P(IN(2),J)
        P(N+NRS,J)=P(N+NRS,J)-P(I,J)
  910 CONTINUE
      IF(P(I,4).LT.P(I,5)) GOTO 640
      KFL(JT)=-KFL(3)
      PMQ(JT)=PMQ(3)
      PX(JT)=-PX(3)
      PY(JT)=-PY(3)
      GAM(JT)=GAM(3)
      IF(IN(3).NE.IN(3*JT+3)) THEN
        DO 920 J=1,4
          P(IN(3*JT+3),J)=P(IN(3),J)
          P(IN(3*JT+3)+1,J)=P(IN(3)+1,J)
  920   CONTINUE
      ENDIF
      DO 930 JQ=1,2
        IN(3*JT+JQ)=IN(JQ)
        P(IN(JQ)+2,3)=P(IN(JQ)+2,3)-P(IN(JQ)+2,4)
        P(IN(JQ)+2,JT)=P(IN(JQ)+2,JT)-JS*(3-2*JQ)*P(IN(JQ)+2,4)
  930 CONTINUE
      GOTO 780

C...Final hadron: side, flavour, hadron, mass.
  940 I=I+1
      K(I,1)=1
      K(I,3)=IE(JR)
      K(I,4)=0
      K(I,5)=0
      CALL LUKFDI(KFL(JR),-KFL(3),KFLDMP,K(I,2))
      IF(K(I,2).EQ.0) GOTO 640
      P(I,5)=ULMASS(K(I,2))
      PR(JR)=P(I,5)**2+(PX(JR)-PX(3))**2+(PY(JR)-PY(3))**2

C...Final two hadrons: find common setup of four-vectors.
      JQ=1
      IF(P(IN(4)+2,3)*P(IN(5)+2,3)*FOUR(IN(4),IN(5)).LT.P(IN(7),3)*
     +P(IN(8),3)*FOUR(IN(7),IN(8))) JQ=2
      DHC12=FOUR(IN(3*JQ+1),IN(3*JQ+2))
      DHR1=FOUR(N+NRS,IN(3*JQ+2))/DHC12
      DHR2=FOUR(N+NRS,IN(3*JQ+1))/DHC12
      IF(IN(4).NE.IN(7).OR.IN(5).NE.IN(8)) THEN
        PX(3-JQ)=-FOUR(N+NRS,IN(3*JQ+3))-PX(JQ)
        PY(3-JQ)=-FOUR(N+NRS,IN(3*JQ+3)+1)-PY(JQ)
        PR(3-JQ)=P(I+(JT+JQ-3)**2-1,5)**2+(PX(3-JQ)+(2*JQ-3)*JS*
     +  PX(3))**2+(PY(3-JQ)+(2*JQ-3)*JS*PY(3))**2
      ENDIF

C...Solve kinematics for final two hadrons, if possible.
      WREM2=WREM2+(PX(1)+PX(2))**2+(PY(1)+PY(2))**2
      FD=(SQRT(PR(1))+SQRT(PR(2)))/SQRT(WREM2)
      IF(MJU(1)+MJU(2).NE.0.AND.I.EQ.ISAV+2.AND.FD.GE.1.) GOTO 200
      IF(FD.GE.1.) GOTO 640
      FA=WREM2+PR(JT)-PR(JR)
      IF(MSTJ(11).NE.2) PREV=0.5*EXP(MAX(-50.,LOG(FD)*PARJ(38)*
     +(PR(1)+PR(2))**2))
      IF(MSTJ(11).EQ.2) PREV=0.5*FD**PARJ(39)
      FB=SIGN(SQRT(MAX(0.,FA**2-4.*WREM2*PR(JT))),JS*(RLU(0)-PREV))
      KFL1A=IABS(KFL(1))
      KFL2A=IABS(KFL(2))
      IF(MAX(MOD(KFL1A,10),MOD(KFL1A/1000,10),MOD(KFL2A,10),
     +MOD(KFL2A/1000,10)).GE.6) FB=SIGN(SQRT(MAX(0.,FA**2-
     +4.*WREM2*PR(JT))),FLOAT(JS))
      DO 950 J=1,4
        P(I-1,J)=(PX(JT)+PX(3))*P(IN(3*JQ+3),J)+(PY(JT)+PY(3))* P(IN(3*
     +  JQ+3)+1,J)+0.5*(DHR1*(FA+FB)*P(IN(3*JQ+1),J)+ DHR2*(FA-FB)*
     +  P(IN(3*JQ+2),J))/WREM2
        P(I,J)=P(N+NRS,J)-P(I-1,J)
  950 CONTINUE
      IF(P(I-1,4).LT.P(I-1,5).OR.P(I,4).LT.P(I,5)) GOTO 640

C...Mark jets as fragmented and give daughter pointers.
      N=I-NRS+1
      DO 960 I=NSAV+1,NSAV+NP
        IM=K(I,3)
        K(IM,1)=K(IM,1)+10
        IF(MSTU(16).NE.2) THEN
          K(IM,4)=NSAV+1
          K(IM,5)=NSAV+1
        ELSE
          K(IM,4)=NSAV+2
          K(IM,5)=N
        ENDIF
  960 CONTINUE

C...Document string system. Move up particles.
      NSAV=NSAV+1
      K(NSAV,1)=11
      K(NSAV,2)=92
      K(NSAV,3)=IP
      K(NSAV,4)=NSAV+1
      K(NSAV,5)=N
      DO 970 J=1,4
        P(NSAV,J)=DPS(J)
        V(NSAV,J)=V(IP,J)
  970 CONTINUE
      P(NSAV,5)=SQRT(MAX(0D0,DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2))
      V(NSAV,5)=0.
      DO 990 I=NSAV+1,N
        DO 980 J=1,5
          K(I,J)=K(I+NRS-1,J)
          P(I,J)=P(I+NRS-1,J)
          V(I,J)=0.
  980   CONTINUE
  990 CONTINUE
      MSTU91=MSTU(90)
      DO 1000 IZ=MSTU90+1,MSTU91
        MSTU9T(IZ)=MSTU(90+IZ)-NRS+1-NSAV+N
        PARU9T(IZ)=PARU(90+IZ)
 1000 CONTINUE
      MSTU(90)=MSTU90

C...Order particles in rank along the chain. Update mother pointer.
      DO 1020 I=NSAV+1,N
        DO 1010 J=1,5
          K(I-NSAV+N,J)=K(I,J)
          P(I-NSAV+N,J)=P(I,J)
 1010   CONTINUE
 1020 CONTINUE
      I1=NSAV
      DO 1050 I=N+1,2*N-NSAV
        IF(K(I,3).NE.IE(1)) GOTO 1050
        I1=I1+1
        DO 1030 J=1,5
          K(I1,J)=K(I,J)
          P(I1,J)=P(I,J)
 1030   CONTINUE
        IF(MSTU(16).NE.2) K(I1,3)=NSAV
        DO 1040 IZ=MSTU90+1,MSTU91
          IF(MSTU9T(IZ).EQ.I) THEN
            MSTU(90)=MSTU(90)+1
            MSTU(90+MSTU(90))=I1
            PARU(90+MSTU(90))=PARU9T(IZ)
          ENDIF
 1040   CONTINUE
 1050 CONTINUE
      DO 1080 I=2*N-NSAV,N+1,-1
        IF(K(I,3).EQ.IE(1)) GOTO 1080
        I1=I1+1
        DO 1060 J=1,5
          K(I1,J)=K(I,J)
          P(I1,J)=P(I,J)
 1060   CONTINUE
        IF(MSTU(16).NE.2) K(I1,3)=NSAV
        DO 1070 IZ=MSTU90+1,MSTU91
          IF(MSTU9T(IZ).EQ.I) THEN
            MSTU(90)=MSTU(90)+1
            MSTU(90+MSTU(90))=I1
            PARU(90+MSTU(90))=PARU9T(IZ)
          ENDIF
 1070   CONTINUE
 1080 CONTINUE

C...Boost back particle system. Set production vertices.
      IF(MBST.EQ.0) THEN
        MSTU(33)=1
        CALL LUDBRB(NSAV+1,N,0.,0.,DPS(1)/DPS(4),DPS(2)/DPS(4),
     +  DPS(3)/DPS(4))
      ELSE
        DO 1090 I=NSAV+1,N
          HHPMT=P(I,1)**2+P(I,2)**2+P(I,5)**2
          IF(P(I,3).GT.0.) THEN
            HHPEZ=(P(I,4)+P(I,3))*HHBZ
            P(I,3)=0.5*(HHPEZ-HHPMT/HHPEZ)
            P(I,4)=0.5*(HHPEZ+HHPMT/HHPEZ)
          ELSE
            HHPEZ=(P(I,4)-P(I,3))/HHBZ
            P(I,3)=-0.5*(HHPEZ-HHPMT/HHPEZ)
            P(I,4)=0.5*(HHPEZ+HHPMT/HHPEZ)
          ENDIF
 1090   CONTINUE
      ENDIF
      DO 1110 I=NSAV+1,N
        DO 1100 J=1,4
          V(I,J)=V(IP,J)
 1100   CONTINUE
 1110 CONTINUE

      RETURN
      END
*CMZ :  1.01/51 24/05/96  10.09.46  by  Piero Zucchelli
*CMZ :  1.01/50 04/04/96  11.50.00  by  Piero Zucchelli
*CMZ :  1.01/40 20/11/95  12.56.04  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.23.54  BY  PIERO ZUCCHELLI
*CMZ :  1.01/12 14/05/95  11.19.11  BY  PIERO ZUCCHELLI
*CMZ :  1.01/10 04/05/95  19.29.38  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.58.58  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 15/08/94  07.10.58  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION LUTOGE(KF)
************************************************************************
*                                                                      *
*       RETURNS IN LUTOGE THE GEANT CODE OF A PARTICLE WITH            *
*       LNUD CODE KF                                                   *
*       IPALUP,IPALUM ESTABLISH A CORRESPONDENCE TABLE BETWEEN         *
*       GEANT AND LUND PARTICLE CODES( FOR "STABLE" PARTICLES).        *
*                                                                      *
************************************************************************
      DIMENSION IPALUP(4232),IPALUM(4232)

*      DATA IPALUP/10*0,3,4,6,4,34,4,0,4,0,0,8,1,10,37,35,39,7,17,
*     + 12*0,16,10,2*0,14,13,19,20,21,22,23,9*0,18,41,11*0,24,4048*0/
*      DATA IPALUM/10*0,2,4,5,4,33,4,0,4,0,0,9,1,10,38,36,40,7,17,
*     + 12*0,16,10,2*0,15,25,27,28,29,30,31,9*0,26,12*0,32,4048*0/


      DATA IPALUP/4232*0/
      DATA IPALUM/4232*0/




* LEPTON SECTOR

* gamma
      IPALUP(22)=1
      IPALUM(22)=1
* electron
      IPALUP(11)=3
      IPALUM(11)=2
* neutrinos
      IPALUP(12)=4
      IPALUM(12)=4
      IPALUP(14)=4
      IPALUM(14)=4
      IPALUP(16)=4
      IPALUM(16)=4
* muons
      IPALUP(13)=6
      IPALUM(13)=5
* taus
      IPALUP(15)=34
      IPALUM(15)=33


* MESON SECTOR

* pi0
      IPALUP(111)=7
      IPALUM(111)=7
* piplus/minus
      IPALUP(211)=8
      IPALUM(211)=9
* K0long
      IPALUP(130)=10
      IPALUM(130)=10
* K0short
      IPALUP(310)=16
      IPALUM(310)=16
* K+
      IPALUP(321)=11
      IPALUM(321)=12
* D+
      IPALUP(411)=35
      IPALUM(411)=36
* D0
      IPALUP(421)=37
      IPALUM(421)=38
* D_s+
      IPALUP(431)=39
      IPALUM(431)=40



* BARYONS SECTOR

* neutron
      IPALUP(2112)=13
      IPALUM(2112)=25
* proton
      IPALUP(2212)=14
      IPALUM(2212)=15
* sigma - and +
      IPALUP(3112)=21
      IPALUM(3112)=29
* lambda 0
      IPALUP(3122)=18
      IPALUM(3122)=26
* sigma + and -
      IPALUP(3222)=19
      IPALUM(3222)=27
* Xi-
      IPALUP(3312)=23
      IPALUM(3312)=31
* Xi0
      IPALUP(3322)=22
      IPALUM(3322)=30
* Lambda_c
      IPALUM(4122)=42
      IPALUP(4122)=41
* Xi_c0
      IPALUM(4132)=45
      IPALUP(4132)=46
* Xi_c+
      IPALUM(4232)=44
      IPALUP(4232)=43


      IF (KF.GT.0) THEN
        LUTOGE=IPALUP(KF)
      ELSE
        LUTOGE=IPALUM(-KF)
      ENDIF



      IF (LUTOGE .EQ. 0.AND.ABS(KF).GT.10) THEN
        WRITE(*,*)' +++ LUTOGE: UNKNOWN LUND CODE',KF,' - GEANTINO '
        LUTOGE=48
        CALL LULIST(3)
      ENDIF
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C ********************************************************************

      SUBROUTINE LWBB(ENU)

C...GIVES ENERGY (ENU) OF A (ANTI-)NEUTRINO CHOSEN FROM A SIMPLE
C...PARAMETRIZATION OF A WIDE BAND BEAM.

      DATA EMEAN,SLOPE,EMIN,EMAX/30.,0.02,12.,300./
      A1=1./(EMEAN-12.)
      A2=EXP(EMEAN*SLOPE)
   10 ENU=300.*RLU(0)
      IF(ENU.LT.EMEAN)THEN
        E=A1*(ENU-12.)
      ELSE
        E=A2*EXP(-ENU*SLOPE)
      ENDIF
      IF(ENU.LT.EMIN.OR.ENU.GT.EMAX) GOTO 10
      IF(E.LT.RLU(0)) GOTO 10
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 26/07/94  18.10.50  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 15/07/94  14.06.01  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LWEITS(LFILE)

C...INTEGRATES THE QCD MATRIX ELEMENTS TO OBTAIN PROBABILITIES FOR
C...QG- AND QQ-EVENTS AS A FUNCTION OF (X,W). ALSO FINDS VARIOUS
C...MAXIMUM VALUES TO BE USED FOR THE QCD SIMULATION. RESULTS STORED
C...IN COMMON LGRID AND OPTIONALLY WRITTEN TO LOGICAL FILE LFILE.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LOPTIM/ OPTX(4),OPTY(4),OPTQ2(4),OPTW2(4),COMFAC
      COMMON /LGRID/ NXX,NWW,XX(20),WW(15),PQG(20,15,3),PQQB(20,15,2),
     +QGMAX(20,15,3),QQBMAX(20,15,2),YCUT(20,15),XTOT(20,15),NP
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      DIMENSION WWI(15,4),XXI(20,4)
      EXTERNAL DSIGMA
      DATA WWI/5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,17.5,20.,22.5,25.,
     +5.,7.5,10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,32.5,35.,40.,45.,
     +5.,10.,20.,30.,50.,75.,100.,125.,150.,175.,200.,225.,250.,300.,
     +350.,5.,10.,20.,35.,60.,100.,150.,225.,350.,500.,700.,1000.,
     +1400.,1900.,2500./
      DATA XXI/
     +.001,.002,.004,.006,.008,.01,.02,.04,.06,.08,
     +          .1,.125,.15,.2,.25,.3,.45,.6,.75,.99,
     +.001,.002,.004,.006,.008,.01,.02,.04,.06,.08,
     +          .1,.125,.15,.2,.25,.3,.45,.6,.75,.99,
     +.0001,.0003,.0006,.001,.0025,.0050,.0075,
     +          .01,.02,.04,.06,.08,.1,.125,.15,.2,.3,.5,.75,.99,
     +.0001,.0003,.0006,.001,.0025,.0050,.0075,
     +          .01,.02,.04,.06,.08,.1,.125,.15,.2,.3,.5,.75,.99/
      DATA NCALL/0/

      NCALL=NCALL+1
      LST2=LST(2)
      LST(2)=-3
      WMAX=SQRT(PARL(21))

      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10000)
     +PARL(11),LST(13),MSTU(112),PARU(112), PARL(8),PARL(9),PARL(12),
     +PARL(13)
      IF(LST(17).EQ.0) THEN
        NP=1
        IPMAX=2
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10100)
      ELSE
        NP=3
        IF(LST(23).EQ.1) NP=2
        IPMAX=3
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10200)
      ENDIF

      IF(LST(19).GE.1.AND.LST(19).LE.4) THEN
C...GRID TAKEN FROM DATA IN ARRAYS WWI, XXI.
        DO 10 IW=1,NWW
   10   WW(IW)=WWI(IW,LST(19))
        DO 20 IX=1,NXX
   20   XX(IX)=XXI(IX,LST(19))
      ELSE
C...GRID SPECIFIED BY USER.
	WRITE(6,*)'  Read next nww,nxx '
        READ(5,*) NWW,NXX
        READ(5,*) (WW(IW),IW=1,NWW)
        READ(5,*) (XX(IX),IX=1,NXX)
        IF(XX(NXX).GT..99) XX(NXX)=.99
      ENDIF
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10300)
     +LST(19),NWW,NXX,WW,XX
      IF(WMAX.GT.WW(NWW)) THEN
        IF(LST(3).GE.1) WRITE(6,10400) WMAX,WW(NWW)
        IF(LST(3).GE.2) THEN
          WRITE(6,10700)
          STOP
        ENDIF
      ENDIF
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10500)

      LW=0
      DO 70  IW=1,NWW
        W=WW(IW)
        IF(LW.GT.0) GOTO 80
        IF(W.GT.WMAX) LW=LW+1
        W2=W**2
        LX=0
        DO 60  IX=1,NXX
          X=XX(IX)
          IF(LX.GT.0) GOTO 70
          IF(X.GT.1.-W2/WMAX**2) LX=LX+1
          CALL LEPTO
          PQCOM=PARI(31)*PQ(17)*COMFAC
*     WRITE(*,*)'PQCOM=',PQCOM
          PARL(25)=ULALPS(Q2)
          PARI(20)=PQ(17)
          XTOT(IX,IW)=PQ(17)
          PARL(27)=MAX(PARL(9)**2/W2,PARL(8))
          YCLOW=PARL(27)
          IYCUT=0
   30     IYCUT=IYCUT+1
          RQG=0.
          RQQB=0.
          XPMIN=DBLE(X)/(1.D0-2.D0*(1.D0-DBLE(X))*DBLE(PARL(27)))
          XPMAX=DBLE(X)/(DBLE(X)+(1.D0-DBLE(X))*DBLE(PARL(27)))
          IF(XPMIN.GE.XPMAX) GOTO 50
C...Y_CUT>0.5 CAN GIVE XPMIN<0
          IF(XPMIN.LE.0.) GOTO 50
          DO 40  IP=1,NP
            IF(LST(17).EQ.0) THEN
              PARI(15)=0.
              PARI(16)=0.
              PARI(18)=0.
              PARI(19)=0.
            ELSE
              PARI(14+IP)=0.
              IF(IP.LE.2) PARI(17+IP)=0.
            ENDIF
            LST(20)=IP
            LST(24)=2
            EPS=PARL(11)
            CALL GADAP(XPMIN,XPMAX,DSIGMA,EPS,RESULT)
*     WRITE(*,*)'AFTER GADAP1 IX,DSIGMA',IX,RESULT
            RQG=RQG+RESULT
            PQG(IX,IW,IP)=RESULT/PARL(25)
            IF(LST(17).EQ.0) THEN
              QGMAX(IX,IW,1)=PARI(15)
              QGMAX(IX,IW,2)=PARI(16)
            ELSE
              PQG(IX,IW,IP)=RESULT*PARI(20)/PARI(23+IP)/PARL(25)
              QGMAX(IX,IW,IP)=PARI(14+IP)
            ENDIF
            IF(IP.EQ.3) GOTO 40
            LST(24)=3
            EPS=PARL(11)
            CALL GADAP(XPMIN,XPMAX,DSIGMA,EPS,RESULT)
*     WRITE(*,*)'AFTER GADAP2 IX,DSIGMA',IX,RESULT
            RQQB=RQQB+RESULT
            PQQB(IX,IW,IP)=RESULT/PARL(25)
            IF(LST(17).EQ.0) THEN
              QQBMAX(IX,IW,1)=PARI(18)
              QQBMAX(IX,IW,2)=PARI(19)
            ELSE
              PQQB(IX,IW,IP)=RESULT*PARI(20)/PARI(23+IP)/PARL(25)
              QGMAX(IX,IW,IP)=PARI(17+IP)
            ENDIF
   40     CONTINUE
   50     CONTINUE
          RQ=1.-RQG-RQQB
*        WRITE(*,*)'RQG,RQQB',RQG,RQBB

          IF(RQ.LT.0.) THEN
C...QCD PROBABILITIES > 1, INCREASE CUTOFF.
            YCLOW=PARL(27)
            POT=SQRT(1./(RQG+RQQB))
            PARL(27)=(1./PARL(12)+0.01)*(PARL(12)*PARL(27))**POT
*        WRITE(*,*)'RQ<=',RQ
            GOTO 30
          ELSEIF(IYCUT.GT.1.AND.RQ.GT.PARL(13)) THEN
C...CUTOFF INCREASED TOO MUCH, TRY LOWER.
            PARL(27)=(PARL(27)+YCLOW)/2.
*        WRITE(*,*)'RQ>=',RQ
            GOTO 30
          ENDIF
          YCUT(IX,IW)=PARL(27)
          IF(LST(39).EQ.-91) THEN
C...INCLUDE 3-JET CROSS SECTION IN DENOMINATOR
            QTOT=1.+RQG+RQQB
            RQG =RQG/QTOT
            RQQB=RQQB/QTOT
            RQ=1.-RQG-RQQB
          ENDIF
          IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,
     +    10600) W,X,Y,Q2,PARL(25),PQCOM,PARL(27),IYCUT, RQ,RQG,RQQB,
     +    (QGMAX(IX,IW,IP),IP=1,IPMAX), (QQBMAX(IX,IW,IP),IP=1,MIN(2,
     +    IPMAX))
   60   CONTINUE
   70 CONTINUE
   80 CONTINUE

      LST(2)=LST2
      IF(LFILE.LT.0) THEN
C...WRITE RESULTS ON LOGICAL FILE NUMBER IABS(LFILE)
        WRITE(IABS(LFILE)) LST,PARL,NXX,NWW,NP,XX,WW
        WRITE(IABS(LFILE))(((PQG(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,NP),
     +  (((PQQB(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,NP),
     +  (((QGMAX(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,IPMAX),
     +  (((QQBMAX(IX,IW,IP),IX=1,NXX),IW=1,NWW),IP=1,MIN(2,IPMAX)),
     +  YCUT
        IF(NP.NE.1) WRITE(IABS(LFILE)) XTOT
        CLOSE(IABS(LFILE))
      ENDIF
      RETURN

10000 FORMAT('1',/,5X,'INTEGRATION OF 1ST ORDER QCD MATRIX ELEMENTS',
     +           /,5X,'============================================',
     +/,' FOR GLUON RADIATION (QG-EVENT) AND BOSON-GLUON FUSION ',
     +'(QQ-EVENT) PROBABILITY.',
     +//,' REQUIRED PRECISION IN INTEGRATION, PARL(11) =',F8.4,
     +//,' HEAVIEST FLAVOUR PRODUCED IN BOSON-GLUON FUSION, LST(13) =',
     +I5,//,' ALPHA-STRONG PARAMETERS: # FLAVOURS, MSTU(112) =',I3,
     +/,25X,' QCD LAMBDA, PARU(112) =',F6.3,' GEV',
     +//,' CUTS ON MATRIX ELEMENTS:',
     +/,' PARL(8), PARL(9), PARL(12), PARL(13) =',4F8.4,/)
10100 FORMAT(' LEPTON ENERGY NOT ALLOWED TO VARY IN SIMULATION.',/)
10200 FORMAT(' LEPTON ENERGY ALLOWED TO VARY IN SIMULATION, ',/,
     +' Y IN TABLE BELOW CALCULATED ASSUMING MAX ENERGY.',/)
10300 FORMAT(' GRID CHOICE, LST(19) =',I3,5X,'# GRID POINTS IN W, X =',
     +2I5,/,' W-VALUES IN ARRAY WW:',/,10F8.1,/,5F8.1,
     +/,' X-VALUES IN ARRAY XX:',/,10F8.4,/,10F8.4,/)
10400 FORMAT(' WARNING: MAX W OUTSIDE GRID, WMAX, GRID-MAX =',2F12.1)
10500 FORMAT(//,6X,'W',7X,'X',7X,'Y',6X,'Q**2',1X,'ALPHA',1X,'DSIGMA',
     +9X,'CUT',' IT',2X,'Q-EVENT',1X,'QG-EVENT',
     +1X,'QQ-EVENT',' MAX OF MATRIX ELEMENTS QG & QQ; L,R OR T,S,I',
     +/,1X,132(1H-),/)
10600 FORMAT(F7.1,2F8.4,1PG10.3,0PF6.2,1PG11.3,0PF8.4,I3,3F9.4,1P,5E9.2)
10700 FORMAT(' EXECUTION STOPPED ',/)
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LXP(XP,IFAIL)

C...CHOOSE VALUE OF XP ACCORDING TO QCD MATRIX ELEMENTS WEIGHTED BY
C...STRUCTURE FUNCTIONS.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      DOUBLE PRECISION DXPMAX

      IFAIL=1
      XPMIN=DBLE(X)/(1.D0-2.D0*(1.D0-DBLE(X))*DBLE(PARL(27)))
      DXPMAX=DBLE(X)/(DBLE(X)+(1.D0-DBLE(X))*DBLE(PARL(27)))
      XPMAX=SNGL(DXPMAX)
      IF(XPMIN.GE.XPMAX) RETURN
      AP=1.-XPMIN
      BP=(1.D0-DXPMAX)/AP
      IF(LST(24).EQ.2) THEN
        QXPMAX=PARI(15)
        IF(LST(17).NE.0) QXPMAX=PARI(24)*PARI(15)+PARI(25)*PARI(16)+
     +  PARI(26)*PARI(17)
      ELSE
        QXPMAX=PARI(18)
        IF(LST(17).NE.0) QXPMAX=PARI(24)*PARI(18)+PARI(25)*PARI(19)
      ENDIF
C...SAFETY FACTOR ON MAX VALUE.
      QXPMAX=QXPMAX*1.05
      LOOP=0
   10 LOOP=LOOP+1
      IF(LOOP.GT.1000) RETURN
      XP=1.-AP*BP**RLU(0)
      XPWEIT=DSIGMA(XP)/QXPMAX
      IF(XPWEIT.LT.RLU(0)) GOTO 10
      IFAIL=0
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.25  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 20/09/94  14.43.37  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LXSECT

C...INTEGRATE DIFFERENTIAL CROSS-SECTION USING GADAP, RIWIAD OR DIVONNE

      COMMON /LINTEG/ NTOT,NPASS
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      DOUBLE PRECISION ACC,VALUE,ERRIW,FLOW,FHIGH
      COMMON /PARAMS/ ACC,NDIM,NSUB,ITER
      COMMON /ANSWER/ VALUE,ERRIW
      COMMON /BNDLMT/ FLOW,FHIGH
      COMMON /SAMPLE/ NPOINT
      DIMENSION XMINUS(2),XPLUS(2)
      EXTERNAL DCROSS,DLOWER,DUPPER,RIWFUN
      COMMON/LINPATCH/NCALLS,NCALL

      NCALL=NCALL+1
      CALL LTIMEX(TI1)
      NTOT=0
      NPASS=0
      SIGMA=0.
      ERREST=0.
      NDIM=2
C...PARAMETERS FOR RIWIAD INTEGRATION.
      ACC=PARL(15)
      NSUB=100
      ITER=100
C...PARAMETERS FOR DIVON INTEGRATION.
      DO 10 I=1,2
        XMINUS(I)=0.
   10 XPLUS(I)=1.
      EPS=PARL(15)
      MAXNUM=50000
      FLOW=-1.D0
      FHIGH=1.D+20
      NPOINT=100
C...ADDITIONAL PARAMETERS FOR DETAILED DIVON INTEGRATION.
      SPRDMX=2.
      MAXPTS=50000
      JDEG=0
      NPT=1000

      SIGMA=0.
      ERREST=0.
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10000)
      IF(LST(10).EQ.1) THEN
C...INTEGRATION USING GADAP.
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10100)
        ACCUR=PARL(15)
        IT=0
   20   IT=IT+1
        ERREST=ACCUR
        CALL GADAP2(XMIN,XMAX,DLOWER,DUPPER,DCROSS,ERREST,SIGMA)
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10200)
     +  IT,NTOT,NPASS,SIGMA
        IF(SIGMA.GT.1.) THEN
          IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,
     +    10300) ACCUR
        ELSE
          IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,
     +    10400) ACCUR,ACCUR/MAX(1.E-22,SIGMA),PARL(15)
          ACCUR=MAX(1.E-22,SIGMA*PARL(15))
          IF(IT.LT.2) GOTO 20
        ENDIF
      ELSEIF(LST(10).EQ.2) THEN
C...INTEGRATION USING RIWIAD. WHEN RIWIAD CANNOT BE LOADED:
C...ACTIVATE NEXT TWO LINES AND DEACTIVATE RIWIAD CALL.
C       WRITE(6,*) ' RIWIAD NOT AVAILABLE, EXECUTION STOPPED.'
C       STOP
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10500)
     +  SNGL(ACC),NSUB,ITER
        CALL RIWIAD(RIWFUN)
        SIGMA=SNGL(VALUE)
        ERREST=SNGL(ERRIW)
      ELSEIF(LST(10).EQ.3) THEN
C...INTEGRATION USING SIMPLE DIVONNE. WHEN DIVONNE CANNOT BE LOADED:
C...ACTIVATE NEXT TWO LINES AND DEACTIVATE DIVONNE CALL.
C       WRITE(6,*) ' DIVONNE NOT AVAILABLE, EXECUTION STOPPED.'
C       STOP
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10600)
     +  EPS,MAXNUM,SNGL(FLOW),SNGL(FHIGH),NPOINT
        CALL DIVON(NDIM,XMINUS,XPLUS,EPS,MAXNUM,SIGMA,ERREST)
      ELSEIF(LST(10).EQ.4) THEN
C...INTEGRATION USING DETAILED DIVONNE. WHEN DIVONNE CANNOT BE LOADED:
C...ACTIVATE NEXT TWO LINES AND DEACTIVATE PARTN AND INTGRL CALLS.
C       WRITE(6,*) ' DIVONNE NOT AVAILABLE, EXECUTION STOPPED.'
C       STOP
        IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)) WRITE(6,10700)
     +  EPS,MAXNUM, SNGL(FLOW),SNGL(FHIGH),NPOINT,SPRDMX,MAXPTS,JDEG,
     +  NPT
        CALL PARTN(NDIM,XMINUS,XPLUS,SPRDMX,MAXPTS)
        CALL INTGRL(NDIM,JDEG,NPT,SIGMA,ERREST)
      ELSE
        IF(LST(3).GE.1) WRITE(6,*) ' WARNING: LST(10) = ',LST(10),
     +  ' NOT ALLOWED.'
      ENDIF
      CALL LTIMEX(TI2)
      IF(LST(3).GE.4.OR.(LST(3).EQ.3.AND.NCALL.EQ.1)
     +.OR.(LST(3).GE.1.AND.NPASS.EQ.0)) THEN
        WRITE(6,10800) SIGMA,ERREST,NTOT,NPASS,TI2-TI1
        IF(LST(3).GE.1.AND.NPASS.EQ.0) WRITE(6,10900)
      ENDIF
      PARL(23)=SIGMA

      RETURN
10000 FORMAT(/,' INTEGRATION OF CROSS SECTION:',/,1X,28('-'))
10100 FORMAT(5X,'USING GADAP = ADAPTIVE GAUSSIAN INTEGRATION')
10200 FORMAT(5X,'ITERATION #',I3,/,
     +10X,'# FUNCTION EVALUATIONS; TOTAL & NON-ZERO =',2I8,/,
     +10X,'SIGMA =',G10.2,' PB')
10300 FORMAT(10X,'REQUIRED  RELATIVE ERROR = ',G10.2)
10400 FORMAT(10X,'EFFECTIVE ABSOLUTE ERROR = ',G10.2,/,
     +       10X,'EFFECTIVE RELATIVE ERROR = ',G10.2,/,
     +       10X,'REQUIRED  RELATIVE ERROR = ',G10.2)
10500 FORMAT(5X,'USING RIWIAD WITH PARAMETERS: REL. ACC. = ',F10.4,
     +/,5X,'# OF SUBVOLUMES = ',I5,5X,'MAX # ITERATIONS = ',I5)
10600 FORMAT(5X,'USING AUTOMATIC DIVONNE WITH PARAMETERS: ',
     +'REL. ACC. = ',F10.4,/,5X,'MAX # FUNCTION CALLS = ',I5,
     +/,5X,'LOWER AND UPPER BOUND ON INTEGRAND =',2E12.4,
     +/,5X,'# SAMPLE POINTS/REGION =',I5)
10700 FORMAT(5X,'USING DETAILED DIVONNE WITH PARAMETERS: ',
     +'REL. ACC. = ',F10.4,/,5X,'MAX # FUNCTION CALLS = ',I5,
     +/,5X,'LOWER AND UPPER BOUND ON INTEGRAND =',2E12.4,
     +/,5X,'# SAMPLE POINTS/REGION =',I5,
     +/,5X,'SPRDMX, MAXPTS, JDEG, NPT =',F5.2,3I10)
10800 FORMAT(/,' ===> CROSS-SECTION =',1P,G12.3,
     +' PB,  ERROR ESTIMATE = ',G12.3,/,
     +6X,'# OF INTEGRAND EVALUATIONS; TOTAL & NON-ZERO =',2I8,/,
     +6X,'CPU TIME FOR INTEGRATION =',G12.3,' SECONDS',/)
10900 FORMAT(' WARNING: INTEGRAND ALWAYS ZERO, PROBABLY NO ALLOWED',
     +' PHASE SPACE DUE TO CUTS',/,
     +10X,'CHECK, IN PARTICULAR, CUT(11) TO CUT(14)')
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE LZP(XP,ZP,IFAIL)

C...CHOOSE VALUE OF ZP ACCORDING TO QCD MATRIX ELEMENTS.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      DATA C1,C2/0.2122066,0.0795775/,DZPMAX,SZP,CP/3*0./
      FQG(DZ,DX,DA,DB,DC)=DA*(DZ**2+DX**2)/(1.-DX)+2.*DA*DX*DZ*(1.-DZ)
     &+2.*DA*(1.-DZ)+4.*DB*DX*DZ*(1.-DZ)+DC*(DZ**2+DX**2)/(1.-DX)+
     &2.*DC*(DX+DZ)*(1.-DZ)
      FQQ(DZ,DX,DA,DB,DC,DD,DE)=DA*DD*(DZ**2+(1.-DZ)**2)+DB*DE*DZ*
     &(1.-DZ)+DC*DD*(2.*DZ-1.)

      IFAIL=1
      IH=1
      IF(LST(30).EQ.1) IH=2
      ZPMIN=(1.-X)*XP/(XP-X)*PARL(27)
      IF(ZPMIN.GE.0.5) RETURN
      ZPMAX=1.-ZPMIN
      I=IABS(LST(25))
      AP=1.-ZPMIN
      BP=ZPMIN/AP
      IF(LST(23).EQ.2) THEN
        A=PARI(24)
        B=PARI(25)
        CSIGN=-LST(30)*ISIGN(1,LST(25))*PARI(26)
      ELSE
        A=(EWQC(1,IH,I)+EWQC(2,IH,I))*PARI(24)
        B=(EWQC(1,IH,I)+EWQC(2,IH,I))*PARI(25)
        C=(EWQC(1,IH,I)-EWQC(2,IH,I))*PARI(26)
        CSIGN=-C*LST(30)*ISIGN(1,LST(25))
      ENDIF
      IF(LST(24).EQ.2) THEN
        DZPMAX=MAX(FQG(ZPMIN,XP,A,B,CSIGN),FQG(ZPMAX,XP,A,B,CSIGN))
        AA=2.*(A+CSIGN)/(1.-XP)-4.*A*XP-8.*B*XP-4.*CSIGN
        IF(ABS(AA).GT.1.E-20) THEN
          BB=2.*A*(XP-1.)+4.*B*XP+2.*CSIGN*(1.-XP)
          Z1=-BB/AA
          IF(Z1.GT.ZPMIN.AND.Z1.LT.ZPMAX) THEN
            DZPMAX=MAX(DZPMAX,FQG(Z1,XP,A,B,CSIGN))
          ENDIF
        ENDIF
        DZPMAX=DZPMAX*C1*1.05
      ELSEIF(LST(24).EQ.3) THEN
        CP=1./BP**2
        D=XP**2+(1.-XP)**2
        E=8.*XP*(1-XP)
        DZPMAX=MAX(FQQ(ZPMIN,XP,A,B,CSIGN,D,E),
     &  FQQ(ZPMAX,XP,A,B,CSIGN,D,E))
        AA=4.*A*D-2.*B*E
        IF(ABS(AA).GT.1.E-20) THEN
          BB=B*E-2.*A*D+2.*CSIGN*D
          Z1=-BB/AA
          IF(Z1.GT.ZPMIN.AND.Z1.LT.ZPMAX) THEN
            DZPMAX=MAX(DZPMAX,FQQ(Z1,XP,A,B,CSIGN,D,E))
          ENDIF
        ENDIF
        DZPMAX=DZPMAX*C2*1.05
      ENDIF
      IPART=LST(24)-1
      LOOP=0
   10 LOOP=LOOP+1
      IF(LOOP.GT.1000) RETURN
      IF(LST(24).EQ.2) THEN
        ZP=1.-AP*BP**RLU(0)
        SZP=1.-ZP
      ELSEIF(LST(24).EQ.3) THEN
        DP=BP*CP**RLU(0)
        ZP=DP/(1.+DP)
        SZP=ZP*(1.-ZP)
      ENDIF
      ZPWEIT=SZP*(A*DQCD(0,IPART,1,XP,ZP,0.)+B*DQCD(0,IPART,2,XP,ZP,0.)
     &+CSIGN*DQCD(0,IPART,3,XP,ZP,0.))/DZPMAX
      IF(ZPWEIT.LT.RLU(0)) GOTO 10
      IFAIL=0
      RETURN
      END
*CMZ :  1.01/45 08/01/96  11.37.02  by  Piero Zucchelli
*CMZ :  1.01/44 05/01/96  18.04.59  by  Piero Zucchelli
*CMZ :  1.01/40 10/11/95  19.03.13  by  Piero Zucchelli
*-- Author :
      SUBROUTINE MZINI
*-----------------------------------------------------*
*                                                     *
*               INITIALIZE BANKS                      *
*                                                     *
*-----------------------------------------------------*
*KEEP,zebra.

      PARAMETER (NNQ=1000000)
*
      DIMENSION LQ(NNQ),IQ(NNQ),Q(NNQ)
      EQUIVALENCE (Q(1),IQ(1),LQ(9),JSTRUC(8))
      COMMON /QUEST/IQUEST(100)
      COMMON /XQSTOR/IXEVT,IFENCE(16),JGEEV,JSTRUC(99),JREFER(100),
     +DIV12(NNQ)
      COMMON /FZLUN/LUNFZ
      COMMON/MZIOALL/IOGENF

*KEEP,info.
         COMMON/INFONEW/IRDATE,IRTIME
*KEND.
*
*---   INITIALISATION OF ZEBRA
*
      CALL MZFORM('GENF','5I 2H',IOGENF)
      RETURN
      END
*CMZ :  1.01/51 23/05/96  16.08.27  by  Piero Zucchelli
*CMZ :  1.01/38 18/10/95  18.27.55  by  Piero Zucchelli
*CMZ :  1.01/37 04/09/95  15.00.03  BY  PIERO ZUCCHELLI
*-- AUTHOR :    PIERO ZUCCHELLI   04/09/95


      SUBROUTINE ORTH(PO,P,PB)

* ASSUMPTION: BEAM ALONG Z!

      REAL*4 BDIR(3),PB2(3),P(3),PB(3),PV(3)


      BDIR(1)=0.
      BDIR(2)=0.
      BDIR(3)=1.

* ORTHOGONALIZE BEAM DIRECTION AND BASE DIRECTION


      BDIRL=SQRT(BDIR(1)**2+BDIR(2)**2+BDIR(3)**2)


* PB DOT BEAM DIRECTION
      PBDB=0

      DO I=1,3
        BDIR(I)=BDIR(I)/BDIRL
        PBDB=PBDB+PB(I)*BDIR(I)
      ENDDO

      DO I=1,3
        PB2(I)=PB(I) - PBDB*BDIR(I)
      ENDDO

* PB2 INOW IS ORTHOGONAL TO BEAM DIRECTION: NORMALIZE IT
      PB2L=SQRT(PB2(1)**2+PB2(2)**2+PB2(3)**2)
      IF (PB2L.EQ.0) THEN
        PO=0
        RETURN
      ENDIF

      DO I=1,3
        PB2(I)=PB2(I)/PB2L
      ENDDO


* CALCULATE  ALPHAS

      A1=0
      A2=0

      DO I=1,3
        A1=A1+P(I)*BDIR(I)
        A2=A2+P(I)*PB2(I)
      ENDDO


      PO=0
      DO I=1,3
        PV(I)=P(I)-A1*BDIR(I)-A2*PB2(I)
        PO=PO+PV(I)**2
      ENDDO

      IF (PO.GT.0) THEN
        PO=SQRT(PO)
      ELSE
        PO=0
      ENDIF

      RETURN
      END
*CMZ :  1.01/08 05/03/95  11.35.13  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 23/09/94  12.17.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :    PIERO ZUCCHELLI   23/09/94

      SUBROUTINE PARUPD

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     +Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)

*KEEP,JETTA.
C--
         PARAMETER (ICENTO=100)
         PARAMETER (ISIZ=93)
         PARAMETER (IOF1=32)
         PARAMETER (IOF2=83)
         PARAMETER (LUX_LEVEL=4)
         INTEGER*4 JTAU(100),JPRI(100),JSTRO(100)
         REAL*4 FTUPLE(ISIZ)
         COMMON/JETTAGL/JTAU,JPRI,JSTRO
         COMMON/NTUPLA/FTUPLE,ISFIRST
         COMMON/BEAM/SPEC(ICENTO)
         COMMON /MAXSPEC/RMAXSPEC,RINTSPEC
         COMMON/SAV/XMINSAV(ICENTO),XMAXSAV(ICENTO),YMINSAV(ICENTO),
     &   YMAXSAV(ICENTO),Q2MINSAV(ICENTO),Q2MAXSAV(ICENTO),
     &   W2MINSAV(ICENTO),W2MAXSAV(ICENTO),PARIMAX(ICENTO),
     &   PPSAVE(ICENTO,3,4,5),PARICOR(ICENTO),INDEX,SIGMASAV(ICENTO),
     &   XMSIGMA,XSECT

*KEND.


      IF (PARI(32).NE.PARICOR(INDEX)) THEN
        PARICOR(INDEX)=PARI(32)
        PARIMAX(INDEX)=PARI(LST(23))
      ENDIF

      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION PHINT(IDUM)
C.----------------------------------------------------------------------
C.
C.    PHINT:   PHOTOS INTERFERENCE
C.
C.    PURPOSE:  CALCULATES INTERFERENCE BETWEEN EMISSION OF PHOTONS FROM
C.              DIFFERENT POSSIBLE CHAGED DAUGHTERS STORED IN
C.              THE  HEP COMMON /PHOEVT/.
C.
C.    INPUT PARAMETER:    COMMONS /PHOEVT/ /PHOMOM/ /PHOPHS/
C.
C.
C.    OUTPUT PARAMETERS:
C.
C.
C.    AUTHOR(S):  Z. WAS,                         CREATED AT:  10/08/93
C.                                                LAST UPDATE:
C.
C.----------------------------------------------------------------------

C--   IMPLICIT NONE
      REAL PHINT
      REAL PHOCHA
      INTEGER IDUM
      INTEGER NMXPHO
      PARAMETER (NMXPHO=2000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     +JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL PNEUTR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
      REAL MPASQR,XX,BETA
      LOGICAL IFINT
      INTEGER K,IDENT
C
      DO  K=JDAPHO(2,1),JDAPHO(1,1),-1
        IF(IDPHO(K).NE.22) THEN
          IDENT=K
          GOTO 10
        ENDIF
      ENDDO
   10 CONTINUE
C CHECK IF THERE IS A PHOTON
      IFINT= NPHO.GT.IDENT
C CHECK IF IT IS TWO BODY + GAMMAS REACTION
      IFINT= IFINT.AND.(IDENT-JDAPHO(1,1)).EQ.1
C CHECK IF TWO BODY WAS PARTICLE ANTIPARTICLE
      IFINT= IFINT.AND.IDPHO(JDAPHO(1,1)).EQ.-IDPHO(IDENT)
C CHECK IF PARTICLES WERE CHARGED
      IFINT= IFINT.AND.PHOCHA(IDENT).NE.0
C CALCULATES INTERFERENCE WEIGHT CONTRIBUTION
      IF(IFINT) THEN
        MPASQR = PPHO(5,1)**2
        XX=4.*MCHSQR/MPASQR*(1.-XPHOTO)/(1.-XPHOTO+(MCHSQR-MNESQR)/
     +     MPASQR)**2
        BETA=SQRT(1.-XX)
        PHINT = 2D0/(1D0+COSTHG**2*BETA**2)
      ELSE
        PHINT = 1D0
      ENDIF
      END
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHLUPA(IPOINT)
C.----------------------------------------------------------------------
C.
C.    PHLUPA:   DEBUGGING TOOL
C.
C.    PURPOSE:  NONE, EVENTUALLY MAY PRINTOUT CONTENT OF THE
C.              /PHOEVT/ COMMON
C.
C.    INPUT PARAMETERS:   COMMON /PHOEVT/ AND /PHNUM/
C.                        LATTER MAY HAVE NUMBER OF THE EVENT.
C.
C.    OUTPUT PARAMETERS:  NONE
C.
C.    AUTHOR(S):  Z. WAS                          CREATED AT:  30/05/93
C.                                                LAST UPDATE: 10/08/93
C.
C.----------------------------------------------------------------------
      INTEGER NMXPHO
      PARAMETER (NMXPHO=2000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     +JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      COMMON /PHNUM/ IEV
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      DIMENSION SUM(5)
      IF (IPOINT.LT.3000) RETURN
      IOUT=56
      IF (IEV.LT.1000) THEN
        DO I=1,5
          SUM(I)=0.0
        ENDDO
        WRITE(PHLUN,*) 'EVENT NR=',IEV, 'WE ARE TESTING /PHOEVT/ AT '
     +  //'IPOINT=',IPOINT
        WRITE(PHLUN,10000)
        I=1
        WRITE(PHLUN,10100) IDPHO(I),PPHO(1,I),PPHO(2,I),PPHO(3,I),
     +  PPHO(4, I),PPHO(5,I),JDAPHO(1,I),JDAPHO(2,I)
        I=2
        WRITE(PHLUN,10100) IDPHO(I),PPHO(1,I),PPHO(2,I),PPHO(3,I),
     +  PPHO(4, I),PPHO(5,I),JDAPHO(1,I),JDAPHO(2,I)
        WRITE(PHLUN,*) ' '
        DO I=3,NPHO
          WRITE(PHLUN,10100) IDPHO(I),PPHO(1,I),PPHO(2,I),PPHO(3,I),
     +    PPHO(4,I),PPHO(5,I),JMOPHO(1,I),JMOPHO(2,I)
          DO J=1,4
            SUM(J)=SUM(J)+PPHO(J,I)
          ENDDO
        ENDDO
        SUM(5)=SQRT(ABS(SUM(4)**2-SUM(1)**2-SUM(2)**2-SUM(3)**2))
        WRITE(PHLUN,10200) SUM
10000 FORMAT(1X,'  ID      ','P_X      ','P_Y      ','P_Z      ',
     +                   'E        ','M        ',
     +                   'ID-MO_DA1','ID-MO DA2' )
10100 FORMAT(1X,I4,5(F9.3),2I9)
10200 FORMAT(1X,' SUM',5(F9.3))
      ENDIF
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION PHOAN1(X,Y)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS CALCULATION OF ANGLE '1'
C.
C.    PURPOSE:  CALCULATE ANGLE FROM X AND Y
C.
C.    INPUT PARAMETERS:  X, Y
C.
C.    OUTPUT PARAMETER:  FUNCTION VALUE
C.
C.    AUTHOR(S):  S. JADACH                       CREATED AT:  01/01/89
C.                B. VAN EIJK                     LAST UPDATE: 02/01/90
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION PHOAN1
      REAL X,Y
      REAL PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      IF (ABS(Y).LT.ABS(X)) THEN
        PHOAN1=ATAN(ABS(Y/X))
        IF (X.LE.0.) PHOAN1=PI-PHOAN1
      ELSE
        PHOAN1=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      IF (Y.LT.0.) PHOAN1=TWOPI-PHOAN1
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION PHOAN2(X,Y)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS CALCULATION OF ANGLE '2'
C.
C.    PURPOSE:  CALCULATE ANGLE FROM X AND Y
C.
C.    INPUT PARAMETERS:  X, Y
C.
C.    OUTPUT PARAMETER:  FUNCTION VALUE
C.
C.    AUTHOR(S):  S. JADACH                       CREATED AT:  01/01/89
C.                B. VAN EIJK                     LAST UPDATE: 02/01/90
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION PHOAN2
      REAL X,Y
      REAL PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      IF (ABS(Y).LT.ABS(X)) THEN
        PHOAN2=ATAN(ABS(Y/X))
        IF (X.LE.0.) PHOAN2=PI-PHOAN2
      ELSE
        PHOAN2=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOBO3(ANGLE,PVEC)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS BOOST ROUTINE '3'
C.
C.    PURPOSE:  BOOST  VECTOR PVEC  ALONG Z-AXIS WHERE ANGLE = EXP(ETA),
C.              ETA IS THE HYPERBOLIC VELOCITY.
C.
C.    INPUT PARAMETERS:  ANGLE, PVEC
C.
C.    OUTPUT PARAMETER:  PVEC
C.
C.    AUTHOR(S):  S. JADACH                       CREATED AT:  01/01/89
C.                B. VAN EIJK                     LAST UPDATE: 02/01/90
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION QPL,QMI,ANGLE
      REAL PVEC(4)
      QPL=(PVEC(4)+PVEC(3))*ANGLE
      QMI=(PVEC(4)-PVEC(3))/ANGLE
      PVEC(3)=(QPL-QMI)/2.
      PVEC(4)=(QPL+QMI)/2.
      RETURN
      END
*CMZ :  1.01/50 23/05/96  10.22.20  by  Piero Zucchelli
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOBOS(IP,PBOOS1,PBOOS2,FIRST,LAST)
C.----------------------------------------------------------------------
C.
C.    PHOBOS:   PHOTON RADIATION IN DECAYS BOOST ROUTINE
C.
C.    PURPOSE:  BOOST PARTICLES  IN  CASCADE DECAY  TO PARENT REST FRAME
C.              AND BOOST BACK WITH MODIFIED BOOST VECTOR.
C.
C.    INPUT PARAMETERS:       IP:  POINTER OF PARTICLE STARTING CHAIN
C.                                 TO BE BOOSTED
C.                        PBOOS1:  BOOST VECTOR TO REST FRAME,
C.                        PBOOS2:  BOOST VECTOR TO MODIFIED FRAME,
C.                        FIRST:   POINTER TO FIRST PARTICLE TO BE BOOS-
C.                                 TED (/HEPEVT/),
C.                        LAST:    POINTER TO LAST  PARTICLE TO BE BOOS-
C.                                 TED (/HEPEVT/).
C.
C.    OUTPUT PARAMETERS:  COMMON /HEPEVT/.
C.
C.    AUTHOR(S):  B. VAN EIJK                     CREATED AT:  13/02/90
C.                Z. WAS                          LAST UPDATE: 16/11/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION BET1(3),BET2(3),GAM1,GAM2,PB,DATA
      INTEGER I,J,FIRST,LAST,MAXSTA,NSTACK,IP
      PARAMETER (MAXSTA=2000)
      INTEGER STACK(MAXSTA)
      REAL PBOOS1(5),PBOOS2(5)
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      DOUBLE PRECISION PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     +JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      IF ((LAST.EQ.0).OR.(LAST.LT.FIRST)) RETURN
      NSTACK=0
      DO 10 J=1,3
        BET1(J)=-PBOOS1(J)/PBOOS1(5)
   10 BET2(J)=PBOOS2(J)/PBOOS2(5)
      GAM1=PBOOS1(4)/PBOOS1(5)
      GAM2=PBOOS2(4)/PBOOS2(5)
C--
C--   BOOST VECTOR TO PARENT REST FRAME...
   20 DO 50 I=FIRST,LAST
        PB=BET1(1)*PHEP(1,I)+BET1(2)*PHEP(2,I)+BET1(3)*PHEP(3,I)
        IF (JMOHEP(1,I).EQ.IP) THEN
          DO 30 J=1,3
   30     PHEP(J,I)=PHEP(J,I)+BET1(J)*(PHEP(4,I)+PB/(GAM1+1.))
          PHEP(4,I)=GAM1*PHEP(4,I)+PB
C--
C--    ...AND BOOST BACK TO MODIFIED PARENT FRAME.
          PB=BET2(1)*PHEP(1,I)+BET2(2)*PHEP(2,I)+BET2(3)*PHEP(3,I)
          DO 40 J=1,3
   40     PHEP(J,I)=PHEP(J,I)+BET2(J)*(PHEP(4,I)+PB/(GAM2+1.))
          PHEP(4,I)=GAM2*PHEP(4,I)+PB
          IF (JDAHEP(1,I).NE.0) THEN
            NSTACK=NSTACK+1
C--
C--    CHECK ON STACK LENGTH...
            IF (NSTACK.GT.MAXSTA) THEN
              DATA=NSTACK
              CALL PHOERR(7,'PHOBOS',DATA)
            ENDIF
            STACK(NSTACK)=I
          ENDIF
        ENDIF
   50 CONTINUE
      IF (NSTACK.NE.0) THEN
C--
C--   NOW GO ONE STEP FURTHER IN THE DECAY TREE...
        FIRST=JDAHEP(1,STACK(NSTACK))
        LAST=JDAHEP(2,STACK(NSTACK))
        IP=STACK(NSTACK)
        NSTACK=NSTACK-1
        GOTO 20
      ENDIF
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION PHOCHA(IDHEP)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS CHARGE DETERMINATION
C.
C.    PURPOSE:  CALCULATE THE CHARGE  OF PARTICLE  WITH CODE IDHEP.  THE
C.              CODE  OF THE  PARTICLE  IS  DEFINED BY THE PARTICLE DATA
C.              GROUP IN PHYS. LETT. B204 (1988) 1.
C.
C.    INPUT PARAMETER:   IDHEP
C.
C.    OUTPUT PARAMETER:  FUNTION VALUE = CHARGE  OF  PARTICLE  WITH CODE
C.                       IDHEP
C.
C.    AUTHOR(S):  E. BARBERIO AND B. VAN EIJK     CREATED AT:  29/11/89
C.                                                LAST UPDATE: 02/01/90
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      REAL PHOCHA
      INTEGER IDHEP,IDABS,Q1,Q2,Q3
C--
C--   ARRAY 'CHARGE' CONTAINS THE CHARGE  OF THE FIRST 101 PARTICLES AC-
C--   CORDING  TO  THE PDG PARTICLE CODE... (0 IS ADDED FOR CONVENIENCE)
      REAL CHARGE(0:100)
      DATA CHARGE/ 0.,
     &-0.3333333333,  0.6666666667, -0.3333333333, 0.6666666667,
     &-0.3333333333,  0.6666666667, -0.3333333333, 0.6666666667,
     & 2*0., -1., 0., -1., 0., -1., 0., -1., 6*0., 1., 12*0., 1., 63*0./
      IDABS=ABS(IDHEP)
      IF (IDABS.LE.100) THEN
C--
C--   CHARGE OF QUARK, LEPTON, BOSON ETC....
        PHOCHA = CHARGE(IDABS)
      ELSE
C--
C--   CHECK ON PARTICLE BUILD OUT OF QUARKS, UNPACK ITS CODE...
        Q3=MOD(IDABS/1000,10)
        Q2=MOD(IDABS/100,10)
        Q1=MOD(IDABS/10,10)
        IF (Q3.EQ.0) THEN
C--
C--   ...MESON...
          IF(MOD(Q2,2).EQ.0) THEN
            PHOCHA=CHARGE(Q2)-CHARGE(Q1)
          ELSE
            PHOCHA=CHARGE(Q1)-CHARGE(Q2)
          ENDIF
        ELSE
C--
C--   ...DIQUARKS OR BARYON.
          PHOCHA=CHARGE(Q1)+CHARGE(Q2)+CHARGE(Q3)
        ENDIF
      ENDIF
C--
C--   FIND THE SIGN OF THE CHARGE...
      IF (IDHEP.LT.0.) PHOCHA=-PHOCHA
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOCHK(JFIRST)
C.----------------------------------------------------------------------
C.
C.    PHOCHK:   CHECKING BRANCH.
C.
C.    PURPOSE:  CHECKS WHETHER PARTICLES IN THE COMMON BLOCK /PHOEVT/
C.              CAN BE SERVED BY PHOMAK.
C.              JFIRST IS THE POSITION IN /HEPEVT/ (!) OF THE FIRST DAUGHTER
C.              OF SUB-BRANCH UNDER ACTION.
C.
C.
C.    AUTHOR(S):  Z. WAS                           CREATED AT: 22/10/92
C.                                                LAST UPDATE: 16/10/93
C.
C.----------------------------------------------------------------------
C     ********************
C--   IMPLICIT NONE
      INTEGER NMXPHO
      PARAMETER (NMXPHO=2000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     +JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      LOGICAL CHKIF
      COMMON/PHOIF/CHKIF(NMXPHO)
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      LOGICAL QEDRAD
      COMMON/PHOQED/QEDRAD(NMXHEP)
      INTEGER JFIRST
      LOGICAL F
      INTEGER IDABS,NLAST,I,IPPAR
      LOGICAL INTERF,ISEC,IFTOP
      REAL FINT,FSEC
      COMMON /PHOKEY/ INTERF,FINT,ISEC,FSEC,IFTOP
      LOGICAL IFRAD
      INTEGER IDENT,K
C THESE ARE OK .... IF YOU DO NOT LIKE SOMEBODY ELSE, ADD HERE.
      F(IDABS)=
     +     ( ((IDABS.GT.9).AND.(IDABS.LE.40)) .OR. (IDABS.GT.100) )
     + .AND.(IDABS.NE.21)
     + .AND.(IDABS.NE.2101).AND.(IDABS.NE.3101).AND.(IDABS.NE.3201)
     + .AND.(IDABS.NE.1103).AND.(IDABS.NE.2103).AND.(IDABS.NE.2203)
     + .AND.(IDABS.NE.3103).AND.(IDABS.NE.3203).AND.(IDABS.NE.3303)
C
      NLAST = NPHO
C
      IPPAR=1
C CHECKING FOR GOOD PARTICLES
      DO 10 I=IPPAR,NLAST
        IDABS = ABS(IDPHO(I))
C POSSIBLY CALL ON PHZODE IS A DEAD (TO BE OMITTED) CODE.
        CHKIF(I)= F(IDABS) .AND.F(ABS(IDPHO(1))) .AND. (IDPHO(2).EQ.0)
        IF(I.GT.2) CHKIF(I)=CHKIF(I).AND.QEDRAD(JFIRST+I-IPPAR-2)
   10 CONTINUE
C--
C NOW WE GO TO SPECIAL CASES, WHERE CHKIF(I) WILL BE OVERWRITTEN
C--
      IF(IFTOP) THEN
C SPECIAL CASE OF TOP PAIR PRODUCTION
        DO  K=JDAPHO(2,1),JDAPHO(1,1),-1
          IF(IDPHO(K).NE.22) THEN
            IDENT=K
            GOTO 20
          ENDIF
        ENDDO
   20   CONTINUE
        IFRAD=((IDPHO(1).EQ.21).AND.(IDPHO(2).EQ.21))
     +  .OR. ((ABS(IDPHO(1)).LE.6).AND.((IDPHO(2)).EQ.(-IDPHO(1))))
        IFRAD=IFRAD
     +        .AND.(ABS(IDPHO(3)).EQ.6).AND.((IDPHO(4)).EQ.(-IDPHO(3)))
     +        .AND.(IDENT.EQ.4)
        IF(IFRAD) THEN
          DO 30 I=IPPAR,NLAST
            CHKIF(I)= .TRUE.
            IF(I.GT.2) CHKIF(I)=CHKIF(I).AND.QEDRAD(JFIRST+I-IPPAR-2)
   30     CONTINUE
        ENDIF
      ENDIF
C--
C--
      IF(IFTOP) THEN
C SPECIAL CASE OF TOP DECAY
        DO  K=JDAPHO(2,1),JDAPHO(1,1),-1
          IF(IDPHO(K).NE.22) THEN
            IDENT=K
            GOTO 40
          ENDIF
        ENDDO
   40   CONTINUE
        IFRAD=((ABS(IDPHO(1)).EQ.6).AND.(IDPHO(2).EQ.0))
        IFRAD=IFRAD
     +        .AND.((ABS(IDPHO(3)).EQ.24).AND.(ABS(IDPHO(4)).EQ.5)
     +        .OR.(ABS(IDPHO(3)).EQ.5).AND.(ABS(IDPHO(4)).EQ.24))
     +        .AND.(IDENT.EQ.4)
        IF(IFRAD) THEN
          DO 50 I=IPPAR,NLAST
            CHKIF(I)= .TRUE.
            IF(I.GT.2) CHKIF(I)=CHKIF(I).AND.QEDRAD(JFIRST+I-IPPAR-2)
   50     CONTINUE
        ENDIF
      ENDIF
C--
C--
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOCIN
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON COMMON INITIALISATION
C.
C.    PURPOSE:  INITIALISATION OF PARAMETERS IN COMMON BLOCKS.
C.
C.    INPUT PARAMETERS:   NONE
C.
C.    OUTPUT PARAMETERS:  COMMONS /PHOLUN/, /PHOPHO/, /PHOCOP/, /PHPICO/
C.                                AND /PHSEED/.
C.
C.    AUTHOR(S):  B. VAN EIJK                     CREATED AT:  26/11/89
C.                Z. WAS                          LAST UPDATE: 10/08/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      LOGICAL QEDRAD
      COMMON/PHOQED/QEDRAD(NMXHEP)
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      REAL ALPHA,XPHCUT
      COMMON/PHOCOP/ALPHA,XPHCUT
      REAL PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      INTEGER ISEED,I97,J97
      REAL URAN,CRAN,CDRAN,CMRAN
      COMMON/PHSEED/ISEED(2),I97,J97,URAN(97),CRAN,CDRAN,CMRAN
      INTEGER PHOMES
      PARAMETER (PHOMES=10)
      INTEGER STATUS
      COMMON/PHOSTA/STATUS(PHOMES)
      LOGICAL INTERF,ISEC,IFTOP
      REAL FINT,FSEC
      COMMON /PHOKEY/ INTERF,FINT,ISEC,FSEC,IFTOP
      INTEGER INIT,I
      SAVE INIT
      DATA INIT/ 0/
C--
C--   RETURN IF ALREADY INITIALIZED...
      IF (INIT.NE.0) RETURN
      INIT=1
C--
C--   PRESET SWITCH  FOR  PHOTON EMISSION TO 'TRUE' FOR EACH PARTICLE IN
C--   /HEPEVT/, THIS INTERFACE IS NEEDED FOR KORALB AND KORALZ...
      DO 10 I=1,NMXHEP
   10 QEDRAD(I)=.TRUE.
C--
C--   LOGICAL OUTPUT UNIT FOR PRINTING OF PHOTOS ERROR MESSAGES
      PHLUN=6
C--
C--   SET CUT PARAMETER FOR PHOTON RADIATION
      XPHCUT=0.01
C--
C--   DEFINE SOME CONSTANTS
      ALPHA=0.00729735039
      PI=3.14159265358979324
      TWOPI=6.28318530717958648
C--
C--   DEFAULT SEEDS MARSAGLIA AND ZAMAN RANDOM NUMBER GENERATOR
      ISEED(1)=1802
      ISEED(2)=9373
C--
C--   IITIALIZATION FOR EXTRA OPTIONS
C--   (1)
C--   INTERFERENCE WEIGHT FOR TWO BODY SYMMETRIC CHANNELS ONLY.
      INTERF=.TRUE.
C--   (2)
C--   SECOND ORDER - DOUBLE PHOTON SWITCH
      ISEC=.TRUE.
C--   (3)
C--   EMISION IN THE HARD PROCESS G G (Q QBAR) --> T TBAR
C--                                 T          --> W B
      IFTOP=.TRUE.
C--
C--   FURTHER INITIALIZATION DONE AUTOMATICALLY
      IF (INTERF) THEN
C--   BEST CHOICE IS IF FINT=2**N WHERE N+1 IS MAXIMAL NUMBER OF CHARGED DAUGHTE
C--   SEE REPORT ON OVERWEIHTED EVENTS
        FINT=2.0
      ELSE
        FINT=1.0
      ENDIF
C--   INITIALISE STATUS COUNTER FOR WARNING MESSAGES
      DO 20 I=1,PHOMES
   20 STATUS(I)=0
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION PHOCOR(MPASQR,MCHREN,ME)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS CORRECTION WEIGHT FROM
C.              MATRIX ELEMENTS
C.
C.    PURPOSE:  CALCULATE  PHOTON  ANGLE.  THE RESHAPING FUNCTIONS  WILL
C.              HAVE  TO  DEPEND  ON THE SPIN S OF THE CHARGED PARTICLE.
C.              WE DEFINE:  ME = 2 * S + 1 !
C.
C.    INPUT PARAMETERS:  MPASQR:  PARENT MASS SQUARED,
C.                       MCHREN:  RENORMALISED MASS OF CHARGED SYSTEM,
C.                       ME:      2 * SPIN + 1 DETERMINES MATRIX ELEMENT
C.
C.    OUTPUT PARAMETER:  FUNCTION VALUE.
C.
C.    AUTHOR(S):  Z. WAS, B. VAN EIJK             CREATED AT:  26/11/89
C.                                                LAST UPDATE: 21/03/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION MPASQR,MCHREN,BETA,XX,YY,DATA
      INTEGER ME
      REAL PHOCOR,PHOFAC,WT1,WT2,WT3
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL PNEUTR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
      INTEGER IREP
      REAL PROBH,CORWT,XF
      COMMON/PHOPRO/IREP,PROBH,CORWT,XF
C--
C--   SHAPING (MODIFIED BY ZW)...
      XX=4.*MCHSQR/MPASQR*(1.-XPHOTO)/(1.-XPHOTO+(MCHSQR-MNESQR)/
     &MPASQR)**2
      IF (ME.EQ.1) THEN
        YY=1.
        WT3=(1.-XPHOTO/XPHMAX)/((1.+(1.-XPHOTO/XPHMAX)**2)/2.)
      ELSEIF (ME.EQ.2) THEN
        YY=0.5*(1.-XPHOTO/XPHMAX+1./(1.-XPHOTO/XPHMAX))
        WT3=1.
      ELSEIF ((ME.EQ.3).OR.(ME.EQ.4).OR.(ME.EQ.5)) THEN
        YY=1.
        WT3=(1.+(1.-XPHOTO/XPHMAX)**2-(XPHOTO/XPHMAX)**3)/(1.+(1.
     &  -XPHOTO/XPHMAX)** 2)
      ELSE
        DATA=(ME-1.)/2.
        CALL PHOERR(6,'PHOCOR',DATA)
        YY=1.
        WT3=1.
      ENDIF
      BETA=SQRT(1.-XX)
      WT1=(1.-COSTHG*SQRT(1.-MCHREN))/(1.-COSTHG*BETA)
      WT2=(1.-XX/YY/(1.-BETA**2*COSTHG**2))*(1.+COSTHG*BETA)/2.
      WT2=WT2*PHOFAC(1)
      PHOCOR=WT1*WT2*WT3
      CORWT=PHOCOR
      IF (PHOCOR.GT.1.) THEN
        DATA=PHOCOR
        CALL PHOERR(3,'PHOCOR',DATA)
      ENDIF
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHODO(IP,NCHARB,NEUDAU)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN  DECAYS DOING OF KINEMATICS
C.
C.    PURPOSE:  STARTING  FROM   THE  CHARGED  PARTICLE ENERGY/MOMENTUM,
C.              PNEUTR, PHOTON  ENERGY  FRACTION AND PHOTON  ANGLE  WITH
C.              RESPECT  TO  THE AXIS FORMED BY CHARGED PARTICLE ENERGY/
C.              MOMENTUM  VECTOR  AND PNEUTR, SCALE THE ENERGY/MOMENTUM,
C.              KEEPING THE ORIGINAL DIRECTION OF THE NEUTRAL SYSTEM  IN
C.              THE LAB. FRAME UNTOUCHED.
C.
C.    INPUT PARAMETERS:   IP:      POINTER  TO   DECAYING  PARTICLE   IN
C.                                 /PHOEVT/  AND   THE   COMMON   ITSELF
C.                        NCHARB:  POINTER TO THE CHARGED RADIATING
C.                                 DAUGHTER IN /PHOEVT/.
C.                        NEUDAU:  POINTER TO THE FIRST NEUTRAL DAUGHTER
C.    OUTPUT PARAMETERS:  COMMON /PHOEVT/, WITH PHOTON ADDED.
C.
C.    AUTHOR(S):  Z. WAS, B. VAN EIJK             CREATED AT:  26/11/89
C.                                                LAST UPDATE: 27/05/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION PHOAN1,PHOAN2,ANGLE,FI1,FI3,FI4,FI5,TH1,TH3,TH4
      DOUBLE PRECISION PARNE,QNEW,QOLD,DATA
      INTEGER IP,FI3DUM,I,J,NEUDAU,FIRST,LAST
      INTEGER NCHARB
      REAL EPHOTO,PMAVIR,PHOTRI
      REAL GNEUT,PHORAN,CCOSTH,SSINTH,PVEC(4)
      INTEGER NMXPHO
      PARAMETER (NMXPHO=2000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     +JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL PNEUTR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
      REAL PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
C--
      EPHOTO=XPHOTO*PPHO(5,IP)/2.
      PMAVIR=SQRT(PPHO(5,IP)*(PPHO(5,IP)-2.*EPHOTO))
C--
C--   RECONSTRUCT  KINEMATICS  OF  CHARGED PARTICLE  AND  NEUTRAL SYSTEM
      FI1=PHOAN1(PNEUTR(1),PNEUTR(2))
C--
C--   CHOOSE AXIS ALONG  Z OF  PNEUTR, CALCULATE  ANGLE  BETWEEN X AND Y
C--   COMPONENTS  AND Z  AND X-Y PLANE AND  PERFORM LORENTZ TRANSFORM...
      TH1=PHOAN2(PNEUTR(3),SQRT(PNEUTR(1)**2+PNEUTR(2)**2))
      CALL PHORO3(-FI1,PNEUTR(1))
      CALL PHORO2(-TH1,PNEUTR(1))
C--
C--   TAKE  AWAY  PHOTON ENERGY FROM CHARGED PARTICLE AND PNEUTR !  THUS
C--   THE ONSHELL CHARGED PARTICLE  DECAYS INTO VIRTUAL CHARGED PARTICLE
C--   AND PHOTON.  THE VIRTUAL CHARGED  PARTICLE MASS BECOMES:
C--   SQRT(PPHO(5,IP)*(PPHO(5,IP)-2*EPHOTO)).  CONSTRUCT  NEW PNEUTR MO-
C--   MENTUM IN THE REST FRAME OF THE PARENT:
C--   1) SCALING PARAMETERS...
      QNEW=PHOTRI(PMAVIR,PNEUTR(5),PPHO(5,NCHARB))
      QOLD=PNEUTR(3)
      GNEUT=(QNEW**2+QOLD**2+MNESQR)/(QNEW*QOLD+SQRT((QNEW**2+MNESQR)*
     +(QOLD**2+MNESQR)))
      IF (GNEUT.LT.1.) THEN
        DATA=0.
        CALL PHOERR(4,'PHOKIN',DATA)
      ENDIF
      PARNE=GNEUT-SQRT(MAX(GNEUT**2-1.0,0.))
C--
C--   2) ...REDUCTIVE BOOST...
      CALL PHOBO3(PARNE,PNEUTR)
C--
C--   ...CALCULATE PHOTON ENERGY IN THE REDUCED SYSTEM...
      NPHO=NPHO+1
      ISTPHO(NPHO)=1
      IDPHO(NPHO) =22
C--   PHOTON MOTHER AND DAUGHTER POINTERS !
      JMOPHO(1,NPHO)=IP
      JMOPHO(2,NPHO)=0
      JDAPHO(1,NPHO)=0
      JDAPHO(2,NPHO)=0
      PPHO(4,NPHO)=EPHOTO*PPHO(5,IP)/PMAVIR
C--
C--   ...AND PHOTON MOMENTA
      CCOSTH=-COSTHG
      SSINTH=SINTHG
      TH3=PHOAN2(CCOSTH,SSINTH)
      FI3=TWOPI*PHORAN(FI3DUM)
      PPHO(1,NPHO)=PPHO(4,NPHO)*SINTHG*COS(FI3)
      PPHO(2,NPHO)=PPHO(4,NPHO)*SINTHG*SIN(FI3)
C--
C--   MINUS SIGN BECAUSE AXIS OPPOSITE DIRECTION OF CHARGED PARTICLE !
      PPHO(3,NPHO)=-PPHO(4,NPHO)*COSTHG
      PPHO(5,NPHO)=0.
C--
C--   ROTATE IN ORDER TO GET PHOTON ALONG Z-AXIS
      CALL PHORO3(-FI3,PNEUTR(1))
      CALL PHORO3(-FI3,PPHO(1,NPHO))
      CALL PHORO2(-TH3,PNEUTR(1))
      CALL PHORO2(-TH3,PPHO(1,NPHO))
      ANGLE=EPHOTO/PPHO(4,NPHO)
C--
C--   BOOST TO THE REST FRAME OF DECAYING PARTICLE
      CALL PHOBO3(ANGLE,PNEUTR(1))
      CALL PHOBO3(ANGLE,PPHO(1,NPHO))
C--
C--   BACK IN THE PARENT REST FRAME BUT PNEUTR NOT YET ORIENTED !
      FI4=PHOAN1(PNEUTR(1),PNEUTR(2))
      TH4=PHOAN2(PNEUTR(3),SQRT(PNEUTR(1)**2+PNEUTR(2)**2))
      CALL PHORO3(FI4,PNEUTR(1))
      CALL PHORO3(FI4,PPHO(1,NPHO))
C--
      DO 10 I=2,4
   10 PVEC(I)=0.
      PVEC(1)=1.
      CALL PHORO3(-FI3,PVEC)
      CALL PHORO2(-TH3,PVEC)
      CALL PHOBO3(ANGLE,PVEC)
      CALL PHORO3(FI4,PVEC)
      CALL PHORO2(-TH4,PNEUTR)
      CALL PHORO2(-TH4,PPHO(1,NPHO))
      CALL PHORO2(-TH4,PVEC)
      FI5=PHOAN1(PVEC(1),PVEC(2))
C--
C--   CHARGED PARTICLE RESTORES ORIGINAL DIRECTION
      CALL PHORO3(-FI5,PNEUTR)
      CALL PHORO3(-FI5,PPHO(1,NPHO))
      CALL PHORO2(TH1,PNEUTR(1))
      CALL PHORO2(TH1,PPHO(1,NPHO))
      CALL PHORO3(FI1,PNEUTR)
      CALL PHORO3(FI1,PPHO(1,NPHO))
C--   SEE WHETHER NEUTRAL SYSTEM HAS MULTIPLICITY LARGER THAN 1...
      IF ((JDAPHO(2,IP)-JDAPHO(1,IP)).GT.1) THEN
C--   FIND POINTERS TO COMPONENTS OF 'NEUTRAL' SYSTEM
C--
        FIRST=NEUDAU
        LAST=JDAPHO(2,IP)
        DO 20 I=FIRST,LAST
          IF (I.NE.NCHARB.AND.(JMOPHO(1,I).EQ.IP)) THEN
C--
C--   RECONSTRUCT KINEMATICS...
            CALL PHORO3(-FI1,PPHO(1,I))
            CALL PHORO2(-TH1,PPHO(1,I))
C--
C--   ...REDUCTIVE BOOST
            CALL PHOBO3(PARNE,PPHO(1,I))
C--
C--   ROTATE IN ORDER TO GET PHOTON ALONG Z-AXIS
            CALL PHORO3(-FI3,PPHO(1,I))
            CALL PHORO2(-TH3,PPHO(1,I))
C--
C--   BOOST TO THE REST FRAME OF DECAYING PARTICLE
            CALL PHOBO3(ANGLE,PPHO(1,I))
C--
C--   BACK IN THE PARENT REST-FRAME BUT PNEUTR NOT YET ORIENTED.
            CALL PHORO3(FI4,PPHO(1,I))
            CALL PHORO2(-TH4,PPHO(1,I))
C--
C--   CHARGED PARTICLE RESTORES ORIGINAL DIRECTION
            CALL PHORO3(-FI5,PPHO(1,I))
            CALL PHORO2(TH1,PPHO(1,I))
            CALL PHORO3(FI1,PPHO(1,I))
          ENDIF
   20   CONTINUE
      ELSE
C--
C--   ...ONLY ONE 'NEUTRAL' PARTICLE IN ADDITION TO PHOTON!
        DO 30 J=1,4
   30   PPHO(J,NEUDAU)=PNEUTR(J)
      ENDIF
C--
C--   ALL 'NEUTRALS' TREATED, FILL /PHOEVT/ FOR CHARGED PARTICLE...
      DO 40 J=1,3
   40 PPHO(J,NCHARB)=-(PPHO(J,NPHO)+PNEUTR(J))
      PPHO(4,NCHARB)=PPHO(5,IP)-(PPHO(4,NPHO)+PNEUTR(4))
C--
      END
*CMZ :  1.01/50 19/04/96  12.03.42  by  Piero Zucchelli
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOENE(MPASQR,MCHREN,BETA,IDENT)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS CALCULATION  OF PHOTON ENERGY
C.              FRACTION
C.
C.    PURPOSE:  SUBROUTINE  RETURNS  PHOTON  ENERGY FRACTION (IN (PARENT
C.              MASS)/2 UNITS) FOR THE DECAY BREMSSTRAHLUNG.
C.
C.    INPUT PARAMETERS:  MPASQR:  MASS OF DECAYING SYSTEM SQUARED,
C.                       XPHCUT:  MINIMUM ENERGY FRACTION OF PHOTON,
C.                       XPHMAX:  MAXIMUM ENERGY FRACTION OF PHOTON.
C.
C.    OUTPUT PARAMETER:  MCHREN:  RENORMALISED MASS SQUARED,
C.                       BETA:    BETA FACTOR DUE TO RENORMALISATION,
C.                       XPHOTO:  PHOTON ENERGY FRACTION,
C.                       XF:      CORRECTION FACTOR FOR PHOFAC.
C.
C.    AUTHOR(S):  S. JADACH, Z. WAS               CREATED AT:  01/01/89
C.                B. VAN EIJK                     LAST UPDATE: 26/03/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE

      DOUBLE PRECISION MPASQR,MCHREN,BIGLOG,BETA,DATA
      INTEGER IWT1,IRN,IWT2
      REAL PRSOFT,PRHARD,PHORAN,PHOFAC
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL PNEUTR
      INTEGER IDENT
      REAL PHOCHA
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
      REAL ALPHA,XPHCUT
      COMMON/PHOCOP/ALPHA,XPHCUT
      REAL PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      INTEGER IREP
      REAL PROBH,CORWT,XF
      COMMON/PHOPRO/IREP,PROBH,CORWT,XF
      LOGICAL INTERF,ISEC,IFTOP
      REAL FINT,FSEC
      COMMON /PHOKEY/ INTERF,FINT,ISEC,FSEC,IFTOP
C--
      IF (XPHMAX.LE.XPHCUT) THEN
        XPHOTO=0.0
        RETURN
      ENDIF
C--   PROBABILITIES FOR HARD AND SOFT BREMSTRAHLUNG...
      MCHREN=4.*MCHSQR/MPASQR/(1.+MCHSQR/MPASQR)**2
      BETA=SQRT(1.-MCHREN)
      BIGLOG=LOG(MPASQR/MCHSQR*(1.+BETA)**2/4.*(1.+MCHSQR/MPASQR)**2)
      PRHARD=ALPHA/PI/BETA*BIGLOG*(LOG(XPHMAX/XPHCUT)-.75+XPHCUT/
     &XPHMAX-.25*XPHCUT**2/XPHMAX**2)
      PRHARD=PRHARD*PHOCHA(IDENT)**2*FINT*FSEC
      IF (IREP.EQ.0) PROBH=0.
      PRHARD=PRHARD*PHOFAC(0)
      PROBH=PRHARD
      PRSOFT=1.-PRHARD
C--
C--   CHECK ON KINEMATICAL BOUNDS
      IF (PRSOFT.LT.0.1) THEN
        DATA=PRSOFT
        CALL PHOERR(2,'PHOENE',DATA)
      ENDIF
      IF (PHORAN(IWT1).LT.PRSOFT) THEN
C--
C--   NO PHOTON... (IE. PHOTON TOO SOFT)
        XPHOTO=0.
      ELSE
C--
C--   HARD  PHOTON... (IE.  PHOTON  HARD ENOUGH).
C--   CALCULATE  ALTARELLI-PARISI KERNEL
   10   XPHOTO=EXP(PHORAN(IRN)*LOG(XPHCUT/XPHMAX))
        XPHOTO=XPHOTO*XPHMAX
        IF (PHORAN(IWT2).GT.((1.+(1.-XPHOTO/XPHMAX)**2)/2.)) GOTO 10
      ENDIF
C--
C--   CALCULATE PARAMETER FOR PHOFAC FUNCTION
      XF=4.*MCHSQR*MPASQR/(MPASQR+MCHSQR-MNESQR)**2
      RETURN
      END
*CMZ :  1.02/02 12/01/97  17.54.22  by  P. Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOERR(IMES,TEXT,DATA)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS ERRROR HANDLING
C.
C.    PURPOSE:  INFORM USER  ABOUT (FATAL) ERRORS AND WARNINGS GENERATED
C.              BY EITHER THE USER OR THE PROGRAM.
C.
C.    INPUT PARAMETERS:   IMES, TEXT, DATA
C.
C.    OUTPUT PARAMETERS:  NONE
C.
C.    AUTHOR(S):  B. VAN EIJK                     CREATED AT:  29/11/89
C.                                                LAST UPDATE: 10/01/92
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION DATA
      INTEGER IMES,IERROR
      REAL SDATA
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      INTEGER PHOMES
      PARAMETER (PHOMES=10)
      INTEGER STATUS
      COMMON/PHOSTA/STATUS(PHOMES)
      CHARACTER TEXT*(*)
      SAVE IERROR
C--   SECURITY STOP SWITCH
      LOGICAL ISEC
      SAVE ISEC
      DATA IERROR/ 0/
      DATA ISEC /.TRUE./
      IF (IMES.LE.PHOMES) STATUS(IMES)=STATUS(IMES)+1
C--
C--   COUNT NUMBER OF NON-FATAL ERRORS...
      IF ((IMES.EQ. 6).AND.(STATUS(IMES).GE.2)) RETURN
      IF ((IMES.EQ.10).AND.(STATUS(IMES).GE.2)) RETURN
      SDATA=DATA
      WRITE(PHLUN,10000)
      WRITE(PHLUN,11100)
      GOTO (10,20,30,40,50,60,70,80,90,100),IMES
      WRITE(PHLUN,11200) IMES
      GOTO 120
   10 WRITE(PHLUN,10100) TEXT,INT(SDATA)
      GOTO 110
   20 WRITE(PHLUN,10200) TEXT,SDATA
      GOTO 110
   30 WRITE(PHLUN,10300) TEXT,SDATA
      GOTO 110
   40 WRITE(PHLUN,10400) TEXT
      GOTO 110
   50 WRITE(PHLUN,10500) TEXT,INT(SDATA)
      GOTO 110
   60 WRITE(PHLUN,10600) TEXT,SDATA
      GOTO 130
   70 WRITE(PHLUN,10700) TEXT,INT(SDATA)
      GOTO 110
   80 WRITE(PHLUN,10800) TEXT,INT(SDATA)
      GOTO 110
   90 WRITE(PHLUN,10900) TEXT,INT(SDATA)
      GOTO 110
  100 WRITE(PHLUN,11000) TEXT,SDATA
      GOTO 130
  110 CONTINUE
      WRITE(PHLUN,11300)
      WRITE(PHLUN,11100)
      WRITE(PHLUN,10000)
      IF (ISEC) THEN
        STOP
      ELSE
        GOTO 130
      ENDIF
  120 IERROR=IERROR+1
      IF (IERROR.GE.10) THEN
        WRITE(PHLUN,11400)
        WRITE(PHLUN,11100)
        WRITE(PHLUN,10000)
        IF (ISEC) THEN
          STOP
        ELSE
          GOTO 130
        ENDIF
      ENDIF
  130 WRITE(PHLUN,11100)
      WRITE(PHLUN,10000)
      RETURN
10000 FORMAT(1H ,80('*'))
10100 FORMAT(1H ,'* ',A,': TOO MANY CHARGED PARTICLES, NCHARG =',I6,T81,
     &'*')
10200 FORMAT(1H ,'* ',A,': TOO MUCH BREMSSTRAHLUNG REQUIRED, PRSOFT = ',
     &F15.6,T81,'*')
10300 FORMAT(1H ,'* ',A,': COMBINED WEIGHT IS EXCEEDING 1., WEIGHT = ',
     &F15.6,T81,'*')
10400 FORMAT(1H ,'* ',A,
     &': ERROR IN RESCALING CHARGED AND NEUTRAL VECTORS',T81,'*')
10500 FORMAT(1H ,'* ',A,
     &': NON MATCHING CHARGED PARTICLE POINTER, NCHARG = ',I5,T81,'*')
10600 FORMAT(1H ,'* ',A,
     &': DO YOU REALLY WORK WITH A PARTICLE OF SPIN: ',F4.1,' ?',T81,
     &'*')
10700 FORMAT(1H ,'* ',A, ': STACK LENGTH EXCEEDED, NSTACK = ',I5 ,T81,
     &'*')
10800 FORMAT(1H ,'* ',A,
     &': RANDOM NUMBER GENERATOR SEED(1) OUT OF RANGE: ',I8,T81,'*')
10900 FORMAT(1H ,'* ',A,
     &': RANDOM NUMBER GENERATOR SEED(2) OUT OF RANGE: ',I8,T81,'*')
11000 FORMAT(1H ,'* ',A,
     &': AVAILABLE PHASE SPACE BELOW CUT-OFF: ',F15.6,' GEV/C^2',T81,
     &'*')
11100 FORMAT(1H ,'*',T81,'*')
11200 FORMAT(1H ,'* FUNNY ERROR MESSAGE: ',I4,' ! WHAT TO DO ?',T81,'*')
11300 FORMAT(1H ,'* FATAL ERROR MESSAGE, I STOP THIS RUN !',T81,'*')
11400 FORMAT(1H ,'* 10 ERROR MESSAGES GENERATED, I STOP THIS RUN !',T81,
     &'*')
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION PHOFAC(MODE)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS CONTROL FACTOR
C.
C.    PURPOSE:  THIS IS THE CONTROL FUNCTION FOR THE PHOTON SPECTRUM AND
C.              FINAL WEIGHTING.  IT IS  CALLED  FROM PHOENE FOR GENERA-
C.              TING THE RAW PHOTON ENERGY SPECTRUM (MODE=0) AND IN PHO-
C.              COR TO SCALE THE FINAL WEIGHT (MODE=1).  THE FACTOR CON-
C.              SISTS OF 3 TERMS.  ADDITION OF  THE FACTOR FF WHICH MUL-
C.              TIPLIES PHOFAC FOR MODE=0 AND DIVIDES PHOFAC FOR MODE=1,
C.              DOES NOT AFFECT  THE RESULTS FOR  THE MC GENERATION.  AN
C.              APPROPRIATE CHOICE  FOR FF CAN SPEED UP THE CALCULATION.
C.              NOTE THAT A TOO SMALL VALUE OF FF MAY CAUSE WEIGHT OVER-
C.              FLOW IN PHOCOR  AND WILL GENERATE A WARNING, HALTING THE
C.              EXECUTION.  PRX  SHOULD  BE  INCLUDED FOR REPEATED CALLS
C.              FOR  THE  SAME EVENT, ALLOWING MORE PARTICLES TO RADIATE
C.              PHOTONS.  AT  THE  FIRST  CALL IREP=0, FOR  MORE  THAN 1
C.              CHARGED  DECAY  PRODUCTS, IREP >= 1.  THUS,  PRSOFT  (NO
C.              PHOTON RADIATION  PROBABILITY  IN  THE  PREVIOUS  CALLS)
C.              APPROPRIATELY SCALES THE STRENGTH OF THE BREMSSTRAHLUNG.
C.
C.    INPUT PARAMETERS:  MODE, PROBH, XF
C.
C.    OUTPUT PARAMETER:  FUNCTION VALUE
C.
C.    AUTHOR(S):  S. JADACH, Z. WAS               CREATED AT:  01/01/89
C.                B. VAN EIJK                     LAST UPDATE: 13/02/90
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      REAL PHOFAC,FF,PRX
      INTEGER MODE
      INTEGER IREP
      REAL PROBH,CORWT,XF
      COMMON/PHOPRO/IREP,PROBH,CORWT,XF
      SAVE PRX,FF
      DATA PRX,FF/ 0., 0./
      IF (MODE.EQ.0) THEN
        IF (IREP.EQ.0) PRX=1.
        PRX=PRX/(1.-PROBH)
        FF=1.
C--
C--   FOLLOWING OPTIONS ARE NOT CONSIDERED FOR THE TIME BEING...
C--   (1) GOOD CHOICE, BUT DOES NOT SAVE VERY MUCH TIME:
C--       FF=(1.0-SQRT(XF)/2.0)/(1.0+SQRT(XF)/2.0)
C--   (2) TAKEN FROM THE BLUE, BUT WORKS WITHOUT WEIGHT OVERFLOWS...
C--       FF=(1.-XF/(1-(1-SQRT(XF))**2))*(1+(1-SQRT(XF))/SQRT(1-XF))/2
        PHOFAC=FF*PRX
      ELSE
        PHOFAC=1./FF
      ENDIF
      END
*CMZ :  1.01/50 23/05/96  10.22.20  by  Piero Zucchelli
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOIN(IP,BOOST,NHEP0)
C.----------------------------------------------------------------------
C.
C.    PHOIN:   PHOTOS INPUT
C.
C.    PURPOSE:  COPIES IP BRANCH OF THE COMMON /HEPEVT/ INTO /PHOEVT/
C.              MOVES BRANCH INTO ITS CMS SYSTEM.
C.
C.    INPUT PARAMETERS:       IP:  POINTER OF PARTICLE STARTING BRANCH
C.                                 TO BE COPIED
C.                        BOOST:   FLAG WHETHER BOOST TO CMS WAS OR WAS
C     .                            NOT PERFORMED.
C.
C.    OUTPUT PARAMETERS:  COMMONS: /PHOEVT/, /PHOCMS/
C.
C.    AUTHOR(S):  Z. WAS                          CREATED AT:  24/05/93
C.                                                LAST UPDATE: 16/11/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      DOUBLE PRECISION PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     +JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      INTEGER NMXPHO
      PARAMETER (NMXPHO=2000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     +JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      INTEGER IP,IP2,I,FIRST,LAST,LL,NA
      LOGICAL BOOST
      INTEGER J,NHEP0
      DOUBLE PRECISION BET(3),GAM,PB
      COMMON /PHOCMS/ BET,GAM
      LOGICAL INTERF,ISEC,IFTOP
      REAL FINT,FSEC
      COMMON /PHOKEY/ INTERF,FINT,ISEC,FSEC,IFTOP
C--
C LET'S CALCULATE SIZE OF THE LITTLE COMMON ENTRY
      FIRST=JDAHEP(1,IP)
      LAST =JDAHEP(2,IP)
      NPHO=3+LAST-FIRST+NHEP-NHEP0
      NEVPHO=NPHO
C LET'S TAKE IN DECAYING PARTICLE
      IDPHO(1)=IDHEP(IP)
      JDAPHO(1,1)=3
      JDAPHO(2,1)=3+LAST-FIRST
      DO I=1,5
        PPHO(I,1)=PHEP(I,IP)
      ENDDO
C LET'S TAKE IN EVENTUAL SECOND MOTHER
      IP2=JMOHEP(2,JDAHEP(1,IP))
      IF((IP2.NE.0).AND.(IP2.NE.IP)) THEN
        IDPHO(2)=IDHEP(IP2)
        JDAPHO(1,2)=3
        JDAPHO(2,2)=3+LAST-FIRST
        DO I=1,5
          PPHO(I,2)=PHEP(I,IP2)
        ENDDO
      ELSE
        IDPHO(2)=0
        DO I=1,5
          PPHO(I,2)=0.0
        ENDDO
      ENDIF
C LET'S TAKE IN DAUGHTERS
      DO LL=0,LAST-FIRST
        IDPHO(3+LL)=IDHEP(FIRST+LL)
        JMOPHO(1,3+LL)=JMOHEP(1,FIRST+LL)
        IF (JMOHEP(1,FIRST+LL).EQ.IP) JMOPHO(1,3+LL)=1
        DO I=1,5
          PPHO(I,3+LL)=PHEP(I,FIRST+LL)
        ENDDO
      ENDDO
      IF (NHEP.GT.NHEP0) THEN
C LET'S TAKE IN ILLEGITIMATE DAUGHTERS
        NA=3+LAST-FIRST
        DO LL=1,NHEP-NHEP0
          IDPHO(NA+LL)=IDHEP(NHEP0+LL)
          JMOPHO(1,NA+LL)=JMOHEP(1,NHEP0+LL)
          IF (JMOHEP(1,NHEP0+LL).EQ.IP) JMOPHO(1,NA+LL)=1
          DO I=1,5
            PPHO(I,NA+LL)=PHEP(I,NHEP0+LL)
          ENDDO
        ENDDO
C--        THERE IS NHEP-NHEP0 DAUGTERS MORE.
        JDAPHO(2,1)=3+LAST-FIRST+NHEP-NHEP0
      ENDIF
      CALL PHLUPA(1)
C SPECIAL CASE OF T TBAR PRODUCTION PROCESS
      IF(IFTOP) CALL PHOTWO(0)
      BOOST=.FALSE.
C--   CHECK WHETHER PARENT IS IN ITS REST FRAME...
      IF (     (ABS(PPHO(4,1)-PPHO(5,1)).GT.PPHO(5,1)*1.E-8)
     +    .AND.(PPHO(5,1).NE.0))                            THEN
        BOOST=.TRUE.
C--
C--   BOOST DAUGHTER PARTICLES TO REST FRAME OF PARENT...
C--   RESULTANT NEUTRAL SYSTEM ALREADY CALCULATED IN REST FRAME !
        DO 10 J=1,3
   10   BET(J)=-PPHO(J,1)/PPHO(5,1)
        GAM=PPHO(4,1)/PPHO(5,1)
        DO 30 I=JDAPHO(1,1),JDAPHO(2,1)
          PB=BET(1)*PPHO(1,I)+BET(2)*PPHO(2,I)+BET(3)*PPHO(3,I)
          DO 20 J=1,3
   20     PPHO(J,I)=PPHO(J,I)+BET(J)*(PPHO(4,I)+PB/(GAM+1.))
   30   PPHO(4,I)=GAM*PPHO(4,I)+PB
C--    FINALLY BOOST MOTHER AS WELL
        I=1
        PB=BET(1)*PPHO(1,I)+BET(2)*PPHO(2,I)+BET(3)*PPHO(3,I)
        DO J=1,3
          PPHO(J,I)=PPHO(J,I)+BET(J)*(PPHO(4,I)+PB/(GAM+1.))
        ENDDO
        PPHO(4,I)=GAM*PPHO(4,I)+PB
      ENDIF
C SPECIAL CASE OF T TBAR PRODUCTION PROCESS
      IF(IFTOP) CALL PHOTWO(1)
      CALL PHLUPA(2)
      END
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOINF
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS GENERAL INFO
C.
C.    PURPOSE:  PRINT PHOTOS INFO
C.
C.    INPUT PARAMETERS:   PHOLUN
C.
C.    OUTPUT PARAMETERS:  PHOVN1, PHOVN2
C.
C.    AUTHOR(S):  B. VAN EIJK                     CREATED AT:  12/04/90
C.                                                LAST UPDATE: 02/10/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      INTEGER IV1,IV2,IV3
      INTEGER PHOVN1,PHOVN2
      COMMON/PHOVER/PHOVN1,PHOVN2
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      LOGICAL INTERF,ISEC,IFTOP
      REAL FINT,FSEC
      COMMON /PHOKEY/ INTERF,FINT,ISEC,FSEC,IFTOP
      REAL ALPHA,XPHCUT
      COMMON/PHOCOP/ALPHA,XPHCUT
C--
C--   PHOTOS VERSION NUMBER AND RELEASE DATE
      PHOVN1=200
      PHOVN2=161193
C--
C--   PRINT INFO
      WRITE(PHLUN,10000)
      WRITE(PHLUN,10200)
      WRITE(PHLUN,10100)
      WRITE(PHLUN,10300)
      IV1=PHOVN1/100
      IV2=PHOVN1-IV1*100
      WRITE(PHLUN,10400) IV1,IV2
      IV1=PHOVN2/10000
      IV2=(PHOVN2-IV1*10000)/100
      IV3=PHOVN2-IV1*10000-IV2*100
      WRITE(PHLUN,10500) IV1,IV2,IV3
      WRITE(PHLUN,10300)
      WRITE(PHLUN,10100)
      WRITE(PHLUN,10600)
      WRITE(PHLUN,10100)
      WRITE(PHLUN,11100)
      WRITE(PHLUN,10100)
      WRITE(PHLUN,10200)
      WRITE(PHLUN,10100)
      WRITE(PHLUN,11000) INTERF,ISEC,IFTOP,ALPHA,XPHCUT
      WRITE(PHLUN,10100)
      IF (INTERF) WRITE(PHLUN,10700)
      IF (ISEC) WRITE(PHLUN,10800)
      IF (IFTOP) WRITE(PHLUN,10900)
      WRITE(PHLUN,10100)
      WRITE(PHLUN,10200)
      RETURN
10000 FORMAT(1H1)
10100 FORMAT(1H ,'*',T81,'*')
10200 FORMAT(1H ,80('*'))
10300 FORMAT(1H ,'*',26X,26('='),T81,'*')
10400 FORMAT(1H ,'*',28X,'PHOTOS, VERSION: ',I2,'.',I2,T81,'*')
10500 FORMAT(1H ,'*',28X,'RELEASED AT:  ',I2,'/',I2,'/',I2,T81,'*')
10600 FORMAT(1H ,'*',18X,'PHOTOS QED CORRECTIONS IN PARTICLE DECAYS',
     &T81,'*')
10700 FORMAT(1H ,'*',18X,'OPTION WITH INTERFERENCE IS ACTIVE       ',
     &T81,'*')
10800 FORMAT(1H ,'*',18X,'OPTION WITH DOUBLE PHOTONS IS ACTIVE     ',
     &T81,'*')
10900 FORMAT(1H ,'*',18X,'EMISION IN T TBAR PRODUCTION IS ACTIVE   ',
     &T81,'*')
11000 FORMAT(1H ,'*',18X,'INTERNAL INPUT PARAMETERS:',T81,'*'
     &,/,    1H ,'*',T81,'*'
     &,/,    1H ,'*',18X,'INTERF=',L2,'  ISEC=',L2,'  IFTOP=',L2,T81,'*'
     &,/,    1H ,'*',18X,'ALPHA_QED=',F8.5,'   XPHCUT=',F8.5,T81,'*')
11100 FORMAT(1H ,'*',9X,'MONTE CARLO PROGRAM - BY E. BARBERIO, B. VAN EI
     &JK AND Z. WAS',T81,'*',/,
     &      1H ,'*',9X,'FROM VERSION 2.0 ON - BY E.B. AND Z.W.',T81,'*')
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTOS CDE'S
C.
C.    PURPOSE:  KEEP DEFINITIONS  FOR PHOTOS QED CORRECTION MONTE CARLO.
C.
C.    INPUT PARAMETERS:   NONE
C.
C.    OUTPUT PARAMETERS:  NONE
C.
C.    AUTHOR(S):  Z. WAS, B. VAN EIJK             CREATED AT:  29/11/89
C.                                                LAST UPDATE: 10/08/93
C.
C. =========================================================
C.    GENERAL STRUCTURE INFORMATION:                       =
C. =========================================================
C:   ROUTINES:
C.             1) INITIALIZATION:
C.                                      PHOCDE
C.                                      PHOINI
C.                                      PHOCIN
C.                                      PHOINF
C.             2) GENERAL INTERFACE:
C.                                      PHOTOS
C.                                      PHOBOS
C.                                      PHOIN
C.                                      PHOTWO (SPECIFIC INTERFACE
C.                                      PHOOUT
C.                                      PHOCHK
C.                                      PHTYPE (SPECIFIC INTERFACE
C.                                      PHOMAK (SPECIFIC INTERFACE
C.             3) QED PHOTON GENERATION:
C.                                      PHINT
C.                                      PHOPRE
C.                                      PHOOMA
C.                                      PHOENE
C.                                      PHOCOR
C.                                      PHOFAC
C.                                      PHODO
C.             4) UTILITIES:
C.                                      PHOTRI
C.                                      PHOAN1
C.                                      PHOAN2
C.                                      PHOBO3
C.                                      PHORO2
C.                                      PHORO3
C.                                      PHORIN
C.                                      PHORAN
C.                                      PHOCHA
C.                                      PHOSPI
C.                                      PHOERR
C.                                      PHOREP
C.                                      PHLUPA
C.   COMMONS:
C.   NAME     USED IN SECT. # OF OCC.     COMMENT
C.   PHOQED   1) 2)            3      FLAGS WHETHER EMISSON TO BE GENER.
C.   PHOLUN   1) 4)            5      OUTPUT DEVICE NUMBER
C.   PHOCOP   1) 3)            4      PHOTON COUPLING & MIN ENERGY
C.   PHPICO   1) 3) 4)         5      PI & 2*PI
C.   PHSEED   1) 4)            3      RN SEED
C.   PHOSTA   1) 4)            3      STATUS INFORMATION
C.   PHOKEY   1) 2) 3)         7      KEYS FOR NONSTANDARD APPLICATION
C.   PHOVER   1)               1      VERSION INFO FOR OUTSIDE
C.   HEPEVT   2)               6      PDG COMMON
C.   PHOEVT   2) 3)            9      PDG BRANCH
C.   PHOIF    2) 3)            2      EMISSION FLAGS FOR PDG BRANCH
C.   PHOMOM   3)               5      PARAM OF CHAR-NEUTR SYSTEM
C.   PHOPHS   3)               5      PHOTON MOMENTUM PARAMETERS
C.   PHOPRO   3)               4      VAR. FOR PHOTON REP. (IN BRANCH)
C.   PHOCMS   2)               3      PARAMETERS OF BOOST TO BRANCH CMS
C.   PHNUM    4)               1      EVENT NUMBER FROM OUTSIDE
C.----------------------------------------------------------------------
      SUBROUTINE PHOINI
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS INITIALISATION
C.
C.    PURPOSE:  INITIALISATION  ROUTINE  FOR  THE  PHOTOS  QED RADIATION
C.              PACKAGE.  SHOULD BE CALLED  AT LEAST ONCE  BEFORE A CALL
C.              TO THE STEERING PROGRAM 'PHOTOS' IS MADE.
C.
C.    INPUT PARAMETERS:   NONE
C.
C.    OUTPUT PARAMETERS:  NONE
C.
C.    AUTHOR(S):  Z. WAS, B. VAN EIJK             CREATED AT:  26/11/89
C.                                                LAST UPDATE: 12/04/90
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      INTEGER INIT
      SAVE INIT
      DATA INIT/ 0/
C--
C--   RETURN IF ALREADY INITIALIZED...
      IF (INIT.NE.0) RETURN
      INIT=1
C--
C--   PRESET PARAMETERS IN PHOTOS COMMONS
      CALL PHOCIN
C--
C--   PRINT INFO
      CALL PHOINF
C--
C--   INITIALIZE MARSAGLIA AND ZAMAN RANDOM NUMBER GENERATOR
      CALL PHORIN
      RETURN
      END
*CMZ :  1.01/50 23/05/96  10.22.20  by  Piero Zucchelli
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOMAK(IPPAR,NHEP0)
C.----------------------------------------------------------------------
C.
C.    PHOMAK:   PHOTOS MAKE
C.
C.    PURPOSE:  SINGLE OR DOUBLE BREMSTRAHLUNG RADIATIVE CORRECTIONS
C.              ARE GENERATED IN  THE DECAY OF THE IPPAR-TH PARTICLE IN
C.              THE  HEP COMMON /HEPEVT/. EXAMPLE OF THE USE OF
C.              GENERAL TOOLS.
C.
C.    INPUT PARAMETER:    IPPAR:  POINTER   TO   DECAYING  PARTICLE  IN
C.                                /HEPEVT/ AND THE COMMON ITSELF
C.
C.    OUTPUT PARAMETERS:  COMMON  /HEPEVT/, EITHER  WITH  OR  WITHOUT A
C.                                PARTICLES ADDED.
C.
C.    AUTHOR(S):  Z. WAS,                         CREATED AT:  26/05/93
C.                                                LAST UPDATE:
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION DATA
      REAL PHORAN
      INTEGER IP,IPPAR,NCHARG
      INTEGER WTDUM,IDUM,NHEP0
      INTEGER NCHARB,NEUDAU
      REAL RN,WT,PHINT
      LOGICAL BOOST
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      DOUBLE PRECISION PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     +JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      LOGICAL INTERF,ISEC,IFTOP
      REAL FINT,FSEC
      COMMON /PHOKEY/ INTERF,FINT,ISEC,FSEC,IFTOP
C--
      IP=IPPAR
      IDUM=1
      NCHARG=0
C--
      CALL PHOIN(IP,BOOST,NHEP0)
      CALL PHOCHK(JDAHEP(1,IP))
      WT=0.0
      CALL PHOPRE(1,WT,NEUDAU,NCHARB)
      IF (WT.EQ.0.0) RETURN
      RN=PHORAN(WTDUM)
C PHODO IS CALLING PHORAN, THUS CHANGE OF SERIES IF IT IS MOVED BEFORE IF.
      CALL PHODO(1,NCHARB,NEUDAU)
      IF (INTERF) WT=WT*PHINT(IDUM)/FINT
      DATA=WT
      IF (WT.GT.1.0) CALL PHOERR(3,'WT_INT',DATA)
      IF (RN.LE.WT) THEN
        CALL PHOOUT(IP,BOOST,NHEP0)
      ENDIF
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOOMA(IFIRST,ILAST,POINTR)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS ORDER MASS VECTOR
C.
C.    PURPOSE:  ORDER  THE  CONTENTS  OF ARRAY 'POINTR' ACCORDING TO THE
C.              DECREASING VALUE IN THE ARRAY 'MASS'.
C.
C.    INPUT PARAMETERS:  IFIRST, ILAST:  POINTERS  TO  THE  VECTOR LOCA-
C.                                       TION BE SORTED,
C.                       POINTR:         UNSORTED ARRAY WITH POINTERS TO
C.                                       /PHOEVT/.
C.
C.    OUTPUT PARAMETER:  POINTR:         SORTED ARRAYS  WITH  RESPECT TO
C.                                       PARTICLE MASS 'PPHO(5,*)'.
C.
C.    AUTHOR(S):  B. VAN EIJK                     CREATED AT:  28/11/89
C.                                                LAST UPDATE: 27/05/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      INTEGER NMXPHO
      PARAMETER (NMXPHO=2000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      INTEGER IFIRST,ILAST,I,J,BUFPOI,POINTR(NMXPHO)
      REAL BUFMAS,MASS(NMXPHO)
      IF (IFIRST.EQ.ILAST) RETURN
C--
C--   COPY PARTICLE MASSES
      DO 10 I=IFIRST,ILAST
   10 MASS(I)=PPHO(5,POINTR(I))
C--
C--   ORDER THE MASSES IN A DECREASING SERIES
      DO 30 I=IFIRST,ILAST-1
        DO 20 J=I+1,ILAST
          IF (MASS(J).LE.MASS(I)) GOTO 20
          BUFPOI=POINTR(J)
          POINTR(J)=POINTR(I)
          POINTR(I)=BUFPOI
          BUFMAS=MASS(J)
          MASS(J)=MASS(I)
          MASS(I)=BUFMAS
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
*CMZ :  1.01/50 23/05/96  10.22.20  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOOUT(IP,BOOST,NHEP0)
C.----------------------------------------------------------------------
C.
C.    PHOOUT:   PHOTOS OUTPUT
C.
C.    PURPOSE:  COPIES BACK IP BRANCH OF THE COMMON /HEPEVT/ FROM /PHOEVT/
C.              MOVES BRANCH BACK FROM ITS CMS SYSTEM.
C.
C.    INPUT PARAMETERS:       IP:  POINTER OF PARTICLE STARTING BRANCH
C.                                 TO BE GIVEN BACK.
C.                        BOOST:   FLAG WHETHER BOOST TO CMS WAS OR WAS
C     .                            NOT PERFORMED.
C.
C.    OUTPUT PARAMETERS:  COMMON /PHOEVT/,
C.
C.    AUTHOR(S):  Z. WAS                          CREATED AT:  24/05/93
C.                                                LAST UPDATE:
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      DOUBLE PRECISION PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     +JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      INTEGER NMXPHO
      PARAMETER (NMXPHO=2000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     +JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      INTEGER IP,LL,FIRST,LAST,I
      LOGICAL BOOST
      INTEGER NN,J,K,NHEP0,NA
      DOUBLE PRECISION BET(3),GAM,PB
      COMMON /PHOCMS/ BET,GAM
      IF(NPHO.EQ.NEVPHO) RETURN
C--   WHEN PARENT WAS NOT IN ITS REST-FRAME, BOOST BACK...
      CALL PHLUPA(10)
      IF (BOOST) THEN
        DO 20  J=JDAPHO(1,1),JDAPHO(2,1)
          PB=-BET(1)*PPHO(1,J)-BET(2)*PPHO(2,J)-BET(3)*PPHO(3,J)
          DO 10  K=1,3
   10     PPHO(K,J)=PPHO(K,J)-BET(K)*(PPHO(4,J)+PB/(GAM+1.))
   20   PPHO(4,J)=GAM*PPHO(4,J)+PB
C--   ...BOOST PHOTON, OR WHATEVER ELSE HAS SHOWN UP
        DO NN=NEVPHO+1,NPHO
          PB=-BET(1)*PPHO(1,NN)-BET(2)*PPHO(2,NN)-BET(3)*PPHO(3,NN)
          DO 30  K=1,3
   30     PPHO(K,NN)=PPHO(K,NN)-BET(K)*(PPHO(4,NN)+PB/(GAM+1.))
          PPHO(4,NN)=GAM*PPHO(4,NN)+PB
        ENDDO
      ENDIF
      FIRST=JDAHEP(1,IP)
      LAST =JDAHEP(2,IP)
C LET'S TAKE IN ORIGINAL DAUGHTERS
      DO LL=0,LAST-FIRST
        IDHEP(FIRST+LL) = IDPHO(3+LL)
        DO I=1,5
          PHEP(I,FIRST+LL) = PPHO(I,3+LL)
        ENDDO
      ENDDO
C LET'S TAKE NEWCOMERS TO THE END OF HEPEVT.
      NA=3+LAST-FIRST
      DO LL=1,NPHO-NA
        IDHEP(NHEP0+LL) = IDPHO(NA+LL)
        ISTHEP(NHEP0+LL)=ISTPHO(NA+LL)
        JMOHEP(1,NHEP0+LL)=IP
        JMOHEP(2,NHEP0+LL)=JMOHEP(2,JDAHEP(1,IP))
        JDAHEP(1,NHEP0+LL)=0
        JDAHEP(2,NHEP0+LL)=0
        DO I=1,5
          PHEP(I,NHEP0+LL) = PPHO(I,NA+LL)
        ENDDO
      ENDDO
      NHEP=NHEP+NPHO-NEVPHO
      CALL PHLUPA(20)
      END
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOPRE(IPARR,WT,NEUDAU,NCHARB)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS
C.
C.    PURPOSE:  ORDER (ALPHA) RADIATIVE CORRECTIONS  ARE  GENERATED  IN
C.              THE DECAY OF THE IPPAR-TH PARTICLE IN THE HEP-LIKE
C.              COMMON /PHOEVT/.  PHOTON RADIATION TAKES PLACE FROM ONE
C.              OF THE CHARGED DAUGHTERS OF THE DECAYING PARTICLE IPPAR
C.              WT IS CALCULATED, EVENTUAL REJECTION WILL BE PERFORMED
C.              LATER AFTER INCLUSION OF INTERFERENCE WEIGHT.
C.
C.    INPUT PARAMETER:    IPPAR:  POINTER   TO   DECAYING  PARTICLE  IN
C.                                /PHOEVT/ AND THE COMMON ITSELF,
C.
C.    OUTPUT PARAMETERS:  COMMON  /PHOEVT/, EITHER  WITH  OR  WITHOUT A
C.                                PHOTON(S) ADDED.
C.                        WT      WEIGHT OF THE CONFIGURATION
C.
C.    AUTHOR(S):  Z. WAS, B. VAN EIJK             CREATED AT:  26/11/89
C.                                                LAST UPDATE: 26/05/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION MINMAS,MPASQR,MCHREN
      DOUBLE PRECISION BETA,EPS,DEL1,DEL2,DATA
      REAL PHOCHA,PHOSPI,PHORAN,PHOCOR,MASSUM
      INTEGER IP,IPARR,IPPAR,I,J,ME,NCHARG,NEUPOI,NLAST,THEDUM
      INTEGER IDABS,IDUM
      INTEGER NCHARB,NEUDAU
      REAL WT
      INTEGER NMXPHO
      PARAMETER (NMXPHO=2000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     +JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      LOGICAL CHKIF
      COMMON/PHOIF/CHKIF(NMXPHO)
      INTEGER CHAPOI(NMXPHO)
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL PNEUTR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
      REAL ALPHA,XPHCUT
      COMMON/PHOCOP/ALPHA,XPHCUT
      INTEGER IREP
      REAL PROBH,CORWT,XF
      COMMON/PHOPRO/IREP,PROBH,CORWT,XF
C--
      IPPAR=IPARR
C--   STORE POINTERS FOR CASCADE TREATEMENT...
      IP=IPPAR
      NLAST=NPHO
      IDUM=1
C--
C--   CHECK DECAY MULTIPLICITY..
      IF (JDAPHO(1,IP).EQ.0) RETURN
C--
C--   LOOP OVER DAUGHTERS, DETERMINE CHARGE MULTIPLICITY
   10 NCHARG=0
      IREP=0
      MINMAS=0.
      MASSUM=0.
      DO 20 I=JDAPHO(1,IP),JDAPHO(2,IP)
C--
C--
C--   EXCLUDE MARKED PARTICLES, QUARKS AND GLUONS ETC...
        IDABS=ABS(IDPHO(I))
        IF (CHKIF(I-JDAPHO(1,IP)+3)) THEN
          IF (PHOCHA(IDPHO(I)).NE.0) THEN
            NCHARG=NCHARG+1
            IF (NCHARG.GT.NMXPHO) THEN
              DATA=NCHARG
              CALL PHOERR(1,'PHOTOS',DATA)
            ENDIF
            CHAPOI(NCHARG)=I
          ENDIF
          MINMAS=MINMAS+PPHO(5,I)**2
        ENDIF
        MASSUM=MASSUM+PPHO(5,I)
   20 CONTINUE
      IF (NCHARG.NE.0) THEN
C--
C--   CHECK THAT SUM OF DAUGHTER MASSES DOES NOT EXCEED PARENT MASS
        IF ((PPHO(5,IP)-MASSUM)/PPHO(5,IP).GT.2.*XPHCUT) THEN
C--
C--   ORDER  CHARGED  PARTICLES  ACCORDING  TO DECREASING MASS, THIS  TO
C--   INCREASE EFFICIENCY (SMALLEST MASS IS TREATED FIRST).
          IF (NCHARG.GT.1) CALL PHOOMA(1,NCHARG,CHAPOI)
C--
   30     CONTINUE
          DO 40 J=1,3
   40     PNEUTR(J)=-PPHO(J,CHAPOI(NCHARG))
          PNEUTR(4)=PPHO(5,IP)-PPHO(4,CHAPOI(NCHARG))
C--
C--   CALCULATE  INVARIANT  MASS OF 'NEUTRAL' ETC. SYSTEMS
          MPASQR=PPHO(5,IP)**2
          MCHSQR=PPHO(5,CHAPOI(NCHARG))**2
          IF ((JDAPHO(2,IP)-JDAPHO(1,IP)).EQ.1) THEN
            NEUPOI=JDAPHO(1,IP)
            IF (NEUPOI.EQ.CHAPOI(NCHARG)) NEUPOI=JDAPHO(2,IP)
            MNESQR=PPHO(5,NEUPOI)**2
            PNEUTR(5)=PPHO(5,NEUPOI)
          ELSE
            MNESQR=PNEUTR(4)**2-PNEUTR(1)**2-PNEUTR(2)**2-PNEUTR(3)**2
            MNESQR=MAX(MNESQR,MINMAS-MCHSQR)
            PNEUTR(5)=SQRT(MNESQR)
          ENDIF
C--
C--   DETERMINE KINEMATICAL LIMIT...
          XPHMAX=(MPASQR-(PNEUTR(5)+PPHO(5,CHAPOI(NCHARG)))**2)/MPASQR
C--
C--   PHOTON ENERGY FRACTION...
          CALL PHOENE(MPASQR,MCHREN,BETA,IDPHO(CHAPOI(NCHARG)))
C--
C--   ENERGY FRACTION NOT TOO LARGE (VERY SELDOM) ? DEFINE ANGLE.
          IF ((XPHOTO.LT.XPHCUT).OR.(XPHOTO.GT.XPHMAX)) THEN
C--
C--   NO RADIATION WAS ACCEPTED, CHECK  FOR MORE DAUGHTERS  THAT MAY RA-
C--   DIATE AND CORRECT RADIATION PROBABILITY...
            NCHARG=NCHARG-1
            IF (NCHARG.GT.0) THEN
              IREP=IREP+1
              GOTO 30
            ENDIF
          ELSE
C--
C--   ANGLE IS GENERATED  IN  THE  FRAME DEFINED  BY  CHARGED VECTOR AND
C--   PNEUTR, DISTRIBUTION IS TAKEN IN THE INFRARED LIMIT...
            EPS=MCHREN/(1.+BETA)
C--
C--   CALCULATE SIN(THETA) AND COS(THETA) FROM INTERVAL VARIABLES
            DEL1=(2.-EPS)*(EPS/(2.-EPS))**PHORAN(THEDUM)
            DEL2=2.-DEL1
            COSTHG=(1.-DEL1)/BETA
            SINTHG=SQRT(DEL1*DEL2-MCHREN)/BETA
C--
C--   DETERMINE SPIN OF  PARTICLE AND CONSTRUCT CODE  FOR MATRIX ELEMENT
            ME=2.*PHOSPI(IDPHO(CHAPOI(NCHARG)))+1.
C--
C--   WEIGHTING PROCEDURE WITH 'EXACT' MATRIX ELEMENT, RECONSTRUCT KINE-
C--   MATICS FOR PHOTON, NEUTRAL AND CHARGED SYSTEM AND UPDATE /PHOEVT/.
C--   FIND POINTER TO THE FIRST COMPONENT OF 'NEUTRAL' SYSTEM
            DO I=JDAPHO(1,IP),JDAPHO(2,IP)
              IF (I.NE.CHAPOI(NCHARG)) THEN
                NEUDAU=I
                GOTO 50
              ENDIF
            ENDDO
C--
C--   POINTER NOT FOUND...
            DATA=NCHARG
            CALL PHOERR(5,'PHOKIN',DATA)
   50       CONTINUE
            NCHARB=CHAPOI(NCHARG)
            NCHARB=NCHARB-JDAPHO(1,IP)+3
            NEUDAU=NEUDAU-JDAPHO(1,IP)+3
            WT=PHOCOR(MPASQR,MCHREN,ME)

          ENDIF
        ELSE
          DATA=PPHO(5,IP)-MASSUM
          CALL PHOERR(10,'PHOTOS',DATA)
        ENDIF
      ENDIF
C--
      RETURN
      END
*CMZ :  1.01/50 23/05/96  12.34.50  by  Piero Zucchelli
*-- Author :    Piero Zucchelli   20/03/96

      REAL*4 FUNCTION PHORAN(IDUMMY)
      CALL RANLUX(RTIM,1)
      PHORAN=RTIM
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOREP
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS RUN SUMMARY REPORT
C.
C.    PURPOSE:  INFORM USER ABOUT SUCCESS AND/OR RESTRICTIONS OF PHOTOS
C.              ENCOUNTERED DURING EXECUTION.
C.
C.    INPUT PARAMETERS:   COMMON /PHOSTA/
C.
C.    OUTPUT PARAMETERS:  NONE
C.
C.    AUTHOR(S):  B. VAN EIJK                     CREATED AT:  10/01/92
C.                                                LAST UPDATE: 10/01/92
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      INTEGER PHOMES
      PARAMETER (PHOMES=10)
      INTEGER STATUS
      COMMON/PHOSTA/STATUS(PHOMES)
      INTEGER I
      LOGICAL ERROR
      ERROR=.FALSE.
      WRITE(PHLUN,10000)
      WRITE(PHLUN,10100)
      WRITE(PHLUN,10200)
      WRITE(PHLUN,10300)
      WRITE(PHLUN,10400)
      WRITE(PHLUN,10300)
      WRITE(PHLUN,10200)
      DO 10 I=1,PHOMES
        IF (STATUS(I).EQ.0) GOTO 10
        IF ((I.EQ.6).OR.(I.EQ.10)) THEN
          WRITE(PHLUN,10500) I,STATUS(I)
        ELSE
          ERROR=.TRUE.
          WRITE(PHLUN,10600) I,STATUS(I)
        ENDIF
   10 CONTINUE
      IF (.NOT.ERROR) WRITE(PHLUN,10700)
      WRITE(PHLUN,10200)
      WRITE(PHLUN,10100)
      RETURN
10000 FORMAT(1H1)
10100 FORMAT(1H ,80('*'))
10200 FORMAT(1H ,'*',T81,'*')
10300 FORMAT(1H ,'*',26X,25('='),T81,'*')
10400 FORMAT(1H ,'*',30X,'PHOTOS RUN SUMMARY',T81,'*')
10500 FORMAT(1H ,'*',22X,'WARNING #',I2,' OCCURED',I6,' TIMES',T81,'*')
10600 FORMAT(1H ,'*',23X,'ERROR #',I2,' OCCURED',I6,' TIMES',T81,'*')
10700 FORMAT(1H ,'*',16X,'PHOTOS EXECUTION HAS SUCCESSFULLY TERMINATED',
     &T81,'*')
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHORIN
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION  IN DECAYS RANDOM NUMBER GENERATOR INIT
C.
C.    PURPOSE:  INITIALSE PHORAN  WITH  THE USER  SPECIFIED SEEDS IN THE
C.              ARRAY ISEED.  FOR DETAILS  SEE ALSO:  F. JAMES  CERN DD-
C.              REPORT NOVEMBER 1988.
C.
C.    INPUT PARAMETERS:   ISEED(*)
C.
C.    OUTPUT PARAMETERS:  URAN, CRAN, CDRAN, CMRAN, I97, J97
C.
C.    AUTHOR(S):  B. VAN EIJK AND F. JAMES        CREATED AT:  27/09/89
C.                                                LAST UPDATE: 22/02/90
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION DATA
      REAL S,T
      INTEGER I,IS1,IS2,IS3,IS4,IS5,J
      INTEGER ISEED,I97,J97
      REAL URAN,CRAN,CDRAN,CMRAN
      COMMON/PHSEED/ISEED(2),I97,J97,URAN(97),CRAN,CDRAN,CMRAN
C--
C--   CHECK VALUE RANGE OF SEEDS
      IF ((ISEED(1).LT.0).OR.(ISEED(1).GE.31328)) THEN
        DATA=ISEED(1)
        CALL PHOERR(8,'PHORIN',DATA)
      ENDIF
      IF ((ISEED(2).LT.0).OR.(ISEED(2).GE.30081)) THEN
        DATA=ISEED(2)
        CALL PHOERR(9,'PHORIN',DATA)
      ENDIF
C--
C--   CALCULATE MARSAGLIA AND ZAMAN SEEDS (BY F. JAMES)
      IS1=MOD(ISEED(1)/177,177)+2
      IS2=MOD(ISEED(1),177)+2
      IS3=MOD(ISEED(2)/169,178)+1
      IS4=MOD(ISEED(2),169)
      DO 20 I=1,97
        S=0.
        T=0.5
        DO 10 J=1,24
          IS5=MOD (MOD(IS1*IS2,179)*IS3,179)
          IS1=IS2
          IS2=IS3
          IS3=IS5
          IS4=MOD(53*IS4+1,169)
          IF (MOD(IS4*IS5,64).GE.32) S=S+T
   10   T=0.5*T
   20 URAN(I)=S
      CRAN=362436./16777216.
      CDRAN=7654321./16777216.
      CMRAN=16777213./16777216.
      I97=97
      J97=33
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHORO2(ANGLE,PVEC)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS ROTATION ROUTINE '2'
C.
C.    PURPOSE:  ROTATE  X AND Z COMPONENTS  OF VECTOR PVEC  AROUND ANGLE
C.              'ANGLE'.
C.
C.    INPUT PARAMETERS:  ANGLE, PVEC
C.
C.    OUTPUT PARAMETER:  PVEC
C.
C.    AUTHOR(S):  S. JADACH                       CREATED AT:  01/01/89
C.                B. VAN EIJK                     LAST UPDATE: 02/01/90
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION CS,SN,ANGLE
      REAL PVEC(4)
      CS=COS(ANGLE)*PVEC(1)+SIN(ANGLE)*PVEC(3)
      SN=-SIN(ANGLE)*PVEC(1)+COS(ANGLE)*PVEC(3)
      PVEC(1)=CS
      PVEC(3)=SN
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHORO3(ANGLE,PVEC)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS ROTATION ROUTINE '3'
C.
C.    PURPOSE:  ROTATE  X AND Y COMPONENTS  OF VECTOR PVEC  AROUND ANGLE
C.              'ANGLE'.
C.
C.    INPUT PARAMETERS:  ANGLE, PVEC
C.
C.    OUTPUT PARAMETER:  PVEC
C.
C.    AUTHOR(S):  S. JADACH                       CREATED AT:  01/01/89
C.                B. VAN EIJK                     LAST UPDATE: 02/01/90
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION CS,SN,ANGLE
      REAL PVEC(4)
      CS=COS(ANGLE)*PVEC(1)-SIN(ANGLE)*PVEC(2)
      SN=SIN(ANGLE)*PVEC(1)+COS(ANGLE)*PVEC(2)
      PVEC(1)=CS
      PVEC(2)=SN
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION PHOSPI(IDHEP)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION  IN DECAYS FUNCTION FOR SPIN DETERMINA-
C.              TION
C.
C.    PURPOSE:  CALCULATE  THE SPIN  OF PARTICLE  WITH  CODE IDHEP.  THE
C.              CODE  OF THE PARTICLE  IS  DEFINED  BY THE PARTICLE DATA
C.              GROUP IN PHYS. LETT. B204 (1988) 1.
C.
C.    INPUT PARAMETER:   IDHEP
C.
C.    OUTPUT PARAMETER:  FUNTION  VALUE = SPIN  OF  PARTICLE  WITH  CODE
C.                       IDHEP
C.
C.    AUTHOR(S):  E. BARBERIO AND B. VAN EIJK     CREATED AT:  29/11/89
C.                                                LAST UPDATE: 02/01/90
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      REAL PHOSPI
      INTEGER IDHEP,IDABS
C--
C--   ARRAY 'SPIN' CONTAINS THE SPIN  OF  THE FIRST 100 PARTICLES ACCOR-
C--   DING TO THE PDG PARTICLE CODE...
      REAL SPIN(100)
      DATA SPIN/ 8*.5, 1., 0., 8*.5, 2*0., 4*1., 76*0./
      IDABS=ABS(IDHEP)
C--
C--   SPIN OF QUARK, LEPTON, BOSON ETC....
      IF (IDABS.LE.100) THEN
        PHOSPI=SPIN(IDABS)
      ELSE
C--
C--   ...OTHER PARTICLES, HOWEVER...
        PHOSPI=(MOD(IDABS,10)-1.)/2.
C--
C--   ...K_SHORT AND K_LONG ARE SPECIAL !!
        PHOSPI=MAX(PHOSPI,0.)
      ENDIF
      RETURN
      END
*CMZ :  1.01/50 23/05/96  10.22.20  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOTOS(IPARR)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   GENERAL SEARCH ROUTINE
C.
C.    PURPOSE:  SEARCH THROUGH THE /HEPEVT/ STANDARD HEP COMMON,  STAR-
C.              TING  FROM  THE IPPAR-TH  PARTICLE.  WHENEVR  BRANCHING
C     .         POINT IS FOUND ROUTINE PHTYPE(IP) IS CALLED.
C.              FINALLY IF CALLS ON PHTYPE(IP) MODIFIED ENTRIES, COMMON
C               /HEPEVT/ IS ORDERED.
C.
C.    INPUT PARAMETER:    IPPAR:  POINTER   TO   DECAYING  PARTICLE  IN
C.                                /HEPEVT/ AND THE COMMON ITSELF,
C.
C.    OUTPUT PARAMETERS:  COMMON  /HEPEVT/, EITHER WITH OR WITHOUT  NEW
C.                                PARTICLES ADDED.
C.
C.    AUTHOR(S):  Z. WAS, B. VAN EIJK             CREATED AT:  26/11/89
C.                                                LAST UPDATE: 30/08/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      REAL PHOTON(5)
      INTEGER IP,IPARR,IPPAR,I,J,K,L,NLAST
      DOUBLE PRECISION DATA
      INTEGER MOTHER,POSPHO
      LOGICAL CASCAD
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      DOUBLE PRECISION PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     +JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      LOGICAL QEDRAD
      COMMON/PHOQED/QEDRAD(NMXHEP)
      INTEGER NMXPHO
      PARAMETER (NMXPHO=2000)
      INTEGER ISTACK(0:NMXPHO),NUMIT,NTRY,KK,LL,II,NA,FIRST,LAST
      INTEGER FIRSTA,LASTA,IPP,IDA1,IDA2,MOTHER2,IDPHO,ISPHO
      REAL PORIG(5,NMXPHO)
C--
      IPPAR=ABS(IPARR)
C--   STORE POINTERS FOR CASCADE TREATEMENT...
      IP=IPPAR
      NLAST=NHEP
      CASCAD=.FALSE.
C--
C--   CHECK DECAY MULTIPLICITY AND MINIMUM OF CORRECTNESS..
      IF ((JDAHEP(1,IP).EQ.0).OR.(JMOHEP(1,JDAHEP(1,IP)).NE.IP)) RETURN
C--
C-- SINGLE BRANCH MODE
C-- WE START LOOKING FOR THE DECAY POINTS IN THE CASCADE
C-- IPPAR IS ORIGINAL POSITION WHERE THE PROGRAM WAS CALLED
      ISTACK(0)=IPPAR
C--   NUMIT DENOTES NUMBER OF SECONDARY DECAY BRANCHES
      NUMIT=0
C--   NTRY DENOTES NUMBER OF SECONDARY BRANCHES ALREADY CHECKED FOR
C--        FOR EXISTENCE OF FURTHER BRANCHES
      NTRY=0
C-- LET'S SEARCH IF IPARR DOES NOT PREVENT SEARCHING.
      IF (IPARR.GT.0)  THEN
   10   CONTINUE
        DO I=JDAHEP(1,IP),JDAHEP(2,IP)
          IF (JDAHEP(1,I).NE.0.AND.JMOHEP(1,JDAHEP(1,I)).EQ.I) THEN
            NUMIT=NUMIT+1
            IF (NUMIT.GT.NMXPHO) THEN
              DATA=NUMIT
              CALL PHOERR(7,'PHOTOS',DATA)
            ENDIF
            ISTACK(NUMIT)=I
          ENDIF
        ENDDO
        IF(NUMIT.GT.NTRY) THEN
          NTRY=NTRY+1
          IP=ISTACK(NTRY)
          GOTO 10
        ENDIF
      ENDIF
C-- LET'S DO GENERATION
      DO 20 KK=0,NUMIT
        NA=NHEP
        FIRST=JDAHEP(1,ISTACK(KK))
        LAST=JDAHEP(2,ISTACK(KK))
        DO II=1,LAST-FIRST+1
          DO LL=1,5
            PORIG(LL,II)=PHEP(LL,FIRST+II-1)
          ENDDO
        ENDDO
C--
        CALL PHTYPE(ISTACK(KK))
C--
C--  CORRECT ENERGY/MOMENTUM OF CASCADE DAUGHTERS
        IF(NHEP.GT.NA) THEN
          DO II=1,LAST-FIRST+1
            IPP=FIRST+II-1
            FIRSTA=JDAHEP(1,IPP)
            LASTA=JDAHEP(2,IPP)
            IF(JMOHEP(1,IPP).EQ.ISTACK(KK)) CALL PHOBOS(IPP,PORIG(1,II)
     +      ,PHEP(1,IPP),FIRSTA,LASTA)
          ENDDO
        ENDIF
   20 CONTINUE
C--
C--   REARRANGE  /HEPEVT/  TO GET CORRECT ORDER..
      IF (NHEP.GT.NLAST) THEN
        DO 100 I=NLAST+1,NHEP
C--
C--   PHOTON MOTHER AND POSITION...
          MOTHER=JMOHEP(1,I)
          POSPHO=JDAHEP(2,MOTHER)+1
C--   INTERMEDIATE SAVE OF PHOTON ENERGY/MOMENTUM AND POINTERS
          DO 30 J=1,5
   30     PHOTON(J)=PHEP(J,I)
          ISPHO =ISTHEP(I)
          IDPHO =IDHEP(I)
          MOTHER2 =JMOHEP(2,I)
          IDA1 =JDAHEP(1,I)
          IDA2 =JDAHEP(2,I)
C--
C--   EXCLUDE PHOTON IN SEQUENCE !
          IF (POSPHO.NE.NHEP) THEN
C--
C--
C--   ORDER /HEPEVT/
            DO 60  K=I,POSPHO+1,-1
              ISTHEP(K)=ISTHEP(K-1)
              QEDRAD(K)=QEDRAD(K-1)
              IDHEP(K)=IDHEP(K-1)
              DO 40  L=1,2
                JMOHEP(L,K)=JMOHEP(L,K-1)
   40         JDAHEP(L,K)=JDAHEP(L,K-1)
              DO 50  L=1,5
   50         PHEP(L,K)=PHEP(L,K-1)
              DO 60  L=1,4
   60       VHEP(L,K)=VHEP(L,K-1)
C--
C--   CORRECT POINTERS ASSUMING MOST DIRTY /HEPEVT/...
            DO 70  K=1,NHEP
              DO 70  L=1,2
                IF ((JMOHEP(L,K).NE.0).AND.(JMOHEP(L,K).GE. POSPHO))
     +          JMOHEP(L,K)=JMOHEP(L,K)+1
                IF ((JDAHEP(L,K).NE.0).AND.(JDAHEP(L,K).GE. POSPHO))
     +          JDAHEP(L,K)=JDAHEP(L,K)+1
   70       CONTINUE
C--
C--   STORE PHOTON ENERGY/MOMENTUM
            DO 80  J=1,5
   80       PHEP(J,POSPHO)=PHOTON(J)
          ENDIF
C--
C--   STORE POINTERS FOR THE PHOTON...
          JDAHEP(2,MOTHER)=POSPHO
          ISTHEP(POSPHO)=ISPHO
          IDHEP(POSPHO)=IDPHO
          JMOHEP(1,POSPHO)=MOTHER
          JMOHEP(2,POSPHO)=MOTHER2
          JDAHEP(1,POSPHO)=IDA1
          JDAHEP(2,POSPHO)=IDA2
C--
C--   GET PHOTON PRODUCTION VERTEX POSITION
          DO 90  J=1,4
   90     VHEP(J,POSPHO)=VHEP(J,POSPHO-1)
  100   CONTINUE
      ENDIF
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION PHOTRI(A,B,C)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOTON RADIATION IN DECAYS CALCULATION OF TRIANGLE FIE
C.
C.    PURPOSE:  CALCULATION OF TRIANGLE FUNCTION FOR PHASE SPACE.
C.
C.    INPUT PARAMETERS:  A, B, C (VIRTUAL) PARTICLE MASSES.
C.
C.    OUTPUT PARAMETER:  FUNCTION VALUE =
C.                       SQRT(LAMBDA(A**2,B**2,C**2))/(2*A)
C.
C.    AUTHOR(S):  B. VAN EIJK                     CREATED AT:  15/11/89
C.                                                LAST UPDATE: 02/01/90
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      DOUBLE PRECISION DA,DB,DC,DAPB,DAMB,DTRIAN
      REAL A,B,C,PHOTRI
      DA=A
      DB=B
      DC=C
      DAPB=DA+DB
      DAMB=DA-DB
      DTRIAN=SQRT((DAMB-DC)*(DAPB+DC)*(DAMB+DC)*(DAPB-DC))
      PHOTRI=DTRIAN/(DA+DA)
      RETURN
      END
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHOTWO(MODE)
C.----------------------------------------------------------------------
C.
C.    PHOTWO:   PHOTOS BUT TWO MOTHERS ALLOWED
C.
C.    PURPOSE:  COMBINES TWO MOTHERS INTO ONE IN /PHOEVT/
C.              NECESSARY EG IN CASE OF G G (Q QBAR) --> T TBAR
C.
C.    INPUT PARAMETERS: COMMON /PHOEVT/ (/PHOCMS/)
C.
C.    OUTPUT PARAMETERS:  COMMON /PHOEVT/, (STORED MOTHERS)
C.
C.    AUTHOR(S):  Z. WAS                          CREATED AT:  5/08/93
C.                                                LAST UPDATE:10/08/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      INTEGER NMXPHO
      PARAMETER (NMXPHO=2000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     +JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      DOUBLE PRECISION BET(3),GAM
      COMMON /PHOCMS/ BET,GAM
      INTEGER I,MODE
      REAL MPASQR
      LOGICAL IFRAD
C LOGICAL IFRAD IS USED TO TAG CASES WHEN TWO MOTHERS MAY BE
C MERGED TO THE SOLE ONE.
C SO FAR USED IN CASE:
C                      1) OF T TBAR PRODUCTION
C
C T TBAR CASE
      IF(MODE.EQ.0) THEN
        IFRAD=(IDPHO(1).EQ.21).AND.(IDPHO(2).EQ.21)
        IFRAD=IFRAD.OR.(IDPHO(1).EQ.-IDPHO(2).AND.ABS(IDPHO(1)).LE.6)
        IFRAD=IFRAD .AND.(ABS(IDPHO(3)).EQ.6).AND.(ABS(IDPHO(4)).EQ.6)
        MPASQR= (PPHO(4,1)+PPHO(4,2))**2-(PPHO(3,1)+PPHO(3,2))**2
     +          -(PPHO(2,1)+PPHO(2,2))**2-(PPHO(1,1)+PPHO(1,2))**2
        IFRAD=IFRAD.AND.(MPASQR.GT.0.0)
        IF(IFRAD) THEN
C.....COMBINING FIRST AND SECOND MOTHER
          DO I=1,4
            PPHO(I,1)=PPHO(I,1)+PPHO(I,2)
          ENDDO
          PPHO(5,1)=SQRT(MPASQR)
C.....REMOVING SECOND MOTHER,
          DO I=1,5
            PPHO(I,2)=0.0
          ENDDO
        ENDIF
      ELSE
C BOOSTING OF THE MOTHERS TO THE REACTION FRAME NOT IMPLEMENTED YET.
C TO DO IT IN MODE 0 ORIGINAL MOTHERS HAVE TO BE STORED IN NEW COMMON (?)
C AND IN MODE 1 BOOSTED TO CMS.
      ENDIF
      END
*CMZ :  1.01/50 23/05/96  10.22.20  by  Piero Zucchelli
*CMZ :  1.00/00 10/08/94  16.28.56  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PHTYPE(ID)
C.----------------------------------------------------------------------
C.
C.    PHTYPE:   CENTRAL MANADGEMENT ROUTINE.
C.
C.    PURPOSE:   DEFINES WHAT KIND OF THE
C.              ACTIONS WILL BE PERFORMED AT POINT ID.
C.
C.    INPUT PARAMETERS:       ID:  POINTER OF PARTICLE STARTING BRANCH
C.                                 IN /HEPEVT/ TO BE TREATED.
C.
C.    OUTPUT PARAMETERS:  COMMON /HEPEVT/.
C.
C.    AUTHOR(S):  Z. WAS                          CREATED AT:  24/05/93
C.                                                LAST UPDATE: 01/10/93
C.
C.----------------------------------------------------------------------
C--   IMPLICIT NONE
      INTEGER NMXHEP
      PARAMETER (NMXHEP=2000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      DOUBLE PRECISION PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      LOGICAL INTERF,ISEC,IFTOP
      REAL FINT,FSEC
      COMMON /PHOKEY/ INTERF,FINT,ISEC,FSEC,IFTOP
      INTEGER ID,NHEP0
      LOGICAL IPAIR
      REAL RN,PHORAN
      INTEGER WTDUM
C--
      IPAIR=.TRUE.
C--   CHECK DECAY MULTIPLICITY..
      IF (JDAHEP(1,ID).EQ.0) RETURN
C      IF (JDAHEP(1,ID).EQ.JDAHEP(2,ID)) RETURN
C--
      NHEP0=NHEP
C--
      IF(ISEC) THEN
C-- DOUBLE PHOTON EMISSION
        FSEC=1.0
        RN=PHORAN(WTDUM)
        IF (RN.GE.0.5) THEN
          CALL PHOMAK(ID,NHEP0)
          CALL PHOMAK(ID,NHEP0)
        ENDIF
      ELSE
C-- SINGLE PHOTON EMISSION
        FSEC=1.0
        CALL PHOMAK(ID,NHEP0)
      ENDIF
C--
C-- ELECTRON POSITRON PAIR (COOMENTED OUT FOR A WHILE
C      IF (IPAIR) CALL PHOPAR(ID,NHEP0)
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE PROD5(P1,P2,P3,PIA)
C ----------------------------------------------------------------------
C EXTERNAL PRODUCT OF P1, P2, P3 4-MOMENTA.
C SIGN IS CHOSEN +/- FOR DECAY OF TAU +/- RESPECTIVELY
C     CALLED BY : DAMPAA, CLNUT
C ----------------------------------------------------------------------
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / IDFC  / IDFF
      REAL PIA(4),P1(4),P2(4),P3(4)
      DET2(I,J)=P1(I)*P2(J)-P2(I)*P1(J)
* -----------------------------------
      IF     (KTOM.EQ.1.OR.KTOM.EQ.-1) THEN
        SIGN= IDFF/ABS(IDFF)
      ELSEIF (KTOM.EQ.2) THEN
        SIGN=-IDFF/ABS(IDFF)
      ELSE
        PRINT *, 'STOP IN PROD5: KTOM=',KTOM
        STOP
      ENDIF
C
C EPSILON( P1(1), P2(2), P3(3), (4) ) = 1
C
      PIA(1)= -P3(3)*DET2(2,4)+P3(4)*DET2(2,3)+P3(2)*DET2(3,4)
      PIA(2)= -P3(4)*DET2(1,3)+P3(3)*DET2(1,4)-P3(1)*DET2(3,4)
      PIA(3)=  P3(4)*DET2(1,2)-P3(2)*DET2(1,4)+P3(1)*DET2(2,4)
      PIA(4)=  P3(3)*DET2(1,2)-P3(2)*DET2(1,3)+P3(1)*DET2(2,3)
C ALL FOUR INDICES ARE UP SO  PIA(3) AND PIA(4) HAVE SAME SIGN
      DO 10 I=1,4
   10 PIA(I)=PIA(I)*SIGN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 25/07/94  19.08.36  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE PYREMM(IPU1,IPU2)

C...ADDS ON TARGET REMNANTS (ONE OR TWO FROM EACH SIDE) AND
C...INCLUDES PRIMORDIAL KT.
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON /PYPROC/ ISUB,KFL(3,2),X(2),SH,TH,UH,Q2,XSEC(0:40)
      DIMENSION KFLCH(2),KFLSP(2),CHI(2),PMS(6),IS(2),ROBO(5)
      DOUBLE PRECISION DBETAX,DBETAZ,DROBO(5)
      DATA IPU,IQ/0,0/,PEI,PE,PZI,PZ,SHS,PZH,PEH/7*0./

C...FIND EVENT TYPE, SET POINTERS
      IF(IPU1.EQ.0.AND.IPU2.EQ.0) RETURN
      ILEP=0
      IF(IPU1.EQ.0) ILEP=1
      IF(IPU2.EQ.0) ILEP=2
      IF(ISUB.EQ.7) ILEP=-1
      IF(ILEP.EQ.1) IQ=21
      IF(ILEP.EQ.2) IQ=23
      IP=MAX(IPU1,IPU2)
      NS=N
C...DEFINE INITIAL PARTONS, INCLUDING PRIMORDIAL KT
   10 DO 30  I=3,4
        IF(I.EQ.3) IPU=IPU1
        IF(I.EQ.4) IPU=IPU2
        K(I,1)=21
        K(I,3)=I-2
        DO 20  J=1,5
   20   P(I,J)=0.
        IF(ISUB.EQ.7) THEN
          K(I,2)=21
          SHS=0.
        ELSEIF(IPU.NE.0) THEN
          K(I,2)=K(IPU,2)
          P(I,5)=P(IPU,5)
          CALL LPRIKT(PARL(3),PTSPL,PHISPL)
          P(I,1)=PTSPL*COS(PHISPL)
          P(I,2)=PTSPL*SIN(PHISPL)
          PMS(I-2)=P(I,5)**2+P(I,1)**2+P(I,2)**2
        ELSE
          K(I,2)=K(IQ,2)
          P(I,5)=-SQRT(Q2)
          PMS(I-2)=-Q2
          SHS=(1.-X(5-I))*Q2/X(5-I)+PYVAR(7-I)**2
        ENDIF
   30 CONTINUE

C...KINEMATICS CONSTRUCTION FOR INITIAL PARTONS
      IF(ILEP.EQ.0) SHS=PYVAR(31)*PYVAR(32)*PYVAR(2)+
     +(P(3,1)+P(4,1))**2+(P(3,2)+P(4,2))**2
      SHR=SQRT(MAX(0.,SHS))
      IF(ILEP.EQ.0) THEN
        IF((SHS-PMS(1)-PMS(2))**2-4.*PMS(1)*PMS(2).LE.0.) GOTO 10
        P(3,4)=0.5*(SHR+(PMS(1)-PMS(2))/SHR)
        P(3,3)=SQRT(MAX(0.,P(3,4)**2-PMS(1)))
        P(4,4)=SHR-P(3,4)
        P(4,3)=-P(3,3)
      ELSEIF(ILEP.EQ.1) THEN
        P(3,4)=P(IQ,4)
        P(3,3)=P(IQ,3)
        P(4,4)=P(IP,4)
        P(4,3)=P(IP,3)
      ELSEIF(ILEP.EQ.2) THEN
        P(3,4)=P(IP,4)
        P(3,3)=P(IP,3)
        P(4,4)=P(IQ,4)
        P(4,3)=P(IQ,3)
      ENDIF

C...TRANSFORM PARTONS TO OVERALL CM-FRAME (NOT FOR LEPTOPRODUCTION)
      IF(ILEP.EQ.0) THEN
        MSTU(1)=3
        MSTU(2)=4
        DROBO(3)=(P(3,1)+P(4,1))/SHR
        DROBO(4)=(P(3,2)+P(4,2))/SHR
* INNOCENT
        CALL LUDBRB(MSTU(1),MSTU(2),0.,0.,-DROBO(3),-DROBO(4),0.D0)
        ROBO(2)=ULANGL(P(3,1),P(3,2))
* INNOCENT
        CALL LUDBRB(MSTU(1),MSTU(2),0.,-ROBO(2),0.D0,0.D0,0.D0)
        ROBO(1)=ULANGL(P(3,3),P(3,1))
* INNOCENT
        CALL LUDBRB(MSTU(1),MSTU(2),-ROBO(1),0.,0.D0,0.D0,0.D0)
        MSTU(2)=MAX(IPY(47),IPU1,IPU2)
* INNOCENT
        CALL LUDBRB(MSTU(1),MSTU(2),
     +  ROBO(1),ROBO(2),DROBO(3),DROBO(4),0.D0)
        DROBO(5)=MAX(-0.999999,MIN(0.999999,(PYVAR(31)-PYVAR(32))/
     +  (PYVAR(31)+PYVAR(32))))
* INNOCENT
        CALL LUDBRB(MSTU(1),MSTU(2),0.,0.,0.D0,0.D0,DROBO(5))
        MSTU(1)=0
        MSTU(2)=0
      ENDIF

C...CHECK INVARIANT MASS OF REMNANT SYSTEM:
C...HADRONIC EVENTS OR LEPTOPRODUCTION
      IF(ILEP.LE.0) THEN
        WRITE(*,*)'ILEP<0!!!!'
        IF(IPY(12).LE.0.OR.ISUB.EQ.7) PYVAR(33)=0.
        IF(IPY(12).LE.0.OR.ISUB.EQ.7) PYVAR(34)=0.
        PEH=P(3,4)+P(4,4)+0.5*PYVAR(1)*(PYVAR(33)+PYVAR(34))
        PZH=P(3,3)+P(4,3)+0.5*PYVAR(1)*(PYVAR(33)-PYVAR(34))
        SHH=(PYVAR(1)-PEH)**2-(P(3,1)+P(4,1))**2-(P(3,2)+P(4,2))**2-
     +  PZH**2
        PMMIN=P(1,5)+P(2,5)+ULMASS(K(3,2))+ULMASS(K(4,2))
        IF(SHR.GE.PYVAR(1).OR.SHH.LE.(PMMIN+PYPAR(12))**2) THEN
          WRITE(*,*)'ERROR 1 IPY(48)'
          IPY(48)=1
          RETURN
        ENDIF
        SHR=SQRT(SHH+(P(3,1)+P(4,1))**2+(P(3,2)+P(4,2))**2)
      ELSE
*        NSAVE=N
*        N=40
*        CALL LULIST(1)
*        N=NSAVE
        PEI=P(IQ,4)+P(IP,4)
        PZI=P(IQ,3)+P(IP,3)
        PMS(ILEP)=MAX(0.,PEI**2-PZI**2+P(5-ILEP,1)**2+P(5-ILEP,2)**2)
        PMMIN=P(3-ILEP,5)+ULMASS(K(5-ILEP,2))+SQRT(PMS(ILEP))
*      WRITE(*,*)'SHR,PMMIN,PYPAR(12),X=',SHR,PMMIN,PYPAR(12),X(1),X(2)
*      WRITE(*,*)'PEI,PZI,PMS,PX,PY',SHR,PMMIN,PYPAR(12),X(1),X(2)
*      WRITE(*,*)'PYVAR1, PYVAR33, PYVAR34',PYVAR(1),PYVAR(33),PYVAR(34)
*      WRITE(*,*)'IQ,IP,IPY(12),ISUB,ILEP=',IQ,IP,IPY(12),ISUB,ILEP
        IF(SHR.LE.PMMIN+PYPAR(12)) THEN
*        WRITE(*,*)'ERROR 2 IPY(48)'
          IPY(48)=1
          RETURN
        ENDIF
      ENDIF

C...SUBDIVIDE REMNANT IF NECESSARY, STORE FIRST PARTON
   40 I=NS-1
      DO 70  JT=1,2
        IF(JT.EQ.ILEP) GOTO 70
        IF(JT.EQ.1) IPU=IPU1
        IF(JT.EQ.2) IPU=IPU2
        CALL PYSPLA(IPY(40+JT),KFL(1,JT),KFLCH(JT),KFLSP(JT))
        I=I+2
        IS(JT)=I
        K(I,1)=3
        K(I,2)=KFLSP(JT)
        K(I,3)=JT
        P(I,5)=ULMASS(K(I,2))
C...FIRST PARTON COLOUR CONNECTIONS AND TRANSVERSE MASS
        K(I+1,1)=-1
        K(I+1,3)=I
        K(I+1,2)=1000
        IF(IPY(34).GE.1) K(I+1,2)=1000+JT
        DO 50  J=1,5
   50   P(I+1,J)=0.
        IF(KFLSP(JT).EQ.21) THEN
          P(I+1,3)=IPU
          P(I+1,4)=IPU
          P(IPU+1,1)=I
          P(IPU+1,2)=I
          K(I,4)=IPU+IPU*MSTU(5)
          K(I,5)=IPU+IPU*MSTU(5)
          K(IPU,4)=MOD(K(IPU,4),MSTU(5))+I*MSTU(5)
          K(IPU,5)=MOD(K(IPU,5),MSTU(5))+I*MSTU(5)
        ELSE
          IFLS=(3-ISIGN(1,KFLSP(JT)*(1102-IABS(KFLSP(JT)))))/2
          P(I+1,IFLS+2)=IPU
          P(IPU+1,3-IFLS)=I
          K(I,IFLS+3)=IPU
          K(IPU,6-IFLS)=MOD(K(IPU,6-IFLS),MSTU(5))+I*MSTU(5)
        ENDIF
        IF(KFLCH(JT).EQ.0) THEN
          P(I,1)=-P(JT+2,1)
          P(I,2)=-P(JT+2,2)
          PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2
        ELSE
C...WHEN EXTRA REMNANT PARTON OR HADRON: FIND RELATIVE PT, STORE
C...PRIMORDIAL KT SPLIT SHARED BETWEEN REMNANTS
          CALL LPRIKT(PARL(14),PTSPL,PHISPL)
C...RELATIVE DISTRIBUTION OF ENERGY; EXTRA PARTON COLOUR CONNECTION
          CALL LREMH(0,KFLSP(JT),KFLCH(JT),CHI(JT))
          P(I,1)=-P(JT+2,1)*(1.-CHI(JT))+PTSPL*COS(PHISPL)
          P(I,2)=-P(JT+2,2)*(1.-CHI(JT))+PTSPL*SIN(PHISPL)
          PMS(JT+2)=P(I,5)**2+P(I,1)**2+P(I,2)**2
          I=I+2
          DO 60  J=1,5
            K(I,J)=0
            K(I+1,J)=0
            P(I,J)=0.
   60     P(I+1,J)=0.
          K(I,1)=1
          K(I,2)=KFLCH(JT)
          K(I,3)=JT
          P(I,5)=ULMASS(K(I,2))
          P(I,1)=-P(JT+2,1)*CHI(JT)-PTSPL*COS(PHISPL)
          P(I,2)=-P(JT+2,2)*CHI(JT)-PTSPL*SIN(PHISPL)
          PMS(JT+4)=P(I,5)**2+P(I,1)**2+P(I,2)**2
C...END OF UPDATE
          PMS(JT)=PMS(JT+4)/CHI(JT)+PMS(JT+2)/(1.-CHI(JT))
          K(I+1,1)=-1
          K(I+1,3)=I
          K(I+1,2)=1000
          IF(IPY(34).GE.1) K(I+1,2)=1000+JT
          IF((IABS(KFLCH(JT)).GE.1.AND.IABS(KFLCH(JT)).LE.8).OR.
     +    IABS(KFLCH(JT)).EQ.21.OR.LUCOMP(IABS(KFLCH(JT))).EQ.90) THEN
            IFLS=(3-ISIGN(1,KFLCH(JT)*(1102-IABS(KFLCH(JT)))))/2
            P(I+1,IFLS+2)=IPU
            P(IPU+1,3-IFLS)=I
            K(I,1)=3
            K(I,IFLS+3)=IPU
            K(IPU,6-IFLS)=MOD(K(IPU,6-IFLS),MSTU(5))+I*MSTU(5)
          ELSE
            IF(IPY(34).GE.1) THEN
              K(I,1)=1
              K(I,3)=JT
            ENDIF
          ENDIF
        ENDIF
   70 CONTINUE
      IF(SHR.LE.SQRT(PMS(1))+SQRT(PMS(2))) GOTO 40
      N=I+1

C...RECONSTRUCT KINEMATICS OF REMNANTS
      DO 80  JT=1,2
        IF(JT.EQ.ILEP) GOTO 80
        PE=0.5*(SHR+(PMS(JT)-PMS(3-JT))/SHR)
        PZ=SQRT(PE**2-PMS(JT))
        IF(KFLCH(JT).EQ.0) THEN
          P(IS(JT),4)=PE
          P(IS(JT),3)=PZ*(-1)**(JT-1)
        ELSE
          PW1=CHI(JT)*(PE+PZ)
          P(IS(JT)+2,4)=0.5*(PW1+PMS(JT+4)/PW1)
          P(IS(JT)+2,3)=0.5*(PW1-PMS(JT+4)/PW1)*(-1)**(JT-1)
          P(IS(JT),4)=PE-P(IS(JT)+2,4)
          P(IS(JT),3)=PZ*(-1)**(JT-1)-P(IS(JT)+2,3)
        ENDIF
   80 CONTINUE

C     CALL GULIST(31,2)
C...HADRONIC EVENTS: BOOST REMNANTS TO CORRECT LONGITUDINAL FRAME
      IF(ILEP.LE.0) THEN
        MSTU(1)=NS+1
* INNOCENT
        CALL LUDBRB(MSTU(1),MSTU(2),
     +  0.,0.,0.D0,0.D0,-DBLE(PZH)/(DBLE(PYVAR(1))-DBLE(PEH)))
        MSTU(1)=0
C...LEPTOPRODUCTION EVENTS: BOOST COLLIDING SUBSYSTEM
      ELSE
        IMIN=21
        IMAX=MAX(IP,IPY(47))
        PEF=SHR-PE
        PZF=PZ*(-1)**(ILEP-1)
        PT2=P(5-ILEP,1)**2+P(5-ILEP,2)**2
        PHIPT=ULANGL(P(5-ILEP,1),P(5-ILEP,2))
        CALL LUDBRB(IMIN,IMAX,0.,-PHIPT,0.D0,0.D0,0.D0)
        RQP=P(IQ,3)*(PT2+PEI**2)-P(IQ,4)*PEI*PZI
        SINTH=P(IQ,4)*SQRT(PT2*(PT2+PEI**2)/(RQP**2+PT2*
     +  P(IQ,4)**2*PZI**2))*SIGN(1.,-RQP)
        CALL LUDBRB(IMIN,IMAX,ASIN(SINTH),0.,0.D0,0.D0,0.D0)
        DBETAX=(-DBLE(PEI)*PZI*SINTH+
     +  SQRT(DBLE(PT2)*(PT2+PEI**2-(PZI*SINTH)**2)))/
     +  (DBLE(PT2)+PEI**2)
        CALL LUDBRB(IMIN,IMAX,0.,0.,DBETAX,0.D0,0.D0)
        CALL LUDBRB(IMIN,IMAX,0.,PHIPT,0.D0,0.D0,0.D0)
        PEM=P(IQ,4)+P(IP,4)
        PZM=P(IQ,3)+P(IP,3)
        DBETAZ=(-DBLE(PEM)*PZM+
     +  PZF*SQRT(DBLE(PZF)**2+PEM**2-PZM**2))/(DBLE(PZF)**2+PEM**2)
        CALL LUDBRB(IMIN,IMAX,0.,0.,0.D0,0.D0,DBETAZ)
        CALL LUDBRB(3,4,ASIN(SINTH),0.,DBETAX,0.D0,0.D0)
        CALL LUDBRB(3,4,0.,PHIPT,0.D0,0.D0,DBETAZ)
      ENDIF

      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE PYSPLA(KPART,KFLIN,KFLCH,KFLSP)

C...IN CASE OF A HADRON REMNANT WHICH IS MORE COMPLICATED THAN JUST A
C...QUARK OR A DIQUARK, SPLIT IT INTO TWO (PARTONS OR HADRON + PARTON).
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP

      IFLIN=KFLIN
      KSIGN=ISIGN(1,KPART)
      IFL=KFLIN*KSIGN
      KFLCH=0
      IDUM=0

      IF(LST(14).EQ.0) THEN
C...IF BARYON PRODUCTION FROM REMNANT EXCLUDED, REMNANT IS ANTIFLAVOUR
        KFLSP=-KFLIN
        IF(KFLIN.EQ.21) KFLSP=21
        RETURN
      ENDIF

      IF(IABS(KPART).EQ.211) THEN
C...DECOMPOSE PI+ (PI-).
        IF(IFL.EQ.2) THEN
C...VALENCE U (UBAR) REMOVED.
          KFLSP=-1*KSIGN
        ELSEIF(IFL.EQ.-1) THEN
C...VALENCE D (DBAR) REMOVED.
          KFLSP=2*KSIGN
        ELSEIF(KFLIN.EQ.21) THEN
C...GLUON REMOVED.
          R=2.*RLU(0)
          IF(R.LT.1.) THEN
            KFLCH=2*KSIGN
            KFLSP=-1*KSIGN
          ELSE
            KFLCH=-1*KSIGN
            KFLSP=2*KSIGN
          ENDIF
        ELSEIF((IFL.GE.1.AND.IFL.LE.8).AND.IFL.NE.2) THEN
C...SEA QUARK (ANTIQUARK) REMOVED.
          CALL LUKFDI(-IFLIN,2*KSIGN,IDUM,KFLCH)
          KFLSP=-1*KSIGN
        ELSEIF((IFL.GE.-8.AND.IFL.LE.-1).AND.IFL.NE.-1) THEN
C...SEA ANTIQUARK (QUARK) REMOVED.
          CALL LUKFDI(-IFLIN,-1*KSIGN,IDUM,KFLCH)
          KFLSP=2*KSIGN
        ENDIF

      ELSEIF(IABS(KPART).EQ.2212) THEN
C...DECOMPOSE PROTON (ANTIPROTON).
        IF(IFL.EQ.2) THEN
C...VALENCE U (UBAR) REMOVED.
          R=4.*RLU(0)
          IF(R.LT.3.) THEN
            KFLSP=2101*KSIGN
          ELSE
            KFLSP=2103*KSIGN
          ENDIF
        ELSEIF(IFL.EQ.1) THEN
C...VALENCE D (DBAR) REMOVED.
          KFLSP=2203*KSIGN
        ELSEIF(KFLIN.EQ.21) THEN
C...GLUON REMOVED.
          R=6.*RLU(0)
          IF(R.LT.3.) THEN
            KFLCH=2*KSIGN
            KFLSP=2101*KSIGN
          ELSEIF(R.LT.4.) THEN
            KFLCH=2*KSIGN
            KFLSP=2103*KSIGN
          ELSE
            KFLCH=1*KSIGN
            KFLSP=2203*KSIGN
          ENDIF
        ELSEIF(IFL.GT.2) THEN
C...SEA QUARK (ANTIQUARK) REMOVED.
          R=6*RLU(0)
          IF(R.LT.3.) THEN
            CALL LUKFDI(-IFLIN,2*KSIGN,IDUM,KFLCH)
            KFLSP=2101*KSIGN
          ELSEIF(R.LT.4.) THEN
            CALL LUKFDI(-IFLIN,2*KSIGN,IDUM,KFLCH)
            KFLSP=2103*KSIGN
          ELSE
            CALL LUKFDI(-IFLIN,1*KSIGN,IDUM,KFLCH)
            KFLSP=2203*KSIGN
          ENDIF
        ELSEIF(IFL.LT.0) THEN
C...SEA ANTIQUARK (QUARK) REMOVED.
   10     R=6*RLU(0)
          IF(R.LT.3.) THEN
            CALL LUKFDI(2101*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=2*KSIGN
          ELSEIF(R.LT.4.) THEN
            CALL LUKFDI(2103*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=2*KSIGN
          ELSE
            CALL LUKFDI(2203*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=1*KSIGN
          ENDIF
          IF(KFLCH.EQ.0) GOTO 10
        ENDIF

      ELSEIF(IABS(KPART).EQ.2112) THEN
C...DECOMPOSE NEUTRON (ANTINEUTRON).
        IF(IFL.EQ.2) THEN
C...VALENCE U (UBAR) REMOVED.
          KFLSP=1103*KSIGN
        ELSEIF(IFL.EQ.1) THEN
C...VALENCE D (DBAR) REMOVED.
          R=4.*RLU(0)
          IF(R.LT.3.) THEN
            KFLSP=2101*KSIGN
          ELSE
            KFLSP=2103*KSIGN
          ENDIF
        ELSEIF(KFLIN.EQ.21) THEN
C...GLUON REMOVED.
          R=6.*RLU(0)
          IF(R.LT.2.) THEN
            KFLCH=2*KSIGN
            KFLSP=1103*KSIGN
          ELSEIF(R.LT.5.) THEN
            KFLCH=1*KSIGN
            KFLSP=2101*KSIGN
          ELSE
            KFLCH=1*KSIGN
            KFLSP=2103*KSIGN
          ENDIF
        ELSEIF(IFL.GT.2) THEN
C...SEA QUARK (ANTIQUARK) REMOVED.
          R=6*RLU(0)
          IF(R.LT.2.) THEN
            CALL LUKFDI(-IFLIN,2*KSIGN,IDUM,KFLCH)
            KFLSP=1103*KSIGN
          ELSEIF(R.LT.5.) THEN
            CALL LUKFDI(-IFLIN,1*KSIGN,IDUM,KFLCH)
            KFLSP=2101*KSIGN
          ELSE
            CALL LUKFDI(-IFLIN,1*KSIGN,IDUM,KFLCH)
            KFLSP=2103*KSIGN
          ENDIF
        ELSEIF(IFL.LT.0) THEN
C...SEA ANTIQUARK (QUARK) REMOVED.
   20     R=6*RLU(0)
          IF(R.LT.2.) THEN
            CALL LUKFDI(1103*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=2*KSIGN
          ELSEIF(R.LT.5.) THEN
            CALL LUKFDI(2101*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=1*KSIGN
          ELSE
            CALL LUKFDI(2103*KSIGN,-IFLIN,IDUM,KFLCH)
            KFLSP=1*KSIGN
          ENDIF
          IF(KFLCH.EQ.0) GOTO 20
        ENDIF
      ENDIF

      RETURN
      END
*CMZ :  1.01/50 22/05/96  12.22.19  by  Piero Zucchelli
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  18.01.14  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C ********************************************************************

      SUBROUTINE PYSSPB(IPU1,IPU2)
*KEEP,FOREFI.
C--
	INTEGER*4 IEVT
         COMMON/FOREFICASS/IEVT


*KEND.
C...NEW X REDEFINITION
C...GENERATES SPACELIKE PARTON SHOWERS
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      COMMON /PYPROC/ ISUB,KFL(3,2),X(2),SH,TH,UH,Q2,XSEC(0:40)
      COMMON /PYINT1/ XQ(2,-6:6)
      DIMENSION IFLS(4),IS(2),XS(2),ZS(2),Q2S(2),TEVS(2),ROBO(5),
     +XFS(2,-6:6),XFA(-6:6),XFB(-6:6),WTAP(-6:6),WTSF(-6:6)
      DOUBLE PRECISION DQ2(3),DSH,DSHZ,DSHR,DPLCM,DPC(3),DPD(4),DMS,
     +DMSMA,DPT2,DPB(4),DBE1(4),DBE2(4),DBEP,DGABEP,DPQ(4),DPQS(2),
     +DM2,DQ2B,DROBO(5),DBEZ
C-GI &DQ23,DPH(4),DM2,DQ2B,DQM2
      DATA IFLA,NQ/0,0/,Z,XE0,XA/3*0./,DSHZ,DMSMA,DPT2,DSHR/4*0.D0/

C...COMMON CONSTANTS, SET UP INITIAL VALUES
      ILEP=0
      IF(IPU1.EQ.0) ILEP=1
      IF(IPU2.EQ.0) ILEP=2
      Q2E=Q2
C-GI  IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.3) Q2E=Q2E/PYPAR(26)
      IF(ISUB.EQ.27) Q2E=PMAS(23,1)**2
      IF(ISUB.EQ.28) Q2E=PMAS(24,1)**2
      TMAX=ALOG(PYPAR(26)*PYPAR(27)*Q2E/PYPAR(21)**2)
      IF(ILEP.GE.1) THEN
        SH=P(25,5)**2
        IF(N.GE.27) SH=P(27,5)**2
        CALL LSCALE(-1,QMAX)
        Q2E=QMAX**2
        Q2E=MAX(PYPAR(21)**2,MIN(Q2E,(0.95/X(3-ILEP)-1.)*Q2-SH,
     +  Q2/2.+SH))
        TMAX=ALOG(Q2E/PYPAR(21)**2)
      ENDIF
      IF(PYPAR(26)*Q2E.LT.MAX(PYPAR(22),2.*PYPAR(21)**2).OR.
     +TMAX.LT.0.2) RETURN
      IF(ILEP.EQ.0) XE0=2.*PYPAR(23)/PYVAR(1)
      B0=(33.-2.*IPY(8))/6.
      NS=N
      MSTU(2)=0
   10 N=NS
      IF(ILEP.GE.1) THEN
        NQ=IPU2-2
        IF(ILEP.EQ.2) NQ=IPU1+2
        DPQS(1)=DBLE(P(NQ,3))
        DPQS(2)=DBLE(P(NQ,4))
        XBMIN=X(3-ILEP)*MAX(0.5,SH/Q2)
        CALL PYSTFU(IPY(43-ILEP),XBMIN,Q2,XFB)
        DO 20  IFL=-6,6
   20   XQ(3-ILEP,IFL)=XFB(IFL)
      ENDIF
      DO 30  JT=1,2
        IFLS(JT)=KFL(2,JT)
        IF(KFL(2,JT).EQ.21) IFLS(JT)=0
        IFLS(JT+2)=IFLS(JT)
        XS(JT)=X(JT)
        ZS(JT)=1.
        IF(ILEP.EQ.0) Q2S(JT)=PYPAR(26)*Q2E
        TEVS(JT)=TMAX
        DO 30  IFL=-6,6
   30 XFS(JT,IFL)=XQ(JT,IFL)
      IF(ILEP.GE.1) THEN
        Q2S(ILEP)=P(NQ,5)**2
        DQ2(ILEP)=Q2S(ILEP)
        Q2S(3-ILEP)=Q2E
      ENDIF
      DSH=SH
      IHFC=0
      IHFX=0

C...PICK UP LEG WITH HIGHEST VIRTUALITY
   40 CONTINUE
      IF(N.GT.MSTU(4)-10) THEN
        WRITE(6,*) ' PYSSPB: NO MORE MEMORY IN LUJETS'
        LST(21)=51
        RETURN
      ENDIF
      DO 50  I=N+1,N+8
        DO 50  J=1,5
          K(I,J)=0
   50 P(I,J)=0.
C     CALL GULIST(21,2)
      N=N+2
      JT=1
      IF((N.GT.NS+2.AND.Q2S(2).GT.Q2S(1).AND.ILEP.EQ.0).OR.ILEP.EQ.1)
     +JT=2
      JR=3-JT
      IFLB=IFLS(JT)
      XB=XS(JT)
      IF(ILEP.GE.1.AND.N.EQ.NS+2) XB=XS(JT)*MAX(SH/Q2,0.5)
      DO 60  IFL=-6,6
   60 XFB(IFL)=XFS(JT,IFL)
      Q2B=Q2S(JT)
      TEVB=TEVS(JT)
      IF(IPY(14).GE.9.AND.N.GT.NS+4) THEN
        Q2B=0.5*(1./ZS(JT)+1.)*Q2S(JT)+0.5*(1./ZS(JT)-1.)*(Q2S(3-JT)-
     +  SNGL(DSH)+SQRT((SNGL(DSH)+Q2S(1)+Q2S(2))**2+8.*Q2S(1)*Q2S(2)*
     +  ZS(JT)/(1.-ZS(JT))))
        TEVB=ALOG(PYPAR(27)*Q2B/PYPAR(21)**2)
      ENDIF
      IF(ILEP.EQ.0) THEN
        DSHR=2.*DSQRT(DSH)
        DSHZ=DSH/DBLE(ZS(JT))
      ELSEIF(ILEP.GE.1) THEN
        DSHZ=DSH
        IF(N.GT.NS+4) DSHZ=(DSH+DQ2(JR)-DQ2(JT))/ZS(JT)-DQ2(JR)+
     +  PYPAR(22)
        DPD(2)=DSHZ+DQ2(JR)+DBLE(PYPAR(22))

        QMASS=ULMASS(IABS(IFLB))
        IF(IABS(IFLB).EQ.0) QMASS=ULMASS(21)
C...CHECK IF QUARK PAIR CREATION ONLY POSSIBILITY
        IF(DQ2(JR).LT.4.*QMASS**2) THEN
          DM2=QMASS**2
          DPC(1)=DQ2(JR)*(DBLE(PYPAR(22))+DM2)**2
          DPC(2)=DPD(2)*(DPD(2)-2D0*PYPAR(22))*(PYPAR(22)+DM2)
          DPC(3)=PYPAR(22)*(DPD(2)-2D0*PYPAR(22))**2
          XE0=1D0-(DPC(2)-DSQRT(DPC(2)**2-4D0*DPC(1)*DPC(3)))/
     +    (2D0*DPC(1))
        ELSE
          XE0=1D0-(DPD(2)-2D0*DBLE(PYPAR(22)))*(DPD(2)-DSQRT(DPD(2)**2-
     +    4D0*DQ2(JR)*DBLE(PYPAR(22))))/(2D0*DQ2(JR)*DBLE(PYPAR(22)))
        ENDIF
      ENDIF
   70 XE=MAX(XE0,XB*(1./(1.-PYPAR(24))-1.))
      IF(XB+XE.GE.0.999) THEN
        Q2B=0.
        GOTO 150
      ENDIF

C...CALCULATE ALTARELLI-PARISI AND STRUCTURE FUNCTION WEIGHTS
      DO 80  IFL=-6,6
        WTAP(IFL)=0.
   80 WTSF(IFL)=0.
      IF(IFLB.EQ.0) THEN
        WTAPQ=16.*(1.-SQRT(XB+XE))/(3.*SQRT(XB))
        DO 90  IFL=-IPY(8),IPY(8)
          IF(IFL.EQ.0) WTAP(IFL)=6.*ALOG((1.-XB)/XE)
   90   IF(IFL.NE.0) WTAP(IFL)=WTAPQ
      ELSE
        WTAP(0)=0.5*XB*(1./(XB+XE)-1.)
        WTAP(IFLB)=8.*ALOG((1.-XB)*(XB+XE)/XE)/3.
      ENDIF
  100 WTSUM=0.
      IF(IHFC.EQ.0) THEN
        DO 110 IFL=-IPY(8),IPY(8)
          WTSF(IFL)=XFB(IFL)/MAX(1E-10,XFB(IFLB))
  110   WTSUM=WTSUM+WTAP(IFL)*WTSF(IFL)
        IF(IABS(IFLB).GE.4.AND.WTSUM.GT.1E3) THEN
          IHFX=1
          DO 120 IFL=-IPY(8),IPY(8)
  120     WTSF(IFL)=WTSF(IFL)*1E3/WTSUM
          WTSUM=1E3
        ENDIF
      ENDIF

C...CHOOSE NEW T AND FLAVOUR
  130 IF(IPY(14).LE.6.OR.IPY(14).GE.9) THEN
        TEVXP=B0/MAX(0.0001,WTSUM)
      ELSE
        TEVXP=B0/MAX(0.0001,5.*WTSUM)
      ENDIF
      TEVB=TEVB*EXP(MAX(-100.,ALOG(RLU(0))*TEVXP))
      Q2REF=PYPAR(21)**2*EXP(TEVB)/PYPAR(27)
      Q2B=Q2REF/PYPAR(27)
      DQ2B=Q2B
      IF(ILEP.GE.1) THEN
        DSHZ=DSH
        IF(N.GT.NS+4) DSHZ=(DSH+DQ2(JR)-DQ2(JT))/DBLE(ZS(JT))-DQ2(JR)+
     +  DQ2B
      ENDIF
      IF(Q2B.LT.PYPAR(22)) THEN
        Q2B=0.
      ELSE
        WTRAN=RLU(0)*WTSUM
        IFLA=-IPY(8)-1
  140   IFLA=IFLA+1
        WTRAN=WTRAN-WTAP(IFLA)*WTSF(IFLA)
        IF(IFLA.LT.IPY(8).AND.WTRAN.GT.0.) GOTO 140

C...CHOOSE Z VALUE AND CORRECTIVE WEIGHT
        IF(IFLB.EQ.0.AND.IFLA.EQ.0) THEN
          Z=1./(1.+((1.-XB)/XB)*(XE/(1.-XB))**RLU(0))
          WTZ=(1.-Z*(1.-Z))**2
        ELSEIF(IFLB.EQ.0) THEN
          Z=XB/(1.-RLU(0)*(1.-SQRT(XB+XE)))**2
          WTZ=0.5*(1.+(1.-Z)**2)*SQRT(Z)
        ELSEIF(IFLA.EQ.0) THEN
          Z=XB*(1.+RLU(0)*(1./(XB+XE)-1.))
          WTZ=1.-2.*Z*(1.-Z)
        ELSE
          Z=1.-(1.-XB)*(XE/((XB+XE)*(1.-XB)))**RLU(0)
          WTZ=0.5*(1.+Z**2)
        ENDIF

C...REWEIGHT FIRST LEG BECAUSE OF MODIFIED XB OR CHECK PHASE SPACE
        IF(ILEP.GE.1.AND.N.EQ.NS+2) THEN
          XBNEW=X(JT)*(1.+(DSH-Q2B)/DQ2(JR))
          IF(XBNEW.GT.MIN(Z,0.999)) GOTO 130
          XB=XBNEW
        ENDIF

C...SUM UP SOFT GLUON EMISSION AS EFFECTIVE Z SHIFT
        IF(IPY(15).GE.1) THEN
          RSOFT=6.
          IF(IFLB.NE.0) RSOFT=8./3.
          Z=Z*(TEVB/TEVS(JT))**(RSOFT*XE/((XB+XE)*B0))
          IF(Z.LE.XB) GOTO 130
        ENDIF

C...CHECK IF HEAVY FLAVOUR BELOW THRESHOLD
        IHFT=0
        IF(ILEP.GE.1.AND.IABS(IFLB).GE.4.AND.(XFB(IFLB).LT.1E-10.OR.
     +    Q2B.LT.5.*ULMASS(IABS(IFLB))**2)) THEN
          IHFT=1
          IFLA=0
        ENDIF

C...FOR LEPTOPRODUCTION, CHECK Z AGAINST NEW LIMIT
        IF(ILEP.GE.1) THEN
          DPD(2)=DSHZ+DQ2(JR)+DQ2B
          DM2=ULMASS(IABS(IFLA-IFLB))**2
          IF(IABS(IFLA-IFLB).EQ.0) DM2=ULMASS(21)**2
          DPC(1)=DQ2(JR)*(DQ2B+DM2)**2
          DPC(2)=DPD(2)*(DPD(2)-2D0*DQ2B)*(DQ2B+DM2)
          DPC(3)=DQ2B*(DPD(2)-2D0*DQ2B)**2
          ZU=(DPC(2)-DSQRT(DPC(2)**2-4D0*DPC(1)*DPC(3)))/(2D0*DPC(1))
          IF(Z.GE.ZU) GOTO 130
        ENDIF

C...OPTION WITH EVOLUTION IN KT2=(1-Z)Q2:
        IF(IPY(14).GE.5.AND.IPY(14).LE.6.AND.N.LE.NS+4) THEN
C...CHECK THAT (Q2)LAST BRANCHING < (Q2)HARD
          IF(Q2B/(1.-Z).GT.PYPAR(26)*Q2) GOTO 130
        ELSEIF(IPY(14).GE.3.AND.IPY(14).LE.6.AND.N.GE.NS+6) THEN
C...CHECK THAT Z,Q2 COMBINATION IS KINEMATICALLY ALLOWED
          Q2MAX=0.5*(1./ZS(JT)+1.)*DQ2(JT)+0.5*(1./ZS(JT)-1.)*
     +    (DQ2(3-JT)-DSH+SQRT((DSH+DQ2(1)+DQ2(2))**2+8.*DQ2(1)*DQ2(2)*
     +    ZS(JT)/(1.-ZS(JT))))
          IF(Q2B/(1.-Z).GE.Q2MAX) GOTO 130

        ELSEIF(IPY(14).EQ.7.OR.IPY(14).EQ.8) THEN
C...OPTION WITH ALPHAS((1-Z)Q2): DEMAND KT2 > CUTOFF, REWEIGHT
          IF((1.-Z)*Q2B.LT.PYPAR(22)) GOTO 130
          ALPRAT=TEVB/(TEVB+ALOG(1.-Z))
          IF(ALPRAT.LT.5.*RLU(0)) GOTO 130
          IF(ALPRAT.GT.5.) WTZ=WTZ*ALPRAT/5.
        ENDIF

C...WEIGHTING WITH NEW STRUCTURE FUNCTIONS
        CALL PYSTFU(IPY(40+JT),XB,Q2REF,XFB)
        XA=XB/Z
        CALL PYSTFU(IPY(40+JT),XA,Q2REF,XFA)
        IF(IHFT.EQ.1.OR.IHFX.EQ.1) THEN
          IF(XFA(IFLA).LT.1E-10) IHFC=1
          GOTO 150
        ELSEIF(XFB(IFLB).LT.1E-20) THEN
          GOTO 10
        ENDIF
        IF(WTZ*XFA(IFLA)/XFB(IFLB).LT.RLU(0)*WTSF(IFLA)) THEN
          IF(ILEP.GE.1.AND.N.EQ.NS+2) GOTO 70
          GOTO 100
        ENDIF
      ENDIF

  150 IF(N.EQ.NS+4-2*MIN(1,ILEP)) THEN
C...DEFINE TWO HARD SCATTERERS IN THEIR CM-FRAME
        DQ2(JT)=Q2B
        IF(IPY(14).GE.3.AND.IPY(14).LE.6) DQ2(JT)=Q2B/(1.-Z)
        IF(ILEP.EQ.0) THEN
          DPLCM=DSQRT((DSH+DQ2(1)+DQ2(2))**2-4.*DQ2(1)*DQ2(2))/DSHR
          DO 160 JR=1,2
            I=NS+2*JR-1
            IPO=19+2*JR
            K(I,1)=14
            K(I,2)=IFLS(JR+2)
            IF(IFLS(JR+2).EQ.0) K(I,2)=21
            K(I,3)=0
            K(I,4)=IPO
            K(I,5)=IPO
            P(I,1)=0.
            P(I,2)=0.
            P(I,3)=DPLCM*(-1)**(JR+1)
            P(I,4)=(DSH+DQ2(3-JR)-DQ2(JR))/DSHR
            P(I,5)=-SQRT(SNGL(DQ2(JR)))
            K(I+1,1)=-1
            K(I+1,2)=K(IPO+1,2)
            K(I+1,3)=I
            K(I+1,4)=0
            K(I+1,5)=0
            P(I+1,1)=0.
            P(I+1,2)=0.
            P(I+1,3)=IPO
            P(I+1,4)=IPO
            P(I+1,5)=0.
            P(IPO+1,1)=I
            P(IPO+1,2)=I
            K(IPO,4)=MOD(K(IPO,4),MSTU(5))+I*MSTU(5)
            K(IPO,5)=MOD(K(IPO,5),MSTU(5))+I*MSTU(5)
  160     CONTINUE
        ELSE
C..LEPTOPRODUCTION EVENTS: BOSON AND HADRON REST FRAME
          I1=NS+2*ILEP-1
          I2=NS-2*ILEP+5
          DO 170 ITEMP=NS+1,NS+4
            DO 170 J=1,5
              K(ITEMP,J)=0
  170     P(ITEMP,J)=0.
          DO 180 J=1,5
  180     P(I1,J)=P(NQ,J)
          K(NS+1,1)=11
          K(NS+3,1)=14
          IF(ILEP.EQ.2) THEN
            K(NS+1,1)=14
            K(NS+3,1)=11
          ENDIF
          K(NS+2,1)=-1
          K(NS+4,1)=-1
          K(NS+1,3)=0
          K(NS+2,3)=NS+1
          K(NS+3,3)=0
          K(NS+4,3)=NS+3
          K(I1,2)=KFL(2,ILEP)
          K(I2,2)=KFL(2,3-ILEP)
          DPD(1)=DSH+DQ2(1)+DQ2(2)
          DPD(3)=(3-2*ILEP)*DSQRT(DPD(1)**2-4D0*DQ2(1)*DQ2(2))
          P(I2,3)=(DPQS(2)*DPD(3)-DPQS(1)*DPD(1))/
     +    (2D0*DQ2(JR))
          P(I2,4)=(DPQS(1)*DPD(3)-DPQS(2)*DPD(1))/
     +    (2D0*DQ2(JR))
          P(I2,5)=-SQRT(SNGL(DQ2(3-ILEP)))
          P(I2+1,3)=MAX(IPU1,IPU2)
          P(I2+1,4)=MAX(IPU1,IPU2)
          K(I2,4)=K(I2,4)-MOD(K(I2,4),MSTU(5))+MAX(IPU1,IPU2)
          K(I2,5)=K(I2,5)-MOD(K(I2,5),MSTU(5))+MAX(IPU1,IPU2)
          P(26-2*ILEP,1)=I2
          P(26-2*ILEP,2)=I2
          K(25-2*ILEP,4)=MOD(K(25-2*ILEP,4),MSTU(5))+I2*MSTU(5)
          K(25-2*ILEP,5)=MOD(K(25-2*ILEP,5),MSTU(5))+I2*MSTU(5)
          N=N+2
        ENDIF

      ELSEIF(N.GT.NS+4) THEN
C...FIND MAXIMUM ALLOWED MASS OF TIMELIKE PARTON
        DQ2(3)=Q2B
        IF(IPY(14).GE.3.AND.IPY(14).LE.6) DQ2(3)=Q2B/(1.-Z)
        IF(IS(1).GE.1.AND.IS(1).LE.MSTU(4)) THEN
          DPC(1)=P(IS(1),4)
          DPC(3)=0.5*(ABS(P(IS(1),3))+ABS(P(IS(2),3)))
        ELSE
C...IS(1) NOT INITIALIZED
          DPC(1)=0.
          DPC(3)=0.5*(       0.      +ABS(P(IS(2),3)))
        ENDIF
        DPC(2)=P(IS(2),4)
        DPD(1)=DSH+DQ2(JR)+DQ2(JT)
        DPD(2)=DSHZ+DQ2(JR)+DQ2(3)
        DPD(3)=DSQRT(DPD(1)**2-4.*DQ2(JR)*DQ2(JT))
        DPD(4)=DSQRT(DPD(2)**2-4.*DQ2(JR)*DQ2(3))
        IKIN=0
        IF((Q2S(JR).GE.0.5*PYPAR(22).AND.DPD(1)-DPD(3).GE.1D-10*DPD(1))
     +  .OR.ILEP.GE.1) IKIN=1
        IF(IKIN.EQ.0) DMSMA=(DQ2(JT)/DBLE(ZS(JT))-DQ2(3))*(DSH/
     +  (DSH+DQ2(JT))-DSH/(DSHZ+DQ2(3)))
        IF(IKIN.EQ.1) DMSMA=(DPD(1)*DPD(2)-DPD(3)*DPD(4))/(2.*
     +  DQ2(JR))-DQ2(JT)-DQ2(3)

C...GENERATE TIMELIKE PARTON SHOWER (IF REQUIRED)
        IT=N-1
        K(IT,1)=3
        K(IT,2)=IFLB-IFLS(JT+2)
        IF(IFLB-IFLS(JT+2).EQ.0) K(IT,2)=21
        P(IT,5)=ULMASS(K(IT,2))
        IF(SNGL(DMSMA).LE.P(IT,5)**2) GOTO 10
        P(IT,2)=0.
        DO 190 J=1,5
          K(IT+1,J)=0
  190   P(IT+1,J)=0.
        K(IT+1,1)=-1
        K(IT+1,2)=K(IS(JT)+1,2)
        K(IT+1,3)=IT
        IF(MOD(IPY(14),2).EQ.0) THEN
          P(IT,1)=0.
          IF(ILEP.EQ.0) P(IT,4)=(DSHZ-DSH-P(IT,5)**2)/DSHR
          IF(ILEP.GE.1) P(IT,4)=0.5*(P(IS(JT),3)*DPD(2)+
     +    DPQS(1)*(DQ2(JT)+DQ2(3)+P(IT,5)**2))/(P(IS(JT),3)*DPQS(2)-
     +    P(IS(JT),4)*DPQS(1))-DPC(JT)
          P(IT,3)=SQRT(MAX(0.,P(IT,4)**2-P(IT,5)**2))
          CALL LUSHOW(IT,0,SQRT(MIN(SNGL(DMSMA),PYPAR(25)*Q2)))
          IF(N.GE.IT+2) P(IT,5)=P(IT+2,5)
          IF(N.GT.MSTU(4)-10) THEN
            WRITE(6,*) ' PYSSPB: NO MORE MEMORY IN LUJETS'
            LST(21)=52
            RETURN
          ENDIF
          DO 200 I=N+1,N+8
            DO 200 J=1,5
              K(I,J)=0
  200     P(I,J)=0.
        ENDIF

C...RECONSTRUCT KINEMATICS OF BRANCHING: TIMELIKE PARTON SHOWER
        DMS=P(IT,5)**2
        IF(IKIN.EQ.0.AND.ILEP.EQ.0) DPT2=(DMSMA-DMS)*(DSHZ+DQ2(3))/
     +  (DSH+DQ2(JT))
        IF(IKIN.EQ.1.AND.ILEP.EQ.0) DPT2=(DMSMA-DMS)*(0.5*DPD(1)*
     +  DPD(2)+0.5*DPD(3)*DPD(4)-DQ2(JR)*(DQ2(JT)+DQ2(3)+DMS))/
     +  (4.*DSH*DPC(3)**2)
        IF(IKIN.EQ.1.AND.ILEP.GE.1) DPT2=(DMSMA-DMS)*(0.5*DPD(1)*
     +  DPD(2)+0.5*DPD(3)*DPD(4)-DQ2(JR)*(DQ2(JT)+DQ2(3)+DMS))/
     +  DPD(3)**2
        IF(DPT2.LT.0.) GOTO 10
        K(IT,3)=N+1
        P(IT,1)=SQRT(SNGL(DPT2))
        IF(ILEP.EQ.0) THEN
          DPB(1)=(0.5*DPD(2)-DPC(JR)*(DSHZ+DQ2(JR)-DQ2(JT)-DMS)/
     +    DSHR)/DPC(3)-DPC(3)
          P(IT,3)=DPB(1)*(-1)**(JT+1)
          P(IT,4)=(DSHZ-DSH-DMS)/DSHR
        ELSE
          DPC(3)=DQ2(JT)+DQ2(3)+DMS
          DPB(2)=DPQS(2)*DBLE(P(IS(JT),3))-DPQS(1)*DPC(JT)
          DPB(1)=0.5D0*(DPC(JT)*DPD(2)+DPQS(2)*DPC(3))/DPB(2)-
     +    DBLE(P(IS(JT),3))
          P(IT,3)=DPB(1)
          P(IT,4)=0.5D0*(DBLE(P(IS(JT),3))*DPD(2)+
     +    DPQS(1)*DPC(3))/DPB(2)-DPC(JT)
        ENDIF
        IF(N.GE.IT+2) THEN
          MSTU(1)=IT+2
          DPB(1)=DSQRT(DPB(1)**2+DPT2)
          DPB(2)=DSQRT(DPB(1)**2+DMS)
          DPB(3)=P(IT+2,3)
          DPB(4)=DSQRT(DPB(3)**2+DMS)
          DBEZ=(DPB(4)*DPB(1)-DPB(3)*DPB(2))/(DPB(4)*DPB(2)-DPB(3)*
     +    DPB(1))
          CALL LUDBRB(MSTU(1),MSTU(2),0.,0.,0.D0,0.D0,DBEZ)
          THE=ULANGL(P(IT,3),P(IT,1))
          CALL LUDBRB(MSTU(1),MSTU(2),THE,0.,0.D0,0.D0,0.D0)
        ENDIF

C...RECONSTRUCT KINEMATICS OF BRANCHING: SPACELIKE PARTON
        K(N+1,1)=14
        K(N+1,2)=IFLB
        IF(IFLB.EQ.0) K(N+1,2)=21
        K(N+1,3)=0
        P(N+1,1)=P(IT,1)
        P(N+1,2)=0.
        P(N+1,3)=P(IT,3)+P(IS(JT),3)
        P(N+1,4)=P(IT,4)+P(IS(JT),4)
        P(N+1,5)=-SQRT(SNGL(DQ2(3)))
        DO 210 J=1,5
          K(N+2,J)=0
  210   P(N+2,J)=0.
        K(N+2,1)=-1
        K(N+2,2)=K(IS(JT)+1,2)
        K(N+2,3)=N+1

C...DEFINE COLOUR FLOW OF BRANCHING
        K(IS(JT),1)=14
        K(IS(JT),3)=N+1
        ID1=IT
        KN1=ISIGN(500+IABS(K(N+1,2)),2*K(N+1,2)+1)
        KD1=ISIGN(500+IABS(K(ID1,2)),2*K(ID1,2)+1)
        IF(K(N+1,2).EQ.21) KN1=500
        IF(K(ID1,2).EQ.21) KD1=500
        IF((KN1.GE.501.AND.KD1.GE.501).OR.(KN1.LT.0.AND.
     +  KD1.EQ.500).OR.(KN1.EQ.500.AND.KD1.EQ.500.AND.
     +  RLU(0).GT.0.5).OR.(KN1.EQ.500.AND.KD1.LT.0))
     +  ID1=IS(JT)
        ID2=IT+IS(JT)-ID1
        P(N+2,3)=ID1
        P(N+2,4)=ID2
        P(ID1+1,1)=N+1
        P(ID1+1,2)=ID2
        P(ID2+1,1)=ID1
        P(ID2+1,2)=N+1
        K(N+1,4)=K(N+1,4)-MOD(K(N+1,4),MSTU(5))+ID1
        K(N+1,5)=K(N+1,5)-MOD(K(N+1,5),MSTU(5))+ID2
        K(ID1,4)=MOD(K(ID1,4),MSTU(5))+(N+1)*MSTU(5)
        K(ID1,5)=MOD(K(ID1,5),MSTU(5))+ID2*MSTU(5)
        K(ID2,4)=MOD(K(ID2,4),MSTU(5))+ID1*MSTU(5)
        K(ID2,5)=MOD(K(ID2,5),MSTU(5))+(N+1)*MSTU(5)
        N=N+2
C     CALL GULIST(22,2)

C...BOOST TO NEW CM-FRAME
        MSTU(1)=NS+1
        IF(ILEP.EQ.0) THEN
          CALL LUDBRB(MSTU(1),MSTU(2),0.,0.,
     +    -DBLE(P(N-1,1)+P(IS(JR),1))/DBLE(P(N-1,4)+P(IS(JR),4)),
     +    0.D0,-DBLE(P(N-1,3)+P(IS(JR),3))/DBLE(P(N-1,4)+P(IS(JR),4)))
          IR=N-1+(JT-1)*(IS(1)-N+1)
          CALL LUDBRB(MSTU(1),MSTU(2),
     +    -ULANGL(P(IR,3),P(IR,1)),PARU(2)*RLU(0),0.D0,0.D0,0.D0)
        ELSE
C...REORIENTATE EVENT WITHOUT CHANGING THE BOSON FOUR MOMENTUM
          DO 220 J=1,4
  220     DPQ(J)=P(NQ,J)
          DBE1(4)=DPQ(4)+DBLE(P(N-1,4))
          DO 230 J=1,3,2
  230     DBE1(J)=-(DPQ(J)+DBLE(P(N-1,J)))/DBE1(4)
          DBE1(4)=1D0/DSQRT(1D0-DBE1(1)**2-DBE1(3)**2)
          DBEP=DBE1(1)*DPQ(1)+DBE1(3)*DPQ(3)
          DGABEP=DBE1(4)*(DBE1(4)*DBEP/(1D0+DBE1(4))+DPQ(4))
          DO 240 J=1,3,2
  240     DPQ(J)=DPQ(J)+DGABEP*DBE1(J)
          DPQ(4)=DBE1(4)*(DPQ(4)+DBEP)
          DPC(1)=DSQRT(DPQ(1)**2+DPQ(3)**2)
          DBE2(4)=-(DPQ(4)*DPC(1)-DPQS(2)*DSQRT(DPQS(2)**2+DPC(1)**2-
     +    DPQ(4)**2))/(DPC(1)**2+DPQS(2)**2)
          THE=ULANGL(SNGL(DPQ(3)),SNGL(DPQ(1)))
          DBE2(1)=DBE2(4)*DSIN(DBLE(THE))
          DBE2(3)=DBE2(4)*DCOS(DBLE(THE))
          DBE2(4)=1D0/(1D0-DBE2(1)**2-DBE2(3)*2)

C...CONSTRUCT THE COMBINED BOOST
          DPB(1)=DBE1(4)**2*DBE2(4)/(1D0+DBE1(4))
          DPB(2)=DBE1(1)*DBE2(1)+DBE1(3)*DBE2(3)
          DPB(3)=DBE1(4)*DBE2(4)*(1D0+DPB(2))
          DO 250 JB=1,3,2
  250     DROBO(JB+2)=(DBE1(4)*DBE2(4)*DBE1(JB)+DBE2(4)*DBE2(JB)+
     +    DPB(1)*DBE1(JB)*DPB(2))/DPB(3)
          CALL LUDBRB(MSTU(1),MSTU(2),0.,0.,DROBO(3),0.D0,DROBO(5))
          IF(ILEP.EQ.1) THE=ULANGL(P(NS+1,3),P(NS+1,1))
          IF(ILEP.EQ.2) THE=PARU(1)+ULANGL(P(NS+3,3),P(NS+3,1))
          CALL LUDBRB(MSTU(1),MSTU(2),-THE,PARU(2)*RLU(0),0D0,0D0,0D0)
        ENDIF
        MSTU(1)=0
      ENDIF

C...SAVE QUANTITIES, LOOP BACK
      IS(JT)=N-1
      IF(ILEP.EQ.2.AND.N.EQ.NS+4) IS(JT)=N-3
      Q2S(JT)=Q2B
      DQ2(JT)=Q2B
      IF(IPY(14).GE.3.AND.IPY(14).LE.6) DQ2(JT)=Q2B/(1.-Z)
      DSH=DSHZ
      IF(Q2B.GE.0.5*PYPAR(22)) THEN
        IFLS(JT+2)=IFLS(JT)
        IFLS(JT)=IFLA
        XS(JT)=XA
        ZS(JT)=Z
        DO 260 IFL=-6,6
  260   XFS(JT,IFL)=XFA(IFL)
        TEVS(JT)=TEVB
      ELSE
        IF(JT.EQ.1) IPU1=N-1
        IF(JT.EQ.2) IPU2=N-1
      ENDIF
      IF(MAX(IABS(1-ILEP)*Q2S(1),MIN(1,2-ILEP)*Q2S(2)).GE.0.5*PYPAR(22)
     +.OR.N.LE.NS+2) GOTO 40
      IF(ILEP.EQ.0) THEN
C...BOOST HARD SCATTERING PARTONS TO FRAME OF SHOWER INITIATORS
        DO 270 J=1,3
  270   DROBO(J+2)=(P(NS+1,J)+P(NS+3,J))/(P(NS+1,4)+P(NS+3,4))
        DO 280 J=1,5
  280   P(N+2,J)=P(NS+1,J)
        MSTU(1)=N+2
        MSTU(2)=N+2
        CALL LUDBRB(N+2,N+2,0.,0.,-DROBO(3),-DROBO(4),-DROBO(5))
        ROBO(2)=ULANGL(P(N+2,1),P(N+2,2))
        ROBO(1)=ULANGL(P(N+2,3),SQRT(P(N+2,1)**2+P(N+2,2)**2))
        MSTU(1)=4
        MSTU(2)=NS
        CALL LUDBRB(4,NS,ROBO(1),ROBO(2),DROBO(3),DROBO(4),DROBO(5))
        MSTU(1)=0
        MSTU(2)=0
      ENDIF

C...STORE USER INFORMATION
      K(21,1)=14
      IF(ILEP.NE.0) K(21,1)=11
      K(23,1)=14
      K(21,3)=NS+1
      K(23,3)=NS+3
      DO 290 JT=1,2
        KFL(1,JT)=IFLS(JT)
        IF(IFLS(JT).EQ.0) KFL(1,JT)=21
  290 PYVAR(30+JT)=XS(JT)

      DO 300 I=NS+1,N
        DO 300 J=1,5
  300 V(I,J)=0.

      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE PYSTFU(KF,X,Q2,XPQ)

C...MODIFIED VERSION OF ROUTINE IN PYTHIA 5.6, COURTESY OF T. SJOSTRAND.
C...GIVES PROTON AND NEUTRON STRUCTURE FUNCTIONS ACCORDING TO A FEW
C...DIFFERENT PARAMETRIZATIONS.
C...NOTE THAT WHAT IS CODED IS X TIMES THE PROBABILITY DISTRIBUTION,
C...I.E. XQ(X,Q2) ETC.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP
      SAVE /LUDAT1/
      DIMENSION XPQ(-6:6),XPPR(-6:6)
      DOUBLE PRECISION XX,QQ,UPV,DNV,SEA,STR,CHM,BOT,TOP,GLU,
     +XPDF(-6:6),XPDG(0:5),VAL(20)
      CHARACTER*20 PARM(20)
      DATA NPDF/0/,UPV,DNV,SEA,STR,CHM,BOT,TOP,GLU/8*0.D0/

C...RESET STRUCTURE FUNCTIONS.
      DO 10  KFL=-6,6
   10 XPQ(KFL)=0.

C...CHECK X AND PARTICLE SPECIES.
      IF(X.LE.0..OR.X.GE.1.) THEN
        WRITE(MSTU(11),10000) X
        RETURN
      ENDIF
      KFA=IABS(KF)
      IF(KFA.NE.2112.AND.KFA.NE.2212) THEN
        WRITE(MSTU(11),10100) KF
        RETURN
      ENDIF

C...CONVERT LST SWITCHES TO MSPT, AND PARL PARAMETERS TO PARP
      MSTP57=1
      IF(LST(15).LT.0) MSTP57=0
      MSTP51=IABS(LST(15))
      MSTP52=LST(16)
      MSTP58=LST(12)
      PARP51=PARL(20)
      PARL(26)=0.

C...PROTON STRUCTURE FUNCTION CALL.
      IF(MSTP52.EQ.1.AND.MSTP51.GE.1.AND.MSTP51.LE.10) THEN
        CALL PYSTPR(X,Q2,XPPR)
        DO 20  KFL=-6,6
   20   XPQ(KFL)=XPPR(KFL)
      ELSEIF(MSTP52.EQ.2) THEN
C...CALL PDFLIB STRUCTURE FUNCTIONS.
        XX=X
        QQ=SQRT(MAX(0.,Q2))
        PARM(1)='MODE'
        VAL(1)=MSTP51
        NPDF=NPDF+1
C!..ENABLE THE NEXT TWO LINES TO USE PDFLIB.
C!        IF(NPDF.EQ.1) CALL PDFSET(PARM,VAL)
C!        CALL STRUCTF(XX,QQ,UPV,DNV,SEA,STR,CHM,BOT,TOP,GLU)
        XPQ(0)=GLU
        XPQ(1)=DNV+SEA
        XPQ(-1)=SEA
        XPQ(2)=UPV+SEA
        XPQ(-2)=SEA
        XPQ(3)=STR
        XPQ(-3)=STR
        XPQ(4)=CHM
        XPQ(-4)=CHM
        XPQ(5)=BOT
        XPQ(-5)=BOT
        XPQ(6)=TOP
        XPQ(-6)=TOP
      ELSEIF(MSTP52.EQ.3) THEN
C...CALL PAKPDF STRUCTURE FUNCTIONS.
        IPARC=(MSTP51+50)/100
        ISETC=MSTP51-100*IPARC
        XX=X
        QQ=Q2
C!..ENABLE THE NEXT LINE TO USE PAKPDF.
C!        CALL PDVAL(IPARC,ISETC,XX,QQ,XPDF,IRETC)
        DO 30  KFL=-6,6
   30   XPQ(KFL)=XPDF(KFL)
      ELSE
        WRITE(MSTU(11),10200) KF,MSTP52,MSTP51
      ENDIF

C...ISOSPIN CONJUGATION FOR NEUTRON.
      IF(KFA.EQ.2112) THEN
        XPS=XPQ(1)
        XPQ(1)=XPQ(2)
        XPQ(2)=XPS
        XPS=XPQ(-1)
        XPQ(-1)=XPQ(-2)
        XPQ(-2)=XPS
      ENDIF

C...CHARGE CONJUGATION FOR ANTIPARTICLE.
      IF(KF.LT.0) THEN
        DO 40  KFL=1,6
          XPS=XPQ(KFL)
          XPQ(KFL)=XPQ(-KFL)
          XPQ(-KFL)=XPS
   40   CONTINUE
      ENDIF

C...CHECK POSITIVITY AND RESET ABOVE MAXIMUM ALLOWED FLAVOUR.
      DO 50  KFL=-6,6
        XPQ(KFL)=MAX(0.,XPQ(KFL))
   50 IF(IABS(KFL).GT.MSTP58) XPQ(KFL)=0.

C...FORMATS FOR ERROR PRINTOUTS.
10000 FORMAT(' ERROR: X VALUE OUTSIDE PHYSICAL RANGE; X =',1P,E12.3)
10100 FORMAT(' ERROR: ILLEGAL PARTICLE CODE FOR STRUCTURE FUNCTION;',
     +' KF =',I5)
10200 FORMAT(' ERROR: UNKNOWN STRUCTURE FUNCTION; KF, LIBRARY, SET =',
     +3I5)

      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.28  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C*********************************************************************

      SUBROUTINE PYSTPR(X,Q2,XPPR)

C...MODIFIED VERSION OF ROUTINE IN PYTHIA 5.6, COURTESY OF T. SJOSTRAND.
C...GIVES PROTON STRUCTURE FUNCTIONS ACCORDING TO A FEW DIFFERENT
C...PARAMETRIZATIONS.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP
      SAVE /LUDAT1/,/LUDAT2/
      DIMENSION XPPR(-6:6),XQ(9),TX(6),TT(6),TS(6),NEHLQ(8,2),
     +CEHLQ(6,6,2,8,2),CDO(3,6,5,2),CMT(0:3,0:2,9,4),EXMT(0:3)
      DATA ALAM,VX/2*0./

C...THE FOLLOWING DATA LINES ARE COEFFICIENTS NEEDED IN THE
C...EICHTEN, HINCHLIFFE, LANE, QUIGG PROTON STRUCTURE FUNCTION
C...PARAMETRIZATIONS, SEE BELOW.
C...POWERS OF 1-X IN DIFFERENT CASES.
      DATA NEHLQ/3,4,7,5,7,7,7,7,3,4,7,6,7,7,7,7/
C...EXPANSION COEFFICIENTS FOR UP VALENCE QUARK DISTRIBUTION.
      DATA (((CEHLQ(IX,IT,NX,1,1),IX=1,6),IT=1,6),NX=1,2)/
     + 7.677E-01,-2.087E-01,-3.303E-01,-2.517E-02,-1.570E-02,-1.000E-04,
     +-5.326E-01,-2.661E-01, 3.201E-01, 1.192E-01, 2.434E-02, 7.620E-03,
     + 2.162E-01, 1.881E-01,-8.375E-02,-6.515E-02,-1.743E-02,-5.040E-03,
     +-9.211E-02,-9.952E-02, 1.373E-02, 2.506E-02, 8.770E-03, 2.550E-03,
     + 3.670E-02, 4.409E-02, 9.600E-04,-7.960E-03,-3.420E-03,-1.050E-03,
     +-1.549E-02,-2.026E-02,-3.060E-03, 2.220E-03, 1.240E-03, 4.100E-04,
     + 2.395E-01, 2.905E-01, 9.778E-02, 2.149E-02, 3.440E-03, 5.000E-04,
     + 1.751E-02,-6.090E-03,-2.687E-02,-1.916E-02,-7.970E-03,-2.750E-03,
     +-5.760E-03,-5.040E-03, 1.080E-03, 2.490E-03, 1.530E-03, 7.500E-04,
     + 1.740E-03, 1.960E-03, 3.000E-04,-3.400E-04,-2.900E-04,-1.800E-04,
     +-5.300E-04,-6.400E-04,-1.700E-04, 4.000E-05, 6.000E-05, 4.000E-05,
     + 1.700E-04, 2.200E-04, 8.000E-05, 1.000E-05,-1.000E-05,-1.000E-05/
      DATA (((CEHLQ(IX,IT,NX,1,2),IX=1,6),IT=1,6),NX=1,2)/
     + 7.237E-01,-2.189E-01,-2.995E-01,-1.909E-02,-1.477E-02, 2.500E-04,
     +-5.314E-01,-2.425E-01, 3.283E-01, 1.119E-01, 2.223E-02, 7.070E-03,
     + 2.289E-01, 1.890E-01,-9.859E-02,-6.900E-02,-1.747E-02,-5.080E-03,
     +-1.041E-01,-1.084E-01, 2.108E-02, 2.975E-02, 9.830E-03, 2.830E-03,
     + 4.394E-02, 5.116E-02,-1.410E-03,-1.055E-02,-4.230E-03,-1.270E-03,
     +-1.991E-02,-2.539E-02,-2.780E-03, 3.430E-03, 1.720E-03, 5.500E-04,
     + 2.410E-01, 2.884E-01, 9.369E-02, 1.900E-02, 2.530E-03, 2.400E-04,
     + 1.765E-02,-9.220E-03,-3.037E-02,-2.085E-02,-8.440E-03,-2.810E-03,
     +-6.450E-03,-5.260E-03, 1.720E-03, 3.110E-03, 1.830E-03, 8.700E-04,
     + 2.120E-03, 2.320E-03, 2.600E-04,-4.900E-04,-3.900E-04,-2.300E-04,
     +-6.900E-04,-8.200E-04,-2.000E-04, 7.000E-05, 9.000E-05, 6.000E-05,
     + 2.400E-04, 3.100E-04, 1.100E-04, 0.000E+00,-2.000E-05,-2.000E-05/
C...EXPANSION COEFFICIENTS FOR DOWN VALENCE QUARK DISTRIBUTION.
      DATA (((CEHLQ(IX,IT,NX,2,1),IX=1,6),IT=1,6),NX=1,2)/
     + 3.813E-01,-8.090E-02,-1.634E-01,-2.185E-02,-8.430E-03,-6.200E-04,
     +-2.948E-01,-1.435E-01, 1.665E-01, 6.638E-02, 1.473E-02, 4.080E-03,
     + 1.252E-01, 1.042E-01,-4.722E-02,-3.683E-02,-1.038E-02,-2.860E-03,
     +-5.478E-02,-5.678E-02, 8.900E-03, 1.484E-02, 5.340E-03, 1.520E-03,
     + 2.220E-02, 2.567E-02,-3.000E-05,-4.970E-03,-2.160E-03,-6.500E-04,
     +-9.530E-03,-1.204E-02,-1.510E-03, 1.510E-03, 8.300E-04, 2.700E-04,
     + 1.261E-01, 1.354E-01, 3.958E-02, 8.240E-03, 1.660E-03, 4.500E-04,
     + 3.890E-03,-1.159E-02,-1.625E-02,-9.610E-03,-3.710E-03,-1.260E-03,
     +-1.910E-03,-5.600E-04, 1.590E-03, 1.590E-03, 8.400E-04, 3.900E-04,
     + 6.400E-04, 4.900E-04,-1.500E-04,-2.900E-04,-1.800E-04,-1.000E-04,
     +-2.000E-04,-1.900E-04, 0.000E+00, 6.000E-05, 4.000E-05, 3.000E-05,
     + 7.000E-05, 8.000E-05, 2.000E-05,-1.000E-05,-1.000E-05,-1.000E-05/
      DATA (((CEHLQ(IX,IT,NX,2,2),IX=1,6),IT=1,6),NX=1,2)/
     + 3.578E-01,-8.622E-02,-1.480E-01,-1.840E-02,-7.820E-03,-4.500E-04,
     +-2.925E-01,-1.304E-01, 1.696E-01, 6.243E-02, 1.353E-02, 3.750E-03,
     + 1.318E-01, 1.041E-01,-5.486E-02,-3.872E-02,-1.038E-02,-2.850E-03,
     +-6.162E-02,-6.143E-02, 1.303E-02, 1.740E-02, 5.940E-03, 1.670E-03,
     + 2.643E-02, 2.957E-02,-1.490E-03,-6.450E-03,-2.630E-03,-7.700E-04,
     +-1.218E-02,-1.497E-02,-1.260E-03, 2.240E-03, 1.120E-03, 3.500E-04,
     + 1.263E-01, 1.334E-01, 3.732E-02, 7.070E-03, 1.260E-03, 3.400E-04,
     + 3.660E-03,-1.357E-02,-1.795E-02,-1.031E-02,-3.880E-03,-1.280E-03,
     +-2.100E-03,-3.600E-04, 2.050E-03, 1.920E-03, 9.800E-04, 4.400E-04,
     + 7.700E-04, 5.400E-04,-2.400E-04,-3.900E-04,-2.400E-04,-1.300E-04,
     +-2.600E-04,-2.300E-04, 2.000E-05, 9.000E-05, 6.000E-05, 4.000E-05,
     + 9.000E-05, 1.000E-04, 2.000E-05,-2.000E-05,-2.000E-05,-1.000E-05/
C...EXPANSION COEFFICIENTS FOR UP AND DOWN SEA QUARK DISTRIBUTIONS.
      DATA (((CEHLQ(IX,IT,NX,3,1),IX=1,6),IT=1,6),NX=1,2)/
     + 6.870E-02,-6.861E-02, 2.973E-02,-5.400E-03, 3.780E-03,-9.700E-04,
     +-1.802E-02, 1.400E-04, 6.490E-03,-8.540E-03, 1.220E-03,-1.750E-03,
     +-4.650E-03, 1.480E-03,-5.930E-03, 6.000E-04,-1.030E-03,-8.000E-05,
     + 6.440E-03, 2.570E-03, 2.830E-03, 1.150E-03, 7.100E-04, 3.300E-04,
     +-3.930E-03,-2.540E-03,-1.160E-03,-7.700E-04,-3.600E-04,-1.900E-04,
     + 2.340E-03, 1.930E-03, 5.300E-04, 3.700E-04, 1.600E-04, 9.000E-05,
     + 1.014E+00,-1.106E+00, 3.374E-01,-7.444E-02, 8.850E-03,-8.700E-04,
     + 9.233E-01,-1.285E+00, 4.475E-01,-9.786E-02, 1.419E-02,-1.120E-03,
     + 4.888E-02,-1.271E-01, 8.606E-02,-2.608E-02, 4.780E-03,-6.000E-04,
     +-2.691E-02, 4.887E-02,-1.771E-02, 1.620E-03, 2.500E-04,-6.000E-05,
     + 7.040E-03,-1.113E-02, 1.590E-03, 7.000E-04,-2.000E-04, 0.000E+00,
     +-1.710E-03, 2.290E-03, 3.800E-04,-3.500E-04, 4.000E-05, 1.000E-05/
      DATA (((CEHLQ(IX,IT,NX,3,2),IX=1,6),IT=1,6),NX=1,2)/
     + 1.008E-01,-7.100E-02, 1.973E-02,-5.710E-03, 2.930E-03,-9.900E-04,
     +-5.271E-02,-1.823E-02, 1.792E-02,-6.580E-03, 1.750E-03,-1.550E-03,
     + 1.220E-02, 1.763E-02,-8.690E-03,-8.800E-04,-1.160E-03,-2.100E-04,
     +-1.190E-03,-7.180E-03, 2.360E-03, 1.890E-03, 7.700E-04, 4.100E-04,
     +-9.100E-04, 2.040E-03,-3.100E-04,-1.050E-03,-4.000E-04,-2.400E-04,
     + 1.190E-03,-1.700E-04,-2.000E-04, 4.200E-04, 1.700E-04, 1.000E-04,
     + 1.081E+00,-1.189E+00, 3.868E-01,-8.617E-02, 1.115E-02,-1.180E-03,
     + 9.917E-01,-1.396E+00, 4.998E-01,-1.159E-01, 1.674E-02,-1.720E-03,
     + 5.099E-02,-1.338E-01, 9.173E-02,-2.885E-02, 5.890E-03,-6.500E-04,
     +-3.178E-02, 5.703E-02,-2.070E-02, 2.440E-03, 1.100E-04,-9.000E-05,
     + 8.970E-03,-1.392E-02, 2.050E-03, 6.500E-04,-2.300E-04, 2.000E-05,
     +-2.340E-03, 3.010E-03, 5.000E-04,-3.900E-04, 6.000E-05, 1.000E-05/
C...EXPANSION COEFFICIENTS FOR GLUON DISTRIBUTION.
      DATA (((CEHLQ(IX,IT,NX,4,1),IX=1,6),IT=1,6),NX=1,2)/
     + 9.482E-01,-9.578E-01, 1.009E-01,-1.051E-01, 3.456E-02,-3.054E-02,
     +-9.627E-01, 5.379E-01, 3.368E-01,-9.525E-02, 1.488E-02,-2.051E-02,
     + 4.300E-01,-8.306E-02,-3.372E-01, 4.902E-02,-9.160E-03, 1.041E-02,
     +-1.925E-01,-1.790E-02, 2.183E-01, 7.490E-03, 4.140E-03,-1.860E-03,
     + 8.183E-02, 1.926E-02,-1.072E-01,-1.944E-02,-2.770E-03,-5.200E-04,
     +-3.884E-02,-1.234E-02, 5.410E-02, 1.879E-02, 3.350E-03, 1.040E-03,
     + 2.948E+01,-3.902E+01, 1.464E+01,-3.335E+00, 5.054E-01,-5.915E-02,
     + 2.559E+01,-3.955E+01, 1.661E+01,-4.299E+00, 6.904E-01,-8.243E-02,
     +-1.663E+00, 1.176E+00, 1.118E+00,-7.099E-01, 1.948E-01,-2.404E-02,
     +-2.168E-01, 8.170E-01,-7.169E-01, 1.851E-01,-1.924E-02,-3.250E-03,
     + 2.088E-01,-4.355E-01, 2.239E-01,-2.446E-02,-3.620E-03, 1.910E-03,
     +-9.097E-02, 1.601E-01,-5.681E-02,-2.500E-03, 2.580E-03,-4.700E-04/
      DATA (((CEHLQ(IX,IT,NX,4,2),IX=1,6),IT=1,6),NX=1,2)/
     + 2.367E+00, 4.453E-01, 3.660E-01, 9.467E-02, 1.341E-01, 1.661E-02,
     +-3.170E+00,-1.795E+00, 3.313E-02,-2.874E-01,-9.827E-02,-7.119E-02,
     + 1.823E+00, 1.457E+00,-2.465E-01, 3.739E-02, 6.090E-03, 1.814E-02,
     +-1.033E+00,-9.827E-01, 2.136E-01, 1.169E-01, 5.001E-02, 1.684E-02,
     + 5.133E-01, 5.259E-01,-1.173E-01,-1.139E-01,-4.988E-02,-2.021E-02,
     +-2.881E-01,-3.145E-01, 5.667E-02, 9.161E-02, 4.568E-02, 1.951E-02,
     + 3.036E+01,-4.062E+01, 1.578E+01,-3.699E+00, 6.020E-01,-7.031E-02,
     + 2.700E+01,-4.167E+01, 1.770E+01,-4.804E+00, 7.862E-01,-1.060E-01,
     +-1.909E+00, 1.357E+00, 1.127E+00,-7.181E-01, 2.232E-01,-2.481E-02,
     +-2.488E-01, 9.781E-01,-8.127E-01, 2.094E-01,-2.997E-02,-4.710E-03,
     + 2.506E-01,-5.427E-01, 2.672E-01,-3.103E-02,-1.800E-03, 2.870E-03,
     +-1.128E-01, 2.087E-01,-6.972E-02,-2.480E-03, 2.630E-03,-8.400E-04/
C...EXPANSION COEFFICIENTS FOR STRANGE SEA QUARK DISTRIBUTION.
      DATA (((CEHLQ(IX,IT,NX,5,1),IX=1,6),IT=1,6),NX=1,2)/
     + 4.968E-02,-4.173E-02, 2.102E-02,-3.270E-03, 3.240E-03,-6.700E-04,
     +-6.150E-03,-1.294E-02, 6.740E-03,-6.890E-03, 9.000E-04,-1.510E-03,
     +-8.580E-03, 5.050E-03,-4.900E-03,-1.600E-04,-9.400E-04,-1.500E-04,
     + 7.840E-03, 1.510E-03, 2.220E-03, 1.400E-03, 7.000E-04, 3.500E-04,
     +-4.410E-03,-2.220E-03,-8.900E-04,-8.500E-04,-3.600E-04,-2.000E-04,
     + 2.520E-03, 1.840E-03, 4.100E-04, 3.900E-04, 1.600E-04, 9.000E-05,
     + 9.235E-01,-1.085E+00, 3.464E-01,-7.210E-02, 9.140E-03,-9.100E-04,
     + 9.315E-01,-1.274E+00, 4.512E-01,-9.775E-02, 1.380E-02,-1.310E-03,
     + 4.739E-02,-1.296E-01, 8.482E-02,-2.642E-02, 4.760E-03,-5.700E-04,
     +-2.653E-02, 4.953E-02,-1.735E-02, 1.750E-03, 2.800E-04,-6.000E-05,
     + 6.940E-03,-1.132E-02, 1.480E-03, 6.500E-04,-2.100E-04, 0.000E+00,
     +-1.680E-03, 2.340E-03, 4.200E-04,-3.400E-04, 5.000E-05, 1.000E-05/
      DATA (((CEHLQ(IX,IT,NX,5,2),IX=1,6),IT=1,6),NX=1,2)/
     + 6.478E-02,-4.537E-02, 1.643E-02,-3.490E-03, 2.710E-03,-6.700E-04,
     +-2.223E-02,-2.126E-02, 1.247E-02,-6.290E-03, 1.120E-03,-1.440E-03,
     +-1.340E-03, 1.362E-02,-6.130E-03,-7.900E-04,-9.000E-04,-2.000E-04,
     + 5.080E-03,-3.610E-03, 1.700E-03, 1.830E-03, 6.800E-04, 4.000E-04,
     +-3.580E-03, 6.000E-05,-2.600E-04,-1.050E-03,-3.800E-04,-2.300E-04,
     + 2.420E-03, 9.300E-04,-1.000E-04, 4.500E-04, 1.700E-04, 1.100E-04,
     + 9.868E-01,-1.171E+00, 3.940E-01,-8.459E-02, 1.124E-02,-1.250E-03,
     + 1.001E+00,-1.383E+00, 5.044E-01,-1.152E-01, 1.658E-02,-1.830E-03,
     + 4.928E-02,-1.368E-01, 9.021E-02,-2.935E-02, 5.800E-03,-6.600E-04,
     +-3.133E-02, 5.785E-02,-2.023E-02, 2.630E-03, 1.600E-04,-8.000E-05,
     + 8.840E-03,-1.416E-02, 1.900E-03, 5.800E-04,-2.500E-04, 1.000E-05,
     +-2.300E-03, 3.080E-03, 5.500E-04,-3.700E-04, 7.000E-05, 1.000E-05/
C...EXPANSION COEFFICIENTS FOR CHARM SEA QUARK DISTRIBUTION.
      DATA (((CEHLQ(IX,IT,NX,6,1),IX=1,6),IT=1,6),NX=1,2)/
     + 9.270E-03,-1.817E-02, 9.590E-03,-6.390E-03, 1.690E-03,-1.540E-03,
     + 5.710E-03,-1.188E-02, 6.090E-03,-4.650E-03, 1.240E-03,-1.310E-03,
     +-3.960E-03, 7.100E-03,-3.590E-03, 1.840E-03,-3.900E-04, 3.400E-04,
     + 1.120E-03,-1.960E-03, 1.120E-03,-4.800E-04, 1.000E-04,-4.000E-05,
     + 4.000E-05,-3.000E-05,-1.800E-04, 9.000E-05,-5.000E-05,-2.000E-05,
     +-4.200E-04, 7.300E-04,-1.600E-04, 5.000E-05, 5.000E-05, 5.000E-05,
     + 8.098E-01,-1.042E+00, 3.398E-01,-6.824E-02, 8.760E-03,-9.000E-04,
     + 8.961E-01,-1.217E+00, 4.339E-01,-9.287E-02, 1.304E-02,-1.290E-03,
     + 3.058E-02,-1.040E-01, 7.604E-02,-2.415E-02, 4.600E-03,-5.000E-04,
     +-2.451E-02, 4.432E-02,-1.651E-02, 1.430E-03, 1.200E-04,-1.000E-04,
     + 1.122E-02,-1.457E-02, 2.680E-03, 5.800E-04,-1.200E-04, 3.000E-05,
     +-7.730E-03, 7.330E-03,-7.600E-04,-2.400E-04, 1.000E-05, 0.000E+00/
      DATA (((CEHLQ(IX,IT,NX,6,2),IX=1,6),IT=1,6),NX=1,2)/
     + 9.980E-03,-1.945E-02, 1.055E-02,-6.870E-03, 1.860E-03,-1.560E-03,
     + 5.700E-03,-1.203E-02, 6.250E-03,-4.860E-03, 1.310E-03,-1.370E-03,
     +-4.490E-03, 7.990E-03,-4.170E-03, 2.050E-03,-4.400E-04, 3.300E-04,
     + 1.470E-03,-2.480E-03, 1.460E-03,-5.700E-04, 1.200E-04,-1.000E-05,
     +-9.000E-05, 1.500E-04,-3.200E-04, 1.200E-04,-6.000E-05,-4.000E-05,
     +-4.200E-04, 7.600E-04,-1.400E-04, 4.000E-05, 7.000E-05, 5.000E-05,
     + 8.698E-01,-1.131E+00, 3.836E-01,-8.111E-02, 1.048E-02,-1.300E-03,
     + 9.626E-01,-1.321E+00, 4.854E-01,-1.091E-01, 1.583E-02,-1.700E-03,
     + 3.057E-02,-1.088E-01, 8.022E-02,-2.676E-02, 5.590E-03,-5.600E-04,
     +-2.845E-02, 5.164E-02,-1.918E-02, 2.210E-03,-4.000E-05,-1.500E-04,
     + 1.311E-02,-1.751E-02, 3.310E-03, 5.100E-04,-1.200E-04, 5.000E-05,
     +-8.590E-03, 8.380E-03,-9.200E-04,-2.600E-04, 1.000E-05,-1.000E-05/
C...EXPANSION COEFFICIENTS FOR BOTTOM SEA QUARK DISTRIBUTION.
      DATA (((CEHLQ(IX,IT,NX,7,1),IX=1,6),IT=1,6),NX=1,2)/
     + 9.010E-03,-1.401E-02, 7.150E-03,-4.130E-03, 1.260E-03,-1.040E-03,
     + 6.280E-03,-9.320E-03, 4.780E-03,-2.890E-03, 9.100E-04,-8.200E-04,
     +-2.930E-03, 4.090E-03,-1.890E-03, 7.600E-04,-2.300E-04, 1.400E-04,
     + 3.900E-04,-1.200E-03, 4.400E-04,-2.500E-04, 2.000E-05,-2.000E-05,
     + 2.600E-04, 1.400E-04,-8.000E-05, 1.000E-04, 1.000E-05, 1.000E-05,
     +-2.600E-04, 3.200E-04, 1.000E-05,-1.000E-05, 1.000E-05,-1.000E-05,
     + 8.029E-01,-1.075E+00, 3.792E-01,-7.843E-02, 1.007E-02,-1.090E-03,
     + 7.903E-01,-1.099E+00, 4.153E-01,-9.301E-02, 1.317E-02,-1.410E-03,
     +-1.704E-02,-1.130E-02, 2.882E-02,-1.341E-02, 3.040E-03,-3.600E-04,
     +-7.200E-04, 7.230E-03,-5.160E-03, 1.080E-03,-5.000E-05,-4.000E-05,
     + 3.050E-03,-4.610E-03, 1.660E-03,-1.300E-04,-1.000E-05, 1.000E-05,
     +-4.360E-03, 5.230E-03,-1.610E-03, 2.000E-04,-2.000E-05, 0.000E+00/
      DATA (((CEHLQ(IX,IT,NX,7,2),IX=1,6),IT=1,6),NX=1,2)/
     + 8.980E-03,-1.459E-02, 7.510E-03,-4.410E-03, 1.310E-03,-1.070E-03,
     + 5.970E-03,-9.440E-03, 4.800E-03,-3.020E-03, 9.100E-04,-8.500E-04,
     +-3.050E-03, 4.440E-03,-2.100E-03, 8.500E-04,-2.400E-04, 1.400E-04,
     + 5.300E-04,-1.300E-03, 5.600E-04,-2.700E-04, 3.000E-05,-2.000E-05,
     + 2.000E-04, 1.400E-04,-1.100E-04, 1.000E-04, 0.000E+00, 0.000E+00,
     +-2.600E-04, 3.200E-04, 0.000E+00,-3.000E-05, 1.000E-05,-1.000E-05,
     + 8.672E-01,-1.174E+00, 4.265E-01,-9.252E-02, 1.244E-02,-1.460E-03,
     + 8.500E-01,-1.194E+00, 4.630E-01,-1.083E-01, 1.614E-02,-1.830E-03,
     +-2.241E-02,-5.630E-03, 2.815E-02,-1.425E-02, 3.520E-03,-4.300E-04,
     +-7.300E-04, 8.030E-03,-5.780E-03, 1.380E-03,-1.300E-04,-4.000E-05,
     + 3.460E-03,-5.380E-03, 1.960E-03,-2.100E-04, 1.000E-05, 1.000E-05,
     +-4.850E-03, 5.950E-03,-1.890E-03, 2.600E-04,-3.000E-05, 0.000E+00/
C...EXPANSION COEFFICIENTS FOR TOP SEA QUARK DISTRIBUTION.
      DATA (((CEHLQ(IX,IT,NX,8,1),IX=1,6),IT=1,6),NX=1,2)/
     + 4.410E-03,-7.480E-03, 3.770E-03,-2.580E-03, 7.300E-04,-7.100E-04,
     + 3.840E-03,-6.050E-03, 3.030E-03,-2.030E-03, 5.800E-04,-5.900E-04,
     +-8.800E-04, 1.660E-03,-7.500E-04, 4.700E-04,-1.000E-04, 1.000E-04,
     +-8.000E-05,-1.500E-04, 1.200E-04,-9.000E-05, 3.000E-05, 0.000E+00,
     + 1.300E-04,-2.200E-04,-2.000E-05,-2.000E-05,-2.000E-05,-2.000E-05,
     +-7.000E-05, 1.900E-04,-4.000E-05, 2.000E-05, 0.000E+00, 0.000E+00,
     + 6.623E-01,-9.248E-01, 3.519E-01,-7.930E-02, 1.110E-02,-1.180E-03,
     + 6.380E-01,-9.062E-01, 3.582E-01,-8.479E-02, 1.265E-02,-1.390E-03,
     +-2.581E-02, 2.125E-02, 4.190E-03,-4.980E-03, 1.490E-03,-2.100E-04,
     + 7.100E-04, 5.300E-04,-1.270E-03, 3.900E-04,-5.000E-05,-1.000E-05,
     + 3.850E-03,-5.060E-03, 1.860E-03,-3.500E-04, 4.000E-05, 0.000E+00,
     +-3.530E-03, 4.460E-03,-1.500E-03, 2.700E-04,-3.000E-05, 0.000E+00/
      DATA (((CEHLQ(IX,IT,NX,8,2),IX=1,6),IT=1,6),NX=1,2)/
     + 4.260E-03,-7.530E-03, 3.830E-03,-2.680E-03, 7.600E-04,-7.300E-04,
     + 3.640E-03,-6.050E-03, 3.030E-03,-2.090E-03, 5.900E-04,-6.000E-04,
     +-9.200E-04, 1.710E-03,-8.200E-04, 5.000E-04,-1.200E-04, 1.000E-04,
     +-5.000E-05,-1.600E-04, 1.300E-04,-9.000E-05, 3.000E-05, 0.000E+00,
     + 1.300E-04,-2.100E-04,-1.000E-05,-2.000E-05,-2.000E-05,-1.000E-05,
     +-8.000E-05, 1.800E-04,-5.000E-05, 2.000E-05, 0.000E+00, 0.000E+00,
     + 7.146E-01,-1.007E+00, 3.932E-01,-9.246E-02, 1.366E-02,-1.540E-03,
     + 6.856E-01,-9.828E-01, 3.977E-01,-9.795E-02, 1.540E-02,-1.790E-03,
     +-3.053E-02, 2.758E-02, 2.150E-03,-4.880E-03, 1.640E-03,-2.500E-04,
     + 9.200E-04, 4.200E-04,-1.340E-03, 4.600E-04,-8.000E-05,-1.000E-05,
     + 4.230E-03,-5.660E-03, 2.140E-03,-4.300E-04, 6.000E-05, 0.000E+00,
     +-3.890E-03, 5.000E-03,-1.740E-03, 3.300E-04,-4.000E-05, 0.000E+00/

C...THE FOLLOWING DATA LINES ARE COEFFICIENTS NEEDED IN THE
C...DUKE, OWENS PROTON STRUCTURE FUNCTION PARAMETRIZATIONS, SEE BELOW.
C...EXPANSION COEFFICIENTS FOR (UP+DOWN) VALENCE QUARK DISTRIBUTION.
      DATA ((CDO(IP,IS,1,1),IS=1,6),IP=1,3)/
     + 4.190E-01, 3.460E+00, 4.400E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     + 4.000E-03, 7.240E-01,-4.860E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     +-7.000E-03,-6.600E-02, 1.330E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,1,2),IS=1,6),IP=1,3)/
     + 3.740E-01, 3.330E+00, 6.030E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     + 1.400E-02, 7.530E-01,-6.220E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     + 0.000E+00,-7.600E-02, 1.560E+00, 0.000E+00, 0.000E+00, 0.000E+00/
C...EXPANSION COEFFICIENTS FOR DOWN VALENCE QUARK DISTRIBUTION.
      DATA ((CDO(IP,IS,2,1),IS=1,6),IP=1,3)/
     + 7.630E-01, 4.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     +-2.370E-01, 6.270E-01,-4.210E-01, 0.000E+00, 0.000E+00, 0.000E+00,
     + 2.600E-02,-1.900E-02, 3.300E-02, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,2,2),IS=1,6),IP=1,3)/
     + 7.610E-01, 3.830E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     +-2.320E-01, 6.270E-01,-4.180E-01, 0.000E+00, 0.000E+00, 0.000E+00,
     + 2.300E-02,-1.900E-02, 3.600E-02, 0.000E+00, 0.000E+00, 0.000E+00/
C...EXPANSION COEFFICIENTS FOR (UP+DOWN+STRANGE) SEA QUARK DISTRIBUTION.
      DATA ((CDO(IP,IS,3,1),IS=1,6),IP=1,3)/
     + 1.265E+00, 0.000E+00, 8.050E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     +-1.132E+00,-3.720E-01, 1.590E+00, 6.310E+00,-1.050E+01, 1.470E+01,
     + 2.930E-01,-2.900E-02,-1.530E-01,-2.730E-01,-3.170E+00, 9.800E+00/
      DATA ((CDO(IP,IS,3,2),IS=1,6),IP=1,3)/
     + 1.670E+00, 0.000E+00, 9.150E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     +-1.920E+00,-2.730E-01, 5.300E-01, 1.570E+01,-1.010E+02, 2.230E+02,
     + 5.820E-01,-1.640E-01,-7.630E-01,-2.830E+00, 4.470E+01,-1.170E+02/
C...EXPANSION COEFFICIENTS FOR CHARM SEA QUARK DISTRIBUTION.
      DATA ((CDO(IP,IS,4,1),IS=1,6),IP=1,3)/
     + 0.000E+00,-3.600E-02, 6.350E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     + 1.350E-01,-2.220E-01, 3.260E+00,-3.030E+00, 1.740E+01,-1.790E+01,
     +-7.500E-02,-5.800E-02,-9.090E-01, 1.500E+00,-1.130E+01, 1.560E+01/
      DATA ((CDO(IP,IS,4,2),IS=1,6),IP=1,3)/ 0.000E+00,-1.200E-01,
     +3.510E+00, 0.000E+00, 0.000E+00, 0.000E+00, 6.700E-02,-2.330E-01,
     +3.660E+00,-4.740E-01, 9.500E+00,-1.660E+01,-3.100E-02,-2.300E-02,
     +-4.530E-01, 3.580E-01,-5.430E+00, 1.550E+01/
C...EXPANSION COEFFICIENTS FOR GLUON DISTRIBUTION.
      DATA ((CDO(IP,IS,5,1),IS=1,6),IP=1,3)/
     + 1.560E+00, 0.000E+00, 6.000E+00, 9.000E+00, 0.000E+00, 0.000E+00,
     +-1.710E+00,-9.490E-01, 1.440E+00,-7.190E+00,-1.650E+01, 1.530E+01,
     + 6.380E-01, 3.250E-01,-1.050E+00, 2.550E-01, 1.090E+01,-1.010E+01/
      DATA ((CDO(IP,IS,5,2),IS=1,6),IP=1,3)/
     + 8.790E-01, 0.000E+00, 4.000E+00, 9.000E+00, 0.000E+00, 0.000E+00,
     +-9.710E-01,-1.160E+00, 1.230E+00,-5.640E+00,-7.540E+00,-5.960E-01,
     + 4.340E-01, 4.760E-01,-2.540E-01,-8.170E-01, 5.500E+00, 1.260E-01/

C...THE FOLLOWING DATA LINES ARE COEFFICIENTS NEEDED IN THE
C...MORFIN AND TUNG STRUCTURE FUNCTION PARAMETRIZATIONS.
C...12 COEFFICIENTS EACH FOR D(VALENCE), U(VALENCE), G, U(SEA),
C...D(SEA), S, C, B AND T, IN THAT ORDER.
C...EXPANSION COEFFICIENTS FOR SET 1 (FIT S1).
      DATA (((CMT(IEX,IPN,IFL,1),IFL=1,9),IPN=0,2),IEX=0,3)/
     +   1.30,  1.64,  1.86, -0.60, -0.45, -1.10, -3.87, -6.14,-12.53,
     +  -0.57, -0.33, -2.76, -1.68, -1.64, -1.66,  0.79,  2.65,  8.13,
     +  -0.08, -0.10,  0.10,  0.08,  0.05,  0.13, -0.70, -1.24, -2.64,
     +   0.18,  0.08, -0.17, -0.19, -0.18, -0.19, -0.03, -0.10, -0.38,
     +   0.16,  0.14, -0.07, -0.16, -0.19, -0.09, -0.17, -0.03,  0.34,
     +  -0.02, -0.01,  0.02,  0.04,  0.06,  0.01,  0.03, -0.02, -0.14,
     +   5.27,  3.74,  7.33,  9.31,  9.36,  9.07,  7.96,  6.90, 16.30,
     +   0.43,  0.54, -0.88, -1.17, -1.01, -1.39,  0.95,  1.52,-13.23,
     +   0.06,  0.03, -0.08,  0.29,  0.20,  0.47, -0.38, -0.50,  4.77,
     +  -1.85, -2.04, -0.88, -1.45, -1.48, -1.26,  0.60,  0.80, -0.57,
     +   1.08,  0.88,  2.47,  1.65,  1.49,  1.96,  0.60,  1.05,  3.58,
     +  -0.03,  0.02, -0.32, -0.20, -0.12, -0.36,  0.08, -0.14, -0.99/
C...EXPANSION COEFFICIENTS FOR SET 2 (FIT B1).
      DATA (((CMT(IEX,IPN,IFL,2),IFL=1,9),IPN=0,2),IEX=0,3)/
     +   1.34,  1.62,  1.88, -0.99, -0.99, -0.99, -3.98, -6.28,-13.08,
     +  -0.57, -0.33, -2.78, -1.54, -1.54, -1.54,  0.72,  2.62,  8.54,
     +  -0.08, -0.10,  0.13,  0.10,  0.10,  0.10, -0.63, -1.18, -2.70,
     +   0.15,  0.11, -0.33, -0.33, -0.33, -0.33, -0.15, -0.18, -0.40,
     +   0.16,  0.14,  0.10,  0.03,  0.03,  0.03, -0.06,  0.02,  0.31,
     +  -0.02, -0.01, -0.04, -0.03, -0.03, -0.03,  0.00, -0.03, -0.12,
     +   5.30,  3.68,  7.52,  8.53,  8.53,  8.53,  7.46,  6.56, 15.35,
     +   0.43,  0.53, -1.13, -1.08, -1.08, -1.08,  0.96,  1.40,-11.83,
     +   0.06,  0.03,  0.04,  0.39,  0.39,  0.39, -0.30, -0.38,  4.16,
     +  -1.96, -1.94, -1.34, -1.55, -1.55, -1.55,  0.35,  0.65, -0.43,
     +   1.08,  0.87,  2.92,  2.02,  2.02,  2.02,  0.89,  1.13,  3.18,
     +  -0.03,  0.02, -0.49, -0.39, -0.39, -0.39, -0.04, -0.16, -0.82/
C...EXPANSION COEFFICIENTS FOR SET 3 (FIT B2).
      DATA (((CMT(IEX,IPN,IFL,3),IFL=1,9),IPN=0,2),IEX=0,3)/
     +   1.38,  1.64,  1.52, -0.85, -0.85, -0.85, -3.74, -6.07,-12.08,
     +  -0.59, -0.33, -2.71, -1.43, -1.43, -1.43,  0.21,  2.33,  7.31,
     +  -0.08, -0.10,  0.15, -0.03, -0.03, -0.03, -0.50, -1.15, -2.35,
     +   0.18,  0.09, -0.72, -0.82, -0.82, -0.82, -0.58, -0.52, -0.73,
     +   0.16,  0.14,  0.45,  0.35,  0.35,  0.35,  0.24,  0.22,  0.54,
     +  -0.02, -0.01, -0.15, -0.09, -0.10, -0.10, -0.07, -0.07, -0.18,
     +   5.40,  3.74,  7.75,  9.19,  9.19,  9.19,  9.63,  8.33, 21.14,
     +   0.42,  0.54, -1.56, -0.92, -0.92, -0.92, -1.13,  0.28,-19.17,
     +   0.06,  0.03,  0.16,  0.12,  0.12,  0.12,  0.25, -0.28,  6.64,
     +  -1.91, -2.02, -2.18, -2.76, -2.76, -2.76, -1.09, -0.52, -1.92,
     +   1.11,  0.88,  3.75,  2.56,  2.56,  2.56,  2.10,  1.91,  4.59,
     +  -0.03,  0.02, -0.76, -0.40, -0.40, -0.40, -0.33, -0.31, -1.25/
C...EXPANSION COEFFICIENTS FOR SET 4 (FIT E1).
      DATA (((CMT(IEX,IPN,IFL,4),IFL=1,9),IPN=0,2),IEX=0,3)/
     +   1.43,  1.69,  2.11, -0.84, -0.84, -0.84, -3.87, -6.09,-12.56,
     +  -0.65, -0.33, -3.01, -1.65, -1.65, -1.65,  0.85,  2.81,  8.69,
     +  -0.08, -0.11,  0.18,  0.12,  0.12,  0.12, -0.73, -1.34, -2.93,
     +   0.16,  0.11, -0.33, -0.32, -0.32, -0.32, -0.15, -0.17, -0.38,
     +   0.16,  0.14,  0.10,  0.02,  0.02,  0.02, -0.07,  0.01,  0.30,
     +  -0.02, -0.01, -0.04, -0.03, -0.03, -0.03,  0.00, -0.03, -0.12,
     +   6.17,  3.69,  7.93,  8.96,  8.96,  8.96,  7.83,  6.75, 14.62,
     +   0.43,  0.54, -1.40, -1.24, -1.24, -1.24,  1.00,  1.74,-11.27,
     +   0.06,  0.03,  0.09,  0.45,  0.45,  0.45, -0.36, -0.56,  4.29,
     +  -1.94, -1.99, -1.51, -1.70, -1.70, -1.70,  0.21,  0.54, -0.41,
     +   1.12,  0.90,  3.14,  2.15,  2.15,  2.15,  0.93,  1.15,  3.19,
     +  -0.02,  0.02, -0.55, -0.43, -0.43, -0.43, -0.03, -0.16, -0.87/

C...EULER'S BETA FUNCTION, REQUIRES ORDINARY GAMMA FUNCTION
      EULBET(X,Y)=GAMMA(X)*GAMMA(Y)/GAMMA(X+Y)

C...CONVERT LST SWITCHES TO MSPT, AND PARL PARAMETERS TO PARP
      MSTP57=1
      IF(LST(15).LT.0) MSTP57=0
      MSTP51=IABS(LST(15))
      MSTP52=LST(16)
      MSTP58=LST(12)
      PARP51=PARL(20)
      PARL(26)=0.

C...RESET OUTPUT ARRAY.
      DO 10  KFL=-6,6
   10 XPPR(KFL)=0.

      IF(MSTP51.EQ.1.OR.MSTP51.EQ.2) THEN
C...PROTON STRUCTURE FUNCTIONS FROM EICHTEN, HINCHLIFFE, LANE, QUIGG.
C...ALLOWED VARIABLE RANGE: 5 GEV^2 < Q^2 < 1E8 GEV^2; 1E-4 < X < 1

C...DETERMINE SET, LAMBDA AND X AND T EXPANSION VARIABLES.
        NSET=MSTP51
        IF(NSET.EQ.1) ALAM=0.2
        IF(NSET.EQ.2) ALAM=0.29
        TMIN=LOG(5./ALAM**2)
        TMAX=LOG(1E8/ALAM**2)
        IF(MSTP57.EQ.0) THEN
          T=TMIN
        ELSE
          T=LOG(MAX(1.,Q2/ALAM**2))
        ENDIF
        VT=MAX(-1.,MIN(1.,(2.*T-TMAX-TMIN)/(TMAX-TMIN)))
        NX=1
        IF(X.LE.0.1) NX=2
        IF(NX.EQ.1) VX=(2.*X-1.1)/0.9
        IF(NX.EQ.2) VX=MAX(-1.,(2.*LOG(X)+11.51293)/6.90776)
        CXS=1.
        IF(X.LT.1E-4.AND.ABS(PARP51-1.).GT.0.01) CXS=
     +  (1E-4/X)**(PARP51-1.)

C...CHEBYSHEV POLYNOMIALS FOR X AND T EXPANSION.
        TX(1)=1.
        TX(2)=VX
        TX(3)=2.*VX**2-1.
        TX(4)=4.*VX**3-3.*VX
        TX(5)=8.*VX**4-8.*VX**2+1.
        TX(6)=16.*VX**5-20.*VX**3+5.*VX
        TT(1)=1.
        TT(2)=VT
        TT(3)=2.*VT**2-1.
        TT(4)=4.*VT**3-3.*VT
        TT(5)=8.*VT**4-8.*VT**2+1.
        TT(6)=16.*VT**5-20.*VT**3+5.*VT

C...CALCULATE STRUCTURE FUNCTIONS.
        DO 30  KFL=1,6
          XQSUM=0.
          DO 20  IT=1,6
            DO 20  IX=1,6
   20     XQSUM=XQSUM+CEHLQ(IX,IT,NX,KFL,NSET)*TX(IX)*TT(IT)
   30   XQ(KFL)=XQSUM*(1.-X)**NEHLQ(KFL,NSET)*CXS

C...PUT INTO OUTPUT ARRAY.
        XPPR(0)=XQ(4)
        XPPR(1)=XQ(2)+XQ(3)
        XPPR(2)=XQ(1)+XQ(3)
        XPPR(3)=XQ(5)
        XPPR(4)=XQ(6)
        XPPR(-1)=XQ(3)
        XPPR(-2)=XQ(3)
        XPPR(-3)=XQ(5)
        XPPR(-4)=XQ(6)

C...SPECIAL EXPANSION FOR BOTTOM (THRESHOLD EFFECTS).
        IF(MSTP58.GE.5) THEN
          IF(NSET.EQ.1) TMIN=8.1905
          IF(NSET.EQ.2) TMIN=7.4474
          IF(T.GT.TMIN) THEN
            VT=MAX(-1.,MIN(1.,(2.*T-TMAX-TMIN)/(TMAX-TMIN)))
            TT(1)=1.
            TT(2)=VT
            TT(3)=2.*VT**2-1.
            TT(4)=4.*VT**3-3.*VT
            TT(5)=8.*VT**4-8.*VT**2+1.
            TT(6)=16.*VT**5-20.*VT**3+5.*VT
            XQSUM=0.
            DO 40  IT=1,6
              DO 40  IX=1,6
   40       XQSUM=XQSUM+CEHLQ(IX,IT,NX,7,NSET)*TX(IX)*TT(IT)
            XPPR(5)=XQSUM*(1.-X)**NEHLQ(7,NSET)*CXS
            XPPR(-5)=XPPR(5)
          ENDIF
        ENDIF

C...SPECIAL EXPANSION FOR TOP (THRESHOLD EFFECTS).
        IF(MSTP58.GE.6) THEN
          IF(NSET.EQ.1) TMIN=11.5528
          IF(NSET.EQ.2) TMIN=10.8097
          TMIN=TMIN+2.*LOG(PMAS(6,1)/30.)
          TMAX=TMAX+2.*LOG(PMAS(6,1)/30.)
          IF(T.GT.TMIN) THEN
            VT=MAX(-1.,MIN(1.,(2.*T-TMAX-TMIN)/(TMAX-TMIN)))
            TT(1)=1.
            TT(2)=VT
            TT(3)=2.*VT**2-1.
            TT(4)=4.*VT**3-3.*VT
            TT(5)=8.*VT**4-8.*VT**2+1.
            TT(6)=16.*VT**5-20.*VT**3+5.*VT
            XQSUM=0.
            DO 50  IT=1,6
              DO 50  IX=1,6
   50       XQSUM=XQSUM+CEHLQ(IX,IT,NX,8,NSET)*TX(IX)*TT(IT)
            XPPR(6)=XQSUM*(1.-X)**NEHLQ(8,NSET)*CXS
            XPPR(-6)=XPPR(6)
          ENDIF
        ENDIF

      ELSEIF(MSTP51.EQ.3.OR.MSTP51.EQ.4) THEN
C...PROTON STRUCTURE FUNCTIONS FROM DUKE, OWENS.
C...ALLOWED VARIABLE RANGE: 4 GEV^2 < Q^2 < APPROX 1E6 GEV^2.

C...DETERMINE SET, LAMBDA AND S EXPANSION PARAMETER.
        NSET=MSTP51-2
        IF(NSET.EQ.1) ALAM=0.2
        IF(NSET.EQ.2) ALAM=0.4
        IF(MSTP57.LE.0) THEN
          SD=0.
        ELSE
          Q2IN=MIN(1E6,MAX(4.,Q2))
          SD=LOG(LOG(Q2IN/ALAM**2)/LOG(4./ALAM**2))
        ENDIF

C...CALCULATE STRUCTURE FUNCTIONS.
        DO 70  KFL=1,5
          DO 60  IS=1,6
   60     TS(IS)=CDO(1,IS,KFL,NSET)+CDO(2,IS,KFL,NSET)*SD+ CDO(3,IS,
     +    KFL,NSET)*SD**2
          IF(KFL.LE.2) THEN
            XQ(KFL)=X**TS(1)*(1.-X)**TS(2)*(1.+TS(3)*X)/(EULBET(TS(1),
     +      TS(2)+1.)*(1.+TS(3)*TS(1)/(TS(1)+TS(2)+1.)))
          ELSE
            XQ(KFL)=TS(1)*X**TS(2)*(1.-X)**TS(3)*(1.+TS(4)*X+TS(5)*X**
     +      2+ TS(6)*X**3)
          ENDIF
   70   CONTINUE

C...PUT INTO OUTPUT ARRAYS.
        XPPR(0)=XQ(5)
        XPPR(1)=XQ(2)+XQ(3)/6.
        XPPR(2)=3.*XQ(1)-XQ(2)+XQ(3)/6.
        XPPR(3)=XQ(3)/6.
        XPPR(4)=XQ(4)
        XPPR(-1)=XQ(3)/6.
        XPPR(-2)=XQ(3)/6.
        XPPR(-3)=XQ(3)/6.
        XPPR(-4)=XQ(4)

      ELSEIF(MSTP51.GE.5.AND.MSTP51.LE.8) THEN
C...PROTON STRUCTURE FUNCTIONS FROM MORFIN AND TUNG.
C...ALLOWED VARIABLE RANGE: 4 GEV^2 < Q^2 < 1E8 GEV^2, 0 < X < 1.

C...CALCULATE EXPANSION PARAMETERS.
        NSET=MSTP51-4
        IF(NSET.EQ.1) ALAM=0.187
        IF(NSET.EQ.2) ALAM=0.212
        IF(NSET.EQ.3) ALAM=0.191
        IF(NSET.EQ.4) ALAM=0.155
        IF(MSTP57.EQ.0) THEN
          SD=0.
        ELSE
          SD=LOG(LOG(MAX(4.,Q2)/ALAM**2)/LOG(4./ALAM**2))
        ENDIF
        XL=LOG(MAX(1E-10,X))
        X1L=LOG(MAX(1E-10,1.-X))
        XLL=LOG(MAX(1E-10,LOG(1.+1./MAX(1E-10,X))))

C...CALCULATE STRUCTURE FUNCTIONS UP TO B.
        DO 90  IP=1,8
          DO 80  IEX=0,3
   80     EXMT(IEX)=CMT(IEX,0,IP,NSET)+CMT(IEX,1,IP,NSET)*SD+ CMT(IEX,
     +    2,IP,NSET)*SD**2
          EXMTS=EXMT(0)+EXMT(1)*XL+EXMT(2)*X1L+EXMT(3)*XLL
          IF(EXMTS.LT.-50.) THEN
            XQ(IP)=0.
          ELSE
            XQ(IP)=EXP(EXMTS)
          ENDIF
   90   CONTINUE
        IF(Q2.LE.2.25) XQ(7)=0.
        IF(Q2.LE.25.0) XQ(8)=0.

C...CALCULATE T STRUCTURE FUNCTION, SHIFTING EFFECTIVE Q SCALE FOR
C...NONDEFAULT T MASS, Q_ACTUAL = Q_NOMINAL * M_T_NOMINAL/M_T_ACTUAL.
        IF(MSTP57.EQ.0.OR.Q2.LE.PMAS(6,1)**2) THEN
          XQ(9)=0.
        ELSE
          SD=LOG(LOG(MAX(4.,Q2)/ALAM**2*(100./PMAS(6,1))**2)/
     +    LOG(4./ALAM**2))
          DO 100 IEX=0,3
  100     EXMT(IEX)=CMT(IEX,0,9,NSET)+CMT(IEX,1,9,NSET)*SD+
     +    CMT(IEX,2,9,NSET)*SD**2
          EXMTS=EXMT(0)+EXMT(1)*XL+EXMT(2)*X1L+EXMT(3)*XLL
          IF(EXMTS.LT.-50.) THEN
            XQ(9)=0.
          ELSE
            XQ(9)=EXP(EXMTS)
          ENDIF
        ENDIF

C...PUT INTO OUTPUT ARRAY.
        XPPR(0)=XQ(3)
        XPPR(1)=XQ(1)+XQ(5)
        XPPR(-1)=XQ(5)
        XPPR(2)=XQ(2)+XQ(4)
        XPPR(-2)=XQ(4)
        XPPR(3)=XQ(6)
        XPPR(-3)=XQ(6)
        XPPR(4)=XQ(7)
        XPPR(-4)=XQ(7)
        XPPR(5)=XQ(8)
        XPPR(-5)=XQ(8)
        XPPR(6)=XQ(9)
        XPPR(-6)=XQ(9)

      ELSEIF(MSTP51.EQ.9) THEN
C...LOWEST ORDER PARAMETRIZATION OF GLUCK, REYA, VOGT.
C...ALLOWED VARIABLE RANGE: 0.2 GEV^2 < Q2 < 1E6 GEV^2; 1E-4 < X < 1;
C...EXTENDED TO 0.2 GEV^2 < Q2 < 1E8 GEV^2; 1E-6 < X < 1
C...AFTER CONSULTATION WITH THE AUTHORS.

C...DETERMINE S AND X.
        ALAM=0.25
        IF(MSTP57.EQ.0) THEN
          SD=0.
        ELSE
          Q2IN=MIN(1E8,MAX(0.2,Q2))
          SD=LOG(LOG(Q2IN/ALAM**2)/LOG(0.2/ALAM**2))
        ENDIF
        XC=MAX(1E-6,X)
        XL=-LOG(XC)

C...CALCULATE STRUCTURE FUNCTIONS.
        XQ(1)=(0.794+0.312*SD)*XC**(0.427-0.011*SD)*
     +  (1.+(6.887-2.227*SD)*XC+(-11.083+2.136*SD)*XC**2+
     +  (3.900+1.079*SD)*XC**3)*(1.-XC)**(1.037+1.077*SD)
        XQ(2)=(0.486+0.139*SD)*XC**(0.434-0.018*SD)*
     +  (1.+(7.716-2.684*SD)*XC+(-12.768+3.122*SD)*XC**2+
     +  (4.564+0.470*SD)*XC**3)*(1.-XC)**(1.889+1.129*SD)
        XQ(3)=(XC**(0.415+0.186*SD)*((0.786+0.942*SD)+
     +  (5.256-5.810*SD)*XC+(-4.599+5.842*SD)*XC**2)+SD**0.592*
     +  EXP(-(0.398+2.135*SD)+SQRT(3.779*SD**1.250*XL)))*
     +  (1.-XC)**(1.622+1.980*SD)
        XQ(4)=SD**0.448*(1.-XC)**(5.540-0.445*SD)*
     +  EXP(-(4.668+1.230*SD)+SQRT((13.173-1.361*SD)*SD**0.442*XL))/
     +  XL**(3.181-0.862*SD)
        XQ(5)=0.
        IF(SD.GT.1.125) XQ(5)=(SD-1.125)*(1.-XC)**(2.038+1.022*SD)*
     +  EXP(-(4.290+1.569*SD)+SQRT((2.982+1.452*SD)*SD**0.5*XL))
        XQ(6)=0.
        IF(SD.GT.1.603) XQ(6)=(SD-1.603)*(1.-XC)**(2.230+1.052*SD)*
     +  EXP(-(4.566+1.559*SD)+SQRT((4.147+1.268*SD)*SD**0.5*XL))

C...PUT INTO OUTPUT ARRAY - SPECIAL FACTOR FOR SMALL X.
        CXS=1.
        IF(X.LT.1E-6.AND.ABS(PARP51-1.).GT.0.01)
     +  CXS=(1E-6/X)**(PARP51-1.)
        XPPR(0)=CXS*XQ(3)
        XPPR(1)=CXS*(XQ(2)+XQ(4))
        XPPR(-1)=CXS*XQ(4)
        XPPR(2)=CXS*(XQ(1)+XQ(4))
        XPPR(-2)=CXS*XQ(4)
        XPPR(3)=CXS*XQ(4)
        XPPR(-3)=CXS*XQ(4)
        XPPR(4)=CXS*XQ(5)
        XPPR(-4)=CXS*XQ(5)
        XPPR(5)=CXS*XQ(6)
        XPPR(-5)=CXS*XQ(6)

      ELSEIF(MSTP51.EQ.10) THEN
C...HIGHER ORDER PARAMETRIZATION OF GLUCK, REYA, VOGT.
C...ALLOWED VARIABLE RANGE: 0.2 GEV^2 < Q2 < 1E6 GEV^2; 1E-4 < X < 1;
C...EXTENDED TO 0.2 GEV^2 < Q2 < 1E8 GEV^2; 1E-6 < X < 1
C...AFTER CONSULTATION WITH THE AUTHORS.

C...DETERMINE S AND X.
        ALAM=0.20
        IF(MSTP57.EQ.0) THEN
          SD=0.
        ELSE
          Q2IN=MIN(1E8,MAX(0.2,Q2))
          SD=LOG(LOG(Q2IN/ALAM**2)/LOG(0.2/ALAM**2))
        ENDIF
        SD2=SD**2
        XC=MAX(1E-6,X)
        XL=-LOG(XC)

C...CALCULATE STRUCTURE FUNCTIONS.
        XQ(1)=(1.364+0.989*SD-0.236*SD2)*XC**(0.593-0.048*SD)*
     +  (1.+(8.912-6.092*SD+0.852*SD2)*XC+(-16.737+7.039*SD)*XC**2+
     +  (10.275+0.806*SD-2.000*SD2)*XC**3)*
     +  (1.-XC)**(2.043+1.408*SD-0.283*SD2)
        XQ(2)=(0.835+0.527*SD-0.144*SD2)*XC**(0.600-0.054*SD)*
     +  (1.+(10.245-7.821*SD+1.325*SD2)*XC+(-19.511+10.940*SD-
     +  1.133*SD2)*XC**2+(12.836-2.570*SD-1.041*SD2)*XC**3)*
     +  (1.-XC)**(3.083+1.382*SD-0.276*SD2)
        XQ(3)=(XC**(0.321-0.135*SD)*((10.51-2.299*SD)+
     +  (-17.28+0.755*SD)*XC+(8.242+2.543*SD)*XC**2)*
     +  XL**(-2.023-0.103*SD)+SD**1.044*
     +  EXP(-(-1.178+2.792*SD)+SQRT(2.318*SD**1.673*XL)))*
     +  (1.-XC)**(3.720+2.337*SD-0.199*SD2)
        XQ(4)=SD**0.761*(1.+(6.078-2.065*SD)*XC)*(1.-XC)**(4.654+
     +  0.603*SD-0.326*SD2)*EXP(-(4.231+1.036*SD)+SQRT(3.419*SD**0.316*
     +  XL))/XL**(0.897-0.618*SD)
        XQ(5)=0.
        IF(SD.GT.0.918) XQ(5)=(SD-0.918)*(1.-XC)**(3.328+0.859*SD)*
     +  EXP(-(3.837+1.504*SD)+SQRT((2.150+1.291*SD)*SD**0.5*XL))
        XQ(6)=0.
        IF(SD.GT.1.353) XQ(6)=(SD-1.353)*(1.-XC)**(3.382+0.909*SD)*
     +  EXP(-(4.130+1.486*SD)+SQRT((2.895+1.240*SD)*SD**0.5*XL))

C...PUT INTO OUTPUT ARRAY - SPECIAL FACTOR FOR SMALL X.
        CXS=1.
        IF(X.LT.1E-6.AND.ABS(PARP51-1.).GT.0.01)
     +  CXS=(1E-6/X)**(PARP51-1.)
        XPPR(0)=CXS*XQ(3)
        XPPR(1)=CXS*(XQ(2)+XQ(4))
        XPPR(-1)=CXS*XQ(4)
        XPPR(2)=CXS*(XQ(1)+XQ(4))
        XPPR(-2)=CXS*XQ(4)
        XPPR(3)=CXS*XQ(4)
        XPPR(-3)=CXS*XQ(4)
        XPPR(4)=CXS*XQ(5)
        XPPR(-4)=CXS*XQ(5)
        XPPR(5)=CXS*XQ(6)
        XPPR(-5)=CXS*XQ(6)
      ENDIF
      PARL(26)=ALAM

      RETURN
      END
*CMZ :  1.01/50 23/05/96  12.34.50  by  Piero Zucchelli
*-- Author :    Piero Zucchelli   20/03/96

      SUBROUTINE RANMAR(RVEC,ISEQ)
      DIMENSION RVEC(*)
      CALL RANLUX(RVEC,ISEQ)
      RETURN
      END
*CMZ :  1.02/11 14/01/97  23.32.59  by  P. Zucchelli
*CMZ :  1.02/04 13/01/97  14.45.27  by  P. Zucchelli
*CMZ :  1.01/51 13/06/96  18.29.39  by  Piero Zucchelli
*CMZ :  1.01/50 26/04/96  14.52.50  by  Piero Zucchelli
*CMZ :  1.01/45 08/01/96  09.37.30  by  Piero Zucchelli
*CMZ :  1.01/43 15/12/95  17.58.41  by  Piero Zucchelli
*CMZ :  1.01/41 14/12/95  14.48.28  by  Piero Zucchelli
*CMZ :  1.01/40 12/12/95  12.30.14  by  Piero Zucchelli
*CMZ :  1.01/39 02/11/95  18.35.07  by  Piero Zucchelli
*CMZ :  1.01/37 21/09/95  11.21.48  BY  PIERO ZUCCHELLI
*CMZ :  1.01/36 26/07/95  18.29.05  BY  PIERO ZUCCHELLI
*CMZ :  1.01/33 03/07/95  17.33.30  BY  PIERO ZUCCHELLI
*CMZ :  1.01/27 02/06/95  14.58.52  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.35.13  BY  PIERO ZUCCHELLI
*CMZ :  1.01/01 14/09/94  16.25.32  BY  PIERO ZUCCHELLI
*CMZ :  1.01/00 01/09/94  17.36.14  BY  PIERO ZUCCHELLI
*-- AUTHOR :    PIERO ZUCCHELLI   01/09/94
      SUBROUTINE READFFKY
*KEEP,KEYS.
        COMMON/CFREAD/SPACE(5000)
 	COMMON/MYKEYS/IEVAR,IF5CC,INEUT,IINTE,IFERM,IFLAT,ICOUN,
     &  REFIX,IMUDO,IKAT1,IKAT2,IKAT3,IKAT4,IKAT5,IKAT6,
     &  INEVT,IJAK1,IJAK2,IITDK,RPTAU,RXK0D,NTGR,IDIMUON,IGLU,
     &  IQDEN,LOME(2),IFILES,ICCHA,ISEED,IDSUBS,IONEM,EHAC,IPSEL,
     &  IHIST


*KEND.
      CALL FFINIT(3000)

      CALL FFKEY('PSEL',IPSEL,1,'INTEGER')
      CALL FFKEY('SEED',ISEED,1,'INTEGER')
      CALL FFKEY('IGLU',IGLU,1,'INTEGER')
      CALL FFKEY('EVAR',IEVAR,1,'INTEGER')
      CALL FFKEY('F5CC',IF5CC,1,'INTEGER')
      CALL FFKEY('NEUT',INEUT,1,'INTEGER')
      CALL FFKEY('INTE',IINTE,1,'INTEGER')
      CALL FFKEY('FERM',IFERM,1,'INTEGER')
      CALL FFKEY('FLAT',IFLAT,1,'INTEGER')
      CALL FFKEY('COUN',ICOUN,1,'INTEGER')
      CALL FFKEY('HIST',IHIST,1,'INTEGER')
      CALL FFKEY('EFIX',REFIX,1,'REAL')
      CALL FFKEY('QDEN',IQDEN,1,'INTEGER')
      CALL FFKEY('MUDO',IMUDO,1,'INTEGER')
      CALL FFKEY('NGTR',NTGR,1,'INTEGER')
      CALL FFKEY('DIMU',IDIMUON,1,'INTEGER')
      CALL FFKEY('CCHA',ICCHA,1,'INTEGER')
      CALL FFKEY('LOME',LOME,2,'INTEGER')
      CALL FFKEY('FILE',IFILES,1,'INTEGER')
      CALL FFKEY('KAT1',IKAT1,1,'INTEGER')
      CALL FFKEY('KAT2',IKAT2,1,'INTEGER')
      CALL FFKEY('KAT3',IKAT3,1,'INTEGER')
      CALL FFKEY('KAT4',IKAT4,1,'INTEGER')
      CALL FFKEY('KAT5',IKAT5,1,'INTEGER')
      CALL FFKEY('KAT6',IKAT6,1,'INTEGER')
      CALL FFKEY('NEVT',INEVT,1,'INTEGER')
      CALL FFKEY('JAK1',IJAK1,1,'INTEGER')
      CALL FFKEY('JAK2',IJAK2,1,'INTEGER')
      CALL FFKEY('ITDK',IITDK,1,'INTEGER')
      CALL FFKEY('DSUBS',IDSUBS,1,'INTEGER')
      CALL FFKEY('PTAU',RPTAU,1,'REAL')
      CALL FFKEY('XK0D',RXK0D,1,'REAL')
      CALL FFKEY('EHAC',EHAC,1,'REAL')

      NINP=17
C     OPEN(NINP,FILE="./jetta.crd",STATUS='OLD')
      OPEN(NINP,FILE="./jetta.crd",STATUS='UNKNOWN')
      CALL FFSET('LINP',NINP)
      CALL FFGO
      CLOSE(NINP)

      IF (NTGR.EQ.0) NTGR=1
      IF (IITDK.EQ.0) IITDK=1
      IF (RXK0D.EQ.0.) RXK0D=0.001


      RETURN
      END
*CMZ :  1.02/03 13/01/97  13.46.08  by  P. Zucchelli
*CMZ :  1.01/50 23/05/96  10.19.16  by  Piero Zucchelli
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE RESLU
C     ****************
C INITIALIZE LUND COMMON
      PARAMETER (NMXHEP=2000)
*KEEP,HEPEVT.
      DOUBLE PRECISION PHEP,VHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      SAVE /HEPEVT/

*KEND.
      NHEP=0
      END
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      DOUBLE PRECISION FUNCTION RIWFUN(V)
      DOUBLE PRECISION V(2)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
     &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      DATA V2MIN,V2MAX/2*0./

      RIWFUN=0.D0
      V1MIN=XMIN
      V1MAX=XMAX
      IF(LST(31).EQ.1) THEN
        V2MIN=Q2MIN
        V2MAX=Q2MAX
      ELSEIF(LST(31).EQ.2) THEN
        V2MIN=YMIN
        V2MAX=YMAX
      ELSEIF(LST(31).EQ.3) THEN
        V2MIN=W2MIN
        V2MAX=W2MAX
      ENDIF
      V1=V1MIN+V(1)*(V1MAX-V1MIN)
      V2=V2MIN+V(2)*(V2MAX-V2MIN)
      RIWFUN=DCROSS(V1,V2)*(V1MAX-V1MIN)*(V2MAX-V2MIN)

      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 04/07/94  15.02.27  BY  PIERO ZUCCHELLI
*-- AUTHOR :
C **********************************************************************

      SUBROUTINE RIWIBD
C   BLOCK DATA SUBSTITUTE FROM RIWIAD
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/STORE/XA(11),XB(11),XC(11),XD(11),MA(11),MB(11),MC(11)
      COMMON/STORE1/R(10000),LR
      COMMON/OPTION/IPRRIW,ICONV,IRESET
      COMMON/RANDOM/NSHOTS
      COMMON/INTERN/FACTOR,ALFA,BETA,GAMMA,DELTA,LEVEL,NMIN
      COMMON /LPFLAG/ LST3
      DATA INIT/0/
      IF(INIT.EQ.1) RETURN
      INIT=1
      MA(1)=0
      LR=10000
      ICONV=1
      IRESET=0
      NSHOTS=2
      FACTOR=1.65
      LEVEL=90
      ALFA=0.3
      BETA=0.3
      GAMMA=0.3
      DELTA=.7
      NMIN=2
C...PRINT FLAG TO BE CHANGED HERE.
      IPRRIW=0
      IF(LST3.GE.4) WRITE(6,10000) IPRRIW
      RETURN
10000 FORMAT(5X,'RIWIAD PRINT FLAG CHANGED: IPRRIW =',I5)
      END
*CMZ :  1.01/50 23/05/96  12.34.50  by  Piero Zucchelli
*-- Author :    Piero Zucchelli   20/03/96

      REAL*4 FUNCTION RLU(IDUMMY)
      CALL RANLUX(RTIM,1)
      RLU=RTIM
      RETURN
      END
*CMZ :  1.01/50 23/05/96  12.34.50  by  Piero Zucchelli
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE ROTOD1(PH1,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     USED BY : KORALZ
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PVEC(4),QVEC(4),RVEC(4)
C
      PHI=PH1
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)=RVEC(1)
      QVEC(2)= CS*RVEC(2)-SN*RVEC(3)
      QVEC(3)= SN*RVEC(2)+CS*RVEC(3)
      QVEC(4)=RVEC(4)
      RETURN
      END
*-- Author :    Piero Zucchelli   20/03/96

      REAL*4 FUNCTION RNDMM(IDUMMY)
      CALL RANLUX(RTIM,1)
      RNDMM=RTIM
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE ROTOD2(PH1,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     USED BY : KORALZ RADKOR
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PVEC(4),QVEC(4),RVEC(4)
C
      PHI=PH1
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)= CS*RVEC(1)+SN*RVEC(3)
      QVEC(2)=RVEC(2)
      QVEC(3)=-SN*RVEC(1)+CS*RVEC(3)
      QVEC(4)=RVEC(4)
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE ROTOD3(PH1,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     USED BY : KORALZ RADKOR
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION PVEC(4),QVEC(4),RVEC(4)
      PHI=PH1
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)= CS*RVEC(1)-SN*RVEC(2)
      QVEC(2)= SN*RVEC(1)+CS*RVEC(2)
      QVEC(3)=RVEC(3)
      QVEC(4)=RVEC(4)
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE ROTOR1(PH1,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     CALLED BY :
C ----------------------------------------------------------------------
      REAL*4 PVEC(4),QVEC(4),RVEC(4)
C
      PHI=PH1
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)=RVEC(1)
      QVEC(2)= CS*RVEC(2)-SN*RVEC(3)
      QVEC(3)= SN*RVEC(2)+CS*RVEC(3)
      QVEC(4)=RVEC(4)
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE ROTOR2(PH1,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     USED BY : TAUOLA
C ----------------------------------------------------------------------
      IMPLICIT REAL*4(A-H,O-Z)
      REAL*4 PVEC(4),QVEC(4),RVEC(4)
C
      PHI=PH1
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)= CS*RVEC(1)+SN*RVEC(3)
      QVEC(2)=RVEC(2)
      QVEC(3)=-SN*RVEC(1)+CS*RVEC(3)
      QVEC(4)=RVEC(4)
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE ROTOR3(PHI,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     USED BY : TAUOLA
C ----------------------------------------------------------------------
      REAL*4 PVEC(4),QVEC(4),RVEC(4)
C
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)= CS*RVEC(1)-SN*RVEC(2)
      QVEC(2)= SN*RVEC(1)+CS*RVEC(2)
      QVEC(3)=RVEC(3)
      QVEC(4)=RVEC(4)
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE ROTPOL(THET,PHI,PP)
C ----------------------------------------------------------------------
C
C     CALLED BY : DADMAA,DPHSAA
C ----------------------------------------------------------------------
      REAL  PP(4)
C
      CALL ROTOR2(THET,PP,PP)
      CALL ROTOR3( PHI,PP,PP)
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE ROTPOX(THET,PHI,PP)
      IMPLICIT REAL*8 (A-H,O-Z)
C ----------------------------------------------------------------------
C
C ----------------------------------------------------------------------
      DIMENSION PP(4)
C
      CALL ROTOD2(THET,PP,PP)
      CALL ROTOD3( PHI,PP,PP)
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION SIGEE(Q2,JNP)
C ----------------------------------------------------------------------
C  E+E- CROSS SECTION IN THE (1.GEV2,AMTAU**2) REGION
C  NORMALISED TO SIG0 = 4/3 PI ALFA2
C  USED IN MATRIX ELEMENT FOR MULTIPION TAU DECAYS
C  CF YS.TSAI        PHYS.REV D4 ,2821(1971)
C     F.GILMAN ET AL PHYS.REV D17,1846(1978)
C     C.KIESLING, TO BE PUB. IN HIGH ENERGY E+E- PHYSICS (1988)
C  DATSIG(*,1) = E+E- -> PI+PI-2PI0
C  DATSIG(*,2) = E+E- -> 2PI+2PI-
C  DATSIG(*,3) = 5-PION CONTRIBUTION (A LA TN.PHAM ET AL)
C                (PHYS LETT 78B,623(1978)
C  DATSIG(*,5) = E+E- -> 6PI
C
C  4- AND 6-PION CROSS SECTIONS FROM DATA
C  5-PION CONTRIBUTION RELATED TO 4-PION CROSS SECTION
C
C     CALLED BY DPHNPI
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      REAL*4 DATSIG(17,6)
C
      DATA DATSIG/
     +  7.40,12.00,16.15,21.25,24.90,29.55,34.15,37.40,37.85,37.40,
     + 36.00,33.25,30.50,27.70,24.50,21.25,18.90,
     +  1.24, 2.50, 3.70, 5.40, 7.45,10.75,14.50,18.20,22.30,28.90,
     + 29.35,25.60,22.30,18.60,14.05,11.60, 9.10,
     + 17*.0,
     + 17*.0,
     + 9*.0,.65,1.25,2.20,3.15,5.00,5.75,7.80,8.25,
     + 17*.0/
      DATA SIG0 / 86.8 /
      DATA PI /3.141592653589793238462643/
      DATA INIT / 0 /
C
      SIGEE = 0.
      JNPI=JNP
      IF(JNP.EQ.4) JNPI=3
      IF(JNP.EQ.3) JNPI=4
      IF(INIT.EQ.0) THEN
        INIT=1
        AMPI2=AMPI**2
        FPI = .943*AMPI
        DO 30  I=1,17
          DATSIG(I,2) = DATSIG(I,2)/2.
          DATSIG(I,1) = DATSIG(I,1) + DATSIG(I,2)
          S = 1.025+(I-1)*.05
          FACT=0.
          S2=S**2
          DO 10  J=1,17
            T= 1.025+(J-1)*.05
            IF(T . GT. S-AMPI ) GO TO 20
            T2=T**2
            FACT=(T2/S2)**2*SQRT((S2-T2-AMPI2)**2-4.*T2*AMPI2)/S2 *2.*
     +      T*.05
            FACT = FACT * (DATSIG(J,1)+DATSIG(J+1,1))
   10     DATSIG(I,3) = DATSIG(I,3) + FACT
   20     DATSIG(I,3) = DATSIG(I,3) /(2*PI*FPI)**2
          DATSIG(I,4) = DATSIG(I,3)
          DATSIG(I,6) = DATSIG(I,5)
   30   CONTINUE
C       WRITE(6,1000) DATSIG
10000   FORMAT(///1X,' EE SIGMA USED IN MULTIPI DECAYS'/
     +        (17F7.2/))
      ENDIF
      Q=SQRT(Q2)
      QMIN=1.
      IF(Q.LT.QMIN) THEN
        SIGEE=DATSIG(1,JNPI)+
     +       (DATSIG(2,JNPI)-DATSIG(1,JNPI))*(Q-1.)/.05
      ELSEIF(Q.LT.1.8) THEN
        DO 40 I=1,16
          QMAX = QMIN + .05
          IF(Q.LT.QMAX) GO TO 50
          QMIN = QMIN + .05
   40   CONTINUE
   50   SIGEE=DATSIG(I,JNPI)+
     +       (DATSIG(I+1,JNPI)-DATSIG(I,JNPI)) * (Q-QMIN)/.05
      ELSEIF(Q.GT.1.8) THEN
        SIGEE=DATSIG(17,JNPI)+
     +       (DATSIG(17,JNPI)-DATSIG(16,JNPI)) * (Q-1.8)/.05
      ENDIF
      IF(SIGEE.LT..0) SIGEE=0.
C
      SIGEE = SIGEE/(6.*PI**2*SIG0)
C
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION SIGOLD(Q2,JNPI)
C ----------------------------------------------------------------------
C  E+E- CROSS SECTION IN THE (1.GEV2,AMTAU**2) REGION
C  NORMALISED TO SIG0 = 4/3 PI ALFA2
C  USED IN MATRIX ELEMENT FOR MULTIPION TAU DECAYS
C  CF YS.TSAI        PHYS.REV D4 ,2821(1971)
C     F.GILMAN ET AL PHYS.REV D17,1846(1978)
C     C.KIESLING, TO BE PUB. IN HIGH ENERGY E+E- PHYSICS (1988)
C  DATSIG(*,1) = E+E- -> PI+PI-2PI0
C  DATSIG(*,2) = E+E- -> 2PI+2PI-
C  DATSIG(*,3) = 5-PION CONTRIBUTION (A LA TN.PHAM ET AL)
C                (PHYS LETT 78B,623(1978)
C  DATSIG(*,4) = E+E- -> 6PI
C
C  4- AND 6-PION CROSS SECTIONS FROM DATA
C  5-PION CONTRIBUTION RELATED TO 4-PION CROSS SECTION
C
C     CALLED BY DPHNPI
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      REAL*4 DATSIG(17,4)
C
      DATA DATSIG/
     +  7.40,12.00,16.15,21.25,24.90,29.55,34.15,37.40,37.85,37.40,
     + 36.00,33.25,30.50,27.70,24.50,21.25,18.90,
     +  1.24, 2.50, 3.70, 5.40, 7.45,10.75,14.50,18.20,22.30,28.90,
     + 29.35,25.60,22.30,18.60,14.05,11.60, 9.10,
     + 17*.0,
     + 9*.0,.65,1.25,2.20,3.15,5.00,5.75,7.80,8.25/
      DATA SIG0 / 86.8 /
      DATA PI /3.141592653589793238462643/
      DATA INIT / 0 /
C
      IF(INIT.EQ.0) THEN
        INIT=1
        AMPI2=AMPI**2
        FPI = .943*AMPI
        DO 30  I=1,17
          DATSIG(I,2) = DATSIG(I,2)/2.
          DATSIG(I,1) = DATSIG(I,1) + DATSIG(I,2)
          S = 1.025+(I-1)*.05
          FACT=0.
          S2=S**2
          DO 10  J=1,17
            T= 1.025+(J-1)*.05
            IF(T . GT. S-AMPI ) GO TO 20
            T2=T**2
            FACT=(T2/S2)**2*SQRT((S2-T2-AMPI2)**2-4.*T2*AMPI2)/S2 *2.*
     +      T*.05
            FACT = FACT * (DATSIG(J,1)+DATSIG(J+1,1))
   10     DATSIG(I,3) = DATSIG(I,3) + FACT
   20     DATSIG(I,3) = DATSIG(I,3) /(2*PI*FPI)**2
   30   CONTINUE
C       WRITE(6,1000) DATSIG
10000   FORMAT(///1X,' EE SIGMA USED IN MULTIPI DECAYS'/
     +        (17F7.2/))
      ENDIF
      Q=SQRT(Q2)
      QMIN=1.
      IF(Q.LT.QMIN) THEN
        SIGOL=DATSIG(1,JNPI)+
     +       (DATSIG(2,JNPI)-DATSIG(1,JNPI))*(Q-1.)/.05
      ELSEIF(Q.LT.1.8) THEN
        DO 40 I=1,16
          QMAX = QMIN + .05
          IF(Q.LT.QMAX) GO TO 50
          QMIN = QMIN + .05
   40   CONTINUE
   50   SIGOL=DATSIG(I,JNPI)+
     +       (DATSIG(I+1,JNPI)-DATSIG(I,JNPI)) * (Q-QMIN)/.05
      ELSEIF(Q.GT.1.8) THEN
        SIGOL=DATSIG(17,JNPI)+
     +       (DATSIG(17,JNPI)-DATSIG(16,JNPI)) * (Q-1.8)/.05
      ENDIF
      IF(SIGOL.LT..0) SIGOL=0.
C
      SIGOL = SIGOL/(6.*PI**2*SIG0)
      SIGOLD=SIGOL
C
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE SPHERA(R,X)
C ----------------------------------------------------------------------
C GENERATES UNIFORMLY THREE-VECTOR X ON SPHERE  OF RADIUS R
C
C     CALLED BY : DPHSXX,DADMPI,DADMKK
C ----------------------------------------------------------------------
      REAL  X(4)
      REAL*4 RRR(2)
      DATA PI /3.141592653589793238462643/
C
      CALL RANMAR(RRR,2)
      COSTH=-1.+2.*RRR(1)
      SINTH=SQRT(1.-COSTH**2)
      X(1)=R*SINTH*COS(2*PI*RRR(2))
      X(2)=R*SINTH*SIN(2*PI*RRR(2))
      X(3)=R*COSTH
      RETURN
      END
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE SPHERD(R,X)
C ----------------------------------------------------------------------
C GENERATES UNIFORMLY THREE-VECTOR X ON SPHERE  OF RADIUS R
C DOUBLE PRECISON VERSION OF SPHERA
C ----------------------------------------------------------------------
      REAL*8  R,X(4),PI,COSTH,SINTH
      REAL*4 RRR(2)
      DATA PI /3.141592653589793238462643D0/
C
      CALL RANMAR(RRR,2)
      COSTH=-1+2*RRR(1)
      SINTH=SQRT(1 -COSTH**2)
      X(1)=R*SINTH*COS(2*PI*RRR(2))
      X(2)=R*SINTH*SIN(2*PI*RRR(2))
      X(3)=R*COSTH
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION SQM2(ITDKRC,QP,XN,XA,XK,AK0,HV)
C
C **********************************************************************
C     REAL PHOTON MATRIX ELEMENT SQUARED                               *
C     PARAMETERS:                                                      *
C     HV- POLARIMETRIC FOUR-VECTOR OF TAU                              *
C     QP,XN,XA,XK - 4-MOMENTA OF ELECTRON (MUON), NU, NUBAR AND PHOTON *
C                   ALL FOUR-VECTORS IN TAU REST FRAME (IN GEV)        *
C     AK0 - INFRARED CUTOFF, MINIMAL ENERGY OF HARD PHOTONS (GEV)      *
C     SQM2 - VALUE FOR S=0                                             *
C     SEE EQS. (2.9)-(2.10) FROM CJK ( NUCL.PHYS.B(1991) )             *
C **********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
      REAL*8    QP(4),XN(4),XA(4),XK(4)
      REAL*8    R(4)
      REAL*8   HV(4)
      REAL*8 S0(3),RXA(3),RXK(3),RQP(3)
      DATA PI /3.141592653589793238462643D0/
C
      TMASS=AMTAU
      GF=GFERMI
      ALPHAI=ALFINV
      TMASS2=TMASS**2
      EMASS2=QP(4)**2-QP(1)**2-QP(2)**2-QP(3)**2
      R(4)=TMASS
C     SCALAR PRODUCTS OF FOUR-MOMENTA
      DO 10 I=1,3
        R(1)=0.D0
        R(2)=0.D0
        R(3)=0.D0
        R(I)=TMASS
        RXA(I)=R(4)*XA(4)-R(1)*XA(1)-R(2)*XA(2)-R(3)*XA(3)
C       RXN(I)=R(4)*XN(4)-R(1)*XN(1)-R(2)*XN(2)-R(3)*XN(3)
        RXK(I)=R(4)*XK(4)-R(1)*XK(1)-R(2)*XK(2)-R(3)*XK(3)
        RQP(I)=R(4)*QP(4)-R(1)*QP(1)-R(2)*QP(2)-R(3)*QP(3)
   10 CONTINUE
      QPXN=QP(4)*XN(4)-QP(1)*XN(1)-QP(2)*XN(2)-QP(3)*XN(3)
      QPXA=QP(4)*XA(4)-QP(1)*XA(1)-QP(2)*XA(2)-QP(3)*XA(3)
      QPXK=QP(4)*XK(4)-QP(1)*XK(1)-QP(2)*XK(2)-QP(3)*XK(3)
C     XNXA=XN(4)*XA(4)-XN(1)*XA(1)-XN(2)*XA(2)-XN(3)*XA(3)
      XNXK=XN(4)*XK(4)-XN(1)*XK(1)-XN(2)*XK(2)-XN(3)*XK(3)
      XAXK=XA(4)*XK(4)-XA(1)*XK(1)-XA(2)*XK(2)-XA(3)*XK(3)
      TXN=TMASS*XN(4)
      TXA=TMASS*XA(4)
      TQP=TMASS*QP(4)
      TXK=TMASS*XK(4)
C
      X= XNXK/QPXN
      Z= TXK/TQP
      A= 1+X
      B= 1+ X*(1+Z)/2+Z/2
      S1= QPXN*TXA*( -EMASS2/QPXK**2*A + 2*TQP/(QPXK*TXK)*B-
     +TMASS2/TXK**2)  +
     +QPXN/TXK**2* ( TMASS2*XAXK - TXA*TXK+ XAXK*TXK) -
     +TXA*TXN/TXK - QPXN/(QPXK*TXK)* (TQP*XAXK-TXK*QPXA)
      CONST4=256*PI/ALPHAI*GF**2
      IF (ITDKRC.EQ.0) CONST4=0D0
      SQM2=S1*CONST4
      DO 20 I=1,3
        S0(I) = QPXN*RXA(I)*(-EMASS2/QPXK**2*A + 2*TQP/(QPXK*TXK)*B-
     +  TMASS2/TXK**2) +
     +  QPXN/TXK**2* (TMASS2*XAXK - TXA*RXK(I)+ XAXK*RXK(I))-
     +  RXA(I)*TXN/TXK - QPXN/(QPXK*TXK)*(RQP(I)*XAXK- RXK(I)*QPXA)
   20 HV(I)=S0(I)/S1-1.D0
      RETURN
      END
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 09/08/94  17.43.59  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE TAUFIL
C     *****************
C SUBSITUTE OF TAU PRODUCTION GENERATOR
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / IDFC  / IDFF
C POSITIONS OF TAUS IN THE LUND COMMON BLOCK
C IT WILL BE USED BY TAUOLA OUTPUT ROUTINES.
      COMMON /TAUPOS / NPA,NPB
      DIMENSION XPB1(4),XPB2(4),AQF1(4),AQF2(4)
C
C --- DEFINING DUMMY EVENTS MOMENTA
      DO 10 K=1,3
        XPB1(K)=0.0
        XPB2(K)=0.0
        AQF1(K)=0.0
        AQF2(K)=0.0
   10 CONTINUE
      AQF1(4)=AMTAU
      AQF2(4)=AMTAU
C --- TAU MOMENTA
      CALL TRALO4(1,AQF1,AQF1,AM)
      CALL TRALO4(2,AQF2,AQF2,AM)
C --- BEAMS MOMENTA AND IDENTIFIERS
      KFB1= 11*IDFF/IABS(IDFF)
      KFB2=-11*IDFF/IABS(IDFF)
      XPB1(4)= AQF1(4)
      XPB1(3)= AQF1(4)
      IF(AQF1(3).NE.0.0) XPB1(3)= AQF1(4)*AQF1(3)/ABS(AQF1(3))
      XPB2(4)= AQF2(4)
      XPB2(3)=-AQF2(4)
      IF(AQF2(3).NE.0.0) XPB2(3)= AQF2(4)*AQF2(3)/ABS(AQF2(3))
C --- POSITION OF FIRST AND SECOND TAU IN LUND COMMON
      NPA=3
      NPB=4
C --- FILL TO LUND COMMON
      CALL FILHEP(  1,3, KFB1,0,0,0,0,XPB1, AMEL,.TRUE.)
      CALL FILHEP(  2,3, KFB2,0,0,0,0,XPB2, AMEL,.TRUE.)
      CALL FILHEP(NPA,1, IDFF,1,2,0,0,AQF1,AMTAU,.TRUE.)
      CALL FILHEP(NPB,1,-IDFF,1,2,0,0,AQF2,AMTAU,.TRUE.)
      END
*CMZ :  1.01/15 14/05/95  11.41.24  BY  PIERO ZUCCHELLI
*CMZ :  1.01/11 13/05/95  18.41.54  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  12.08.29  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE TAUINIT
C     **************
C NOTE THAT THE ROUTINES ARE NOT LIKE IN CPC DECK THIS IS HISTORICAL !!
C=======================================================================
C====================== DECTES    : TEST OF TAU DECAY LIBRARY===========
C====================== KTORY = 1 : INTERFACE OF KORAL-Z TYPE ==========
C====================== KTORY = 2 : INTERFACE OF KORAL-B TYPE =========
C=======================================================================
      COMMON  / / BLAN(10000)
      COMMON / INOUT / INUT,IOUT
      INUT=5
      IOUT=6
      KTORY=1
      CALL DECTES(KTORY)
C      KTORY=2
C      CALL DECTES(KTORY)
      END
*CMZ :  1.01/14 14/05/95  11.26.28  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.27  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 09/08/94  17.43.59  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE TAURDF(KTO)
C THIS ROUTINE CAN BE CALLED BEFORE ANY TAU+ OR TAU- EVENT IS GENERATED
C IT CAN BE USED TO GENERATE TAU+ AND TAU- SAMPLES OF DIFFERENT
C CONTENTS
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      COMMON / TAUBRA / GAMPRT(30),JLIST(30),NCHAN
      IF (KTO.EQ.1) THEN
C     ==================
C LIST OF BRANCHING RATIOS
        NCHAN = 19
        DO 10 I = 1,30
          IF (I.LE.NCHAN) THEN
            JLIST(I) = I
            IF(I.EQ. 1) GAMPRT(I) = .0000
            IF(I.EQ. 2) GAMPRT(I) = .0000
            IF(I.EQ. 3) GAMPRT(I) = .0000
            IF(I.EQ. 4) GAMPRT(I) = .0000
            IF(I.EQ. 5) GAMPRT(I) = .0000
            IF(I.EQ. 6) GAMPRT(I) = .0000
            IF(I.EQ. 7) GAMPRT(I) = .0000
            IF(I.EQ. 8) GAMPRT(I) = 1.0000
            IF(I.EQ. 9) GAMPRT(I) = 1.0000
            IF(I.EQ.10) GAMPRT(I) = 1.0000
            IF(I.EQ.11) GAMPRT(I) = 1.0000
            IF(I.EQ.12) GAMPRT(I) = 1.0000
            IF(I.EQ.13) GAMPRT(I) = 1.0000
            IF(I.EQ.14) GAMPRT(I) = 1.0000
            IF(I.EQ.15) GAMPRT(I) = 1.0000
            IF(I.EQ.16) GAMPRT(I) = 1.0000
            IF(I.EQ.17) GAMPRT(I) = 1.0000
            IF(I.EQ.18) GAMPRT(I) = 1.0000
            IF(I.EQ.19) GAMPRT(I) = 1.0000
          ELSE
            JLIST(I) = 0
            GAMPRT(I) = 0.
          ENDIF
   10   CONTINUE
C --- COEFFICIENTS TO FIX RATIO OF:
C --- A1 3CHARGED/ A1 1CHARGED 2 NEUTRALS MATRIX ELEMENTS (MASLESS LIM.)
C --- PROBABILITY OF K0 TO BE KS
C --- PROBABILITY OF K0B TO BE KS
C --- RATIO OF COEFFICIENTS FOR K*--> K0 PI-
C --- ALL COEFFICENTS SHOULD BE IN THE RANGE (0.0,1.0)
C --- THEY MEANING IS PROBABILITY OF THE FIRST CHOICE ONLY IF ONE
C --- NEGLECTS MASS-PHASE SPACE EFFECTS
        BRA1=0.5
        BRK0=0.5
        BRK0B=0.5
        BRKS=0.6667
      ELSE
C     ====
C LIST OF BRANCHING RATIOS
        NCHAN = 19
        DO 20 I = 1,30
          IF (I.LE.NCHAN) THEN
            JLIST(I) = I
            IF(I.EQ. 1) GAMPRT(I) = .0000
            IF(I.EQ. 2) GAMPRT(I) = .0000
            IF(I.EQ. 3) GAMPRT(I) = .0000
            IF(I.EQ. 4) GAMPRT(I) = .0000
            IF(I.EQ. 5) GAMPRT(I) = .0000
            IF(I.EQ. 6) GAMPRT(I) = .0000
            IF(I.EQ. 7) GAMPRT(I) = .0000
            IF(I.EQ. 8) GAMPRT(I) = 1.0000
            IF(I.EQ. 9) GAMPRT(I) = 1.0000
            IF(I.EQ.10) GAMPRT(I) = 1.0000
            IF(I.EQ.11) GAMPRT(I) = 1.0000
            IF(I.EQ.12) GAMPRT(I) = 1.0000
            IF(I.EQ.13) GAMPRT(I) = 1.0000
            IF(I.EQ.14) GAMPRT(I) = 1.0000
            IF(I.EQ.15) GAMPRT(I) = 1.0000
            IF(I.EQ.16) GAMPRT(I) = 1.0000
            IF(I.EQ.17) GAMPRT(I) = 1.0000
            IF(I.EQ.18) GAMPRT(I) = 1.0000
            IF(I.EQ.19) GAMPRT(I) = 1.0000
          ELSE
            JLIST(I) = 0
            GAMPRT(I) = 0.
          ENDIF
   20   CONTINUE
C --- COEFFICIENTS TO FIX RATIO OF:
C --- A1 3CHARGED/ A1 1CHARGED 2 NEUTRALS MATRIX ELEMENTS (MASLESS LIM.)
C --- PROBABILITY OF K0 TO BE KS
C --- PROBABILITY OF K0B TO BE KS
C --- RATIO OF COEFFICIENTS FOR K*--> K0 PI-
C --- ALL COEFFICENTS SHOULD BE IN THE RANGE (0.0,1.0)
C --- THEY MEANING IS PROBABILITY OF THE FIRST CHOICE ONLY IF ONE
C --- NEGLECTS MASS-PHASE SPACE EFFECTS
        BRA1=0.5
        BRK0=0.5
        BRK0B=0.5
        BRKS=0.6667
      ENDIF
C     =====
      END
*CMZ :  1.01/14 14/05/95  11.26.27  BY  PIERO ZUCCHELLI
*CMZ :  1.01/08 05/03/95  11.39.26  BY  PIERO ZUCCHELLI
*CMZ :  1.00/00 10/08/94  16.29.40  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION THB(ITDKRC,QP,XN,XA,AK0,HV)
C
C **********************************************************************
C     BORN +VIRTUAL+SOFT PHOTON MATRIX ELEMENT**2  O(ALPHA)            *
C     PARAMETERS:                                                      *
C     HV- POLARIMETRIC FOUR-VECTOR OF TAU                              *
C     QP,XN,XA - FOUR-MOMENTA OF ELECTRON (MUON), NU AND NUBAR IN GEV  *
C     ALL FOUR-VECTORS IN TAU REST FRAME                               *
C     AK0 - INFRARED CUTOFF, MINIMAL ENERGY OF HARD PHOTONS            *
C     THB - VALUE FOR S=0                                              *
C     SEE EQS. (2.2),(2.4)-(2.5) FROM CJK (NUCL.PHYS.B351(1991)70      *
C     AND (C.2) FROM JK (NUCL.PHYS.B320(1991)20 )                      *
C **********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     +                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     +                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
      DIMENSION QP(4),XN(4),XA(4)
      REAL*8 HV(4)
      DIMENSION R(4)
      REAL*8 RXA(3),RXN(3),RQP(3)
      REAL*8 BORNPL(3),AM3POL(3),XM3POL(3)
      DATA PI /3.141592653589793238462643D0/
C
      TMASS=AMTAU
      GF=GFERMI
      ALPHAI=ALFINV
C
      TMASS2=TMASS**2
      R(4)=TMASS
      DO 10 I=1,3
        R(1)=0.D0
        R(2)=0.D0
        R(3)=0.D0
        R(I)=TMASS
        RXA(I)=R(4)*XA(4)-R(1)*XA(1)-R(2)*XA(2)-R(3)*XA(3)
        RXN(I)=R(4)*XN(4)-R(1)*XN(1)-R(2)*XN(2)-R(3)*XN(3)
C       RXK(I)=R(4)*XK(4)-R(1)*XK(1)-R(2)*XK(2)-R(3)*XK(3)
        RQP(I)=R(4)*QP(4)-R(1)*QP(1)-R(2)*QP(2)-R(3)*QP(3)
   10 CONTINUE
C     QUASI TWO-BODY VARIABLES
      U0=QP(4)/TMASS
      U3=SQRT(QP(1)**2+QP(2)**2+QP(3)**2)/TMASS
      W3=U3
      W0=(XN(4)+XA(4))/TMASS
      UP=U0+U3
      UM=U0-U3
      WP=W0+W3
      WM=W0-W3
      YU=LOG(UP/UM)/2
      YW=LOG(WP/WM)/2
      EPS2=U0**2-U3**2
      EPS=SQRT(EPS2)
      Y=W0**2-W3**2
      AL=AK0/TMASS
C     FORMFACTORS
      F0=2*U0/U3*(  DILOG(1-(UM*WM/(UP*WP)))- DILOG(1-WM/WP) +
     +DILOG(1-UM/UP) -2*YU+ 2*LOG(UP)*(YW+YU) ) +
     +1/Y* ( 2*U3*YU + (1-EPS2- 2*Y)*LOG(EPS) ) +
     + 2 - 4*(U0/U3*YU -1)* LOG(2*AL)
      FP= YU/(2*U3)*(1 + (1-EPS2)/Y ) + LOG(EPS)/Y
      FM= YU/(2*U3)*(1 - (1-EPS2)/Y ) - LOG(EPS)/Y
      F3= EPS2*(FP+FM)/2
C     SCALAR PRODUCTS OF FOUR-MOMENTA
      QPXN=QP(4)*XN(4)-QP(1)*XN(1)-QP(2)*XN(2)-QP(3)*XN(3)
      QPXA=QP(4)*XA(4)-QP(1)*XA(1)-QP(2)*XA(2)-QP(3)*XA(3)
      XNXA=XN(4)*XA(4)-XN(1)*XA(1)-XN(2)*XA(2)-XN(3)*XA(3)
      TXN=TMASS*XN(4)
      TXA=TMASS*XA(4)
      TQP=TMASS*QP(4)
C     DECAY DIFFERENTIAL WIDTH WITHOUT AND WITH POLARIZATION
      CONST3=1/(2*ALPHAI*PI)*64*GF**2
      IF (ITDKRC.EQ.0) CONST3=0D0
      XM3= -( F0* QPXN*TXA +  FP*EPS2* TXN*TXA +
     +FM* QPXN*QPXA + F3* TMASS2*XNXA )
      AM3=XM3*CONST3
C V-A  AND  V+A COUPLINGS, BUT IN THE BORN PART ONLY
      BRAK= (GV+GA)**2*TQP*XNXA+(GV-GA)**2*TXA*QPXN
     +     -(GV**2-GA**2)*TMASS*AMNUTA*QPXA
      BORN= 32*(GFERMI**2/2.)*BRAK
      DO 20 I=1,3
        XM3POL(I)= -( F0* QPXN*RXA(I) +  FP*EPS2* TXN*RXA(I) +
     +  FM* QPXN* (QPXA + (RXA(I)*TQP-TXA*RQP(I))/TMASS2 ) +
     +  F3* (TMASS2*XNXA +TXN*RXA(I) -RXN(I)*TXA)  )
        AM3POL(I)=XM3POL(I)*CONST3
C V-A  AND  V+A COUPLINGS, BUT IN THE BORN PART ONLY
        BORNPL(I)=BORN+(
     +            (GV+GA)**2*TMASS*XNXA*QP(I)
     +           -(GV-GA)**2*TMASS*QPXN*XA(I)
     +           +(GV**2-GA**2)*AMNUTA*TXA*QP(I)
     +           -(GV**2-GA**2)*AMNUTA*TQP*XA(I) )*
     +                                             32*(GFERMI**2/2.)
   20 HV(I)=(BORNPL(I)+AM3POL(I))/(BORN+AM3)-1.D0
      THB=BORN+AM3
      IF (THB/BORN.LT.0.1D0) THEN
        PRINT *, 'ERROR IN THB, THB/BORN=',THB/BORN
        STOP
      ENDIF
      RETURN
      END
*CMZ :  1.00/00 09/08/94  17.43.59  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      SUBROUTINE TRALO4(KTO,P,Q,AM)
C     **************************
C SUBSITUTE OF TRALO4
      REAL  P(4),Q(4)
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON /PTAU/ PTAU
      AM=AMAS4(P)
      ETAU=SQRT(PTAU**2+AMTAU**2)
      EXE=(ETAU+PTAU)/AMTAU
      IF(KTO.EQ.2) EXE=(ETAU-PTAU)/AMTAU
      CALL BOSTR3(EXE,P,Q)
C ======================================================================
C         END OF THE TEST JOB
C ======================================================================
      END
*CMZ :  1.00/00 10/08/94  16.29.32  BY  PIERO ZUCCHELLI
*-- AUTHOR :
      FUNCTION WIGFOR(S,XM,XGAM)
      COMPLEX WIGFOR,WIGNOR
      WIGNOR=CMPLX(-XM**2,XM*XGAM)
      WIGFOR=WIGNOR/CMPLX(S-XM**2,XM*XGAM)
      END
