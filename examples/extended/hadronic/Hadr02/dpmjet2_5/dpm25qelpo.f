      SUBROUTINE QEL_POL(ENU, LTYP,P21,P22,P23,P24,P25)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/req/k2
C     COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/POL/POLARX(4),PMODUL
      PARAMETER (NMXHKK= 89998)
c     PARAMETER (NMXHKK=25000)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
C
C     REAL*4 P,V
      COMMON /CBAD/  LBAD, NBAD
      COMMON /CNUC/ XMN,XMN2,PFERMI,EFERMI,EBIND,EB2,C0
      COMMON /CCONLUN/ LUNA
      DATA INIPRI/0/
C     WRITE(*,*)'Input Enu, Nu-typ, N-ev '
C      READ(*,*) ENU,LTYP,NUM
      CALL MASS_INI
C     num=1000
C     ltyp=5        !select neutrino tau (anutau=6, numu=3, anumu = 4)
C     DO j=5,30     !loop for different neutrino energies
C       enu=j
        enup=enu
C       DO K2=1,NUM
C         CALL GEN_QEL(ENU,LTYP)
          CALL GEN_QEL(ENU, LTYP,P21,P22,P23,P24,P25)
	  IF(INIPRI.LT.10)THEN
	    INIPRI=INIPRI+1
            WRITE(6,23) enup,(POLARX(J1),J1=1,3)
            WRITE(6,*) enup,p(4,4),p(5,4)
	    WRITE(6,*)P21,P22,P23,P24,P25
	  ENDIF
C         WHKK(1,4)=POLARX(1)
C         WHKK(2,4)=POLARX(2)
C         WHKK(3,4)=POLARX(3)
C         WHKK(4,4)=POLARX(4)
   23     FORMAT(4(1x,f10.6))
C       END DO
C     END DO
      RETURN
      END


C==================================================================
C.  Generation of  a Quasi-Elastic neutrino scattering
C==================================================================

C     SUBROUTINE GEN_QEL (ENU, LTYP)
      SUBROUTINE GEN_QEL (ENU, LTYP,P21,P22,P23,P24,P25)
C...Generate a quasi-elastic   neutrino/antineutrino
C.  Interaction on a nuclear target
C.  INPUT  : LTYP = neutrino type (1,...,6)
C.           ENU (GeV) = neutrino energy
C----------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/req/k2
C     COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      DIMENSION ROT(3,3),PI(3),PO(3)
C     REAL*4 P,V
      COMMON /CBAD/  LBAD, NBAD
      COMMON /CNUC/ XMN,XMN2,PFERMI,EFERMI,EBIND,EB2,C0
      COMMON /CCONLUN/ LUNA
CJR+
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
      COMMON /TAUTAU/Q(4,5),
     +        ETL,PXL,PYL,PZL,
     +        ETB,PXB,PYB,PZB,ETN,PXN,PYN,PZN
      COMMON /NEUTYY/NEUTYP,NEUDEC
      COMMON /NUCROS/DSIGSU,DSIGMC,NDSIG
      COMMON /QHELP/EEENU,KKTYP
      DATA ININU/0/
CJR-
C     REAL*8 DBETA(3)
C     REAL*8 MN(2), ML0(6), ML, ML2, MI, MI2, MF, MF2
      DIMENSION DBETA(3),DBETB(3),AMN(2),AML0(6)
      DATA AMN  /0.93827231D0, 0.93956563D0/
      DATA AML0 /2*0.51100D-03,2*0.105659D0, 2*1.777D0/
C     DATA PFERMI/0.22D0/
      DATA INIPRI/0/
CGB+...Binding Energy
      DATA EBIND/0.008D0/
CGB-...
      ININU=ININU+1
      IF(ININU.EQ.1)NDSIG=0
      LBAD = 0
      enu0=enu
c      write(*,*) enu0
C...Lepton mass
      AML = AML0(LTYP)       !  massa leptoni
      AML2 = AML**2          !  massa leptoni **2
C...Particle labels (LUND)
      N = 5
      K(1,1) = 21
      K(2,1) = 21
      K(3,1) = 21
      K(3,3) = 1
      K(4,1) = 1
      K(4,3) = 1
      K(5,1) = 1
      K(5,3) = 2
      K0 = (LTYP-1)/2          !  2
      K1 = LTYP/2              !  2
      KA = 12 + 2*K0           !  16
      IS = -1 + 2*LTYP - 4*K1  !  -1 +10 -8 = 1
      K(1,2) = IS*KA
      K(4,2) = IS*(KA-1)
      K(3,2) = IS*24
      LNU = 2 - LTYP + 2*K1    !  2 - 5 + 2 = - 1
      IF (LNU .EQ. 2)  THEN
        K(2,2) = 2212
        K(5,2) = 2112
        AMI = AMN(1)
        AMF = AMN(2)
CJR+
	PFERMI=TAMFEN
CJR-
      ELSE
        K(2,2) = 2112
        K(5,2) = 2212
        AMI = AMN(2)
        AMF = AMN(1)
CJR+
	PFERMI=TAMFEP
CJR-
      ENDIF
      AMI2 = AMI**2
      AMF2 = AMF**2

      DO IGB=1,5
        P(3,IGB) = 0.
        P(4,IGB) = 0.
        P(5,IGB) = 0.
      END DO

      NTRY = 0
CGB+...
      EFMAX  = SQRT(PFERMI**2 + AMI2) -AMI             ! max. Fermi Energy
      ENWELL = EFMAX + EBIND ! depth of nuclear potential well
CGB-...

  100 CONTINUE

C...4-momentum initial lepton
      P(1,5) = 0.     ! massa
      P(1,4) = ENU0    ! energia
      P(1,1) = 0.     ! px
      P(1,2) = 0.     ! py
      P(1,3) = ENU0    ! pz

C     PF = PFERMI*PYR(0)**(1./3.)
c       write(23,*) PYR(0)
c      write(*,*) 'Pfermi=',PF
c      PF = 0.
      NTRY=NTRY+1
      IF(ntry.GT.2) WRITE(*,*) ntry,enu0,k2
      IF (NTRY .GT. 500)  THEN
        LBAD = 1
        LBAD = 1
        WRITE (   6,1001)  NBAD, ENU
        RETURN
      ENDIF
C     CT = -1. + 2.*PYR(0)
c      CT = -1.
C     ST =  SQRT(1.-CT*CT)
C     F = 2.*3.1415926*PYR(0)
c      F = 0.

C     P(2,4) = SQRT(PF*PF + MI2) - EBIND  ! energia
C     P(2,1) = PF*ST*COS(F)               ! px
C     P(2,2) = PF*ST*SIN(F)               ! py
C     P(2,3) = PF*CT                      ! pz
C     P(2,5) = SQRT(P(2,4)**2-PF*PF)      ! massa
       P(2,1) = P21
       P(2,2) = P22
       P(2,3) = P23
       P(2,4) = P24
       P(2,5) = P25
c      call lulist(1)
      beta1=-p(2,1)/p(2,4)
      beta2=-p(2,2)/p(2,4)
      beta3=-p(2,3)/p(2,4)
      N=2
C      WRITE(6,*)' before transforming into target rest frame'
C      call PYlist(1)
      CALL PYROBO (0,0,0.,0.,BETA1,BETA2,BETA3)
C      print*,' nucl. rest fram ( fermi incl.) prima della rotazione'
C      call PYlist(1)
      N=5

      phi11=atan(p(1,2)/p(1,3))
      pi(1)=p(1,1)
      pi(2)=p(1,2)
      pi(3)=p(1,3)

      CALL TESTROT1(PI,Po,PHI11)
      DO ll=1,3
        IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
      END DO
c        WRITE(*,*) po
      p(1,1)=po(1)
      p(1,2)=po(2)
      p(1,3)=po(3)
      phi12=atan(p(1,1)/p(1,3))

      pi(1)=p(1,1)
      pi(2)=p(1,2)
      pi(3)=p(1,3)
      CALL TESTROT2(Pi,Po,PHI12)
      DO ll=1,3
        IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
      END DO
c        WRITE(*,*) po
      p(1,1)=po(1)
      p(1,2)=po(2)
      p(1,3)=po(3)


C      call PYlist(1)


      enu=p(1,4)

C...Kinematical limits in Q**2
c      S = P(2,5)**2 + 2.*ENU*(P(2,4)-P(2,3)) !            ????
      S = P(2,5)**2 + 2.*ENU*P(2,5)
      SQS = SQRT(S)                          ! E centro massa
      IF (SQS .LT. (AML + AMF + 3.E-03)) GOTO 100
      ELF = (S-AMF2+AML2)/(2.*SQS)           ! energia leptone finale p
      PSTAR = (S-P(2,5)**2)/(2.*SQS)       ! p* neutrino nel c.m.
      PLF = SQRT(ELF**2-AML2)               ! 3-momento leptone finale
      Q2MIN = -AML2 + 2.*PSTAR*(ELF-PLF)    ! + o -
      Q2MAX = -AML2 + 2.*PSTAR*(ELF+PLF)    ! according con cos(theta)
      IF (Q2MIN .LT. 0.)   Q2MIN = 0.      ! ??? non fisico

C...Generate Q**2
      DSIGMAX = DSQEL_Q2 (LTYP,ENU, Q2MIN)
  200 Q2 = Q2MIN + (Q2MAX-Q2MIN)*PYR(0)
      DSIG = DSQEL_Q2 (LTYP,ENU, Q2)
      IF (DSIG .LT.  DSIGMAX*PYR(0)) GOTO 200
      CALL QGAUS(Q2MIN,Q2MAX,DSIGEV,ENU,LTYP)
      NDSIG=NDSIG+1
C     WRITE(6,*)' Q2,Q2min,Q2MAX,DSIGEV',
C    &Q2,Q2min,Q2MAX,DSIGEV 


C...c.m. frame. Neutrino along z axis
      DETOT = (P(1,4)) + (P(2,4)) ! e totale
      DBETA(1) = ((P(1,1)) + (P(2,1)))/DETOT ! px1+px2/etot = beta_x
      DBETA(2) = ((P(1,2)) + (P(2,2)))/DETOT !
      DBETA(3) = ((P(1,3)) + (P(2,3)))/DETOT !
c      WRITE(*,*)
c      WRITE(*,*)
C      WRITE(*,*) 'Input values laboratory frame'
C      CALL PYLIST(1)
      N=2
      CALL PYROBO (0,0,0.,0.,-DBETA(1),-DBETA(2),-DBETA(3))
      N=5
c      STHETA = ULANGL(P(1,3),P(1,1))
c      write(*,*) 'stheta' ,stheta
c      stheta=0.
c      CALL PYROBO (0,0,-STHETA,0.,0.D0,0.D0,0.D0)
c      WRITE(*,*)
c      WRITE(*,*)
C      WRITE(*,*) 'Output values cm frame'
C      CALL PYLIST(1)
C...Kinematic in c.m. frame
      CTSTAR = ELF/PLF - (Q2 + AML2)/(2.*PSTAR*PLF) ! cos(theta) cm
      STSTAR = SQRT(1.-CTSTAR**2)
      PHI = 6.28319*PYR(0) ! random phi tra 0 e 2*pi
      P(4,5) = AML                  ! massa leptone
      P(4,4) = ELF                 ! e leptone
      P(4,3) = PLF*CTSTAR          ! px
      P(4,1) = PLF*STSTAR*COS(PHI) ! py
      P(4,2) = PLF*STSTAR*SIN(PHI) ! pz


      P(5,5) = AMF                  ! barione
      P(5,4) = (S+AMF2-AML2)/(2.*SQS)! e barione
      P(5,3) = -P(4,3)             ! px
      P(5,1) = -P(4,1)             ! py
      P(5,2) = -P(4,2)             ! pz


      P(3,5) = -Q2
      P(3,1) = P(1,1)-P(4,1)
      P(3,2) = P(1,2)-P(4,2)
      P(3,3) = P(1,3)-P(4,3)
      P(3,4) = P(1,4)-P(4,4)

C...Transform back to laboratory  frame
c      WRITE(*,*)
c      WRITE(*,*)
C      WRITE(*,*) 'before going back to nucl rest frame'
C      CALL PYLIST(1)
c      CALL PYROBO (0,0,STHETA,0.,0.D0,0.D0,0.D0)
      N=5
      CALL PYROBO (0,0,0.,0.,DBETA(1),DBETA(2),DBETA(3))
c      WRITE(*,*)
c      WRITE(*,*)
C      WRITE(*,*) 'Now back in nucl rest frame'
C      CALL PYLIST(1)
      IF(LTYP.GE.3) CALL PREPOLA(Q2,LTYP,ENU)

c********************************************

      DO kw=1,5
        pi(1)=p(kw,1)
        pi(2)=p(kw,2)
        pi(3)=p(kw,3)
        CALL TESTROT3(Pi,Po,PHI12)
        DO ll=1,3
          IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
        END DO
        p(kw,1)=po(1)
        p(kw,2)=po(2)
        p(kw,3)=po(3)
      END DO
c********************************************

C       call PYlist(1)




c********************************************

      DO kw=1,5
        pi(1)=p(kw,1)
        pi(2)=p(kw,2)
        pi(3)=p(kw,3)
        CALL TESTROT4(Pi,Po,PHI11)
        DO ll=1,3
          IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
        END DO
        p(kw,1)=po(1)
        p(kw,2)=po(2)
        p(kw,3)=po(3)
      END DO

c********************************************

C      WRITE(*,*) 'Now back in lab frame'
C       call PYlist(1)
      CALL PYROBO (1,5,0.,0.,-BETA1,-BETA2,-BETA3)
CGB+...
C...test (on final momentum of nucleon) if Fermi-blocking  
C...is operating
      ENUCL = SQRT(P(5,1)**2 + P(5,2)**2 + P(5,3)**2 + P(5,5)**2) 
     &  - P(5,5)
      IF (ENUCL.LT. EFMAX) THEN
	IF(INIPRI.LT.10)THEN
	  INIPRI=INIPRI+1
          WRITE(6,*)' qel: Pauli ENUCL.LT.EFMAX ', ENUCL,EFMAX
C...the interaction is not possible due to Pauli-Blocking and 
C...it must be resampled
	ENDIF
        GOTO 100
      ELSE IF (ENUCL.LT.ENWELL.and.ENUCL.GE.EFMAX) THEN
        IF(INIPRI.LT.10)THEN
          INIPRI=INIPRI+1
      WRITE(6,*)' qel: inside ENUCL.LT.ENWELL ', ENUCL,ENWELL
        ENDIF
C                      Reject (J:R) here all these events
C                      are otherwise rejected in dpmjet
        GOTO 100
C...the interaction is possible, but the nucleon remains inside
C...the nucleus. The nucleus is therefore left excited.
C...We treat this case as a nucleon with 0 kinetic energy.
C       P(5,5) = AMF
C       P(5,4) = AMF
C       P(5,1) = 0.
C       P(5,2) = 0.
C       P(5,3) = 0.
      ELSE IF (ENUCL.GE.ENWELL) THEN
C     WRITE(6,*)' qel ENUCL.GE.ENWELL ',ENUCL,ENWELL
C...the interaction is possible, the nucleon can exit the nucleus
C...but the nuclear well depth must be subtracted. The nucleus could be
C...left in an excited state.
        Pstart = SQRT(P(5,1)**2 + P(5,2)**2 + P(5,3)**2)
C       P(5,4) = ENUCL-ENWELL + AMF
        Pnucl = SQRT(P(5,4)**2-AMF**2)
C...The 3-momentum is scaled assuming that the direction remains
C...unaffected
        P(5,1) = P(5,1) * Pnucl/Pstart
        P(5,2) = P(5,2) * Pnucl/Pstart
        P(5,3) = P(5,3) * Pnucl/Pstart
C     WRITE(6,*)' qel new P(5,4) ',P(5,4)
      ENDIF
CGB-...
      DSIGSU=DSIGSU+DSIGEV

C        gbsumx = pxl + pxb + pxn
C        gbsumy = pyl + pyb + pyn
C        gbsumz = pzl + pzb + pzn
C        gbsume = etl + etb + etn
C        write(*,*) 'GB output to LEPDCY 1',ETL,PXL,PYL,PZL
C        write(*,*) 'GB output to LEPDCY 2',ETb,PXb,PYb,PZb
C        write(*,*) 'GB output to LEPDCY 3',ETn,PXn,PYn,PZn
C        write(*,*) 'GB output to LEPDCY e',gbsumx,gbsumy,gbsumz,gbsume
C WRITE(6,*)' P(4,1-5) ',P(4,1),P(4,2),P(4,3),P(4,4),P(4,5)
C WRITE(6,*)' Q(4,1-5) ',Q(4,1),Q(4,2),Q(4,3),Q(4,4),Q(4,5)
	 GA=P(4,4)/P(4,5)
	 BGX=P(4,1)/P(4,5)
	 BGY=P(4,2)/P(4,5)
	 BGZ=P(4,3)/P(4,5)
C WRITE(6,*)' g,bg1-3 ',GA,BGX,BGY,BGZ
C CALL DALTRA(GA,BGX,BGY,BGZ,PXL,PYL,PZL,ETL,
C    &   PPL,PPXL,PPYL,PPZL,EETL)
C CALL DALTRA(GA,BGX,BGY,BGZ,PXB,PYB,PZB,ETB,
C    &   PPB,PPXB,PPYB,PPZB,EETB)
C CALL DALTRA(GA,BGX,BGY,BGZ,PXN,PYN,PZN,ETN,
C    &   PPN,PPXN,PPYN,PPZN,EETN)
C        rrsumx = ppxl + ppxb + ppxn
C        rrsumy = ppyl + ppyb + ppyn
C        rrsumz = ppzl + ppzb + ppzn
C        rrsume = eetl + eetb + eetn
C        write(*,*) 'rr output to LEPDCY 1',EETL,PPXL,PPYL,PPZL
C        write(*,*) 'rr output to LEPDCY 2',EETb,PPXb,PPYb,PPZb
C        write(*,*) 'rr output to LEPDCY 3',EETn,PPXn,PPYn,PPZn
C        write(*,*) 'rr output to LEPDCY e',rrsumx,rrsumy,rrsumz,rrsume
         DBETB(1)=BGX/GA
         DBETB(2)=BGY/GA
         DBETB(3)=BGZ/GA
	 IF(NEUDEC.EQ.1.OR.NEUDEC.EQ.2) THEN
            CALL PYROBO (6,8,0.,0.,DBETB(1),DBETB(2),DBETB(3))
	 ENDIF
c
C      PRINT*,' FINE   EVENTO '
C      call PYlist(1)
      enu=enu0
      RETURN

 1001 FORMAT(2X, 'GEN_QEL   : event rejected ', I5,  G10.3)
      END


C====================================================================
C.  Masses
C====================================================================

      SUBROUTINE MASS_INI
C...Initialize  the kinematics for the quasi-elastic cross section
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CMAQEL/ EML(6), EMLSQ(6), EMN1(6), EMN2(6), ETQE(6)
     &  ,EMN1SQ(6),EMN2SQ(6),EMPROT,EMNEUT,EMN,EMPROTSQ,EMNEUTSQ,EMNSQ
      EML(1) = 0.51100D-03   ! e-
      EML(2) = EML(1)        ! e+
      EML(3) = 0.105659D0      ! mu-
      EML(4) = EML(3)        ! mu+
      EML(5) = 1.7777D0        ! tau-
      EML(6) = EML(5)        ! tau+
      EMPROT = 0.93827231D0    ! p
      EMNEUT = 0.93956563D0    ! n
      EMPROTSQ = EMPROT**2
      EMNEUTSQ = EMNEUT**2
      EMN = (EMPROT + EMNEUT)/2.
      EMNSQ = EMN**2
      DO J=1,3
        J0 = 2*(J-1)
        EMN1(J0+1) = EMNEUT
        EMN1(J0+2) = EMPROT
        EMN2(J0+1) = EMPROT
        EMN2(J0+2) = EMNEUT
      ENDDO
      DO J=1,6
        EMLSQ(J) = EML(J)**2
        ETQE(J)  = ((EMN2(J)+ EML(J))**2-EMN1(J)**2)/(2.*EMN1(J))
      ENDDO
      RETURN
      END

      FUNCTION DSQEL_Q2 (JTYP,ENU, Q2)
C...differential cross section for  Quasi-Elastic scattering
C.       nu + N -> l + N'
C.  From Llewellin Smith  Phys.Rep.  3C, 261, (1971).
C.
C.  INPUT :  JTYP = 1,...,6    nu_e, ...., nubar_tau
C.           ENU (GeV) =  Neutrino energy
C.           Q2  (GeV**2) =  (Transfer momentum)**2
C.
C.  OUTPUT : DSQEL_Q2  = differential  cross section :
C.                       dsigma/dq**2  (10**-38 cm+2/GeV**2)
C------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CAXIAL/ FA0, AXIAL2
      COMMON /CMAQEL/ EML(6), EMLSQ(6), EMN1(6), EMN2(6), ETQE(6)
     &  ,EMN1SQ(6),EMN2SQ(6),EMPROT,EMNEUT,EMN,EMPROTSQ,EMNEUTSQ,EMNSQ
      DIMENSION SS(6)
      DATA C0 /0.17590D0 /  ! G_F**2 cos(theta_c)**2 M**2 /(8 pi) 10**-38 cm+2
      DATA SS /1.D0, -1.D0, 1.D0, -1.D0, 1.D0, -1.D0/
      DATA AXIAL2 /1.03D0/  ! to be checked
      FA0=-1.253D0
      CSI = 3.71D0                   !  ???
      GVE = 1.D0/ (1.D0 + Q2/0.84D0**2)**2   ! G_e(q**2)
      GVM = (1.D0+CSI)*GVE           ! G_m (q**2)
      X = Q2/(EMN*EMN)     ! emn=massa barione
      XA = X/4.D0
      FV1 = 1.D0/(1.D0+XA)*(GVE+XA*GVM)
      FV2 = 1.D0/(1.D0+XA)*(GVM-GVE)
      FA = FA0/(1.D0 + Q2/AXIAL2)**2
      FFA = FA*FA
      FFV1 = FV1*FV1
      FFV2 = FV2*FV2
      RM = EMLSQ(JTYP)/(EMN*EMN)            ! emlsq(jtyp)
      A1 = (4.D0+X)*FFA - (4.D0-X)*FFV1 + X*FFV2*(1.D0-XA)+4*X*FV1*FV2
      A2 = -RM * ((FV1 + FV2)**2 +  FFA)
      AA = (XA+0.25D0*RM)*(A1 + A2)
      BB = -X*FA*(FV1 + FV2)
      CC = 0.25D0*(FFA + FFV1 + XA*FFV2)
      SU = (4.D0*ENU*EMN - Q2 - EMLSQ(JTYP))/(EMN*EMN)
      DSQEL_Q2 = C0*(AA + SS(JTYP)*BB*SU + CC*SU*SU) / (ENU*ENU)  !
      IF(DSQEL_Q2 .LT. 0.D0) DSQEL_Q2 = 0.D0
      RETURN
      END

      SUBROUTINE PREPOLA(Q2,JTYP,ENU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c By G. Battistoni and E. Scapparone (sept. 1997)
c According to:
c     Albright & Jarlskog, Nucl Phys B84 (1975) 467
c
c
C    COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C     REAL*4 P,V
      COMMON /CAXIAL/ FA0, AXIAL2
      COMMON /CMAQEL/ EML(6), EMLSQ(6), EMN1(6), EMN2(6), ETQE(6)
     &  ,EMN1SQ(6),EMN2SQ(6),EMPROT,EMNEUT,EMN,EMPROTSQ,EMNEUTSQ,EMNSQ
      REAL*8 MUMASS,POL(4,4),BB2(3),NU
      COMMON/POL/POLARX(4),PMODUL
      COMMON /TAUTAU/Q(4,5),
     +        ETL,PXL,PYL,PZL,
     +        ETB,PXB,PYB,PZB,ETN,PXN,PYN,PZN
      COMMON /NEUTYY/NEUTYP,NEUDEC
      DIMENSION SS(6)
      DATA C0 /0.17590D0 /  ! G_F**2 cos(theta_c)**2 M**2 /(8 pi) 10**-38 cm+2
      DATA SS /1.D0, -1.D0, 1.D0, -1.D0, 1.D0, -1.D0/
C     DATA AXIAL2 /1.03D0/  ! to be checked
      RML=P(4,5)
      RMM=0.93960D+00
      FM2 = RMM**2
      MPI = 0.135D+00
      OLDQ2=Q2
      FA0=-1.253D+00
      CSI = 3.71D+00                      !
      GVE = 1.D0/ (1.D0 + Q2/(0.84D+00)**2)**2   ! G_e(q**2)
      GVM = (1.D0+CSI)*GVE           ! G_m (q**2)
      X = Q2/(EMN*EMN)     ! emn=massa barione
      XA = X/4.D0
      FV1 = 1.D0/(1.D0+XA)*(GVE+XA*GVM)
      FV2 = 1.D0/(1.D0+XA)*(GVM-GVE)
      FA = FA0/(1.D0 + Q2/AXIAL2**2)**2
      FFA = FA*FA
      FFV1 = FV1*FV1
      FFV2 = FV2*FV2
      FP=2.D0*FA*RMM/(MPI**2 + Q2)
      RM = EMLSQ(JTYP)/(EMN*EMN)            ! emlsq(jtyp)
      A1 = (4.D0+X)*FFA-(4.D0-X)*FFV1+X*FFV2*(1.D0-XA)+4.D0*X*FV1*FV2
      A2 = -RM * ((FV1 + FV2)**2 +  FFA)
      AA = (XA+0.25D+00*RM)*(A1 + A2)
      BB = -X*FA*(FV1 + FV2)
      CC = 0.25D+00*(FFA + FFV1 + XA*FFV2)
      SU = (4.D+00*ENU*EMN - Q2 - EMLSQ(JTYP))/(EMN*EMN)

      OMEGA1=FFA+XA*(FFA+(FV1+FV2)**2   )  ! articolo di ll...-smith
      OMEGA2=4.D+00*CC
      OMEGA3=2.D+00*FA*(FV1+FV2)
      OMEGA4P=(-(FV1+FV2)**2-(FA+2*FP)**2+(4.0D+00+
     1     (Q2/FM2))*FP**2)
      OMEGA5=OMEGA2
      OMEGA4=(OMEGA4P-OMEGA2+2*OMEGA5)/4.D+00
      WW1=2.D+00*OMEGA1*EMN**2
      WW2=2.D+00*OMEGA2*EMN**2
      WW3=2.D+00*OMEGA3*EMN**2
      WW4=2.D+00*OMEGA4*EMN**2
      WW5=2.D+00*OMEGA5*EMN**2

      DO I=1,3
        BB2(I)=-P(4,I)/P(4,4)
      END DO
c      WRITE(*,*)
c      WRITE(*,*)
c      WRITE(*,*) 'Prepola: ready to transform to lepton rest frame'
c      CALL LULIST(1)
      N=5
      CALL PYROBO (0,0,0.,0.,BB2(1),BB2(2),BB2(3) )
* NOW PARTICLES ARE IN THE SCATTERED LEPTON  REST FRAME
c      WRITE(*,*)
c      WRITE(*,*)
c      WRITE(*,*) 'Prepola: now in lepton rest frame'
c      CALL LULIST(1)
      EE=ENU
      QM2=Q2+RML**2
      U=Q2/(2.*RMM)
      FRAC=QM2*WW1 + (2.D+00*EE*(EE-U) - 0.5D+00*QM2)*WW2 - SS(JTYP)*
     +     (0.5D+00/(RMM**2))*(2.D+00*RMM*EE*Q2 - U*QM2)*WW3 +
     +     ((RML**2)/(2.D+00*FM2))*(QM2*WW4-2.D+00*RMM*EE*WW5) !<=FM2 inv di RMM!!

      FACTK=2.D+00*WW1 -WW2 -SS(JTYP)*(EE/RMM)*WW3 +((EE-U)/RMM)*WW5
     +     - ((RML**2)/FM2)*WW4                        !<=FM2 inv di RMM!!

      FACTP=2.D+00*EE/RMM*WW2 - (QM2/(2.D+00*RMM**2))*(SS(JTYP)*WW3+WW5)

      DO I=1,3
        POL(4,I)=RML*SS(JTYP)*(FACTK*P(1,I)+FACTP*P(2,I))/FRAC
        POLARX(I)=POL(4,I)
      END DO


      PMODUL=0.D0
      DO I=1,3
        PMODUL=PMODUL+POL(4,I)**2
      END DO

      IF(JTYP.GT.4.AND.NEUDEC.GT.0) THEN
         IF(NEUDEC.EQ.1) THEN
            CALL LEPDCYP(EML(JTYP),EML(JTYP-2),POLARX(3),
     +        ETL,PXL,PYL,PZL,
     +        ETB,PXB,PYB,PZB,ETN,PXN,PYN,PZN)
c
c     Tau has decayed in muon
c
         ENDIF
         IF(NEUDEC.EQ.2) THEN
            CALL LEPDCYP(EML(JTYP),EML(JTYP-4),POLARX(3),
     +        ETL,PXL,PYL,PZL,
     +        ETB,PXB,PYB,PZB,ETN,PXN,PYN,PZN)
c
c     Tau has decayed in electron
c
         ENDIF
         K(4,1)=15
         K(4,4) = 6
         K(4,5) = 8
         N=N+3
c
c     fill common for muon(electron) 
c
         P(6,1)=PXL
         P(6,2)=PYL
         P(6,3)=PZL
         P(6,4)=ETL
         K(6,1)=1
         IF(JTYP.EQ.5) THEN
            IF(NEUDEC.EQ.1) THEN
               P(6,5)=EML(JTYP-2)
               K(6,2)=13
            ELSEIF(NEUDEC.EQ.2) THEN
               P(6,5)=EML(JTYP-4)
               K(6,2)=11
            ENDIF
         ELSEIF(JTYP.EQ.6) THEN
            IF(NEUDEC.EQ.1) THEN
               K(6,2)=-13
            ELSEIF(NEUDEC.EQ.2) THEN
               K(6,2)=-11
            ENDIF
         END IF
         K(6,3)=4
         K(6,4)=0
         K(6,5)=0
c
c     fill common for tau_(anti)neutrino 
c
         P(7,1)=PXB
         P(7,2)=PYB
         P(7,3)=PZB
         P(7,4)=ETB
         P(7,5)=0.
         K(7,1)=1
         IF(JTYP.EQ.5) THEN
            K(7,2)=16
         ELSEIF(JTYP.EQ.6) THEN
            K(7,2)=-16
         END IF
         K(7,3)=4
         K(7,4)=0
         K(7,5)=0
c
c     Fill common for muon(electron)_(anti)neutrino
c
         P(8,1)=PXN
         P(8,2)=PYN
         P(8,3)=PZN
         P(8,4)=ETN
         P(8,5)=0.
         K(8,1)=1
         IF(JTYP.EQ.5) THEN
            IF(NEUDEC.EQ.1) THEN
               K(8,2)=-14
            ELSEIF(NEUDEC.EQ.2) THEN
               K(8,2)=-12
            ENDIF
         ELSEIF(JTYP.EQ.6) THEN
            IF(NEUDEC.EQ.1) THEN
               K(8,2)=14
            ELSEIF(NEUDEC.EQ.2) THEN
               K(8,2)=12
            ENDIF
         END IF
         K(8,3)=4
         K(8,4)=0
         K(8,5)=0
      ENDIF
c      WRITE(*,*)
c      WRITE(*,*)

c      IF(PMODUL.GE.1.D+00) THEN
c        WRITE(*,*) 'Pol',(POLARX(I),I=1,3)
c        write(*,*) pmodul
c        DO I=1,3
c          POL(4,I)=POL(4,I)/PMODUL
c          POLARX(I)=POL(4,I)
c        END DO
c        PMODUL=0.
c        DO I=1,3
c          PMODUL=PMODUL+POL(4,I)**2
c        END DO
c        WRITE(*,*) 'Pol',(POLARX(I),I=1,3)
c
c      ENDIF

c      WRITE(*,*) 'PMODUL = ',PMODUL

c      WRITE(*,*)
c      WRITE(*,*)
c      WRITE(*,*) 'prepola: Now back to nucl rest frame'
      CALL PYROBO (1,5,0.,0.,-BB2(1),-BB2(2),-BB2(3))
c      call lulist(1)
      XDC = V(4,1)+V(4,5)*P(4,1)/P(4,5)
      YDC = V(4,2)+V(4,5)*P(4,2)/P(4,5)
      ZDC = V(4,3)+V(4,5)*P(4,3)/P(4,5)
      DO NDC =6,8
         V(NDC,1) = XDC
         V(NDC,2) = YDC
         V(NDC,3) = ZDC
      END DO

      RETURN
      END

      SUBROUTINE TESTROT1(PI,PO,PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ROT(3,3),PI(3),PO(3)
      ROT(1,1)=1.D0
      ROT(1,2)=0.D0
      ROT(1,3)=0.D0
      ROT(2,1)=0.D0
      ROT(2,2)=cos(phi)
      ROT(2,3)=-sin(phi)
      ROT(3,1)=0.D0
      ROT(3,2)=SIN(phi)
      ROT(3,3)=COS(phi)
      DO 140 J=1,3
        PO(J)=ROT(J,1)*PI(1)+ROT(J,2)*PI(2)+ROT(J,3)*PI(3)
  140 CONTINUE
      RETURN
      END


      SUBROUTINE TESTROT2(PI,PO,PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ROT(3,3),PI(3),PO(3)
      ROT(1,1)=0.D0
      ROT(1,2)=1.D0
      ROT(1,3)=0.D0
      ROT(2,1)=cos(phi)
      ROT(2,2)=0.D0
      ROT(2,3)=-sin(phi)
      ROT(3,1)=sin(phi)
      ROT(3,2)=0.D0
      ROT(3,3)=COS(phi)
      DO 140 J=1,3
        PO(J)=ROT(J,1)*PI(1)+ROT(J,2)*PI(2)+ROT(J,3)*PI(3)
  140 CONTINUE
      RETURN
      END

      SUBROUTINE TESTROT3(PI,PO,PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ROT(3,3),PI(3),PO(3)
      ROT(1,1)=0.D0
      ROT(2,1)=1.D0
      ROT(3,1)=0.D0
      ROT(1,2)=cos(phi)
      ROT(2,2)=0.D0
      ROT(3,2)=-sin(phi)
      ROT(1,3)=sin(phi)
      ROT(2,3)=0.D0
      ROT(3,3)=COS(phi)
      DO 140 J=1,3
        PO(J)=ROT(J,1)*PI(1)+ROT(J,2)*PI(2)+ROT(J,3)*PI(3)
  140 CONTINUE
      RETURN
      END

      SUBROUTINE TESTROT4(PI,PO,PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ROT(3,3),PI(3),PO(3)
      ROT(1,1)=1.D0
      ROT(2,1)=0.D0
      ROT(3,1)=0.D0
      ROT(1,2)=0.D0
      ROT(2,2)=cos(phi)
      ROT(3,2)=-sin(phi)
      ROT(1,3)=0.D0
      ROT(2,3)=SIN(phi)
      ROT(3,3)=COS(phi)
      DO 140 J=1,3
        PO(J)=ROT(J,1)*PI(1)+ROT(J,2)*PI(2)+ROT(J,3)*PI(3)
  140 CONTINUE
      RETURN
      END

      SUBROUTINE LEPDCYP (AMA,AML,POL,ETL,PXL,PYL,PZL,
     &  ETB,PXB,PYB,PZB,ETN,PXN,PYN,PZN)
C
C-----------------------------------------------------------------
C 
C   Author   :- G. Battistoni         10-NOV-1995
C
C=================================================================
C
C   Purpose   : performs decay of polarized lepton in
C               its rest frame: a => b + l + anti-nu
C               (Example: mu- => nu-mu + e- + anti-nu-e)
C               Polarization is assumed along Z-axis
C               WARNING: 
C               1) b AND anti-nu ARE ASSUMED TO BE NEUTRINOS
C                  OF NEGLIGIBLE MASS
C               2) RADIATIVE CORRECTIONS ARE NOT CONSIDERED
C                  IN THIS VERSION
C
C   Method    : modifies phase space distribution obtained
C               by routine EXPLOD using a rejection against the
C               matrix element for unpolarized lepton decay
C
C   Inputs    : Mass of a :  AMA
C               Mass of l :  AML
C               Polar. of a: POL
C               (Example: fully polar. mu- decay: AMA=AMMUON, AML=AMELCT,
C                                                 POL = -1)
C
C   Outputs   : kinematic variables in the rest frame of decaying lepton
C               ETL,PXL,PYL,PZL 4-moment of l
C               ETB,PXB,PYB,PZB 4-moment of b
C               ETN,PXN,PYN,PZN 4-moment of anti-nu 
C
C============================================================
C +
C Declarations.
C -
*
*=== dblprc ==========================================================*
*                                                                     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      PARAMETER ( PIPIPI = 3.1415926535897932270 D+00 )
      PARAMETER ( ENEPER = 2.7182818284590452354 D+00 )
      PARAMETER ( SQRENT = 1.6487212707001281468 D+00 )
      PARAMETER ( CLIGHT = 2.99792458         D+10 )
      PARAMETER ( AVOGAD = 6.0221367          D+23 )
      PARAMETER ( AMELGR = 9.1093897          D-28 )
      PARAMETER ( PLCKBR = 1.05457266         D-27 )
      PARAMETER ( ELCCGS = 4.8032068          D-10 )
      PARAMETER ( ELCMKS = 1.60217733         D-19 )
      PARAMETER ( AMUGRM = 1.6605402          D-24 )
      PARAMETER ( AMMUMU = 0.113428913        D+00 )
      PARAMETER ( ALPFSC = 7.2973530791728595 D-03 )
      PARAMETER ( FSCTO2 = 5.3251361962113614 D-05 )
      PARAMETER ( FSCTO3 = 3.8859399018437826 D-07 )
      PARAMETER ( FSCTO4 = 2.8357075508200407 D-09 )
      PARAMETER ( PLABRC = 0.197327053        D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMMUON = 0.105658389        D+00 )
      PARAMETER ( RCLSEL = 2.8179409183694872 D-13 )
      PARAMETER ( GEVMEV = 1.0                D+03 )
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( ALGVMV = 6.90775527898214   D+00 )
      PARAMETER ( RADDEG = 180.D+00 / PIPIPI )
      PARAMETER ( DEGRAD = PIPIPI / 180.D+00 )
C +
C    variables for EXPLOD
C - 
      PARAMETER ( KPMX = 10 )
      DIMENSION AMEXPL (KPMX), PXEXPL (KPMX), PYEXPL (KPMX),
     &          PZEXPL (KPMX), PTEXPL (KPMX), ETEXPL (KPMX),
     &          AMHELP (KPMX)
C +
C      test variables
C -
      COMMON /GBATNU/ ELERAT,NTRY
C +
C     Initializes test variables
C - 
      NTRY = 0
      ELERAT = 0.D+00
C +
C     Maximum value for matrix element
C - 
      ELEMAX = ( AMA**2 + AML**2 )**2 / AMA**2 * ( AMA**2 - AML**2 +
     &  SQRT( AMA**4 + AML**4 - 3.D+00 * AMA**2 * AML**2 ) )
C + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
C     Inputs for EXPLOD
C part. no. 1 is l       (e- in mu- decay)
C part. no. 2 is b       (nu-mu in mu- decay)
C part. no. 3 is anti-nu (anti-nu-e in mu- decay)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      NPEXPL = 3
      ETOTEX = AMA
      AMEXPL(1) = AML
      AMEXPL(2) = 0.D+00
      AMEXPL(3) = 0.D+00
C +
C     phase space distribution
C -
  100 CONTINUE
      NTRY = NTRY + 1
      CALL EXPLOD ( NPEXPL, AMEXPL, ETOTEX, ETEXPL, PXEXPL,
     &                 PYEXPL, PZEXPL )
C +
C  Calculates matrix element: 
C  64*GF**2{[P(a)-ama*S(a)]*P(anti-nu)}{P(l)*P(b)}
C  Here CTH is the cosine of the angle between anti-nu and Z axis
C -
      CTH = PZEXPL(3) / SQRT ( PXEXPL(3)**2 + PYEXPL(3)**2 + 
     &  PZEXPL(3)**2 ) 
      PROD1 = ETEXPL(3) * AMA * (1.D+00 - POL * CTH)
      PROD2 = ETEXPL(1) * ETEXPL(2) - PXEXPL(1)*PXEXPL(2) -
     &     PYEXPL(1)*PYEXPL(2) - PZEXPL(1)*PZEXPL(2)
      ELEMAT = 16.D+00 * PROD1 * PROD2 
      IF(ELEMAT.GT.ELEMAX) THEN
        WRITE(6,*) 'Problems in LEPDCY',ELEMAX,ELEMAT
        STOP
      ENDIF
C +
C     Here performs the rejection
C - 
      TEST = RNDM(ETOTEX) * ELEMAX
      IF ( TEST .GT. ELEMAT ) GO TO 100
C +
C     final assignment of variables
C -
      ELERAT = ELEMAT/ELEMAX
      ETL = ETEXPL(1)
      PXL = PXEXPL(1)
      PYL = PYEXPL(1)
      PZL = PZEXPL(1)
      ETB = ETEXPL(2)
      PXB = PXEXPL(2)
      PYB = PYEXPL(2)
      PZB = PZEXPL(2)
      ETN = ETEXPL(3)
      PXN = PXEXPL(3)
      PYN = PYEXPL(3)
      PZN = PZEXPL(3)
  999 RETURN
      END

C==================================================================
C.  Generation of  Delta resonance events
C==================================================================

      SUBROUTINE GEN_DELTA(ENU,LLEP,LTARG,JINT,P21,P22,P23,P24,P25)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...Generate a Delta-production neutrino/antineutrino
C.  CC-interaction on a nucleon
C
C.  INPUT  ENU (GeV) = Neutrino Energy
C.         LLEP = neutrino type
C.         LTARG = nucleon target type 1=p, 2=n.
C.         JINT = 1:CC, 2::NC
C.
C.  OUTPUT PPL(4)  4-monentum of final lepton
C----------------------------------------------------
C     COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON /KINQEL/ PPL(5), PPN(5)
      COMMON /CBAD/  LBAD, NBAD
      COMMON /CCONLUN/ LUNA
      DIMENSION ROT(3,3),PI(3),PO(3)
C     REAL*4 AMD0, AMD, AMN(2), AML0(6), AML, AML2, AMDMIN
      DIMENSION AML0(6),AMN(2)
      DATA AMD0 /1.231/, GAMD /0.12/, DELD/0.169/, AMDMIN/1.084/
      DATA AMN  /0.93827231, 0.93956563/
      DATA AML0 /2*0.51100E-03,2*0.105659, 2*1.777/

c     WRITE(6,*)' GEN_DEL',ENU,LLEP,LTARG,JINT,P21,P22,P23,P24,P25
      LBAD = 0
C...Final lepton mass
      IF (JINT.EQ.1) THEN
	AML = AML0(LLEP)
      ELSE
	AML = 0.
      ENDIF
      AML2 = AML**2

C...Particle labels (LUND)
      N = 5
      K(1,1) = 21
      K(2,1) = 21
      K(3,1) = 21
      K(4,1) = 1
      K(3,3) = 1
      K(4,3) = 1
      IF (LTARG .EQ. 1)  THEN
	 K(2,2) = 2212
      ELSE
	 K(2,2) = 2112
      ENDIF
      K0 = (LLEP-1)/2
      K1 = LLEP/2
      KA = 12 + 2*K0
      IS = -1 + 2*LLEP - 4*K1
      LNU = 2 - LLEP + 2*K1
      K(1,2) = IS*KA
      K(5,1) = 1
      K(5,3) = 2
      IF (JINT .EQ. 1)  THEN                    ! CC interactions
	 K(3,2) = IS*24
	 K(4,2) = IS*(KA-1)
	IF(LNU.EQ.1) THEN
	  IF (LTARG .EQ. 1)  THEN
	      K(5,2) = 2224
	  ELSE
	      K(5,2) = 2214
	  ENDIF
	ELSE
	  IF (LTARG .EQ. 1)  THEN
	      K(5,2) = 2114
	  ELSE
	      K(5,2) = 1114
	  ENDIF
	ENDIF
      ELSE
	 K(3,2) = 23                           ! NC (Z0) interactions
	 K(4,2) = K(1,2)
	 IF (LTARG .EQ. 1)  THEN
	     K(5,2) = 2114
	 ELSE
	     K(5,2) = 2214
	 ENDIF
      ENDIF

C...4-momentum initial lepton
      P(1,5) = 0.
      P(1,4) = ENU
      P(1,1) = 0.
      P(1,2) = 0.
      P(1,3) = ENU
C...4-momentum initial nucleon
      P(2,5) = AMN(LTARG)
C     P(2,4) = P(2,5)
C     P(2,1) = 0.
C     P(2,2) = 0.
C     P(2,3) = 0.
       P(2,1) = P21
       P(2,2) = P22
       P(2,3) = P23
       P(2,4) = P24
       P(2,5) = P25
      N=2
C     call PYlist(1)
      beta1=-p(2,1)/p(2,4)
      beta2=-p(2,2)/p(2,4)
      beta3=-p(2,3)/p(2,4)
      N=2
      CALL PYROBO (0,0,0.,0.,BETA1,BETA2,BETA3)

C     print*,' nucl. rest fram ( fermi incl.) prima della rotazione'
C     call PYlist(1)

      phi11=atan(p(1,2)/p(1,3))
      pi(1)=p(1,1)
      pi(2)=p(1,2)
      pi(3)=p(1,3)

      CALL TESTROT1(PI,Po,PHI11)
      DO ll=1,3
       IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
      END DO
      p(1,1)=po(1)
      p(1,2)=po(2)
      p(1,3)=po(3)
      phi12=atan(p(1,1)/p(1,3))

      pi(1)=p(1,1)
      pi(2)=p(1,2)
      pi(3)=p(1,3)
      CALL TESTROT2(Pi,Po,PHI12)
      DO ll=1,3
        IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
      END DO
      p(1,1)=po(1)
      p(1,2)=po(2)
      p(1,3)=po(3)
C     call PYlist(1)

      ENUU=P(1,4)

C...Generate the Mass of the Delta
      NTRY = 0
100   R = PYR(0)
      AMD=AMD0+0.5*GAMD*TAN((2.*R-1.)*ATAN(2.*DELD/GAMD))
      NTRY = NTRY + 1
      IF (NTRY .GT. 1000)  THEN
	 LBAD = 1
	 WRITE (   6,1001)  NBAD, ENUU,AMD,AMDMIN,AMD0,GAMD,ET
	 WRITE (LUNA,1001)  NBAD, ENUU
	 RETURN
      ENDIF
      IF (AMD .LT. AMDMIN)  GOTO 100
      ET = ((AMD+AML)**2 - AMN(LTARG)**2)/(2.*AMN(LTARG))
      IF (ENUU .LT. ET) GOTO 100

C...Kinematical  limits in Q**2
      S = AMN(LTARG)**2 + 2.*AMN(LTARG)*ENUU
      SQS = SQRT(S)
      PSTAR = (S - AMN(LTARG)**2)/(2.*SQS)
      ELF = (S - AMD**2 + AML2)/(2.*SQS)
      PLF = SQRT(ELF**2 - AML2)
      Q2MIN = -AML2 + 2.*PSTAR*(ELF-PLF)
      Q2MAX = -AML2 + 2.*PSTAR*(ELF+PLF)
      IF (Q2MIN .LT. 0.)   Q2MIN = 0.

      DSIGMAX = DSIGMA_DELTA(LNU,-Q2MIN, S, AML, AMD)
200   Q2 = Q2MIN + (Q2MAX-Q2MIN)*PYR(0)
      DSIG = DSIGMA_DELTA(LNU,-Q2, S, AML, AMD)
      IF (DSIG .LT.  DSIGMAX*PYR(0)) GOTO 200


C...Generate the kinematics of the final particles
      EISTAR = (S + AMN(LTARG)**2)/(2.*SQS)
      GAM = EISTAR/AMN(LTARG)
      BET = PSTAR/EISTAR
      CTSTAR = ELF/PLF - (Q2 + AML2)/(2.*PSTAR*PLF)
      EL  = GAM*(ELF + BET*PLF*CTSTAR)
      PLZ = GAM*(PLF*CTSTAR + BET*ELF)
      PL  = SQRT(EL**2 - AML2)
      PLT = SQRT(MAX(1.E-06,(PL*PL - PLZ*PLZ)))
      PHI = 6.28319*PYR(0)
      P(4,1) = PLT*COS(PHI)
      P(4,2) = PLT*SIN(PHI)
      P(4,3) = PLZ
      P(4,4) = EL
      P(4,5) = AML

C...4-momentum of Delta
      P(5,1) = -P(4,1)
      P(5,2) = -P(4,2)
      P(5,3) = ENUU-P(4,3)
      P(5,4) = ENUU+AMN(LTARG)-P(4,4)
      P(5,5) = AMD

C...4-momentum  of intermediate boson
      P(3,5) = -Q2
      P(3,4) = P(1,4)-P(4,4)
      P(3,1) = P(1,1)-P(4,1)
      P(3,2) = P(1,2)-P(4,2)
      P(3,3) = P(1,3)-P(4,3)
      N=5
C      call PYlist(1)

      DO kw=1,5
        pi(1)=p(kw,1)
        pi(2)=p(kw,2)
        pi(3)=p(kw,3)
        CALL TESTROT3(Pi,Po,PHI12)
        DO ll=1,3
          IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
        END DO
        p(kw,1)=po(1)
        p(kw,2)=po(2)
        p(kw,3)=po(3)
      END DO

c********************************************
C      call PYlist(1)

c********************************************

        DO kw=1,5
          pi(1)=p(kw,1)
          pi(2)=p(kw,2)
          pi(3)=p(kw,3)
          CALL TESTROT4(Pi,Po,PHI11)
          DO ll=1,3
            IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
          END DO
          p(kw,1)=po(1)
          p(kw,2)=po(2)
          p(kw,3)=po(3)
       END DO
c********************************************
C      call PYlist(1)

C         transform back into Lab.
      CALL PYROBO (0,0,0.,0.,-BETA1,-BETA2,-BETA3)

C     WRITE(6,*)' Lab fram ( fermi incl.) '
C     call PYlist(1)
      N=5
      CALL PYEXEC
C     call PYlist(1)


      RETURN
1001  FORMAT(2X, 'GEN_DELTA : event rejected ', I5,  6G10.3)
      END

      FUNCTION DSIGMA_DELTA (LNU, QQ, S, AML, MD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...Reaction nu + N -> lepton + Delta
C.  returns the  cross section
C.  dsigma/dt
C.  INPUT  LNU = 1, 2  (neutrino-antineutrino)
C.         QQ = t (always negative)  GeV**2
C.         S  = (c.m energy)**2      GeV**2
C.  OUTPUT =  10**-38 cm+2/GeV**2
C-----------------------------------------------------
      REAL*8 MN, MN2, MN4, MD,MD2, MD4
      DATA MN /0.938/
      DATA PI /3.1415926/
      GF = (1.1664 * 1.97)
      GF2 = GF*GF
      MN2 = MN*MN
      MN4 = MN2*MN2
      MD2 = MD*MD
      MD4 = MD2*MD2
      AML2 = AML*AML
      AML4 = AML2*AML2
      VQ  = (MN2 - MD2 - QQ)/2.
      VPI = (MN2 + MD2 - QQ)/2.
      VK  = (S + QQ - MN2 - AML2)/2.
      PIK = (S - MN2)/2.
      QK = (AML2 - QQ)/2.
      PIQ = (QQ + MN2 - MD2)/2.
      Q = SQRT(-QQ)
      C3V = 2.07*SQRT(EXP(-6.3*Q)*(1.+9*Q))
      C3 = SQRT(3.)*C3V/MN
      C4 = -C3/MD             ! attenzione al segno
      C5A = 1.18/(1.-QQ/0.4225)**2
      C32 = C3**2
      C42 = C4**2
      C5A2 = C5A**2

      IF (LNU .EQ. 1)  THEN
      ANS3=-MD2*VPI*QK*QQ*C32+MD2*VPI*QK*C5A2+2.*MD2*VQ*
     . PIK*QK*C32+2.*MD2*VQ*QK*PIQ*C32+MD4*VPI*QK*QQ*C42-
     . 2.*VK**2*VPI*QQ*C32+2.*VK**2*VPI*C5A2+4.*VK*VPI*VQ*
     . QK*C32+2.*VK*VPI*VQ*C5A2+2.*VPI*VQ**2*QK*C32
      ANS2=2.*MN*MD*MD2*VK**2*QQ*C42-4.*MN*MD*MD2*VK*VQ*QK
     . *C42-2.*MN*MD*MD2*VQ**2*QK*C42-2.*MN*MD*MD2*QK**2*
     . C32-3.*MN*MD*MD2*QK*QQ*C32+MN*MD*MD2*QK*C5A2-MN*MD*
     . MD4*QK*QQ*C42+2.*MN*MD*VK**2*C5A2+2.*MN*MD*VK*VQ*
     . C5A2+4.*MN*C3*C4*MD2*VK**2*QQ-8.*MN*C3*C4*MD2*VK*VQ
     . *QK-4.*MN*C3*C4*MD2*VQ**2*QK-2.*MN*C3*C4*MD4*QK*QQ-
     . 4.*MN*C3*C5A*MD2*VK*QQ+4.*MN*C3*C5A*MD2*VQ*QK-2.*MD*
     . C3*C4*MD2*VK*PIK*QQ+2.*MD*C3*C4*MD2*VK*QK*PIQ+2.*MD
     . *C3*C4*MD2*VPI*QK*QQ+2.*MD*C3*C4*MD2*VQ*PIK*QK+2.*
     . MD*C3*C4*MD2*VQ*QK*PIQ-2.*MD*C3*C4*VK**2*VPI*QQ+4.*
     . MD*C3*C4*VK*VPI*VQ*QK+2.*MD*C3*C4*VPI*VQ**2*QK-MD*
     . C3*C5A*MD2*PIK*QQ+MD*C3*C5A*MD2*QK*PIQ-3.*MD*C3*C5A
     . *VK*VPI*QQ+MD*C3*C5A*VK*VQ*PIQ+3.*MD*C3*C5A*VPI*VQ*
     . QK-MD*C3*C5A*VQ**2*PIK+C4*C5A*MD2*VK*VPI*QQ+C4*C5A*
     . MD2*VK*VQ*PIQ-C4*C5A*MD2*VPI*VQ*QK-C4*C5A*MD2*VQ**2
     . *PIK-C4*C5A*MD4*PIK*QQ+C4*C5A*MD4*QK*PIQ-2.*MD2*VK
     . **2*VPI*QQ*C42+4.*MD2*VK*VPI*VQ*QK*C42-2.*MD2*VK*
     . PIK*QQ*C32+2.*MD2*VK*QK*PIQ*C32+2.*MD2*VPI*VQ**2*QK
     . *C42-2.*MD2*VPI*QK**2*C32+ANS3
      ELSE
      ANS3=-MD2*VPI*QK*QQ*C32+MD2*VPI*QK*C5A2+2.*MD2*VQ*
     . PIK*QK*C32+2.*MD2*VQ*QK*PIQ*C32+MD4*VPI*QK*QQ*C42-
     . 2.*VK**2*VPI*QQ*C32+2.*VK**2*VPI*C5A2+4.*VK*VPI*VQ*
     . QK*C32+2.*VK*VPI*VQ*C5A2+2.*VPI*VQ**2*QK*C32
      ANS2=2.*MN*MD*MD2*VK**2*QQ*C42-4.*MN*MD*MD2*VK*VQ*QK
     . *C42-2.*MN*MD*MD2*VQ**2*QK*C42-2.*MN*MD*MD2*QK**2*
     . C32-3.*MN*MD*MD2*QK*QQ*C32+MN*MD*MD2*QK*C5A2-MN*MD*
     . MD4*QK*QQ*C42+2.*MN*MD*VK**2*C5A2+2.*MN*MD*VK*VQ*
     . C5A2+4.*MN*C3*C4*MD2*VK**2*QQ-8.*MN*C3*C4*MD2*VK*VQ
     . *QK-4.*MN*C3*C4*MD2*VQ**2*QK-2.*MN*C3*C4*MD4*QK*QQ+
     . 4.*MN*C3*C5A*MD2*VK*QQ-4.*MN*C3*C5A*MD2*VQ*QK-2.*MD*
     . C3*C4*MD2*VK*PIK*QQ+2.*MD*C3*C4*MD2*VK*QK*PIQ+2.*MD
     . *C3*C4*MD2*VPI*QK*QQ+2.*MD*C3*C4*MD2*VQ*PIK*QK+2.*
     . MD*C3*C4*MD2*VQ*QK*PIQ-2.*MD*C3*C4*VK**2*VPI*QQ+4.*
     . MD*C3*C4*VK*VPI*VQ*QK+2.*MD*C3*C4*VPI*VQ**2*QK+MD*
     . C3*C5A*MD2*PIK*QQ-MD*C3*C5A*MD2*QK*PIQ+3.*MD*C3*C5A
     . *VK*VPI*QQ-MD*C3*C5A*VK*VQ*PIQ-3.*MD*C3*C5A*VPI*VQ*
     . QK+MD*C3*C5A*VQ**2*PIK-C4*C5A*MD2*VK*VPI*QQ-C4*C5A*
     . MD2*VK*VQ*PIQ+C4*C5A*MD2*VPI*VQ*QK+C4*C5A*MD2*VQ**2
     . *PIK+C4*C5A*MD4*PIK*QQ-C4*C5A*MD4*QK*PIQ-2.*MD2*VK
     . **2*VPI*QQ*C42+4.*MD2*VK*VPI*VQ*QK*C42-2.*MD2*VK*
     . PIK*QQ*C32+2.*MD2*VK*QK*PIQ*C32+2.*MD2*VPI*VQ**2*QK
     . *C42-2.*MD2*VPI*QK**2*C32+ANS3
      ENDIF
      ANS1=32.*ANS2
      ANS=ANS1/(3.*MD2)
      P1CM = (S-MN2)/(2.*SQRT(S))
      DSIGMA_DELTA  = GF2/2. * ANS/(64.*PI*S*P1CM**2)
      RETURN
      END

      SUBROUTINE QGAUS(A,B,SS,ENU,LTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(5),W(5)
      DATA X/.1488743389D0,.4333953941D0,
     & .6794095682D0,.8650633666D0,.9739065285D0
     */
      DATA W/.2955242247D0,.2692667193D0,
     & .2190863625D0,.1494513491D0,.0666713443D0
     */
      XM=0.5D0*(B+A)
      XR=0.5D0*(B-A)
      SS=0
      DO 11 J=1,5
        DX=XR*X(J)
        SS=SS+W(J)*(DSQEL_Q2(LTYP,ENU,XM+DX)+
     *	DSQEL_Q2(LTYP,ENU,XM-DX))
11    CONTINUE
      SS=XR*SS
      RETURN
      END
