C********************************************************************* 
C********************************************************************* 
C*                                                                  ** 
C*                                                 December 1993    ** 
C*                                                                  ** 
C*   The Lund Monte Carlo for Jet Fragmentation and e+e- Physics    ** 
C*                                                                  ** 
C*                        JETSET version 7.4                        ** 
C*                                                                  ** 
C*                        Torbjorn Sjostrand                        ** 
C*                Department of theoretical physics 2               ** 
C*                        University of Lund                        ** 
C*               Solvegatan 14A, S-223 62 Lund, Sweden              **
C*                    E-mail torbjorn@thep.lu.se                    ** 
C*                    phone +46 - 46 - 222 48 16                    ** 
C*                                                                  ** 
C*          LUSHOW is written together with Mats Bengtsson          ** 
C*                                                                  **
C*   The latest program version and documentation is found on WWW   **
C*         http://thep.lu.se/tf2/staff/torbjorn/Welcome.html        **
C*                                                                  ** 
C*        Copyright Torbjorn Sjostrand and CERN, Geneva 1993        ** 
C*                                                                  ** 
C********************************************************************* 
C********************************************************************* 
C                                                                    * 
C  List of subprograms in order of appearance, with main purpose     * 
C  (S = subroutine, F = function, B = block data)                    * 
C                                                                    * 
C  S   LU1ENT   to fill one entry (= parton or particle)             * 
C  S   LU2ENT   to fill two entries                                  * 
C  S   LU3ENT   to fill three entries                                * 
C  S   LU4ENT   to fill four entries                                 * 
C  S   LUJOIN   to connect entries with colour flow information      * 
C  S   LUGIVE   to fill (or query) commonblock variables             * 
C  S   LUEXEC   to administrate fragmentation and decay chain        * 
C  S   LUPREP   to rearrange showered partons along strings          * 
C  S   LUSTRF   to do string fragmentation of jet system             * 
C  S   LUINDF   to do independent fragmentation of one or many jets  * 
C  S   LUDECY   to do the decay of a particle                        * 
C  S   LUKFDI   to select parton and hadron flavours in fragm        * 
C  S   LUPTDI   to select transverse momenta in fragm                * 
C  S   LUZDIS   to select longitudinal scaling variable in fragm     * 
C  S   LUSHOW   to do timelike parton shower evolution               * 
C  S   LUBOEI   to include Bose-Einstein effects (crudely)           * 
C  F   ULMASS   to give the mass of a particle or parton             * 
C  S   LUNAME   to give the name of a particle or parton             * 
C  F   LUCHGE   to give three times the electric charge              * 
C  F   LUCOMP   to compress standard KF flavour code to internal KC  * 
C  S   LUERRM   to write error messages and abort faulty run         * 
C  F   ULALEM   to give the alpha_electromagnetic value              * 
C  F   ULALPS   to give the alpha_strong value                       * 
C  F   ULANGL   to give the angle from known x and y components      * 
C  F   RLU      to provide a random number generator                 * 
C  S   RLUGET   to save the state of the random number generator     * 
C  S   RLUSET   to set the state of the random number generator      * 
C  S   LUROBO   to rotate and/or boost an event                      * 
C  S   LUEDIT   to remove unwanted entries from record               * 
C  S   LULIST   to list event record or particle data                * 
C  S   LULOGO   to write a logo for JETSET and PYTHIA                * 
C  S   LUUPDA   to update particle data                              * 
C  F   KLU      to provide integer-valued event information          * 
C  F   PLU      to provide real-valued event information             * 
C  S   LUSPHE   to perform sphericity analysis                       * 
C  S   LUTHRU   to perform thrust analysis                           * 
C  S   LUCLUS   to perform three-dimensional cluster analysis        * 
C  S   LUCELL   to perform cluster analysis in (eta, phi, E_T)       * 
C  S   LUJMAS   to give high and low jet mass of event               * 
C  S   LUFOWO   to give Fox-Wolfram moments                          * 
C  S   LUTABU   to analyze events, with tabular output               * 
C                                                                    * 
C  S   LUEEVT   to administrate the generation of an e+e- event      * 
C  S   LUXTOT   to give the total cross-section at given CM energy   * 
C  S   LURADK   to generate initial state photon radiation           * 
C  S   LUXKFL   to select flavour of primary qqbar pair              * 
C  S   LUXJET   to select (matrix element) jet multiplicity          * 
C  S   LUX3JT   to select kinematics of three-jet event              * 
C  S   LUX4JT   to select kinematics of four-jet event               * 
C  S   LUXDIF   to select angular orientation of event               * 
C  S   LUONIA   to perform generation of onium decay to gluons       * 
C                                                                    * 
C  S   LUHEPC   to convert between /LUJETS/ and /HEPEVT/ records     * 
C  S   LUTEST   to test the proper functioning of the package        * 
C  B   LUDATA   to contain default values and particle data          * 
C                                                                    * 
C********************************************************************* 
 
      SUBROUTINE LU1ENT(IP,KF,PE,THE,PHI) 
 
C...Purpose: to store one parton/particle in commonblock LUJETS. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Standard checks. 
      MSTU(28)=0 
      IF(MSTU(12).GE.1) CALL LULIST(0) 
      IPA=MAX(1,IABS(IP)) 
      IF(IPA.GT.MSTU(4)) CALL LUERRM(21, 
     &'(LU1ENT:) writing outside LUJETS memory') 
      KC=LUCOMP(KF) 
      IF(KC.EQ.0) CALL LUERRM(12,'(LU1ENT:) unknown flavour code') 
 
C...Find mass. Reset K, P and V vectors. 
      PM=0. 
      IF(MSTU(10).EQ.1) PM=P(IPA,5) 
      IF(MSTU(10).GE.2) PM=ULMASS(KF) 
      DO 100 J=1,5 
      K(IPA,J)=0 
      P(IPA,J)=0. 
      V(IPA,J)=0. 
  100 CONTINUE 
 
C...Store parton/particle in K and P vectors. 
      K(IPA,1)=1 
      IF(IP.LT.0) K(IPA,1)=2 
      K(IPA,2)=KF 
      P(IPA,5)=PM 
      P(IPA,4)=MAX(PE,PM) 
      PA=SQRT(P(IPA,4)**2-P(IPA,5)**2) 
      P(IPA,1)=PA*SIN(THE)*COS(PHI) 
      P(IPA,2)=PA*SIN(THE)*SIN(PHI) 
      P(IPA,3)=PA*COS(THE) 
 
C...Set N. Optionally fragment/decay. 
      N=IPA 
      IF(IP.EQ.0) CALL LUEXEC 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LU2ENT(IP,KF1,KF2,PECM) 
 
C...Purpose: to store two partons/particles in their CM frame, 
C...with the first along the +z axis. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Standard checks. 
      MSTU(28)=0 
      IF(MSTU(12).GE.1) CALL LULIST(0) 
      IPA=MAX(1,IABS(IP)) 
      IF(IPA.GT.MSTU(4)-1) CALL LUERRM(21, 
     &'(LU2ENT:) writing outside LUJETS memory') 
      KC1=LUCOMP(KF1) 
      KC2=LUCOMP(KF2) 
      IF(KC1.EQ.0.OR.KC2.EQ.0) CALL LUERRM(12, 
     &'(LU2ENT:) unknown flavour code') 
 
C...Find masses. Reset K, P and V vectors. 
      PM1=0. 
      IF(MSTU(10).EQ.1) PM1=P(IPA,5) 
      IF(MSTU(10).GE.2) PM1=ULMASS(KF1) 
      PM2=0. 
      IF(MSTU(10).EQ.1) PM2=P(IPA+1,5) 
      IF(MSTU(10).GE.2) PM2=ULMASS(KF2) 
      DO 110 I=IPA,IPA+1 
      DO 100 J=1,5 
      K(I,J)=0 
      P(I,J)=0. 
      V(I,J)=0. 
  100 CONTINUE 
  110 CONTINUE 
 
C...Check flavours. 
      KQ1=KCHG(KC1,2)*ISIGN(1,KF1) 
      KQ2=KCHG(KC2,2)*ISIGN(1,KF2) 
      IF(MSTU(19).EQ.1) THEN 
        MSTU(19)=0 
      ELSE 
        IF(KQ1+KQ2.NE.0.AND.KQ1+KQ2.NE.4) CALL LUERRM(2, 
     &  '(LU2ENT:) unphysical flavour combination') 
      ENDIF 
      K(IPA,2)=KF1 
      K(IPA+1,2)=KF2 
 
C...Store partons/particles in K vectors for normal case. 
      IF(IP.GE.0) THEN 
        K(IPA,1)=1 
        IF(KQ1.NE.0.AND.KQ2.NE.0) K(IPA,1)=2 
        K(IPA+1,1)=1 
 
C...Store partons in K vectors for parton shower evolution. 
      ELSE 
        K(IPA,1)=3 
        K(IPA+1,1)=3 
        K(IPA,4)=MSTU(5)*(IPA+1) 
        K(IPA,5)=K(IPA,4) 
        K(IPA+1,4)=MSTU(5)*IPA 
        K(IPA+1,5)=K(IPA+1,4) 
      ENDIF 
 
C...Check kinematics and store partons/particles in P vectors. 
      IF(PECM.LE.PM1+PM2) CALL LUERRM(13, 
     &'(LU2ENT:) energy smaller than sum of masses') 
      PA=SQRT(MAX(0.,(PECM**2-PM1**2-PM2**2)**2-(2.*PM1*PM2)**2))/ 
     &(2.*PECM) 
      P(IPA,3)=PA 
      P(IPA,4)=SQRT(PM1**2+PA**2) 
      P(IPA,5)=PM1 
      P(IPA+1,3)=-PA 
      P(IPA+1,4)=SQRT(PM2**2+PA**2) 
      P(IPA+1,5)=PM2 
 
C...Set N. Optionally fragment/decay. 
      N=IPA+1 
      IF(IP.EQ.0) CALL LUEXEC 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LU3ENT(IP,KF1,KF2,KF3,PECM,X1,X3) 
 
C...Purpose: to store three partons or particles in their CM frame, 
C...with the first along the +z axis and the third in the (x,z) 
C...plane with x > 0. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Standard checks. 
      MSTU(28)=0 
      IF(MSTU(12).GE.1) CALL LULIST(0) 
      IPA=MAX(1,IABS(IP)) 
      IF(IPA.GT.MSTU(4)-2) CALL LUERRM(21, 
     &'(LU3ENT:) writing outside LUJETS memory') 
      KC1=LUCOMP(KF1) 
      KC2=LUCOMP(KF2) 
      KC3=LUCOMP(KF3) 
      IF(KC1.EQ.0.OR.KC2.EQ.0.OR.KC3.EQ.0) CALL LUERRM(12, 
     &'(LU3ENT:) unknown flavour code') 
 
C...Find masses. Reset K, P and V vectors. 
      PM1=0. 
      IF(MSTU(10).EQ.1) PM1=P(IPA,5) 
      IF(MSTU(10).GE.2) PM1=ULMASS(KF1) 
      PM2=0. 
      IF(MSTU(10).EQ.1) PM2=P(IPA+1,5) 
      IF(MSTU(10).GE.2) PM2=ULMASS(KF2) 
      PM3=0. 
      IF(MSTU(10).EQ.1) PM3=P(IPA+2,5) 
      IF(MSTU(10).GE.2) PM3=ULMASS(KF3) 
      DO 110 I=IPA,IPA+2 
      DO 100 J=1,5 
      K(I,J)=0 
      P(I,J)=0. 
      V(I,J)=0. 
  100 CONTINUE 
  110 CONTINUE 
 
C...Check flavours. 
      KQ1=KCHG(KC1,2)*ISIGN(1,KF1) 
      KQ2=KCHG(KC2,2)*ISIGN(1,KF2) 
      KQ3=KCHG(KC3,2)*ISIGN(1,KF3) 
      IF(MSTU(19).EQ.1) THEN 
        MSTU(19)=0 
      ELSEIF(KQ1.EQ.0.AND.KQ2.EQ.0.AND.KQ3.EQ.0) THEN 
      ELSEIF(KQ1.NE.0.AND.KQ2.EQ.2.AND.(KQ1+KQ3.EQ.0.OR. 
     &KQ1+KQ3.EQ.4)) THEN 
      ELSE 
        CALL LUERRM(2,'(LU3ENT:) unphysical flavour combination') 
      ENDIF 
      K(IPA,2)=KF1 
      K(IPA+1,2)=KF2 
      K(IPA+2,2)=KF3 
 
C...Store partons/particles in K vectors for normal case. 
      IF(IP.GE.0) THEN 
        K(IPA,1)=1 
        IF(KQ1.NE.0.AND.(KQ2.NE.0.OR.KQ3.NE.0)) K(IPA,1)=2 
        K(IPA+1,1)=1 
        IF(KQ2.NE.0.AND.KQ3.NE.0) K(IPA+1,1)=2 
        K(IPA+2,1)=1 
 
C...Store partons in K vectors for parton shower evolution. 
      ELSE 
        K(IPA,1)=3 
        K(IPA+1,1)=3 
        K(IPA+2,1)=3 
        KCS=4 
        IF(KQ1.EQ.-1) KCS=5 
        K(IPA,KCS)=MSTU(5)*(IPA+1) 
        K(IPA,9-KCS)=MSTU(5)*(IPA+2) 
        K(IPA+1,KCS)=MSTU(5)*(IPA+2) 
        K(IPA+1,9-KCS)=MSTU(5)*IPA 
        K(IPA+2,KCS)=MSTU(5)*IPA 
        K(IPA+2,9-KCS)=MSTU(5)*(IPA+1) 
      ENDIF 
 
C...Check kinematics. 
      MKERR=0 
      IF(0.5*X1*PECM.LE.PM1.OR.0.5*(2.-X1-X3)*PECM.LE.PM2.OR. 
     &0.5*X3*PECM.LE.PM3) MKERR=1 
      PA1=SQRT(MAX(1E-10,(0.5*X1*PECM)**2-PM1**2)) 
      PA2=SQRT(MAX(1E-10,(0.5*(2.-X1-X3)*PECM)**2-PM2**2)) 
      PA3=SQRT(MAX(1E-10,(0.5*X3*PECM)**2-PM3**2)) 
      CTHE2=(PA3**2-PA1**2-PA2**2)/(2.*PA1*PA2) 
      CTHE3=(PA2**2-PA1**2-PA3**2)/(2.*PA1*PA3) 
      IF(ABS(CTHE2).GE.1.001.OR.ABS(CTHE3).GE.1.001) MKERR=1 
      CTHE3=MAX(-1.,MIN(1.,CTHE3)) 
      IF(MKERR.NE.0) CALL LUERRM(13, 
     &'(LU3ENT:) unphysical kinematical variable setup') 
 
C...Store partons/particles in P vectors. 
      P(IPA,3)=PA1 
      P(IPA,4)=SQRT(PA1**2+PM1**2) 
      P(IPA,5)=PM1 
      P(IPA+2,1)=PA3*SQRT(1.-CTHE3**2) 
      P(IPA+2,3)=PA3*CTHE3 
      P(IPA+2,4)=SQRT(PA3**2+PM3**2) 
      P(IPA+2,5)=PM3 
      P(IPA+1,1)=-P(IPA+2,1) 
      P(IPA+1,3)=-P(IPA,3)-P(IPA+2,3) 
      P(IPA+1,4)=SQRT(P(IPA+1,1)**2+P(IPA+1,3)**2+PM2**2) 
      P(IPA+1,5)=PM2 
 
C...Set N. Optionally fragment/decay. 
      N=IPA+2 
      IF(IP.EQ.0) CALL LUEXEC 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LU4ENT(IP,KF1,KF2,KF3,KF4,PECM,X1,X2,X4,X12,X14) 
 
C...Purpose: to store four partons or particles in their CM frame, with 
C...the first along the +z axis, the last in the xz plane with x > 0 
C...and the second having y < 0 and y > 0 with equal probability. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Standard checks. 
      MSTU(28)=0 
      IF(MSTU(12).GE.1) CALL LULIST(0) 
      IPA=MAX(1,IABS(IP)) 
      IF(IPA.GT.MSTU(4)-3) CALL LUERRM(21, 
     &'(LU4ENT:) writing outside LUJETS momory') 
      KC1=LUCOMP(KF1) 
      KC2=LUCOMP(KF2) 
      KC3=LUCOMP(KF3) 
      KC4=LUCOMP(KF4) 
      IF(KC1.EQ.0.OR.KC2.EQ.0.OR.KC3.EQ.0.OR.KC4.EQ.0) CALL LUERRM(12, 
     &'(LU4ENT:) unknown flavour code') 
 
C...Find masses. Reset K, P and V vectors. 
      PM1=0. 
      IF(MSTU(10).EQ.1) PM1=P(IPA,5) 
      IF(MSTU(10).GE.2) PM1=ULMASS(KF1) 
      PM2=0. 
      IF(MSTU(10).EQ.1) PM2=P(IPA+1,5) 
      IF(MSTU(10).GE.2) PM2=ULMASS(KF2) 
      PM3=0. 
      IF(MSTU(10).EQ.1) PM3=P(IPA+2,5) 
      IF(MSTU(10).GE.2) PM3=ULMASS(KF3) 
      PM4=0. 
      IF(MSTU(10).EQ.1) PM4=P(IPA+3,5) 
      IF(MSTU(10).GE.2) PM4=ULMASS(KF4) 
      DO 110 I=IPA,IPA+3 
      DO 100 J=1,5 
      K(I,J)=0 
      P(I,J)=0. 
      V(I,J)=0. 
  100 CONTINUE 
  110 CONTINUE 
 
C...Check flavours. 
      KQ1=KCHG(KC1,2)*ISIGN(1,KF1) 
      KQ2=KCHG(KC2,2)*ISIGN(1,KF2) 
      KQ3=KCHG(KC3,2)*ISIGN(1,KF3) 
      KQ4=KCHG(KC4,2)*ISIGN(1,KF4) 
      IF(MSTU(19).EQ.1) THEN 
        MSTU(19)=0 
      ELSEIF(KQ1.EQ.0.AND.KQ2.EQ.0.AND.KQ3.EQ.0.AND.KQ4.EQ.0) THEN 
      ELSEIF(KQ1.NE.0.AND.KQ2.EQ.2.AND.KQ3.EQ.2.AND.(KQ1+KQ4.EQ.0.OR. 
     &KQ1+KQ4.EQ.4)) THEN 
      ELSEIF(KQ1.NE.0.AND.KQ1+KQ2.EQ.0.AND.KQ3.NE.0.AND.KQ3+KQ4.EQ.0.) 
     &THEN 
      ELSE 
        CALL LUERRM(2,'(LU4ENT:) unphysical flavour combination') 
      ENDIF 
      K(IPA,2)=KF1 
      K(IPA+1,2)=KF2 
      K(IPA+2,2)=KF3 
      K(IPA+3,2)=KF4 
 
C...Store partons/particles in K vectors for normal case. 
      IF(IP.GE.0) THEN 
        K(IPA,1)=1 
        IF(KQ1.NE.0.AND.(KQ2.NE.0.OR.KQ3.NE.0.OR.KQ4.NE.0)) K(IPA,1)=2 
        K(IPA+1,1)=1 
        IF(KQ2.NE.0.AND.KQ1+KQ2.NE.0.AND.(KQ3.NE.0.OR.KQ4.NE.0)) 
     &  K(IPA+1,1)=2 
        K(IPA+2,1)=1 
        IF(KQ3.NE.0.AND.KQ4.NE.0) K(IPA+2,1)=2 
        K(IPA+3,1)=1 
 
C...Store partons for parton shower evolution from q-g-g-qbar or 
C...g-g-g-g event. 
      ELSEIF(KQ1+KQ2.NE.0) THEN 
        K(IPA,1)=3 
        K(IPA+1,1)=3 
        K(IPA+2,1)=3 
        K(IPA+3,1)=3 
        KCS=4 
        IF(KQ1.EQ.-1) KCS=5 
        K(IPA,KCS)=MSTU(5)*(IPA+1) 
        K(IPA,9-KCS)=MSTU(5)*(IPA+3) 
        K(IPA+1,KCS)=MSTU(5)*(IPA+2) 
        K(IPA+1,9-KCS)=MSTU(5)*IPA 
        K(IPA+2,KCS)=MSTU(5)*(IPA+3) 
        K(IPA+2,9-KCS)=MSTU(5)*(IPA+1) 
        K(IPA+3,KCS)=MSTU(5)*IPA 
        K(IPA+3,9-KCS)=MSTU(5)*(IPA+2) 
 
C...Store partons for parton shower evolution from q-qbar-q-qbar event. 
      ELSE 
        K(IPA,1)=3 
        K(IPA+1,1)=3 
        K(IPA+2,1)=3 
        K(IPA+3,1)=3 
        K(IPA,4)=MSTU(5)*(IPA+1) 
        K(IPA,5)=K(IPA,4) 
        K(IPA+1,4)=MSTU(5)*IPA 
        K(IPA+1,5)=K(IPA+1,4) 
        K(IPA+2,4)=MSTU(5)*(IPA+3) 
        K(IPA+2,5)=K(IPA+2,4) 
        K(IPA+3,4)=MSTU(5)*(IPA+2) 
        K(IPA+3,5)=K(IPA+3,4) 
      ENDIF 
 
C...Check kinematics. 
      MKERR=0 
      IF(0.5*X1*PECM.LE.PM1.OR.0.5*X2*PECM.LE.PM2.OR.0.5*(2.-X1-X2-X4)* 
     &PECM.LE.PM3.OR.0.5*X4*PECM.LE.PM4) MKERR=1 
      PA1=SQRT(MAX(1E-10,(0.5*X1*PECM)**2-PM1**2)) 
      PA2=SQRT(MAX(1E-10,(0.5*X2*PECM)**2-PM2**2)) 
      PA4=SQRT(MAX(1E-10,(0.5*X4*PECM)**2-PM4**2)) 
      X24=X1+X2+X4-1.-X12-X14+(PM3**2-PM1**2-PM2**2-PM4**2)/PECM**2 
      CTHE4=(X1*X4-2.*X14)*PECM**2/(4.*PA1*PA4) 
      IF(ABS(CTHE4).GE.1.002) MKERR=1 
      CTHE4=MAX(-1.,MIN(1.,CTHE4)) 
      STHE4=SQRT(1.-CTHE4**2) 
      CTHE2=(X1*X2-2.*X12)*PECM**2/(4.*PA1*PA2) 
      IF(ABS(CTHE2).GE.1.002) MKERR=1 
      CTHE2=MAX(-1.,MIN(1.,CTHE2)) 
      STHE2=SQRT(1.-CTHE2**2) 
      CPHI2=((X2*X4-2.*X24)*PECM**2-4.*PA2*CTHE2*PA4*CTHE4)/ 
     &MAX(1E-8*PECM**2,4.*PA2*STHE2*PA4*STHE4) 
      IF(ABS(CPHI2).GE.1.05) MKERR=1 
      CPHI2=MAX(-1.,MIN(1.,CPHI2)) 
      IF(MKERR.EQ.1) CALL LUERRM(13, 
     &'(LU4ENT:) unphysical kinematical variable setup') 
 
C...Store partons/particles in P vectors. 
      P(IPA,3)=PA1 
      P(IPA,4)=SQRT(PA1**2+PM1**2) 
      P(IPA,5)=PM1 
      P(IPA+3,1)=PA4*STHE4 
      P(IPA+3,3)=PA4*CTHE4 
      P(IPA+3,4)=SQRT(PA4**2+PM4**2) 
      P(IPA+3,5)=PM4 
      P(IPA+1,1)=PA2*STHE2*CPHI2 
      P(IPA+1,2)=PA2*STHE2*SQRT(1.-CPHI2**2)*(-1.)**INT(RLU(0)+0.5) 
      P(IPA+1,3)=PA2*CTHE2 
      P(IPA+1,4)=SQRT(PA2**2+PM2**2) 
      P(IPA+1,5)=PM2 
      P(IPA+2,1)=-P(IPA+1,1)-P(IPA+3,1) 
      P(IPA+2,2)=-P(IPA+1,2) 
      P(IPA+2,3)=-P(IPA,3)-P(IPA+1,3)-P(IPA+3,3) 
      P(IPA+2,4)=SQRT(P(IPA+2,1)**2+P(IPA+2,2)**2+P(IPA+2,3)**2+PM3**2) 
      P(IPA+2,5)=PM3 
 
C...Set N. Optionally fragment/decay. 
      N=IPA+3 
      IF(IP.EQ.0) CALL LUEXEC 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUJOIN(NJOIN,IJOIN) 
 
C...Purpose: to connect a sequence of partons with colour flow indices, 
C...as required for subsequent shower evolution (or other operations). 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION IJOIN(*) 
 
C...Check that partons are of right types to be connected. 
      IF(NJOIN.LT.2) GOTO 120 
      KQSUM=0 
      DO 100 IJN=1,NJOIN 
      I=IJOIN(IJN) 
      IF(I.LE.0.OR.I.GT.N) GOTO 120 
      IF(K(I,1).LT.1.OR.K(I,1).GT.3) GOTO 120 
      KC=LUCOMP(K(I,2)) 
      IF(KC.EQ.0) GOTO 120 
      KQ=KCHG(KC,2)*ISIGN(1,K(I,2)) 
      IF(KQ.EQ.0) GOTO 120 
      IF(IJN.NE.1.AND.IJN.NE.NJOIN.AND.KQ.NE.2) GOTO 120 
      IF(KQ.NE.2) KQSUM=KQSUM+KQ 
      IF(IJN.EQ.1) KQS=KQ 
  100 CONTINUE 
      IF(KQSUM.NE.0) GOTO 120 
 
C...Connect the partons sequentially (closing for gluon loop). 
      KCS=(9-KQS)/2 
      IF(KQS.EQ.2) KCS=INT(4.5+RLU(0)) 
      DO 110 IJN=1,NJOIN 
      I=IJOIN(IJN) 
      K(I,1)=3 
      IF(IJN.NE.1) IP=IJOIN(IJN-1) 
      IF(IJN.EQ.1) IP=IJOIN(NJOIN) 
      IF(IJN.NE.NJOIN) IN=IJOIN(IJN+1) 
      IF(IJN.EQ.NJOIN) IN=IJOIN(1) 
      K(I,KCS)=MSTU(5)*IN 
      K(I,9-KCS)=MSTU(5)*IP 
      IF(IJN.EQ.1.AND.KQS.NE.2) K(I,9-KCS)=0 
      IF(IJN.EQ.NJOIN.AND.KQS.NE.2) K(I,KCS)=0 
  110 CONTINUE 
 
C...Error exit: no action taken. 
      RETURN 
  120 CALL LUERRM(12, 
     &'(LUJOIN:) given entries can not be joined by one string') 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUGIVE(CHIN) 
 
C...Purpose: to set values of commonblock variables (also in PYTHIA!). 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      COMMON/LUDAT4/CHAF(500) 
      CHARACTER CHAF*8 
      COMMON/LUDATR/MRLU(6),RRLU(100) 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      COMMON/PYINT1/MINT(400),VINT(400) 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000) 
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      COMMON/PYINT6/PROC(0:200) 
      COMMON/PYINT7/SIGT(0:6,0:6,0:5) 
      CHARACTER PROC*28 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/,/LUDAT4/,/LUDATR/ 
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT3/,/PYINT4/, 
     &/PYINT5/,/PYINT6/,/PYINT7/ 
      CHARACTER CHIN*(*),CHFIX*104,CHBIT*104,CHOLD*8,CHNEW*8,CHOLD2*28, 
     &CHNEW2*28,CHNAM*4,CHVAR(43)*4,CHALP(2)*26,CHIND*8,CHINI*10, 
     &CHINR*16 
      DIMENSION MSVAR(43,8) 
 
C...For each variable to be translated give: name, 
C...integer/real/character, no. of indices, lower&upper index bounds. 
      DATA CHVAR/'N','K','P','V','MSTU','PARU','MSTJ','PARJ','KCHG', 
     &'PMAS','PARF','VCKM','MDCY','MDME','BRAT','KFDP','CHAF','MRLU', 
     &'RRLU','MSEL','MSUB','KFIN','CKIN','MSTP','PARP','MSTI','PARI', 
     &'MINT','VINT','ISET','KFPR','COEF','ICOL','XSFX','ISIG','SIGH', 
     &'WIDP','WIDE','WIDS','NGEN','XSEC','PROC','SIGT'/ 
      DATA ((MSVAR(I,J),J=1,8),I=1,43)/ 1,7*0,  1,2,1,4000,1,5,2*0, 
     & 2,2,1,4000,1,5,2*0,  2,2,1,4000,1,5,2*0,  1,1,1,200,4*0, 
     & 2,1,1,200,4*0,  1,1,1,200,4*0,  2,1,1,200,4*0, 
     & 1,2,1,500,1,3,2*0,  2,2,1,500,1,4,2*0,  2,1,1,2000,4*0, 
     & 2,2,1,4,1,4,2*0,  1,2,1,500,1,3,2*0,  1,2,1,2000,1,2,2*0, 
     & 2,1,1,2000,4*0,  1,2,1,2000,1,5,2*0,  3,1,1,500,4*0, 
     & 1,1,1,6,4*0,  2,1,1,100,4*0, 
     & 1,7*0,  1,1,1,200,4*0,  1,2,1,2,-40,40,2*0,  2,1,1,200,4*0, 
     & 1,1,1,200,4*0,  2,1,1,200,4*0,  1,1,1,200,4*0,  2,1,1,200,4*0, 
     & 1,1,1,400,4*0,  2,1,1,400,4*0,  1,1,1,200,4*0, 
     & 1,2,1,200,1,2,2*0,  2,2,1,200,1,20,2*0,  1,3,1,40,1,4,1,2, 
     & 2,2,1,2,-40,40,2*0,  1,2,1,1000,1,3,2*0,  2,1,1,1000,4*0, 
     & 2,2,21,40,0,40,2*0,  2,2,21,40,0,40,2*0,  2,2,21,40,1,3,2*0, 
     & 1,2,0,200,1,3,2*0,  2,2,0,200,1,3,2*0,  4,1,0,200,4*0, 
     & 2,3,0,6,0,6,0,5/ 
      DATA CHALP/'abcdefghijklmnopqrstuvwxyz', 
     &'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/ 
 
C...Length of character variable. Subdivide it into instructions. 
      IF(MSTU(12).GE.1) CALL LULIST(0) 
      CHBIT=CHIN//' ' 
      LBIT=101 
  100 LBIT=LBIT-1 
      IF(CHBIT(LBIT:LBIT).EQ.' ') GOTO 100 
      LTOT=0 
      DO 110 LCOM=1,LBIT 
      IF(CHBIT(LCOM:LCOM).EQ.' ') GOTO 110 
      LTOT=LTOT+1 
      CHFIX(LTOT:LTOT)=CHBIT(LCOM:LCOM) 
  110 CONTINUE 
      LLOW=0 
  120 LHIG=LLOW+1 
  130 LHIG=LHIG+1 
      IF(LHIG.LE.LTOT.AND.CHFIX(LHIG:LHIG).NE.';') GOTO 130 
      LBIT=LHIG-LLOW-1 
      CHBIT(1:LBIT)=CHFIX(LLOW+1:LHIG-1) 
 
C...Identify commonblock variable. 
      LNAM=1 
  140 LNAM=LNAM+1 
      IF(CHBIT(LNAM:LNAM).NE.'('.AND.CHBIT(LNAM:LNAM).NE.'='.AND. 
     &LNAM.LE.4) GOTO 140 
      CHNAM=CHBIT(1:LNAM-1)//' ' 
      DO 160 LCOM=1,LNAM-1 
      DO 150 LALP=1,26 
      IF(CHNAM(LCOM:LCOM).EQ.CHALP(1)(LALP:LALP)) CHNAM(LCOM:LCOM)= 
     &CHALP(2)(LALP:LALP) 
  150 CONTINUE 
  160 CONTINUE 
      IVAR=0 
      DO 170 IV=1,43 
      IF(CHNAM.EQ.CHVAR(IV)) IVAR=IV 
  170 CONTINUE 
      IF(IVAR.EQ.0) THEN 
        CALL LUERRM(18,'(LUGIVE:) do not recognize variable '//CHNAM) 
        LLOW=LHIG 
        IF(LLOW.LT.LTOT) GOTO 120 
        RETURN 
      ENDIF 
 
C...Identify any indices. 
      I1=0 
      I2=0 
      I3=0 
      NINDX=0 
      IF(CHBIT(LNAM:LNAM).EQ.'(') THEN 
        LIND=LNAM 
  180   LIND=LIND+1 
        IF(CHBIT(LIND:LIND).NE.')'.AND.CHBIT(LIND:LIND).NE.',') GOTO 180 
        CHIND=' ' 
        IF((CHBIT(LNAM+1:LNAM+1).EQ.'C'.OR.CHBIT(LNAM+1:LNAM+1).EQ.'c'). 
     &  AND.(IVAR.EQ.9.OR.IVAR.EQ.10.OR.IVAR.EQ.13.OR.IVAR.EQ.17)) THEN 
          CHIND(LNAM-LIND+11:8)=CHBIT(LNAM+2:LIND-1) 
          READ(CHIND,'(I8)') KF 
          I1=LUCOMP(KF) 
        ELSEIF(CHBIT(LNAM+1:LNAM+1).EQ.'C'.OR.CHBIT(LNAM+1:LNAM+1).EQ. 
     &  'c') THEN 
          CALL LUERRM(18,'(LUGIVE:) not allowed to use C index for '// 
     &    CHNAM) 
          LLOW=LHIG 
          IF(LLOW.LT.LTOT) GOTO 120 
          RETURN 
        ELSE 
          CHIND(LNAM-LIND+10:8)=CHBIT(LNAM+1:LIND-1) 
          READ(CHIND,'(I8)') I1 
        ENDIF 
        LNAM=LIND 
        IF(CHBIT(LNAM:LNAM).EQ.')') LNAM=LNAM+1 
        NINDX=1 
      ENDIF 
      IF(CHBIT(LNAM:LNAM).EQ.',') THEN 
        LIND=LNAM 
  190   LIND=LIND+1 
        IF(CHBIT(LIND:LIND).NE.')'.AND.CHBIT(LIND:LIND).NE.',') GOTO 190 
        CHIND=' ' 
        CHIND(LNAM-LIND+10:8)=CHBIT(LNAM+1:LIND-1) 
        READ(CHIND,'(I8)') I2 
        LNAM=LIND 
        IF(CHBIT(LNAM:LNAM).EQ.')') LNAM=LNAM+1 
        NINDX=2 
      ENDIF 
      IF(CHBIT(LNAM:LNAM).EQ.',') THEN 
        LIND=LNAM 
  200   LIND=LIND+1 
        IF(CHBIT(LIND:LIND).NE.')'.AND.CHBIT(LIND:LIND).NE.',') GOTO 200 
        CHIND=' ' 
        CHIND(LNAM-LIND+10:8)=CHBIT(LNAM+1:LIND-1) 
        READ(CHIND,'(I8)') I3 
        LNAM=LIND+1 
        NINDX=3 
      ENDIF 
 
C...Check that indices allowed. 
      IERR=0 
      IF(NINDX.NE.MSVAR(IVAR,2)) IERR=1 
      IF(NINDX.GE.1.AND.(I1.LT.MSVAR(IVAR,3).OR.I1.GT.MSVAR(IVAR,4))) 
     &IERR=2 
      IF(NINDX.GE.2.AND.(I2.LT.MSVAR(IVAR,5).OR.I2.GT.MSVAR(IVAR,6))) 
     &IERR=3 
      IF(NINDX.EQ.3.AND.(I3.LT.MSVAR(IVAR,7).OR.I3.GT.MSVAR(IVAR,8))) 
     &IERR=4 
      IF(CHBIT(LNAM:LNAM).NE.'=') IERR=5 
      IF(IERR.GE.1) THEN 
        CALL LUERRM(18,'(LUGIVE:) unallowed indices for '// 
     &  CHBIT(1:LNAM-1)) 
        LLOW=LHIG 
        IF(LLOW.LT.LTOT) GOTO 120 
        RETURN 
      ENDIF 
 
C...Save old value of variable. 
      IF(IVAR.EQ.1) THEN 
        IOLD=N 
      ELSEIF(IVAR.EQ.2) THEN 
        IOLD=K(I1,I2) 
      ELSEIF(IVAR.EQ.3) THEN 
        ROLD=P(I1,I2) 
      ELSEIF(IVAR.EQ.4) THEN 
        ROLD=V(I1,I2) 
      ELSEIF(IVAR.EQ.5) THEN 
        IOLD=MSTU(I1) 
      ELSEIF(IVAR.EQ.6) THEN 
        ROLD=PARU(I1) 
      ELSEIF(IVAR.EQ.7) THEN 
        IOLD=MSTJ(I1) 
      ELSEIF(IVAR.EQ.8) THEN 
        ROLD=PARJ(I1) 
      ELSEIF(IVAR.EQ.9) THEN 
        IOLD=KCHG(I1,I2) 
      ELSEIF(IVAR.EQ.10) THEN 
        ROLD=PMAS(I1,I2) 
      ELSEIF(IVAR.EQ.11) THEN 
        ROLD=PARF(I1) 
      ELSEIF(IVAR.EQ.12) THEN 
        ROLD=VCKM(I1,I2) 
      ELSEIF(IVAR.EQ.13) THEN 
        IOLD=MDCY(I1,I2) 
      ELSEIF(IVAR.EQ.14) THEN 
        IOLD=MDME(I1,I2) 
      ELSEIF(IVAR.EQ.15) THEN 
        ROLD=BRAT(I1) 
      ELSEIF(IVAR.EQ.16) THEN 
        IOLD=KFDP(I1,I2) 
      ELSEIF(IVAR.EQ.17) THEN 
        CHOLD=CHAF(I1) 
      ELSEIF(IVAR.EQ.18) THEN 
        IOLD=MRLU(I1) 
      ELSEIF(IVAR.EQ.19) THEN 
        ROLD=RRLU(I1) 
      ELSEIF(IVAR.EQ.20) THEN 
        IOLD=MSEL 
      ELSEIF(IVAR.EQ.21) THEN 
        IOLD=MSUB(I1) 
      ELSEIF(IVAR.EQ.22) THEN 
        IOLD=KFIN(I1,I2) 
      ELSEIF(IVAR.EQ.23) THEN 
        ROLD=CKIN(I1) 
      ELSEIF(IVAR.EQ.24) THEN 
        IOLD=MSTP(I1) 
      ELSEIF(IVAR.EQ.25) THEN 
        ROLD=PARP(I1) 
      ELSEIF(IVAR.EQ.26) THEN 
        IOLD=MSTI(I1) 
      ELSEIF(IVAR.EQ.27) THEN 
        ROLD=PARI(I1) 
      ELSEIF(IVAR.EQ.28) THEN 
        IOLD=MINT(I1) 
      ELSEIF(IVAR.EQ.29) THEN 
        ROLD=VINT(I1) 
      ELSEIF(IVAR.EQ.30) THEN 
        IOLD=ISET(I1) 
      ELSEIF(IVAR.EQ.31) THEN 
        IOLD=KFPR(I1,I2) 
      ELSEIF(IVAR.EQ.32) THEN 
        ROLD=COEF(I1,I2) 
      ELSEIF(IVAR.EQ.33) THEN 
        IOLD=ICOL(I1,I2,I3) 
      ELSEIF(IVAR.EQ.34) THEN 
        ROLD=XSFX(I1,I2) 
      ELSEIF(IVAR.EQ.35) THEN 
        IOLD=ISIG(I1,I2) 
      ELSEIF(IVAR.EQ.36) THEN 
        ROLD=SIGH(I1) 
      ELSEIF(IVAR.EQ.37) THEN 
        ROLD=WIDP(I1,I2) 
      ELSEIF(IVAR.EQ.38) THEN 
        ROLD=WIDE(I1,I2) 
      ELSEIF(IVAR.EQ.39) THEN 
        ROLD=WIDS(I1,I2) 
      ELSEIF(IVAR.EQ.40) THEN 
        IOLD=NGEN(I1,I2) 
      ELSEIF(IVAR.EQ.41) THEN 
        ROLD=XSEC(I1,I2) 
      ELSEIF(IVAR.EQ.42) THEN 
        CHOLD2=PROC(I1) 
      ELSEIF(IVAR.EQ.43) THEN 
        ROLD=SIGT(I1,I2,I3) 
      ENDIF 
 
C...Print current value of variable. Loop back. 
      IF(LNAM.GE.LBIT) THEN 
        CHBIT(LNAM:14)=' ' 
        CHBIT(15:60)=' has the value                                ' 
        IF(MSVAR(IVAR,1).EQ.1) THEN 
          WRITE(CHBIT(51:60),'(I10)') IOLD 
        ELSEIF(MSVAR(IVAR,1).EQ.2) THEN 
          WRITE(CHBIT(47:60),'(F14.5)') ROLD 
        ELSEIF(MSVAR(IVAR,1).EQ.3) THEN 
          CHBIT(53:60)=CHOLD 
        ELSE 
          CHBIT(33:60)=CHOLD 
        ENDIF 
        IF(MSTU(13).GE.1) WRITE(MSTU(11),5000) CHBIT(1:60) 
        LLOW=LHIG 
        IF(LLOW.LT.LTOT) GOTO 120 
        RETURN 
      ENDIF 
 
C...Read in new variable value. 
      IF(MSVAR(IVAR,1).EQ.1) THEN 
        CHINI=' ' 
        CHINI(LNAM-LBIT+11:10)=CHBIT(LNAM+1:LBIT) 
        READ(CHINI,'(I10)') INEW 
      ELSEIF(MSVAR(IVAR,1).EQ.2) THEN 
        CHINR=' ' 
        CHINR(LNAM-LBIT+17:16)=CHBIT(LNAM+1:LBIT) 
        READ(CHINR,'(F16.2)') RNEW 
      ELSEIF(MSVAR(IVAR,1).EQ.3) THEN 
        CHNEW=CHBIT(LNAM+1:LBIT)//' ' 
      ELSE 
        CHNEW2=CHBIT(LNAM+1:LBIT)//' ' 
      ENDIF 
 
C...Store new variable value. 
      IF(IVAR.EQ.1) THEN 
        N=INEW 
      ELSEIF(IVAR.EQ.2) THEN 
        K(I1,I2)=INEW 
      ELSEIF(IVAR.EQ.3) THEN 
        P(I1,I2)=RNEW 
      ELSEIF(IVAR.EQ.4) THEN 
        V(I1,I2)=RNEW 
      ELSEIF(IVAR.EQ.5) THEN 
        MSTU(I1)=INEW 
      ELSEIF(IVAR.EQ.6) THEN 
        PARU(I1)=RNEW 
      ELSEIF(IVAR.EQ.7) THEN 
        MSTJ(I1)=INEW 
      ELSEIF(IVAR.EQ.8) THEN 
        PARJ(I1)=RNEW 
      ELSEIF(IVAR.EQ.9) THEN 
        KCHG(I1,I2)=INEW 
      ELSEIF(IVAR.EQ.10) THEN 
        PMAS(I1,I2)=RNEW 
      ELSEIF(IVAR.EQ.11) THEN 
        PARF(I1)=RNEW 
      ELSEIF(IVAR.EQ.12) THEN 
        VCKM(I1,I2)=RNEW 
      ELSEIF(IVAR.EQ.13) THEN 
        MDCY(I1,I2)=INEW 
      ELSEIF(IVAR.EQ.14) THEN 
        MDME(I1,I2)=INEW 
      ELSEIF(IVAR.EQ.15) THEN 
        BRAT(I1)=RNEW 
      ELSEIF(IVAR.EQ.16) THEN 
        KFDP(I1,I2)=INEW 
      ELSEIF(IVAR.EQ.17) THEN 
        CHAF(I1)=CHNEW 
      ELSEIF(IVAR.EQ.18) THEN 
        MRLU(I1)=INEW 
      ELSEIF(IVAR.EQ.19) THEN 
        RRLU(I1)=RNEW 
      ELSEIF(IVAR.EQ.20) THEN 
        MSEL=INEW 
      ELSEIF(IVAR.EQ.21) THEN 
        MSUB(I1)=INEW 
      ELSEIF(IVAR.EQ.22) THEN 
        KFIN(I1,I2)=INEW 
      ELSEIF(IVAR.EQ.23) THEN 
        CKIN(I1)=RNEW 
      ELSEIF(IVAR.EQ.24) THEN 
        MSTP(I1)=INEW 
      ELSEIF(IVAR.EQ.25) THEN 
        PARP(I1)=RNEW 
      ELSEIF(IVAR.EQ.26) THEN 
        MSTI(I1)=INEW 
      ELSEIF(IVAR.EQ.27) THEN 
        PARI(I1)=RNEW 
      ELSEIF(IVAR.EQ.28) THEN 
        MINT(I1)=INEW 
      ELSEIF(IVAR.EQ.29) THEN 
        VINT(I1)=RNEW 
      ELSEIF(IVAR.EQ.30) THEN 
        ISET(I1)=INEW 
      ELSEIF(IVAR.EQ.31) THEN 
        KFPR(I1,I2)=INEW 
      ELSEIF(IVAR.EQ.32) THEN 
        COEF(I1,I2)=RNEW 
      ELSEIF(IVAR.EQ.33) THEN 
        ICOL(I1,I2,I3)=INEW 
      ELSEIF(IVAR.EQ.34) THEN 
        XSFX(I1,I2)=RNEW 
      ELSEIF(IVAR.EQ.35) THEN 
        ISIG(I1,I2)=INEW 
      ELSEIF(IVAR.EQ.36) THEN 
        SIGH(I1)=RNEW 
      ELSEIF(IVAR.EQ.37) THEN 
        WIDP(I1,I2)=RNEW 
      ELSEIF(IVAR.EQ.38) THEN 
        WIDE(I1,I2)=RNEW 
      ELSEIF(IVAR.EQ.39) THEN 
        WIDS(I1,I2)=RNEW 
      ELSEIF(IVAR.EQ.40) THEN 
        NGEN(I1,I2)=INEW 
      ELSEIF(IVAR.EQ.41) THEN 
        XSEC(I1,I2)=RNEW 
      ELSEIF(IVAR.EQ.42) THEN 
        PROC(I1)=CHNEW2 
      ELSEIF(IVAR.EQ.43) THEN 
        SIGT(I1,I2,I3)=RNEW 
      ENDIF 
 
C...Write old and new value. Loop back. 
      CHBIT(LNAM:14)=' ' 
      CHBIT(15:60)=' changed from                to               ' 
      IF(MSVAR(IVAR,1).EQ.1) THEN 
        WRITE(CHBIT(33:42),'(I10)') IOLD 
        WRITE(CHBIT(51:60),'(I10)') INEW 
        IF(MSTU(13).GE.1) WRITE(MSTU(11),5000) CHBIT(1:60) 
      ELSEIF(MSVAR(IVAR,1).EQ.2) THEN 
        WRITE(CHBIT(29:42),'(F14.5)') ROLD 
        WRITE(CHBIT(47:60),'(F14.5)') RNEW 
        IF(MSTU(13).GE.1) WRITE(MSTU(11),5000) CHBIT(1:60) 
      ELSEIF(MSVAR(IVAR,1).EQ.3) THEN 
        CHBIT(35:42)=CHOLD 
        CHBIT(53:60)=CHNEW 
        IF(MSTU(13).GE.1) WRITE(MSTU(11),5000) CHBIT(1:60) 
      ELSE 
        CHBIT(15:88)=' changed from '//CHOLD2//' to '//CHNEW2 
        IF(MSTU(13).GE.1) WRITE(MSTU(11),5100) CHBIT(1:88) 
      ENDIF 
      LLOW=LHIG 
      IF(LLOW.LT.LTOT) GOTO 120 
 
C...Format statement for output on unit MSTU(11) (by default 6). 
 5000 FORMAT(5X,A60) 
 5100 FORMAT(5X,A88) 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUEXEC 
 
C...Purpose: to administrate the fragmentation and decay chain. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/ 
      DIMENSION PS(2,6) 
 
C...Initialize and reset. 
      MSTU(24)=0 
      IF(MSTU(12).GE.1) CALL LULIST(0) 
      MSTU(31)=MSTU(31)+1 
      MSTU(1)=0 
      MSTU(2)=0 
      MSTU(3)=0 
      IF(MSTU(17).LE.0) MSTU(90)=0 
      MCONS=1 
 
C...Sum up momentum, energy and charge for starting entries. 
      NSAV=N 
      DO 110 I=1,2 
      DO 100 J=1,6 
      PS(I,J)=0. 
  100 CONTINUE 
  110 CONTINUE 
      DO 130 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 130 
      DO 120 J=1,4 
      PS(1,J)=PS(1,J)+P(I,J) 
  120 CONTINUE 
      PS(1,6)=PS(1,6)+LUCHGE(K(I,2)) 
  130 CONTINUE 
      PARU(21)=PS(1,4) 
 
C...Prepare system for subsequent fragmentation/decay. 
      CALL LUPREP(0) 
 
C...Loop through jet fragmentation and particle decays. 
      MBE=0 
  140 MBE=MBE+1 
      IP=0 
  150 IP=IP+1 
      KC=0 
      IF(K(IP,1).GT.0.AND.K(IP,1).LE.10) KC=LUCOMP(K(IP,2)) 
      IF(KC.EQ.0) THEN 
 
C...Particle decay if unstable and allowed. Save long-lived particle 
C...decays until second pass after Bose-Einstein effects. 
      ELSEIF(KCHG(KC,2).EQ.0) THEN 
        IF(MSTJ(21).GE.1.AND.MDCY(KC,1).GE.1.AND.(MSTJ(51).LE.0.OR.MBE 
     &  .EQ.2.OR.PMAS(KC,2).GE.PARJ(91).OR.IABS(K(IP,2)).EQ.311)) 
     &  CALL LUDECY(IP) 
 
C...Decay products may develop a shower. 
        IF(MSTJ(92).GT.0) THEN 
          IP1=MSTJ(92) 
          QMAX=SQRT(MAX(0.,(P(IP1,4)+P(IP1+1,4))**2-(P(IP1,1)+P(IP1+1, 
     &    1))**2-(P(IP1,2)+P(IP1+1,2))**2-(P(IP1,3)+P(IP1+1,3))**2)) 
          CALL LUSHOW(IP1,IP1+1,QMAX) 
          CALL LUPREP(IP1) 
          MSTJ(92)=0 
        ELSEIF(MSTJ(92).LT.0) THEN 
          IP1=-MSTJ(92) 
          CALL LUSHOW(IP1,-3,P(IP,5)) 
          CALL LUPREP(IP1) 
          MSTJ(92)=0 
        ENDIF 
 
C...Jet fragmentation: string or independent fragmentation. 
      ELSEIF(K(IP,1).EQ.1.OR.K(IP,1).EQ.2) THEN 
        MFRAG=MSTJ(1) 
        IF(MFRAG.GE.1.AND.K(IP,1).EQ.1) MFRAG=2 
        IF(MSTJ(21).GE.2.AND.K(IP,1).EQ.2.AND.N.GT.IP) THEN 
          IF(K(IP+1,1).EQ.1.AND.K(IP+1,3).EQ.K(IP,3).AND. 
     &    K(IP,3).GT.0.AND.K(IP,3).LT.IP) THEN 
            IF(KCHG(LUCOMP(K(K(IP,3),2)),2).EQ.0) MFRAG=MIN(1,MFRAG) 
          ENDIF 
        ENDIF 
        IF(MFRAG.EQ.1) CALL LUSTRF(IP) 
        IF(MFRAG.EQ.2) CALL LUINDF(IP) 
        IF(MFRAG.EQ.2.AND.K(IP,1).EQ.1) MCONS=0 
        IF(MFRAG.EQ.2.AND.(MSTJ(3).LE.0.OR.MOD(MSTJ(3),5).EQ.0)) MCONS=0 
      ENDIF 
 
C...Loop back if enough space left in LUJETS and no error abort. 
      IF(MSTU(24).NE.0.AND.MSTU(21).GE.2) THEN 
      ELSEIF(IP.LT.N.AND.N.LT.MSTU(4)-20-MSTU(32)) THEN 
        GOTO 150 
      ELSEIF(IP.LT.N) THEN 
        CALL LUERRM(11,'(LUEXEC:) no more memory left in LUJETS') 
      ENDIF 
 
C...Include simple Bose-Einstein effect parametrization if desired. 
      IF(MBE.EQ.1.AND.MSTJ(51).GE.1) THEN 
        CALL LUBOEI(NSAV) 
        GOTO 140 
      ENDIF 
 
C...Check that momentum, energy and charge were conserved. 
      DO 170 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 170 
      DO 160 J=1,4 
      PS(2,J)=PS(2,J)+P(I,J) 
  160 CONTINUE 
      PS(2,6)=PS(2,6)+LUCHGE(K(I,2)) 
  170 CONTINUE 
      PDEV=(ABS(PS(2,1)-PS(1,1))+ABS(PS(2,2)-PS(1,2))+ABS(PS(2,3)- 
     &PS(1,3))+ABS(PS(2,4)-PS(1,4)))/(1.+ABS(PS(2,4))+ABS(PS(1,4))) 
      IF(MCONS.EQ.1.AND.PDEV.GT.PARU(11)) CALL LUERRM(15, 
     &'(LUEXEC:) four-momentum was not conserved') 
      IF(MCONS.EQ.1.AND.ABS(PS(2,6)-PS(1,6)).GT.0.1) CALL LUERRM(15, 
     &'(LUEXEC:) charge was not conserved') 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUPREP(IP) 
 
C...Purpose: to rearrange partons along strings, to allow small systems 
C...to collapse into one or two particles and to check flavours. 
      IMPLICIT DOUBLE PRECISION(D) 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/ 
      DIMENSION DPS(5),DPC(5),UE(3) 
 
C...Rearrange parton shower product listing along strings: begin loop. 
      I1=N 
      DO 130 MQGST=1,2 
      DO 120 I=MAX(1,IP),N 
      IF(K(I,1).NE.3) GOTO 120 
      KC=LUCOMP(K(I,2)) 
      IF(KC.EQ.0) GOTO 120 
      KQ=KCHG(KC,2) 
      IF(KQ.EQ.0.OR.(MQGST.EQ.1.AND.KQ.EQ.2)) GOTO 120 
 
C...Pick up loose string end. 
      KCS=4 
      IF(KQ*ISIGN(1,K(I,2)).LT.0) KCS=5 
      IA=I 
      NSTP=0 
  100 NSTP=NSTP+1 
      IF(NSTP.GT.4*N) THEN 
        CALL LUERRM(14,'(LUPREP:) caught in infinite loop') 
        RETURN 
      ENDIF 
 
C...Copy undecayed parton. 
      IF(K(IA,1).EQ.3) THEN 
        IF(I1.GE.MSTU(4)-MSTU(32)-5) THEN 
          CALL LUERRM(11,'(LUPREP:) no more memory left in LUJETS') 
          RETURN 
        ENDIF 
        I1=I1+1 
        K(I1,1)=2 
        IF(NSTP.GE.2.AND.IABS(K(IA,2)).NE.21) K(I1,1)=1 
        K(I1,2)=K(IA,2) 
        K(I1,3)=IA 
        K(I1,4)=0 
        K(I1,5)=0 
        DO 110 J=1,5 
        P(I1,J)=P(IA,J) 
        V(I1,J)=V(IA,J) 
  110   CONTINUE 
        K(IA,1)=K(IA,1)+10 
        IF(K(I1,1).EQ.1) GOTO 120 
      ENDIF 
 
C...Go to next parton in colour space. 
      IB=IA 
      IF(MOD(K(IB,KCS)/MSTU(5)**2,2).EQ.0.AND.MOD(K(IB,KCS),MSTU(5)) 
     &.NE.0) THEN 
        IA=MOD(K(IB,KCS),MSTU(5)) 
        K(IB,KCS)=K(IB,KCS)+MSTU(5)**2 
        MREV=0 
      ELSE 
        IF(K(IB,KCS).GE.2*MSTU(5)**2.OR.MOD(K(IB,KCS)/MSTU(5),MSTU(5)) 
     &  .EQ.0) KCS=9-KCS 
        IA=MOD(K(IB,KCS)/MSTU(5),MSTU(5)) 
        K(IB,KCS)=K(IB,KCS)+2*MSTU(5)**2 
        MREV=1 
      ENDIF 
      IF(IA.LE.0.OR.IA.GT.N) THEN 
        CALL LUERRM(12,'(LUPREP:) colour rearrangement failed') 
        RETURN 
      ENDIF 
      IF(MOD(K(IA,4)/MSTU(5),MSTU(5)).EQ.IB.OR.MOD(K(IA,5)/MSTU(5), 
     &MSTU(5)).EQ.IB) THEN 
        IF(MREV.EQ.1) KCS=9-KCS 
        IF(MOD(K(IA,KCS)/MSTU(5),MSTU(5)).NE.IB) KCS=9-KCS 
        K(IA,KCS)=K(IA,KCS)+2*MSTU(5)**2 
      ELSE 
        IF(MREV.EQ.0) KCS=9-KCS 
        IF(MOD(K(IA,KCS),MSTU(5)).NE.IB) KCS=9-KCS 
        K(IA,KCS)=K(IA,KCS)+MSTU(5)**2 
      ENDIF 
      IF(IA.NE.I) GOTO 100 
      K(I1,1)=1 
  120 CONTINUE 
  130 CONTINUE 
      N=I1 
      IF(MSTJ(14).LT.0) RETURN 
 
C...Find lowest-mass colour singlet jet system, OK if above threshold. 
      IF(MSTJ(14).EQ.0) GOTO 320 
      NS=N 
  140 NSIN=N-NS 
      PDM=1.+PARJ(32) 
      IC=0 
      DO 190 I=MAX(1,IP),NS 
      IF(K(I,1).NE.1.AND.K(I,1).NE.2) THEN 
      ELSEIF(K(I,1).EQ.2.AND.IC.EQ.0) THEN 
        NSIN=NSIN+1 
        IC=I 
        DO 150 J=1,4 
        DPS(J)=P(I,J) 
  150   CONTINUE 
        MSTJ(93)=1 
        DPS(5)=ULMASS(K(I,2)) 
      ELSEIF(K(I,1).EQ.2) THEN 
        DO 160 J=1,4 
        DPS(J)=DPS(J)+P(I,J) 
  160   CONTINUE 
      ELSEIF(IC.NE.0.AND.KCHG(LUCOMP(K(I,2)),2).NE.0) THEN 
        DO 170 J=1,4 
        DPS(J)=DPS(J)+P(I,J) 
  170   CONTINUE 
        MSTJ(93)=1 
        DPS(5)=DPS(5)+ULMASS(K(I,2)) 
        PD=SQRT(MAX(0D0,DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2))-DPS(5) 
        IF(PD.LT.PDM) THEN 
          PDM=PD 
          DO 180 J=1,5 
          DPC(J)=DPS(J) 
  180     CONTINUE 
          IC1=IC 
          IC2=I 
        ENDIF 
        IC=0 
      ELSE 
        NSIN=NSIN+1 
      ENDIF 
  190 CONTINUE 
      IF(PDM.GE.PARJ(32)) GOTO 320 
 
C...Fill small-mass system as cluster. 
      NSAV=N 
      PECM=SQRT(MAX(0D0,DPC(4)**2-DPC(1)**2-DPC(2)**2-DPC(3)**2)) 
      K(N+1,1)=11 
      K(N+1,2)=91 
      K(N+1,3)=IC1 
      K(N+1,4)=N+2 
      K(N+1,5)=N+3 
      P(N+1,1)=DPC(1) 
      P(N+1,2)=DPC(2) 
      P(N+1,3)=DPC(3) 
      P(N+1,4)=DPC(4) 
      P(N+1,5)=PECM 
 
C...Form two particles from flavours of lowest-mass system, if feasible. 
      K(N+2,1)=1 
      K(N+3,1)=1 
      IF(MSTU(16).NE.2) THEN 
        K(N+2,3)=N+1 
        K(N+3,3)=N+1 
      ELSE 
        K(N+2,3)=IC1 
        K(N+3,3)=IC2 
      ENDIF 
      K(N+2,4)=0 
      K(N+3,4)=0 
      K(N+2,5)=0 
      K(N+3,5)=0 
      IF(IABS(K(IC1,2)).NE.21) THEN 
        KC1=LUCOMP(K(IC1,2)) 
        KC2=LUCOMP(K(IC2,2)) 
        IF(KC1.EQ.0.OR.KC2.EQ.0) GOTO 320 
        KQ1=KCHG(KC1,2)*ISIGN(1,K(IC1,2)) 
        KQ2=KCHG(KC2,2)*ISIGN(1,K(IC2,2)) 
        IF(KQ1+KQ2.NE.0) GOTO 320 
  200   CALL LUKFDI(K(IC1,2),0,KFLN,K(N+2,2)) 
        CALL LUKFDI(K(IC2,2),-KFLN,KFLDMP,K(N+3,2)) 
        IF(K(N+2,2).EQ.0.OR.K(N+3,2).EQ.0) GOTO 200 
      ELSE 
        IF(IABS(K(IC2,2)).NE.21) GOTO 320 
  210   CALL LUKFDI(1+INT((2.+PARJ(2))*RLU(0)),0,KFLN,KFDMP) 
        CALL LUKFDI(KFLN,0,KFLM,K(N+2,2)) 
        CALL LUKFDI(-KFLN,-KFLM,KFLDMP,K(N+3,2)) 
        IF(K(N+2,2).EQ.0.OR.K(N+3,2).EQ.0) GOTO 210 
      ENDIF 
      P(N+2,5)=ULMASS(K(N+2,2)) 
      P(N+3,5)=ULMASS(K(N+3,2)) 
      IF(P(N+2,5)+P(N+3,5)+PARJ(64).GE.PECM.AND.NSIN.EQ.1) GOTO 320 
      IF(P(N+2,5)+P(N+3,5)+PARJ(64).GE.PECM) GOTO 260 
 
C...Perform two-particle decay of jet system, if possible. 
      IF(PECM.GE.0.02*DPC(4)) THEN 
        PA=SQRT((PECM**2-(P(N+2,5)+P(N+3,5))**2)*(PECM**2- 
     &  (P(N+2,5)-P(N+3,5))**2))/(2.*PECM) 
        UE(3)=2.*RLU(0)-1. 
        PHI=PARU(2)*RLU(0) 
        UE(1)=SQRT(1.-UE(3)**2)*COS(PHI) 
        UE(2)=SQRT(1.-UE(3)**2)*SIN(PHI) 
        DO 220 J=1,3 
        P(N+2,J)=PA*UE(J) 
        P(N+3,J)=-PA*UE(J) 
  220   CONTINUE 
        P(N+2,4)=SQRT(PA**2+P(N+2,5)**2) 
        P(N+3,4)=SQRT(PA**2+P(N+3,5)**2) 
        MSTU(33)=1 
        CALL LUDBRB(N+2,N+3,0.,0.,DPC(1)/DPC(4),DPC(2)/DPC(4), 
     &  DPC(3)/DPC(4)) 
      ELSE 
        NP=0 
        DO 230 I=IC1,IC2 
        IF(K(I,1).EQ.1.OR.K(I,1).EQ.2) NP=NP+1 
  230   CONTINUE 
        HA=P(IC1,4)*P(IC2,4)-P(IC1,1)*P(IC2,1)-P(IC1,2)*P(IC2,2)- 
     &  P(IC1,3)*P(IC2,3) 
        IF(NP.GE.3.OR.HA.LE.1.25*P(IC1,5)*P(IC2,5)) GOTO 260 
        HD1=0.5*(P(N+2,5)**2-P(IC1,5)**2) 
        HD2=0.5*(P(N+3,5)**2-P(IC2,5)**2) 
        HR=SQRT(MAX(0.,((HA-HD1-HD2)**2-(P(N+2,5)*P(N+3,5))**2)/ 
     &  (HA**2-(P(IC1,5)*P(IC2,5))**2)))-1. 
        HC=P(IC1,5)**2+2.*HA+P(IC2,5)**2 
        HK1=((P(IC2,5)**2+HA)*HR+HD1-HD2)/HC 
        HK2=((P(IC1,5)**2+HA)*HR+HD2-HD1)/HC 
        DO 240 J=1,4 
        P(N+2,J)=(1.+HK1)*P(IC1,J)-HK2*P(IC2,J) 
        P(N+3,J)=(1.+HK2)*P(IC2,J)-HK1*P(IC1,J) 
  240   CONTINUE 
      ENDIF 
      DO 250 J=1,4 
      V(N+1,J)=V(IC1,J) 
      V(N+2,J)=V(IC1,J) 
      V(N+3,J)=V(IC2,J) 
  250 CONTINUE 
      V(N+1,5)=0. 
      V(N+2,5)=0. 
      V(N+3,5)=0. 
      N=N+3 
      GOTO 300 
 
C...Else form one particle from the flavours available, if possible. 
  260 K(N+1,5)=N+2 
      IF(IABS(K(IC1,2)).GT.100.AND.IABS(K(IC2,2)).GT.100) THEN 
        GOTO 320 
      ELSEIF(IABS(K(IC1,2)).NE.21) THEN 
        CALL LUKFDI(K(IC1,2),K(IC2,2),KFLDMP,K(N+2,2)) 
      ELSE 
        KFLN=1+INT((2.+PARJ(2))*RLU(0)) 
        CALL LUKFDI(KFLN,-KFLN,KFLDMP,K(N+2,2)) 
      ENDIF 
      IF(K(N+2,2).EQ.0) GOTO 260 
      P(N+2,5)=ULMASS(K(N+2,2)) 
 
C...Find parton/particle which combines to largest extra mass. 
      IR=0 
      HA=0. 
      HSM=0. 
      DO 280 MCOMB=1,3 
      IF(IR.NE.0) GOTO 280 
      DO 270 I=MAX(1,IP),N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10.OR.(I.GE.IC1.AND.I.LE.IC2 
     &.AND.K(I,1).GE.1.AND.K(I,1).LE.2)) GOTO 270 
      IF(MCOMB.EQ.1) KCI=LUCOMP(K(I,2)) 
      IF(MCOMB.EQ.1.AND.KCI.EQ.0) GOTO 270 
      IF(MCOMB.EQ.1.AND.KCHG(KCI,2).EQ.0.AND.I.LE.NS) GOTO 270 
      IF(MCOMB.EQ.2.AND.IABS(K(I,2)).GT.10.AND.IABS(K(I,2)).LE.100) 
     &GOTO 270 
      HCR=DPC(4)*P(I,4)-DPC(1)*P(I,1)-DPC(2)*P(I,2)-DPC(3)*P(I,3) 
      HSR=2.*HCR+PECM**2-P(N+2,5)**2-2.*P(N+2,5)*P(I,5) 
      IF(HSR.GT.HSM) THEN 
        IR=I 
        HA=HCR 
        HSM=HSR 
      ENDIF 
  270 CONTINUE 
  280 CONTINUE 
 
C...Shuffle energy and momentum to put new particle on mass shell. 
      IF(IR.NE.0) THEN 
        HB=PECM**2+HA 
        HC=P(N+2,5)**2+HA 
        HD=P(IR,5)**2+HA 
        HK2=0.5*(HB*SQRT(MAX(0.,((HB+HC)**2-4.*(HB+HD)*P(N+2,5)**2)/ 
     &  (HA**2-(PECM*P(IR,5))**2)))-(HB+HC))/(HB+HD) 
        HK1=(0.5*(P(N+2,5)**2-PECM**2)+HD*HK2)/HB 
        DO 290 J=1,4 
        P(N+2,J)=(1.+HK1)*DPC(J)-HK2*P(IR,J) 
        P(IR,J)=(1.+HK2)*P(IR,J)-HK1*DPC(J) 
        V(N+1,J)=V(IC1,J) 
        V(N+2,J)=V(IC1,J) 
  290   CONTINUE 
        V(N+1,5)=0. 
        V(N+2,5)=0. 
        N=N+2 
      ELSE 
        CALL LUERRM(3,'(LUPREP:) no match for collapsing cluster') 
        RETURN 
      ENDIF 
 
C...Mark collapsed system and store daughter pointers. Iterate. 
  300 DO 310 I=IC1,IC2 
      IF((K(I,1).EQ.1.OR.K(I,1).EQ.2).AND.KCHG(LUCOMP(K(I,2)),2).NE.0) 
     &THEN 
        K(I,1)=K(I,1)+10 
        IF(MSTU(16).NE.2) THEN 
          K(I,4)=NSAV+1 
          K(I,5)=NSAV+1 
        ELSE 
          K(I,4)=NSAV+2 
          K(I,5)=N 
        ENDIF 
      ENDIF 
  310 CONTINUE 
      IF(N.LT.MSTU(4)-MSTU(32)-5) GOTO 140 
 
C...Check flavours and invariant masses in parton systems. 
  320 NP=0 
      KFN=0 
      KQS=0 
      DO 330 J=1,5 
      DPS(J)=0. 
  330 CONTINUE 
      DO 360 I=MAX(1,IP),N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 360 
      KC=LUCOMP(K(I,2)) 
      IF(KC.EQ.0) GOTO 360 
      KQ=KCHG(KC,2)*ISIGN(1,K(I,2)) 
      IF(KQ.EQ.0) GOTO 360 
      NP=NP+1 
      IF(KQ.NE.2) THEN 
        KFN=KFN+1 
        KQS=KQS+KQ 
        MSTJ(93)=1 
        DPS(5)=DPS(5)+ULMASS(K(I,2)) 
      ENDIF 
      DO 340 J=1,4 
      DPS(J)=DPS(J)+P(I,J) 
  340 CONTINUE 
      IF(K(I,1).EQ.1) THEN 
        IF(NP.NE.1.AND.(KFN.EQ.1.OR.KFN.GE.3.OR.KQS.NE.0)) CALL 
     &  LUERRM(2,'(LUPREP:) unphysical flavour combination') 
        IF(NP.NE.1.AND.DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2.LT. 
     &  (0.9*PARJ(32)+DPS(5))**2) CALL LUERRM(3, 
     &  '(LUPREP:) too small mass in jet system') 
        NP=0 
        KFN=0 
        KQS=0 
        DO 350 J=1,5 
        DPS(J)=0. 
  350   CONTINUE 
      ENDIF 
  360 CONTINUE 
 
      RETURN 
      END 
 
C********************************************************************* 
 
C********************************************************************* 
 
      SUBROUTINE LUINDF(IP) 
 
C...Purpose: to handle the fragmentation of a jet system (or a single 
C...jet) according to independent fragmentation models. 
      IMPLICIT DOUBLE PRECISION(D) 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION DPS(5),PSI(4),NFI(3),NFL(3),IFET(3),KFLF(3), 
     &KFLO(2),PXO(2),PYO(2),WO(2) 
 
C...Reset counters. Identify parton system and take copy. Check flavour. 
      NSAV=N 
      MSTU90=MSTU(90) 
      NJET=0 
      KQSUM=0 
      DO 100 J=1,5 
      DPS(J)=0. 
  100 CONTINUE 
      I=IP-1 
  110 I=I+1 
      IF(I.GT.MIN(N,MSTU(4)-MSTU(32))) THEN 
        CALL LUERRM(12,'(LUINDF:) failed to reconstruct jet system') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      IF(K(I,1).NE.1.AND.K(I,1).NE.2) GOTO 110 
      KC=LUCOMP(K(I,2)) 
      IF(KC.EQ.0) GOTO 110 
      KQ=KCHG(KC,2)*ISIGN(1,K(I,2)) 
      IF(KQ.EQ.0) GOTO 110 
      NJET=NJET+1 
      IF(KQ.NE.2) KQSUM=KQSUM+KQ 
      DO 120 J=1,5 
      K(NSAV+NJET,J)=K(I,J) 
      P(NSAV+NJET,J)=P(I,J) 
      DPS(J)=DPS(J)+P(I,J) 
  120 CONTINUE 
      K(NSAV+NJET,3)=I 
      IF(K(I,1).EQ.2.OR.(MSTJ(3).LE.5.AND.N.GT.I.AND. 
     &K(I+1,1).EQ.2)) GOTO 110 
      IF(NJET.NE.1.AND.KQSUM.NE.0) THEN 
        CALL LUERRM(12,'(LUINDF:) unphysical flavour combination') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
 
C...Boost copied system to CM frame. Find CM energy and sum flavours. 
      IF(NJET.NE.1) THEN 
        MSTU(33)=1 
        CALL LUDBRB(NSAV+1,NSAV+NJET,0.,0.,-DPS(1)/DPS(4), 
     &  -DPS(2)/DPS(4),-DPS(3)/DPS(4)) 
      ENDIF 
      PECM=0. 
      DO 130 J=1,3 
      NFI(J)=0 
  130 CONTINUE 
      DO 140 I=NSAV+1,NSAV+NJET 
      PECM=PECM+P(I,4) 
      KFA=IABS(K(I,2)) 
      IF(KFA.LE.3) THEN 
        NFI(KFA)=NFI(KFA)+ISIGN(1,K(I,2)) 
      ELSEIF(KFA.GT.1000) THEN 
        KFLA=MOD(KFA/1000,10) 
        KFLB=MOD(KFA/100,10) 
        IF(KFLA.LE.3) NFI(KFLA)=NFI(KFLA)+ISIGN(1,K(I,2)) 
        IF(KFLB.LE.3) NFI(KFLB)=NFI(KFLB)+ISIGN(1,K(I,2)) 
      ENDIF 
  140 CONTINUE 
 
C...Loop over attempts made. Reset counters. 
      NTRY=0 
  150 NTRY=NTRY+1 
      IF(NTRY.GT.200) THEN 
        CALL LUERRM(14,'(LUINDF:) caught in infinite loop') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      N=NSAV+NJET 
      MSTU(90)=MSTU90 
      DO 160 J=1,3 
      NFL(J)=NFI(J) 
      IFET(J)=0 
      KFLF(J)=0 
  160 CONTINUE 
 
C...Loop over jets to be fragmented. 
      DO 230 IP1=NSAV+1,NSAV+NJET 
      MSTJ(91)=0 
      NSAV1=N 
      MSTU91=MSTU(90) 
 
C...Initial flavour and momentum values. Jet along +z axis. 
      KFLH=IABS(K(IP1,2)) 
      IF(KFLH.GT.10) KFLH=MOD(KFLH/1000,10) 
      KFLO(2)=0 
      WF=P(IP1,4)+SQRT(P(IP1,1)**2+P(IP1,2)**2+P(IP1,3)**2) 
 
C...Initial values for quark or diquark jet. 
  170 IF(IABS(K(IP1,2)).NE.21) THEN 
        NSTR=1 
        KFLO(1)=K(IP1,2) 
        CALL LUPTDI(0,PXO(1),PYO(1)) 
        WO(1)=WF 
 
C...Initial values for gluon treated like random quark jet. 
      ELSEIF(MSTJ(2).LE.2) THEN 
        NSTR=1 
        IF(MSTJ(2).EQ.2) MSTJ(91)=1 
        KFLO(1)=INT(1.+(2.+PARJ(2))*RLU(0))*(-1)**INT(RLU(0)+0.5) 
        CALL LUPTDI(0,PXO(1),PYO(1)) 
        WO(1)=WF 
 
C...Initial values for gluon treated like quark-antiquark jet pair, 
C...sharing energy according to Altarelli-Parisi splitting function. 
      ELSE 
        NSTR=2 
        IF(MSTJ(2).EQ.4) MSTJ(91)=1 
        KFLO(1)=INT(1.+(2.+PARJ(2))*RLU(0))*(-1)**INT(RLU(0)+0.5) 
        KFLO(2)=-KFLO(1) 
        CALL LUPTDI(0,PXO(1),PYO(1)) 
        PXO(2)=-PXO(1) 
        PYO(2)=-PYO(1) 
        WO(1)=WF*RLU(0)**(1./3.) 
        WO(2)=WF-WO(1) 
      ENDIF 
 
C...Initial values for rank, flavour, pT and W+. 
      DO 220 ISTR=1,NSTR 
  180 I=N 
      MSTU(90)=MSTU91 
      IRANK=0 
      KFL1=KFLO(ISTR) 
      PX1=PXO(ISTR) 
      PY1=PYO(ISTR) 
      W=WO(ISTR) 
 
C...New hadron. Generate flavour and hadron species. 
  190 I=I+1 
      IF(I.GE.MSTU(4)-MSTU(32)-NJET-5) THEN 
        CALL LUERRM(11,'(LUINDF:) no more memory left in LUJETS') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      IRANK=IRANK+1 
      K(I,1)=1 
      K(I,3)=IP1 
      K(I,4)=0 
      K(I,5)=0 
  200 CALL LUKFDI(KFL1,0,KFL2,K(I,2)) 
      IF(K(I,2).EQ.0) GOTO 180 
      IF(MSTJ(12).GE.3.AND.IRANK.EQ.1.AND.IABS(KFL1).LE.10.AND. 
     &IABS(KFL2).GT.10) THEN 
        IF(RLU(0).GT.PARJ(19)) GOTO 200 
      ENDIF 
 
C...Find hadron mass. Generate four-momentum. 
      P(I,5)=ULMASS(K(I,2)) 
      CALL LUPTDI(KFL1,PX2,PY2) 
      P(I,1)=PX1+PX2 
      P(I,2)=PY1+PY2 
      PR=P(I,5)**2+P(I,1)**2+P(I,2)**2 
      CALL LUZDIS(KFL1,KFL2,PR,Z) 
      MZSAV=0 
      IF(IABS(KFL1).GE.4.AND.IABS(KFL1).LE.8.AND.MSTU(90).LT.8) THEN 
        MZSAV=1 
        MSTU(90)=MSTU(90)+1 
        MSTU(90+MSTU(90))=I 
        PARU(90+MSTU(90))=Z 
      ENDIF 
      P(I,3)=0.5*(Z*W-PR/MAX(1E-4,Z*W)) 
      P(I,4)=0.5*(Z*W+PR/MAX(1E-4,Z*W)) 
      IF(MSTJ(3).GE.1.AND.IRANK.EQ.1.AND.KFLH.GE.4.AND. 
     &P(I,3).LE.0.001) THEN 
        IF(W.GE.P(I,5)+0.5*PARJ(32)) GOTO 180 
        P(I,3)=0.0001 
        P(I,4)=SQRT(PR) 
        Z=P(I,4)/W 
      ENDIF 
 
C...Remaining flavour and momentum. 
      KFL1=-KFL2 
      PX1=-PX2 
      PY1=-PY2 
      W=(1.-Z)*W 
      DO 210 J=1,5 
      V(I,J)=0. 
  210 CONTINUE 
 
C...Check if pL acceptable. Go back for new hadron if enough energy. 
      IF(MSTJ(3).GE.0.AND.P(I,3).LT.0.) THEN 
        I=I-1 
        IF(MZSAV.EQ.1) MSTU(90)=MSTU(90)-1 
      ENDIF 
      IF(W.GT.PARJ(31)) GOTO 190 
      N=I 
  220 CONTINUE 
      IF(MOD(MSTJ(3),5).EQ.4.AND.N.EQ.NSAV1) WF=WF+0.1*PARJ(32) 
      IF(MOD(MSTJ(3),5).EQ.4.AND.N.EQ.NSAV1) GOTO 170 
 
C...Rotate jet to new direction. 
      THE=ULANGL(P(IP1,3),SQRT(P(IP1,1)**2+P(IP1,2)**2)) 
      PHI=ULANGL(P(IP1,1),P(IP1,2)) 
      MSTU(33)=1 
      CALL LUDBRB(NSAV1+1,N,THE,PHI,0D0,0D0,0D0) 
      K(K(IP1,3),4)=NSAV1+1 
      K(K(IP1,3),5)=N 
 
C...End of jet generation loop. Skip conservation in some cases. 
  230 CONTINUE 
      IF(NJET.EQ.1.OR.MSTJ(3).LE.0) GOTO 490 
      IF(MOD(MSTJ(3),5).NE.0.AND.N-NSAV-NJET.LT.2) GOTO 150 
 
C...Subtract off produced hadron flavours, finished if zero. 
      DO 240 I=NSAV+NJET+1,N 
      KFA=IABS(K(I,2)) 
      KFLA=MOD(KFA/1000,10) 
      KFLB=MOD(KFA/100,10) 
      KFLC=MOD(KFA/10,10) 
      IF(KFLA.EQ.0) THEN 
        IF(KFLB.LE.3) NFL(KFLB)=NFL(KFLB)-ISIGN(1,K(I,2))*(-1)**KFLB 
        IF(KFLC.LE.3) NFL(KFLC)=NFL(KFLC)+ISIGN(1,K(I,2))*(-1)**KFLB 
      ELSE 
        IF(KFLA.LE.3) NFL(KFLA)=NFL(KFLA)-ISIGN(1,K(I,2)) 
        IF(KFLB.LE.3) NFL(KFLB)=NFL(KFLB)-ISIGN(1,K(I,2)) 
        IF(KFLC.LE.3) NFL(KFLC)=NFL(KFLC)-ISIGN(1,K(I,2)) 
      ENDIF 
  240 CONTINUE 
      NREQ=(IABS(NFL(1))+IABS(NFL(2))+IABS(NFL(3))-IABS(NFL(1)+ 
     &NFL(2)+NFL(3)))/2+IABS(NFL(1)+NFL(2)+NFL(3))/3 
      IF(NREQ.EQ.0) GOTO 320 
 
C...Take away flavour of low-momentum particles until enough freedom. 
      NREM=0 
  250 IREM=0 
      P2MIN=PECM**2 
      DO 260 I=NSAV+NJET+1,N 
      P2=P(I,1)**2+P(I,2)**2+P(I,3)**2 
      IF(K(I,1).EQ.1.AND.P2.LT.P2MIN) IREM=I 
      IF(K(I,1).EQ.1.AND.P2.LT.P2MIN) P2MIN=P2 
  260 CONTINUE 
      IF(IREM.EQ.0) GOTO 150 
      K(IREM,1)=7 
      KFA=IABS(K(IREM,2)) 
      KFLA=MOD(KFA/1000,10) 
      KFLB=MOD(KFA/100,10) 
      KFLC=MOD(KFA/10,10) 
      IF(KFLA.GE.4.OR.KFLB.GE.4) K(IREM,1)=8 
      IF(K(IREM,1).EQ.8) GOTO 250 
      IF(KFLA.EQ.0) THEN 
        ISGN=ISIGN(1,K(IREM,2))*(-1)**KFLB 
        IF(KFLB.LE.3) NFL(KFLB)=NFL(KFLB)+ISGN 
        IF(KFLC.LE.3) NFL(KFLC)=NFL(KFLC)-ISGN 
      ELSE 
        IF(KFLA.LE.3) NFL(KFLA)=NFL(KFLA)+ISIGN(1,K(IREM,2)) 
        IF(KFLB.LE.3) NFL(KFLB)=NFL(KFLB)+ISIGN(1,K(IREM,2)) 
        IF(KFLC.LE.3) NFL(KFLC)=NFL(KFLC)+ISIGN(1,K(IREM,2)) 
      ENDIF 
      NREM=NREM+1 
      NREQ=(IABS(NFL(1))+IABS(NFL(2))+IABS(NFL(3))-IABS(NFL(1)+ 
     &NFL(2)+NFL(3)))/2+IABS(NFL(1)+NFL(2)+NFL(3))/3 
      IF(NREQ.GT.NREM) GOTO 250 
      DO 270 I=NSAV+NJET+1,N 
      IF(K(I,1).EQ.8) K(I,1)=1 
  270 CONTINUE 
 
C...Find combination of existing and new flavours for hadron. 
  280 NFET=2 
      IF(NFL(1)+NFL(2)+NFL(3).NE.0) NFET=3 
      IF(NREQ.LT.NREM) NFET=1 
      IF(IABS(NFL(1))+IABS(NFL(2))+IABS(NFL(3)).EQ.0) NFET=0 
      DO 290 J=1,NFET 
      IFET(J)=1+(IABS(NFL(1))+IABS(NFL(2))+IABS(NFL(3)))*RLU(0) 
      KFLF(J)=ISIGN(1,NFL(1)) 
      IF(IFET(J).GT.IABS(NFL(1))) KFLF(J)=ISIGN(2,NFL(2)) 
      IF(IFET(J).GT.IABS(NFL(1))+IABS(NFL(2))) KFLF(J)=ISIGN(3,NFL(3)) 
  290 CONTINUE 
      IF(NFET.EQ.2.AND.(IFET(1).EQ.IFET(2).OR.KFLF(1)*KFLF(2).GT.0)) 
     &GOTO 280 
      IF(NFET.EQ.3.AND.(IFET(1).EQ.IFET(2).OR.IFET(1).EQ.IFET(3).OR. 
     &IFET(2).EQ.IFET(3).OR.KFLF(1)*KFLF(2).LT.0.OR.KFLF(1)*KFLF(3) 
     &.LT.0.OR.KFLF(1)*(NFL(1)+NFL(2)+NFL(3)).LT.0)) GOTO 280 
      IF(NFET.EQ.0) KFLF(1)=1+INT((2.+PARJ(2))*RLU(0)) 
      IF(NFET.EQ.0) KFLF(2)=-KFLF(1) 
      IF(NFET.EQ.1) KFLF(2)=ISIGN(1+INT((2.+PARJ(2))*RLU(0)),-KFLF(1)) 
      IF(NFET.LE.2) KFLF(3)=0 
      IF(KFLF(3).NE.0) THEN 
        KFLFC=ISIGN(1000*MAX(IABS(KFLF(1)),IABS(KFLF(3)))+ 
     &  100*MIN(IABS(KFLF(1)),IABS(KFLF(3)))+1,KFLF(1)) 
        IF(KFLF(1).EQ.KFLF(3).OR.(1.+3.*PARJ(4))*RLU(0).GT.1.) 
     &  KFLFC=KFLFC+ISIGN(2,KFLFC) 
      ELSE 
        KFLFC=KFLF(1) 
      ENDIF 
      CALL LUKFDI(KFLFC,KFLF(2),KFLDMP,KF) 
      IF(KF.EQ.0) GOTO 280 
      DO 300 J=1,MAX(2,NFET) 
      NFL(IABS(KFLF(J)))=NFL(IABS(KFLF(J)))-ISIGN(1,KFLF(J)) 
  300 CONTINUE 
 
C...Store hadron at random among free positions. 
      NPOS=MIN(1+INT(RLU(0)*NREM),NREM) 
      DO 310 I=NSAV+NJET+1,N 
      IF(K(I,1).EQ.7) NPOS=NPOS-1 
      IF(K(I,1).EQ.1.OR.NPOS.NE.0) GOTO 310 
      K(I,1)=1 
      K(I,2)=KF 
      P(I,5)=ULMASS(K(I,2)) 
      P(I,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2) 
  310 CONTINUE 
      NREM=NREM-1 
      NREQ=(IABS(NFL(1))+IABS(NFL(2))+IABS(NFL(3))-IABS(NFL(1)+ 
     &NFL(2)+NFL(3)))/2+IABS(NFL(1)+NFL(2)+NFL(3))/3 
      IF(NREM.GT.0) GOTO 280 
 
C...Compensate for missing momentum in global scheme (3 options). 
  320 IF(MOD(MSTJ(3),5).NE.0.AND.MOD(MSTJ(3),5).NE.4) THEN 
        DO 340 J=1,3 
        PSI(J)=0. 
        DO 330 I=NSAV+NJET+1,N 
        PSI(J)=PSI(J)+P(I,J) 
  330   CONTINUE 
  340   CONTINUE 
        PSI(4)=PSI(1)**2+PSI(2)**2+PSI(3)**2 
        PWS=0. 
        DO 350 I=NSAV+NJET+1,N 
        IF(MOD(MSTJ(3),5).EQ.1) PWS=PWS+P(I,4) 
        IF(MOD(MSTJ(3),5).EQ.2) PWS=PWS+SQRT(P(I,5)**2+(PSI(1)*P(I,1)+ 
     &  PSI(2)*P(I,2)+PSI(3)*P(I,3))**2/PSI(4)) 
        IF(MOD(MSTJ(3),5).EQ.3) PWS=PWS+1. 
  350   CONTINUE 
        DO 370 I=NSAV+NJET+1,N 
        IF(MOD(MSTJ(3),5).EQ.1) PW=P(I,4) 
        IF(MOD(MSTJ(3),5).EQ.2) PW=SQRT(P(I,5)**2+(PSI(1)*P(I,1)+ 
     &  PSI(2)*P(I,2)+PSI(3)*P(I,3))**2/PSI(4)) 
        IF(MOD(MSTJ(3),5).EQ.3) PW=1. 
        DO 360 J=1,3 
        P(I,J)=P(I,J)-PSI(J)*PW/PWS 
  360   CONTINUE 
        P(I,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2) 
  370   CONTINUE 
 
C...Compensate for missing momentum withing each jet separately. 
      ELSEIF(MOD(MSTJ(3),5).EQ.4) THEN 
        DO 390 I=N+1,N+NJET 
        K(I,1)=0 
        DO 380 J=1,5 
        P(I,J)=0. 
  380   CONTINUE 
  390   CONTINUE 
        DO 410 I=NSAV+NJET+1,N 
        IR1=K(I,3) 
        IR2=N+IR1-NSAV 
        K(IR2,1)=K(IR2,1)+1 
        PLS=(P(I,1)*P(IR1,1)+P(I,2)*P(IR1,2)+P(I,3)*P(IR1,3))/ 
     &  (P(IR1,1)**2+P(IR1,2)**2+P(IR1,3)**2) 
        DO 400 J=1,3 
        P(IR2,J)=P(IR2,J)+P(I,J)-PLS*P(IR1,J) 
  400   CONTINUE 
        P(IR2,4)=P(IR2,4)+P(I,4) 
        P(IR2,5)=P(IR2,5)+PLS 
  410   CONTINUE 
        PSS=0. 
        DO 420 I=N+1,N+NJET 
        IF(K(I,1).NE.0) PSS=PSS+P(I,4)/(PECM*(0.8*P(I,5)+0.2)) 
  420   CONTINUE 
        DO 440 I=NSAV+NJET+1,N 
        IR1=K(I,3) 
        IR2=N+IR1-NSAV 
        PLS=(P(I,1)*P(IR1,1)+P(I,2)*P(IR1,2)+P(I,3)*P(IR1,3))/ 
     &  (P(IR1,1)**2+P(IR1,2)**2+P(IR1,3)**2) 
        DO 430 J=1,3 
        P(I,J)=P(I,J)-P(IR2,J)/K(IR2,1)+(1./(P(IR2,5)*PSS)-1.)*PLS* 
     &  P(IR1,J) 
  430   CONTINUE 
        P(I,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2) 
  440   CONTINUE 
      ENDIF 
 
C...Scale momenta for energy conservation. 
      IF(MOD(MSTJ(3),5).NE.0) THEN 
        PMS=0. 
        PES=0. 
        PQS=0. 
        DO 450 I=NSAV+NJET+1,N 
        PMS=PMS+P(I,5) 
        PES=PES+P(I,4) 
        PQS=PQS+P(I,5)**2/P(I,4) 
  450   CONTINUE 
        IF(PMS.GE.PECM) GOTO 150 
        NECO=0 
  460   NECO=NECO+1 
        PFAC=(PECM-PQS)/(PES-PQS) 
        PES=0. 
        PQS=0. 
        DO 480 I=NSAV+NJET+1,N 
        DO 470 J=1,3 
        P(I,J)=PFAC*P(I,J) 
  470   CONTINUE 
        P(I,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2) 
        PES=PES+P(I,4) 
        PQS=PQS+P(I,5)**2/P(I,4) 
  480   CONTINUE 
        IF(NECO.LT.10.AND.ABS(PECM-PES).GT.2E-6*PECM) GOTO 460 
      ENDIF 
 
C...Origin of produced particles and parton daughter pointers. 
  490 DO 500 I=NSAV+NJET+1,N 
      IF(MSTU(16).NE.2) K(I,3)=NSAV+1 
      IF(MSTU(16).EQ.2) K(I,3)=K(K(I,3),3) 
  500 CONTINUE 
      DO 510 I=NSAV+1,NSAV+NJET 
      I1=K(I,3) 
      K(I1,1)=K(I1,1)+10 
      IF(MSTU(16).NE.2) THEN 
        K(I1,4)=NSAV+1 
        K(I1,5)=NSAV+1 
      ELSE 
        K(I1,4)=K(I1,4)-NJET+1 
        K(I1,5)=K(I1,5)-NJET+1 
        IF(K(I1,5).LT.K(I1,4)) THEN 
          K(I1,4)=0 
          K(I1,5)=0 
        ENDIF 
      ENDIF 
  510 CONTINUE 
 
C...Document independent fragmentation system. Remove copy of jets. 
      NSAV=NSAV+1 
      K(NSAV,1)=11 
      K(NSAV,2)=93 
      K(NSAV,3)=IP 
      K(NSAV,4)=NSAV+1 
      K(NSAV,5)=N-NJET+1 
      DO 520 J=1,4 
      P(NSAV,J)=DPS(J) 
      V(NSAV,J)=V(IP,J) 
  520 CONTINUE 
      P(NSAV,5)=SQRT(MAX(0D0,DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2)) 
      V(NSAV,5)=0. 
      DO 540 I=NSAV+NJET,N 
      DO 530 J=1,5 
      K(I-NJET+1,J)=K(I,J) 
      P(I-NJET+1,J)=P(I,J) 
      V(I-NJET+1,J)=V(I,J) 
  530 CONTINUE 
  540 CONTINUE 
      N=N-NJET+1 
      DO 550 IZ=MSTU90+1,MSTU(90) 
      MSTU(90+IZ)=MSTU(90+IZ)-NJET+1 
  550 CONTINUE 
 
C...Boost back particle system. Set production vertices. 
      IF(NJET.NE.1) CALL LUDBRB(NSAV+1,N,0.,0.,DPS(1)/DPS(4), 
     &DPS(2)/DPS(4),DPS(3)/DPS(4)) 
      DO 570 I=NSAV+1,N 
      DO 560 J=1,4 
      V(I,J)=V(IP,J) 
  560 CONTINUE 
  570 CONTINUE 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUDECY(IP) 
 
C...Purpose: to handle the decay of unstable particles. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/ 
      DIMENSION VDCY(4),KFLO(4),KFL1(4),PV(10,5),RORD(10),UE(3),BE(3), 
     &WTCOR(10),PTAU(4),PCMTAU(4) 
      DOUBLE PRECISION DBETAU(3) 
      DATA WTCOR/2.,5.,15.,60.,250.,1500.,1.2E4,1.2E5,150.,16./ 
 
C...Functions: momentum in two-particle decays, four-product and 
C...matrix element times phase space in weak decays. 
      PAWT(A,B,C)=SQRT((A**2-(B+C)**2)*(A**2-(B-C)**2))/(2.*A) 
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3) 
      HMEPS(HA)=((1.-HRQ-HA)**2+3.*HA*(1.+HRQ-HA))* 
     &SQRT((1.-HRQ-HA)**2-4.*HRQ*HA) 
 
C...Initial values. 
      NTRY=0 
      NSAV=N 
      KFA=IABS(K(IP,2)) 
      KFS=ISIGN(1,K(IP,2)) 
      KC=LUCOMP(KFA) 
      MSTJ(92)=0 
 
C...Choose lifetime and determine decay vertex. 
      IF(K(IP,1).EQ.5) THEN 
        V(IP,5)=0. 
      ELSEIF(K(IP,1).NE.4) THEN 
        V(IP,5)=-PMAS(KC,4)*LOG(RLU(0)) 
      ENDIF 
      DO 100 J=1,4 
      VDCY(J)=V(IP,J)+V(IP,5)*P(IP,J)/P(IP,5) 
  100 CONTINUE 
 
C...Determine whether decay allowed or not. 
      MOUT=0 
      IF(MSTJ(22).EQ.2) THEN 
        IF(PMAS(KC,4).GT.PARJ(71)) MOUT=1 
      ELSEIF(MSTJ(22).EQ.3) THEN 
        IF(VDCY(1)**2+VDCY(2)**2+VDCY(3)**2.GT.PARJ(72)**2) MOUT=1 
      ELSEIF(MSTJ(22).EQ.4) THEN 
        IF(VDCY(1)**2+VDCY(2)**2.GT.PARJ(73)**2) MOUT=1 
        IF(ABS(VDCY(3)).GT.PARJ(74)) MOUT=1 
      ENDIF 
      IF(MOUT.EQ.1.AND.K(IP,1).NE.5) THEN 
        K(IP,1)=4 
        RETURN 
      ENDIF 
 
C...Interface to external tau decay library (for tau polarization). 
      IF(KFA.EQ.15.AND.MSTJ(28).GE.1) THEN 
 
C...Starting values for pointers and momenta. 
        ITAU=IP 
        DO 110 J=1,4 
        PTAU(J)=P(ITAU,J) 
        PCMTAU(J)=P(ITAU,J) 
  110   CONTINUE 
 
C...Iterate to find position and code of mother of tau. 
        IMTAU=ITAU 
  120   IMTAU=K(IMTAU,3) 
 
        IF(IMTAU.EQ.0) THEN 
C...If no known origin then impossible to do anything further. 
          KFORIG=0 
          IORIG=0 
 
        ELSEIF(K(IMTAU,2).EQ.K(ITAU,2)) THEN 
C...If tau -> tau + gamma then add gamma energy and loop. 
          IF(K(K(IMTAU,4),2).EQ.22) THEN 
            DO 130 J=1,4 
            PCMTAU(J)=PCMTAU(J)+P(K(IMTAU,4),J) 
  130       CONTINUE 
          ELSEIF(K(K(IMTAU,5),2).EQ.22) THEN 
            DO 140 J=1,4 
            PCMTAU(J)=PCMTAU(J)+P(K(IMTAU,5),J) 
  140       CONTINUE 
          ENDIF 
          GOTO 120 
 
        ELSEIF(IABS(K(IMTAU,2)).GT.100) THEN 
C...If coming from weak decay of hadron then W is not stored in record, 
C...but can be reconstructed by adding neutrino momentum. 
          KFORIG=-ISIGN(24,K(ITAU,2)) 
          IORIG=0 
          DO 160 II=K(IMTAU,4),K(IMTAU,5) 
          IF(K(II,2)*ISIGN(1,K(ITAU,2)).EQ.-16) THEN 
            DO 150 J=1,4 
            PCMTAU(J)=PCMTAU(J)+P(II,J) 
  150       CONTINUE 
          ENDIF 
  160     CONTINUE 
 
        ELSE 
C...If coming from resonance decay then find latest copy of this 
C...resonance (may not completely agree). 
          KFORIG=K(IMTAU,2) 
          IORIG=IMTAU 
          DO 170 II=IMTAU+1,IP-1 
          IF(K(II,2).EQ.KFORIG.AND.K(II,3).EQ.IORIG.AND. 
     &    ABS(P(II,5)-P(IORIG,5)).LT.1E-5*P(IORIG,5)) IORIG=II 
  170     CONTINUE 
          DO 180 J=1,4 
          PCMTAU(J)=P(IORIG,J) 
  180     CONTINUE 
        ENDIF 
 
C...Boost tau to rest frame of production process (where known) 
C...and rotate it to sit along +z axis. 
        DO 190 J=1,3 
        DBETAU(J)=PCMTAU(J)/PCMTAU(4) 
  190   CONTINUE 
        IF(KFORIG.NE.0) CALL LUDBRB(ITAU,ITAU,0.,0.,-DBETAU(1), 
     &  -DBETAU(2),-DBETAU(3)) 
        PHITAU=ULANGL(P(ITAU,1),P(ITAU,2)) 
        CALL LUDBRB(ITAU,ITAU,0.,-PHITAU,0D0,0D0,0D0) 
        THETAU=ULANGL(P(ITAU,3),P(ITAU,1)) 
        CALL LUDBRB(ITAU,ITAU,-THETAU,0.,0D0,0D0,0D0) 
 
C...Call tau decay routine (if meaningful) and fill extra info. 
        IF(KFORIG.NE.0.OR.MSTJ(28).EQ.2) THEN 
          CALL LUTAUD(ITAU,IORIG,KFORIG,NDECAY) 
          DO 200 II=NSAV+1,NSAV+NDECAY 
          K(II,1)=1 
          K(II,3)=IP 
          K(II,4)=0 
          K(II,5)=0 
  200     CONTINUE 
          N=NSAV+NDECAY 
        ENDIF 
 
C...Boost back decay tau and decay products. 
        DO 210 J=1,4 
        P(ITAU,J)=PTAU(J) 
  210   CONTINUE 
        IF(KFORIG.NE.0.OR.MSTJ(28).EQ.2) THEN 
          CALL LUDBRB(NSAV+1,N,THETAU,PHITAU,0D0,0D0,0D0) 
          IF(KFORIG.NE.0) CALL LUDBRB(NSAV+1,N,0.,0.,DBETAU(1), 
     &    DBETAU(2),DBETAU(3)) 
 
C...Skip past ordinary tau decay treatment. 
          MMAT=0 
          MBST=0 
          ND=0 
          GOTO 660 
        ENDIF 
      ENDIF 
 
C...B-B~ mixing: flip sign of meson appropriately. 
      MMIX=0 
      IF((KFA.EQ.511.OR.KFA.EQ.531).AND.MSTJ(26).GE.1) THEN 
        XBBMIX=PARJ(76) 
        IF(KFA.EQ.531) XBBMIX=PARJ(77) 
        IF(SIN(0.5*XBBMIX*V(IP,5)/PMAS(KC,4))**2.GT.RLU(0)) MMIX=1 
        IF(MMIX.EQ.1) KFS=-KFS 
      ENDIF 
 
C...Check existence of decay channels. Particle/antiparticle rules. 
      KCA=KC 
      IF(MDCY(KC,2).GT.0) THEN 
        MDMDCY=MDME(MDCY(KC,2),2) 
        IF(MDMDCY.GT.80.AND.MDMDCY.LE.90) KCA=MDMDCY 
      ENDIF 
      IF(MDCY(KCA,2).LE.0.OR.MDCY(KCA,3).LE.0) THEN 
        CALL LUERRM(9,'(LUDECY:) no decay channel defined') 
        RETURN 
      ENDIF 
      IF(MOD(KFA/1000,10).EQ.0.AND.(KCA.EQ.85.OR.KCA.EQ.87)) KFS=-KFS 
      IF(KCHG(KC,3).EQ.0) THEN 
        KFSP=1 
        KFSN=0 
        IF(RLU(0).GT.0.5) KFS=-KFS 
      ELSEIF(KFS.GT.0) THEN 
        KFSP=1 
        KFSN=0 
      ELSE 
        KFSP=0 
        KFSN=1 
      ENDIF 
 
C...Sum branching ratios of allowed decay channels. 
  220 NOPE=0 
      BRSU=0. 
      DO 230 IDL=MDCY(KCA,2),MDCY(KCA,2)+MDCY(KCA,3)-1 
      IF(MDME(IDL,1).NE.1.AND.KFSP*MDME(IDL,1).NE.2.AND. 
     &KFSN*MDME(IDL,1).NE.3) GOTO 230 
      IF(MDME(IDL,2).GT.100) GOTO 230 
      NOPE=NOPE+1 
      BRSU=BRSU+BRAT(IDL) 
  230 CONTINUE 
      IF(NOPE.EQ.0) THEN 
        CALL LUERRM(2,'(LUDECY:) all decay channels closed by user') 
        RETURN 
      ENDIF 
 
C...Select decay channel among allowed ones. 
  240 RBR=BRSU*RLU(0) 
      IDL=MDCY(KCA,2)-1 
  250 IDL=IDL+1 
      IF(MDME(IDL,1).NE.1.AND.KFSP*MDME(IDL,1).NE.2.AND. 
     &KFSN*MDME(IDL,1).NE.3) THEN 
        IF(IDL.LT.MDCY(KCA,2)+MDCY(KCA,3)-1) GOTO 250 
      ELSEIF(MDME(IDL,2).GT.100) THEN 
        IF(IDL.LT.MDCY(KCA,2)+MDCY(KCA,3)-1) GOTO 250 
      ELSE 
        IDC=IDL 
        RBR=RBR-BRAT(IDL) 
        IF(IDL.LT.MDCY(KCA,2)+MDCY(KCA,3)-1.AND.RBR.GT.0.) GOTO 250 
      ENDIF 
 
C...Start readout of decay channel: matrix element, reset counters. 
      MMAT=MDME(IDC,2) 
  260 NTRY=NTRY+1 
      IF(NTRY.GT.1000) THEN 
        CALL LUERRM(14,'(LUDECY:) caught in infinite loop') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      I=N 
      NP=0 
      NQ=0 
      MBST=0 
      IF(MMAT.GE.11.AND.MMAT.NE.46.AND.P(IP,4).GT.20.*P(IP,5)) MBST=1 
      DO 270 J=1,4 
      PV(1,J)=0. 
      IF(MBST.EQ.0) PV(1,J)=P(IP,J) 
  270 CONTINUE 
      IF(MBST.EQ.1) PV(1,4)=P(IP,5) 
      PV(1,5)=P(IP,5) 
      PS=0. 
      PSQ=0. 
      MREM=0 
      MHADDY=0 
      IF(KFA.GT.80) MHADDY=1 
 
C...Read out decay products. Convert to standard flavour code. 
      JTMAX=5 
      IF(MDME(IDC+1,2).EQ.101) JTMAX=10 
      DO 280 JT=1,JTMAX 
      IF(JT.LE.5) KP=KFDP(IDC,JT) 
      IF(JT.GE.6) KP=KFDP(IDC+1,JT-5) 
      IF(KP.EQ.0) GOTO 280 
      KPA=IABS(KP) 
      KCP=LUCOMP(KPA) 
      IF(KPA.GT.80) MHADDY=1 
      IF(KCHG(KCP,3).EQ.0.AND.KPA.NE.81.AND.KPA.NE.82) THEN 
        KFP=KP 
      ELSEIF(KPA.NE.81.AND.KPA.NE.82) THEN 
        KFP=KFS*KP 
      ELSEIF(KPA.EQ.81.AND.MOD(KFA/1000,10).EQ.0) THEN 
        KFP=-KFS*MOD(KFA/10,10) 
      ELSEIF(KPA.EQ.81.AND.MOD(KFA/100,10).GE.MOD(KFA/10,10)) THEN 
        KFP=KFS*(100*MOD(KFA/10,100)+3) 
      ELSEIF(KPA.EQ.81) THEN 
        KFP=KFS*(1000*MOD(KFA/10,10)+100*MOD(KFA/100,10)+1) 
      ELSEIF(KP.EQ.82) THEN 
        CALL LUKFDI(-KFS*INT(1.+(2.+PARJ(2))*RLU(0)),0,KFP,KDUMP) 
        IF(KFP.EQ.0) GOTO 260 
        MSTJ(93)=1 
        IF(PV(1,5).LT.PARJ(32)+2.*ULMASS(KFP)) GOTO 260 
      ELSEIF(KP.EQ.-82) THEN 
        KFP=-KFP 
        IF(IABS(KFP).GT.10) KFP=KFP+ISIGN(10000,KFP) 
      ENDIF 
      IF(KPA.EQ.81.OR.KPA.EQ.82) KCP=LUCOMP(KFP) 
 
C...Add decay product to event record or to quark flavour list. 
      KFPA=IABS(KFP) 
      KQP=KCHG(KCP,2) 
      IF(MMAT.GE.11.AND.MMAT.LE.30.AND.KQP.NE.0) THEN 
        NQ=NQ+1 
        KFLO(NQ)=KFP 
        MSTJ(93)=2 
        PSQ=PSQ+ULMASS(KFLO(NQ)) 
      ELSEIF((MMAT.EQ.42.OR.MMAT.EQ.43.OR.MMAT.EQ.48).AND.NP.EQ.3.AND. 
     &MOD(NQ,2).EQ.1) THEN 
        NQ=NQ-1 
        PS=PS-P(I,5) 
        K(I,1)=1 
        KFI=K(I,2) 
        CALL LUKFDI(KFP,KFI,KFLDMP,K(I,2)) 
        IF(K(I,2).EQ.0) GOTO 260 
        MSTJ(93)=1 
        P(I,5)=ULMASS(K(I,2)) 
        PS=PS+P(I,5) 
      ELSE 
        I=I+1 
        NP=NP+1 
        IF(MMAT.NE.33.AND.KQP.NE.0) NQ=NQ+1 
        IF(MMAT.EQ.33.AND.KQP.NE.0.AND.KQP.NE.2) NQ=NQ+1 
        K(I,1)=1+MOD(NQ,2) 
        IF(MMAT.EQ.4.AND.JT.LE.2.AND.KFP.EQ.21) K(I,1)=2 
        IF(MMAT.EQ.4.AND.JT.EQ.3) K(I,1)=1 
        K(I,2)=KFP 
        K(I,3)=IP 
        K(I,4)=0 
        K(I,5)=0 
        P(I,5)=ULMASS(KFP) 
        IF(MMAT.EQ.45.AND.KFPA.EQ.89) P(I,5)=PARJ(32) 
        PS=PS+P(I,5) 
      ENDIF 
  280 CONTINUE 
 
C...Check masses for resonance decays. 
      IF(MHADDY.EQ.0) THEN 
        IF(PS+PARJ(64).GT.PV(1,5)) GOTO 240 
      ENDIF 
 
C...Choose decay multiplicity in phase space model. 
  290 IF(MMAT.GE.11.AND.MMAT.LE.30) THEN 
        PSP=PS 
        CNDE=PARJ(61)*LOG(MAX((PV(1,5)-PS-PSQ)/PARJ(62),1.1)) 
        IF(MMAT.EQ.12) CNDE=CNDE+PARJ(63) 
  300   NTRY=NTRY+1 
        IF(NTRY.GT.1000) THEN 
          CALL LUERRM(14,'(LUDECY:) caught in infinite loop') 
          IF(MSTU(21).GE.1) RETURN 
        ENDIF 
        IF(MMAT.LE.20) THEN 
          GAUSS=SQRT(-2.*CNDE*LOG(MAX(1E-10,RLU(0))))* 
     &    SIN(PARU(2)*RLU(0)) 
          ND=0.5+0.5*NP+0.25*NQ+CNDE+GAUSS 
          IF(ND.LT.NP+NQ/2.OR.ND.LT.2.OR.ND.GT.10) GOTO 300 
          IF(MMAT.EQ.13.AND.ND.EQ.2) GOTO 300 
          IF(MMAT.EQ.14.AND.ND.LE.3) GOTO 300 
          IF(MMAT.EQ.15.AND.ND.LE.4) GOTO 300 
        ELSE 
          ND=MMAT-20 
        ENDIF 
 
C...Form hadrons from flavour content. 
        DO 310 JT=1,4 
        KFL1(JT)=KFLO(JT) 
  310   CONTINUE 
        IF(ND.EQ.NP+NQ/2) GOTO 330 
        DO 320 I=N+NP+1,N+ND-NQ/2 
        JT=1+INT((NQ-1)*RLU(0)) 
        CALL LUKFDI(KFL1(JT),0,KFL2,K(I,2)) 
        IF(K(I,2).EQ.0) GOTO 300 
        KFL1(JT)=-KFL2 
  320   CONTINUE 
  330   JT=2 
        JT2=3 
        JT3=4 
        IF(NQ.EQ.4.AND.RLU(0).LT.PARJ(66)) JT=4 
        IF(JT.EQ.4.AND.ISIGN(1,KFL1(1)*(10-IABS(KFL1(1))))* 
     &  ISIGN(1,KFL1(JT)*(10-IABS(KFL1(JT)))).GT.0) JT=3 
        IF(JT.EQ.3) JT2=2 
        IF(JT.EQ.4) JT3=2 
        CALL LUKFDI(KFL1(1),KFL1(JT),KFLDMP,K(N+ND-NQ/2+1,2)) 
        IF(K(N+ND-NQ/2+1,2).EQ.0) GOTO 300 
        IF(NQ.EQ.4) CALL LUKFDI(KFL1(JT2),KFL1(JT3),KFLDMP,K(N+ND,2)) 
        IF(NQ.EQ.4.AND.K(N+ND,2).EQ.0) GOTO 300 
 
C...Check that sum of decay product masses not too large. 
        PS=PSP 
        DO 340 I=N+NP+1,N+ND 
        K(I,1)=1 
        K(I,3)=IP 
        K(I,4)=0 
        K(I,5)=0 
        P(I,5)=ULMASS(K(I,2)) 
        PS=PS+P(I,5) 
  340   CONTINUE 
        IF(PS+PARJ(64).GT.PV(1,5)) GOTO 300 
 
C...Rescale energy to subtract off spectator quark mass. 
      ELSEIF((MMAT.EQ.31.OR.MMAT.EQ.33.OR.MMAT.EQ.44.OR.MMAT.EQ.45) 
     &.AND.NP.GE.3) THEN 
        PS=PS-P(N+NP,5) 
        PQT=(P(N+NP,5)+PARJ(65))/PV(1,5) 
        DO 350 J=1,5 
        P(N+NP,J)=PQT*PV(1,J) 
        PV(1,J)=(1.-PQT)*PV(1,J) 
  350   CONTINUE 
        IF(PS+PARJ(64).GT.PV(1,5)) GOTO 260 
        ND=NP-1 
        MREM=1 
 
C...Phase space factors imposed in W decay. 
      ELSEIF(MMAT.EQ.46) THEN 
        MSTJ(93)=1 
        PSMC=ULMASS(K(N+1,2)) 
        MSTJ(93)=1 
        PSMC=PSMC+ULMASS(K(N+2,2)) 
        IF(MAX(PS,PSMC)+PARJ(32).GT.PV(1,5)) GOTO 240 
        HR1=(P(N+1,5)/PV(1,5))**2 
        HR2=(P(N+2,5)/PV(1,5))**2 
        IF((1.-HR1-HR2)*(2.+HR1+HR2)*SQRT((1.-HR1-HR2)**2-4.*HR1*HR2) 
     &  .LT.2.*RLU(0)) GOTO 240 
        ND=NP 
 
C...Fully specified final state: check mass broadening effects. 
      ELSE 
        IF(NP.GE.2.AND.PS+PARJ(64).GT.PV(1,5)) GOTO 260 
        ND=NP 
      ENDIF 
 
C...Select W mass in decay Q -> W + q, without W propagator. 
      IF(MMAT.EQ.45.AND.MSTJ(25).LE.0) THEN 
        HLQ=(PARJ(32)/PV(1,5))**2 
        HUQ=(1.-(P(N+2,5)+PARJ(64))/PV(1,5))**2 
        HRQ=(P(N+2,5)/PV(1,5))**2 
  360   HW=HLQ+RLU(0)*(HUQ-HLQ) 
        IF(HMEPS(HW).LT.RLU(0)) GOTO 360 
        P(N+1,5)=PV(1,5)*SQRT(HW) 
 
C...Ditto, including W propagator. Divide mass range into three regions. 
      ELSEIF(MMAT.EQ.45) THEN 
        HQW=(PV(1,5)/PMAS(24,1))**2 
        HLW=(PARJ(32)/PMAS(24,1))**2 
        HUW=((PV(1,5)-P(N+2,5)-PARJ(64))/PMAS(24,1))**2 
        HRQ=(P(N+2,5)/PV(1,5))**2 
        HG=PMAS(24,2)/PMAS(24,1) 
        HATL=ATAN((HLW-1.)/HG) 
        HM=MIN(1.,HUW-0.001) 
        HMV1=HMEPS(HM/HQW)/((HM-1.)**2+HG**2) 
  370   HM=HM-HG 
        HMV2=HMEPS(HM/HQW)/((HM-1.)**2+HG**2) 
        IF(HMV2.GT.HMV1.AND.HM-HG.GT.HLW) THEN 
          HMV1=HMV2 
          GOTO 370 
        ENDIF 
        HMV=MIN(2.*HMV1,HMEPS(HM/HQW)/HG**2) 
        HM1=1.-SQRT(1./HMV-HG**2) 
        IF(HM1.GT.HLW.AND.HM1.LT.HM) THEN 
          HM=HM1 
        ELSEIF(HMV2.LE.HMV1) THEN 
          HM=MAX(HLW,HM-MIN(0.1,1.-HM)) 
        ENDIF 
        HATM=ATAN((HM-1.)/HG) 
        HWT1=(HATM-HATL)/HG 
        HWT2=HMV*(MIN(1.,HUW)-HM) 
        HWT3=0. 
        IF(HUW.GT.1.) THEN 
          HATU=ATAN((HUW-1.)/HG) 
          HMP1=HMEPS(1./HQW) 
          HWT3=HMP1*HATU/HG 
        ENDIF 
 
C...Select mass region and W mass there. Accept according to weight. 
  380   HREG=RLU(0)*(HWT1+HWT2+HWT3) 
        IF(HREG.LE.HWT1) THEN 
          HW=1.+HG*TAN(HATL+RLU(0)*(HATM-HATL)) 
          HACC=HMEPS(HW/HQW) 
        ELSEIF(HREG.LE.HWT1+HWT2) THEN 
          HW=HM+RLU(0)*(MIN(1.,HUW)-HM) 
          HACC=HMEPS(HW/HQW)/((HW-1.)**2+HG**2)/HMV 
        ELSE 
          HW=1.+HG*TAN(RLU(0)*HATU) 
          HACC=HMEPS(HW/HQW)/HMP1 
        ENDIF 
        IF(HACC.LT.RLU(0)) GOTO 380 
        P(N+1,5)=PMAS(24,1)*SQRT(HW) 
      ENDIF 
 
C...Determine position of grandmother, number of sisters, Q -> W sign. 
      NM=0 
      KFAS=0 
      MSGN=0 
      IF(MMAT.EQ.3.OR.MMAT.EQ.46) THEN 
        IM=K(IP,3) 
        IF(IM.LT.0.OR.IM.GE.IP) IM=0 
        IF(MMAT.EQ.46.AND.MSTJ(27).EQ.1) THEN 
          IM=0 
        ELSEIF(MMAT.EQ.46.AND.MSTJ(27).GE.2.AND.IM.NE.0) THEN 
          IF(K(IM,2).EQ.94) THEN 
            IM=K(K(IM,3),3) 
            IF(IM.LT.0.OR.IM.GE.IP) IM=0 
          ENDIF 
        ENDIF 
        IF(IM.NE.0) KFAM=IABS(K(IM,2)) 
        IF(IM.NE.0.AND.MMAT.EQ.3) THEN 
          DO 390 IL=MAX(IP-2,IM+1),MIN(IP+2,N) 
          IF(K(IL,3).EQ.IM) NM=NM+1 
          IF(K(IL,3).EQ.IM.AND.IL.NE.IP) ISIS=IL 
  390     CONTINUE 
          IF(NM.NE.2.OR.KFAM.LE.100.OR.MOD(KFAM,10).NE.1.OR. 
     &    MOD(KFAM/1000,10).NE.0) NM=0 
          IF(NM.EQ.2) THEN 
            KFAS=IABS(K(ISIS,2)) 
            IF((KFAS.LE.100.OR.MOD(KFAS,10).NE.1.OR. 
     &      MOD(KFAS/1000,10).NE.0).AND.KFAS.NE.22) NM=0 
          ENDIF 
        ELSEIF(IM.NE.0.AND.MMAT.EQ.46) THEN 
          MSGN=ISIGN(1,K(IM,2)*K(IP,2)) 
          IF(KFAM.GT.100.AND.MOD(KFAM/1000,10).EQ.0) MSGN= 
     &    MSGN*(-1)**MOD(KFAM/100,10) 
        ENDIF 
      ENDIF 
 
C...Kinematics of one-particle decays. 
      IF(ND.EQ.1) THEN 
        DO 400 J=1,4 
        P(N+1,J)=P(IP,J) 
  400   CONTINUE 
        GOTO 660 
      ENDIF 
 
C...Calculate maximum weight ND-particle decay. 
      PV(ND,5)=P(N+ND,5) 
      IF(ND.GE.3) THEN 
        WTMAX=1./WTCOR(ND-2) 
        PMAX=PV(1,5)-PS+P(N+ND,5) 
        PMIN=0. 
        DO 410 IL=ND-1,1,-1 
        PMAX=PMAX+P(N+IL,5) 
        PMIN=PMIN+P(N+IL+1,5) 
        WTMAX=WTMAX*PAWT(PMAX,PMIN,P(N+IL,5)) 
  410   CONTINUE 
      ENDIF 
 
C...Find virtual gamma mass in Dalitz decay. 
  420 IF(ND.EQ.2) THEN 
      ELSEIF(MMAT.EQ.2) THEN 
        PMES=4.*PMAS(11,1)**2 
        PMRHO2=PMAS(131,1)**2 
        PGRHO2=PMAS(131,2)**2 
  430   PMST=PMES*(P(IP,5)**2/PMES)**RLU(0) 
        WT=(1+0.5*PMES/PMST)*SQRT(MAX(0.,1.-PMES/PMST))* 
     &  (1.-PMST/P(IP,5)**2)**3*(1.+PGRHO2/PMRHO2)/ 
     &  ((1.-PMST/PMRHO2)**2+PGRHO2/PMRHO2) 
        IF(WT.LT.RLU(0)) GOTO 430 
        PV(2,5)=MAX(2.00001*PMAS(11,1),SQRT(PMST)) 
 
C...M-generator gives weight. If rejected, try again. 
      ELSE 
  440   RORD(1)=1. 
        DO 470 IL1=2,ND-1 
        RSAV=RLU(0) 
        DO 450 IL2=IL1-1,1,-1 
        IF(RSAV.LE.RORD(IL2)) GOTO 460 
        RORD(IL2+1)=RORD(IL2) 
  450   CONTINUE 
  460   RORD(IL2+1)=RSAV 
  470   CONTINUE 
        RORD(ND)=0. 
        WT=1. 
        DO 480 IL=ND-1,1,-1 
        PV(IL,5)=PV(IL+1,5)+P(N+IL,5)+(RORD(IL)-RORD(IL+1))*(PV(1,5)-PS) 
        WT=WT*PAWT(PV(IL,5),PV(IL+1,5),P(N+IL,5)) 
  480   CONTINUE 
        IF(WT.LT.RLU(0)*WTMAX) GOTO 440 
      ENDIF 
 
C...Perform two-particle decays in respective CM frame. 
  490 DO 510 IL=1,ND-1 
      PA=PAWT(PV(IL,5),PV(IL+1,5),P(N+IL,5)) 
      UE(3)=2.*RLU(0)-1. 
      PHI=PARU(2)*RLU(0) 
      UE(1)=SQRT(1.-UE(3)**2)*COS(PHI) 
      UE(2)=SQRT(1.-UE(3)**2)*SIN(PHI) 
      DO 500 J=1,3 
      P(N+IL,J)=PA*UE(J) 
      PV(IL+1,J)=-PA*UE(J) 
  500 CONTINUE 
      P(N+IL,4)=SQRT(PA**2+P(N+IL,5)**2) 
      PV(IL+1,4)=SQRT(PA**2+PV(IL+1,5)**2) 
  510 CONTINUE 
 
C...Lorentz transform decay products to lab frame. 
      DO 520 J=1,4 
      P(N+ND,J)=PV(ND,J) 
  520 CONTINUE 
      DO 560 IL=ND-1,1,-1 
      DO 530 J=1,3 
      BE(J)=PV(IL,J)/PV(IL,4) 
  530 CONTINUE 
      GA=PV(IL,4)/PV(IL,5) 
      DO 550 I=N+IL,N+ND 
      BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3) 
      DO 540 J=1,3 
      P(I,J)=P(I,J)+GA*(GA*BEP/(1.+GA)+P(I,4))*BE(J) 
  540 CONTINUE 
      P(I,4)=GA*(P(I,4)+BEP) 
  550 CONTINUE 
  560 CONTINUE 
 
C...Check that no infinite loop in matrix element weight. 
      NTRY=NTRY+1 
      IF(NTRY.GT.800) GOTO 590 
 
C...Matrix elements for omega and phi decays. 
      IF(MMAT.EQ.1) THEN 
        WT=(P(N+1,5)*P(N+2,5)*P(N+3,5))**2-(P(N+1,5)*FOUR(N+2,N+3))**2 
     &  -(P(N+2,5)*FOUR(N+1,N+3))**2-(P(N+3,5)*FOUR(N+1,N+2))**2 
     &  +2.*FOUR(N+1,N+2)*FOUR(N+1,N+3)*FOUR(N+2,N+3) 
        IF(MAX(WT*WTCOR(9)/P(IP,5)**6,0.001).LT.RLU(0)) GOTO 420 
 
C...Matrix elements for pi0 or eta Dalitz decay to gamma e+ e-. 
      ELSEIF(MMAT.EQ.2) THEN 
        FOUR12=FOUR(N+1,N+2) 
        FOUR13=FOUR(N+1,N+3) 
        WT=(PMST-0.5*PMES)*(FOUR12**2+FOUR13**2)+ 
     &  PMES*(FOUR12*FOUR13+FOUR12**2+FOUR13**2) 
        IF(WT.LT.RLU(0)*0.25*PMST*(P(IP,5)**2-PMST)**2) GOTO 490 
 
C...Matrix element for S0 -> S1 + V1 -> S1 + S2 + S3 (S scalar, 
C...V vector), of form cos**2(theta02) in V1 rest frame, and for 
C...S0 -> gamma + V1 -> gamma + S2 + S3, of form sin**2(theta02). 
      ELSEIF(MMAT.EQ.3.AND.NM.EQ.2) THEN 
        FOUR10=FOUR(IP,IM) 
        FOUR12=FOUR(IP,N+1) 
        FOUR02=FOUR(IM,N+1) 
        PMS1=P(IP,5)**2 
        PMS0=P(IM,5)**2 
        PMS2=P(N+1,5)**2 
        IF(KFAS.NE.22) HNUM=(FOUR10*FOUR12-PMS1*FOUR02)**2 
        IF(KFAS.EQ.22) HNUM=PMS1*(2.*FOUR10*FOUR12*FOUR02- 
     &  PMS1*FOUR02**2-PMS0*FOUR12**2-PMS2*FOUR10**2+PMS1*PMS0*PMS2) 
        HNUM=MAX(1E-6*PMS1**2*PMS0*PMS2,HNUM) 
        HDEN=(FOUR10**2-PMS1*PMS0)*(FOUR12**2-PMS1*PMS2) 
        IF(HNUM.LT.RLU(0)*HDEN) GOTO 490 
 
C...Matrix element for "onium" -> g + g + g or gamma + g + g. 
      ELSEIF(MMAT.EQ.4) THEN 
        HX1=2.*FOUR(IP,N+1)/P(IP,5)**2 
        HX2=2.*FOUR(IP,N+2)/P(IP,5)**2 
        HX3=2.*FOUR(IP,N+3)/P(IP,5)**2 
        WT=((1.-HX1)/(HX2*HX3))**2+((1.-HX2)/(HX1*HX3))**2+ 
     &  ((1.-HX3)/(HX1*HX2))**2 
        IF(WT.LT.2.*RLU(0)) GOTO 420 
        IF(K(IP+1,2).EQ.22.AND.(1.-HX1)*P(IP,5)**2.LT.4.*PARJ(32)**2) 
     &  GOTO 420 
 
C...Effective matrix element for nu spectrum in tau -> nu + hadrons. 
      ELSEIF(MMAT.EQ.41) THEN 
        HX1=2.*FOUR(IP,N+1)/P(IP,5)**2 
        HXM=MIN(0.75,2.*(1.-PS/P(IP,5))) 
        IF(HX1*(3.-2.*HX1).LT.RLU(0)*HXM*(3.-2.*HXM)) GOTO 420 
 
C...Matrix elements for weak decays (only semileptonic for c and b) 
      ELSEIF((MMAT.EQ.42.OR.MMAT.EQ.43.OR.MMAT.EQ.44.OR.MMAT.EQ.48) 
     &.AND.ND.EQ.3) THEN 
        IF(MBST.EQ.0) WT=FOUR(IP,N+1)*FOUR(N+2,N+3) 
        IF(MBST.EQ.1) WT=P(IP,5)*P(N+1,4)*FOUR(N+2,N+3) 
        IF(WT.LT.RLU(0)*P(IP,5)*PV(1,5)**3/WTCOR(10)) GOTO 420 
      ELSEIF(MMAT.EQ.42.OR.MMAT.EQ.43.OR.MMAT.EQ.44.OR.MMAT.EQ.48) THEN 
        DO 580 J=1,4 
        P(N+NP+1,J)=0. 
        DO 570 IS=N+3,N+NP 
        P(N+NP+1,J)=P(N+NP+1,J)+P(IS,J) 
  570   CONTINUE 
  580   CONTINUE 
        IF(MBST.EQ.0) WT=FOUR(IP,N+1)*FOUR(N+2,N+NP+1) 
        IF(MBST.EQ.1) WT=P(IP,5)*P(N+1,4)*FOUR(N+2,N+NP+1) 
        IF(WT.LT.RLU(0)*P(IP,5)*PV(1,5)**3/WTCOR(10)) GOTO 420 
 
C...Angular distribution in W decay. 
      ELSEIF(MMAT.EQ.46.AND.MSGN.NE.0) THEN 
        IF(MSGN.GT.0) WT=FOUR(IM,N+1)*FOUR(N+2,IP+1) 
        IF(MSGN.LT.0) WT=FOUR(IM,N+2)*FOUR(N+1,IP+1) 
        IF(WT.LT.RLU(0)*P(IM,5)**4/WTCOR(10)) GOTO 490 
      ENDIF 
 
C...Scale back energy and reattach spectator. 
  590 IF(MREM.EQ.1) THEN 
        DO 600 J=1,5 
        PV(1,J)=PV(1,J)/(1.-PQT) 
  600   CONTINUE 
        ND=ND+1 
        MREM=0 
      ENDIF 
 
C...Low invariant mass for system with spectator quark gives particle, 
C...not two jets. Readjust momenta accordingly. 
      IF((MMAT.EQ.31.OR.MMAT.EQ.45).AND.ND.EQ.3) THEN 
        MSTJ(93)=1 
        PM2=ULMASS(K(N+2,2)) 
        MSTJ(93)=1 
        PM3=ULMASS(K(N+3,2)) 
        IF(P(N+2,5)**2+P(N+3,5)**2+2.*FOUR(N+2,N+3).GE. 
     &  (PARJ(32)+PM2+PM3)**2) GOTO 660 
        K(N+2,1)=1 
        KFTEMP=K(N+2,2) 
        CALL LUKFDI(KFTEMP,K(N+3,2),KFLDMP,K(N+2,2)) 
        IF(K(N+2,2).EQ.0) GOTO 260 
        P(N+2,5)=ULMASS(K(N+2,2)) 
        PS=P(N+1,5)+P(N+2,5) 
        PV(2,5)=P(N+2,5) 
        MMAT=0 
        ND=2 
        GOTO 490 
      ELSEIF(MMAT.EQ.44) THEN 
        MSTJ(93)=1 
        PM3=ULMASS(K(N+3,2)) 
        MSTJ(93)=1 
        PM4=ULMASS(K(N+4,2)) 
        IF(P(N+3,5)**2+P(N+4,5)**2+2.*FOUR(N+3,N+4).GE. 
     &  (PARJ(32)+PM3+PM4)**2) GOTO 630 
        K(N+3,1)=1 
        KFTEMP=K(N+3,2) 
        CALL LUKFDI(KFTEMP,K(N+4,2),KFLDMP,K(N+3,2)) 
        IF(K(N+3,2).EQ.0) GOTO 260 
        P(N+3,5)=ULMASS(K(N+3,2)) 
        DO 610 J=1,3 
        P(N+3,J)=P(N+3,J)+P(N+4,J) 
  610   CONTINUE 
        P(N+3,4)=SQRT(P(N+3,1)**2+P(N+3,2)**2+P(N+3,3)**2+P(N+3,5)**2) 
        HA=P(N+1,4)**2-P(N+2,4)**2 
        HB=HA-(P(N+1,5)**2-P(N+2,5)**2) 
        HC=(P(N+1,1)-P(N+2,1))**2+(P(N+1,2)-P(N+2,2))**2+ 
     &  (P(N+1,3)-P(N+2,3))**2 
        HD=(PV(1,4)-P(N+3,4))**2 
        HE=HA**2-2.*HD*(P(N+1,4)**2+P(N+2,4)**2)+HD**2 
        HF=HD*HC-HB**2 
        HG=HD*HC-HA*HB 
        HH=(SQRT(HG**2+HE*HF)-HG)/(2.*HF) 
        DO 620 J=1,3 
        PCOR=HH*(P(N+1,J)-P(N+2,J)) 
        P(N+1,J)=P(N+1,J)+PCOR 
        P(N+2,J)=P(N+2,J)-PCOR 
  620   CONTINUE 
        P(N+1,4)=SQRT(P(N+1,1)**2+P(N+1,2)**2+P(N+1,3)**2+P(N+1,5)**2) 
        P(N+2,4)=SQRT(P(N+2,1)**2+P(N+2,2)**2+P(N+2,3)**2+P(N+2,5)**2) 
        ND=ND-1 
      ENDIF 
 
C...Check invariant mass of W jets. May give one particle or start over. 
  630 IF((MMAT.EQ.42.OR.MMAT.EQ.43.OR.MMAT.EQ.44.OR.MMAT.EQ.48) 
     &.AND.IABS(K(N+1,2)).LT.10) THEN 
        PMR=SQRT(MAX(0.,P(N+1,5)**2+P(N+2,5)**2+2.*FOUR(N+1,N+2))) 
        MSTJ(93)=1 
        PM1=ULMASS(K(N+1,2)) 
        MSTJ(93)=1 
        PM2=ULMASS(K(N+2,2)) 
        IF(PMR.GT.PARJ(32)+PM1+PM2) GOTO 640 
        KFLDUM=INT(1.5+RLU(0)) 
        CALL LUKFDI(K(N+1,2),-ISIGN(KFLDUM,K(N+1,2)),KFLDMP,KF1) 
        CALL LUKFDI(K(N+2,2),-ISIGN(KFLDUM,K(N+2,2)),KFLDMP,KF2) 
        IF(KF1.EQ.0.OR.KF2.EQ.0) GOTO 260 
        PSM=ULMASS(KF1)+ULMASS(KF2) 
        IF((MMAT.EQ.42.OR.MMAT.EQ.48).AND.PMR.GT.PARJ(64)+PSM) GOTO 640 
        IF(MMAT.GE.43.AND.PMR.GT.0.2*PARJ(32)+PSM) GOTO 640 
        IF(MMAT.EQ.48) GOTO 420 
        IF(ND.EQ.4.OR.KFA.EQ.15) GOTO 260 
        K(N+1,1)=1 
        KFTEMP=K(N+1,2) 
        CALL LUKFDI(KFTEMP,K(N+2,2),KFLDMP,K(N+1,2)) 
        IF(K(N+1,2).EQ.0) GOTO 260 
        P(N+1,5)=ULMASS(K(N+1,2)) 
        K(N+2,2)=K(N+3,2) 
        P(N+2,5)=P(N+3,5) 
        PS=P(N+1,5)+P(N+2,5) 
        IF(PS+PARJ(64).GT.PV(1,5)) GOTO 260 
        PV(2,5)=P(N+3,5) 
        MMAT=0 
        ND=2 
        GOTO 490 
      ENDIF 
 
C...Phase space decay of partons from W decay. 
  640 IF((MMAT.EQ.42.OR.MMAT.EQ.48).AND.IABS(K(N+1,2)).LT.10) THEN 
        KFLO(1)=K(N+1,2) 
        KFLO(2)=K(N+2,2) 
        K(N+1,1)=K(N+3,1) 
        K(N+1,2)=K(N+3,2) 
        DO 650 J=1,5 
        PV(1,J)=P(N+1,J)+P(N+2,J) 
        P(N+1,J)=P(N+3,J) 
  650   CONTINUE 
        PV(1,5)=PMR 
        N=N+1 
        NP=0 
        NQ=2 
        PS=0. 
        MSTJ(93)=2 
        PSQ=ULMASS(KFLO(1)) 
        MSTJ(93)=2 
        PSQ=PSQ+ULMASS(KFLO(2)) 
        MMAT=11 
        GOTO 290 
      ENDIF 
 
C...Boost back for rapidly moving particle. 
  660 N=N+ND 
      IF(MBST.EQ.1) THEN 
        DO 670 J=1,3 
        BE(J)=P(IP,J)/P(IP,4) 
  670   CONTINUE 
        GA=P(IP,4)/P(IP,5) 
        DO 690 I=NSAV+1,N 
        BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3) 
        DO 680 J=1,3 
        P(I,J)=P(I,J)+GA*(GA*BEP/(1.+GA)+P(I,4))*BE(J) 
  680   CONTINUE 
        P(I,4)=GA*(P(I,4)+BEP) 
  690   CONTINUE 
      ENDIF 
 
C...Fill in position of decay vertex. 
      DO 710 I=NSAV+1,N 
      DO 700 J=1,4 
      V(I,J)=VDCY(J) 
  700 CONTINUE 
      V(I,5)=0. 
  710 CONTINUE 
 
C...Set up for parton shower evolution from jets. 
      IF(MSTJ(23).GE.1.AND.MMAT.EQ.4.AND.K(NSAV+1,2).EQ.21) THEN 
        K(NSAV+1,1)=3 
        K(NSAV+2,1)=3 
        K(NSAV+3,1)=3 
        K(NSAV+1,4)=MSTU(5)*(NSAV+2) 
        K(NSAV+1,5)=MSTU(5)*(NSAV+3) 
        K(NSAV+2,4)=MSTU(5)*(NSAV+3) 
        K(NSAV+2,5)=MSTU(5)*(NSAV+1) 
        K(NSAV+3,4)=MSTU(5)*(NSAV+1) 
        K(NSAV+3,5)=MSTU(5)*(NSAV+2) 
        MSTJ(92)=-(NSAV+1) 
      ELSEIF(MSTJ(23).GE.1.AND.MMAT.EQ.4) THEN 
        K(NSAV+2,1)=3 
        K(NSAV+3,1)=3 
        K(NSAV+2,4)=MSTU(5)*(NSAV+3) 
        K(NSAV+2,5)=MSTU(5)*(NSAV+3) 
        K(NSAV+3,4)=MSTU(5)*(NSAV+2) 
        K(NSAV+3,5)=MSTU(5)*(NSAV+2) 
        MSTJ(92)=NSAV+2 
      ELSEIF(MSTJ(23).GE.1.AND.(MMAT.EQ.32.OR.MMAT.EQ.44.OR.MMAT.EQ.46) 
     &.AND.IABS(K(NSAV+1,2)).LE.10.AND.IABS(K(NSAV+2,2)).LE.10) THEN 
        K(NSAV+1,1)=3 
        K(NSAV+2,1)=3 
        K(NSAV+1,4)=MSTU(5)*(NSAV+2) 
        K(NSAV+1,5)=MSTU(5)*(NSAV+2) 
        K(NSAV+2,4)=MSTU(5)*(NSAV+1) 
        K(NSAV+2,5)=MSTU(5)*(NSAV+1) 
        MSTJ(92)=NSAV+1 
      ELSEIF(MSTJ(23).GE.1.AND.(MMAT.EQ.32.OR.MMAT.EQ.44.OR.MMAT.EQ.46) 
     &.AND.IABS(K(NSAV+1,2)).LE.20.AND.IABS(K(NSAV+2,2)).LE.20) THEN 
        MSTJ(92)=NSAV+1 
      ELSEIF(MSTJ(23).GE.1.AND.MMAT.EQ.33.AND.IABS(K(NSAV+2,2)).EQ.21) 
     &THEN 
        K(NSAV+1,1)=3 
        K(NSAV+2,1)=3 
        K(NSAV+3,1)=3 
        KCP=LUCOMP(K(NSAV+1,2)) 
        KQP=KCHG(KCP,2)*ISIGN(1,K(NSAV+1,2)) 
        JCON=4 
        IF(KQP.LT.0) JCON=5 
        K(NSAV+1,JCON)=MSTU(5)*(NSAV+2) 
        K(NSAV+2,9-JCON)=MSTU(5)*(NSAV+1) 
        K(NSAV+2,JCON)=MSTU(5)*(NSAV+3) 
        K(NSAV+3,9-JCON)=MSTU(5)*(NSAV+2) 
        MSTJ(92)=NSAV+1 
      ELSEIF(MSTJ(23).GE.1.AND.MMAT.EQ.33) THEN 
        K(NSAV+1,1)=3 
        K(NSAV+3,1)=3 
        K(NSAV+1,4)=MSTU(5)*(NSAV+3) 
        K(NSAV+1,5)=MSTU(5)*(NSAV+3) 
        K(NSAV+3,4)=MSTU(5)*(NSAV+1) 
        K(NSAV+3,5)=MSTU(5)*(NSAV+1) 
        MSTJ(92)=NSAV+1 
 
C...Set up for parton shower evolution in t -> W + b. 
      ELSEIF(MSTJ(27).GE.1.AND.MMAT.EQ.45.AND.ND.EQ.3) THEN 
        K(NSAV+2,1)=3 
        K(NSAV+3,1)=3 
        K(NSAV+2,4)=MSTU(5)*(NSAV+3) 
        K(NSAV+2,5)=MSTU(5)*(NSAV+3) 
        K(NSAV+3,4)=MSTU(5)*(NSAV+2) 
        K(NSAV+3,5)=MSTU(5)*(NSAV+2) 
        MSTJ(92)=NSAV+1 
      ENDIF 
 
C...Mark decayed particle; special option for B-B~ mixing. 
      IF(K(IP,1).EQ.5) K(IP,1)=15 
      IF(K(IP,1).LE.10) K(IP,1)=11 
      IF(MMIX.EQ.1.AND.MSTJ(26).EQ.2.AND.K(IP,1).EQ.11) K(IP,1)=12 
      K(IP,4)=NSAV+1 
      K(IP,5)=N 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUKFDI(KFL1,KFL2,KFL3,KF) 
 
C...Purpose: to generate a new flavour pair and combine off a hadron. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT1/,/LUDAT2/ 
 
C...Default flavour values. Input consistency checks. 
      KF1A=IABS(KFL1) 
      KF2A=IABS(KFL2) 
      KFL3=0 
      KF=0 
      IF(KF1A.EQ.0) RETURN 
      IF(KF2A.NE.0) THEN 
        IF(KF1A.LE.10.AND.KF2A.LE.10.AND.KFL1*KFL2.GT.0) RETURN 
        IF(KF1A.GT.10.AND.KF2A.GT.10) RETURN 
        IF((KF1A.GT.10.OR.KF2A.GT.10).AND.KFL1*KFL2.LT.0) RETURN 
      ENDIF 
 
C...Check if tabulated flavour probabilities are to be used. 
      IF(MSTJ(15).EQ.1) THEN 
        KTAB1=-1 
        IF(KF1A.GE.1.AND.KF1A.LE.6) KTAB1=KF1A 
        KFL1A=MOD(KF1A/1000,10) 
        KFL1B=MOD(KF1A/100,10) 
        KFL1S=MOD(KF1A,10) 
        IF(KFL1A.GE.1.AND.KFL1A.LE.4.AND.KFL1B.GE.1.AND.KFL1B.LE.4) 
     &  KTAB1=6+KFL1A*(KFL1A-2)+2*KFL1B+(KFL1S-1)/2 
        IF(KFL1A.GE.1.AND.KFL1A.LE.4.AND.KFL1A.EQ.KFL1B) KTAB1=KTAB1-1 
        IF(KF1A.GE.1.AND.KF1A.LE.6) KFL1A=KF1A 
        KTAB2=0 
        IF(KF2A.NE.0) THEN 
          KTAB2=-1 
          IF(KF2A.GE.1.AND.KF2A.LE.6) KTAB2=KF2A 
          KFL2A=MOD(KF2A/1000,10) 
          KFL2B=MOD(KF2A/100,10) 
          KFL2S=MOD(KF2A,10) 
          IF(KFL2A.GE.1.AND.KFL2A.LE.4.AND.KFL2B.GE.1.AND.KFL2B.LE.4) 
     &    KTAB2=6+KFL2A*(KFL2A-2)+2*KFL2B+(KFL2S-1)/2 
          IF(KFL2A.GE.1.AND.KFL2A.LE.4.AND.KFL2A.EQ.KFL2B) KTAB2=KTAB2-1 
        ENDIF 
        IF(KTAB1.GE.0.AND.KTAB2.GE.0) GOTO 150 
      ENDIF 
 
C...Parameters and breaking diquark parameter combinations. 
  100 PAR2=PARJ(2) 
      PAR3=PARJ(3) 
      PAR4=3.*PARJ(4) 
      IF(MSTJ(12).GE.2) THEN 
        PAR3M=SQRT(PARJ(3)) 
        PAR4M=1./(3.*SQRT(PARJ(4))) 
        PARDM=PARJ(7)/(PARJ(7)+PAR3M*PARJ(6)) 
        PARS0=PARJ(5)*(2.+(1.+PAR2*PAR3M*PARJ(7))*(1.+PAR4M)) 
        PARS1=PARJ(7)*PARS0/(2.*PAR3M)+PARJ(5)*(PARJ(6)*(1.+PAR4M)+ 
     &  PAR2*PAR3M*PARJ(6)*PARJ(7)) 
        PARS2=PARJ(5)*2.*PARJ(6)*PARJ(7)*(PAR2*PARJ(7)+(1.+PAR4M)/PAR3M) 
        PARSM=MAX(PARS0,PARS1,PARS2) 
        PAR4=PAR4*(1.+PARSM)/(1.+PARSM/(3.*PAR4M)) 
      ENDIF 
 
C...Choice of whether to generate meson or baryon. 
  110 MBARY=0 
      KFDA=0 
      IF(KF1A.LE.10) THEN 
        IF(KF2A.EQ.0.AND.MSTJ(12).GE.1.AND.(1.+PARJ(1))*RLU(0).GT.1.) 
     &  MBARY=1 
        IF(KF2A.GT.10) MBARY=2 
        IF(KF2A.GT.10.AND.KF2A.LE.10000) KFDA=KF2A 
      ELSE 
        MBARY=2 
        IF(KF1A.LE.10000) KFDA=KF1A 
      ENDIF 
 
C...Possibility of process diquark -> meson + new diquark. 
      IF(KFDA.NE.0.AND.MSTJ(12).GE.2) THEN 
        KFLDA=MOD(KFDA/1000,10) 
        KFLDB=MOD(KFDA/100,10) 
        KFLDS=MOD(KFDA,10) 
        WTDQ=PARS0 
        IF(MAX(KFLDA,KFLDB).EQ.3) WTDQ=PARS1 
        IF(MIN(KFLDA,KFLDB).EQ.3) WTDQ=PARS2 
        IF(KFLDS.EQ.1) WTDQ=WTDQ/(3.*PAR4M) 
        IF((1.+WTDQ)*RLU(0).GT.1.) MBARY=-1 
        IF(MBARY.EQ.-1.AND.KF2A.NE.0) RETURN 
      ENDIF 
 
C...Flavour for meson, possibly with new flavour. 
      IF(MBARY.LE.0) THEN 
        KFS=ISIGN(1,KFL1) 
        IF(MBARY.EQ.0) THEN 
          IF(KF2A.EQ.0) KFL3=ISIGN(1+INT((2.+PAR2)*RLU(0)),-KFL1) 
          KFLA=MAX(KF1A,KF2A+IABS(KFL3)) 
          KFLB=MIN(KF1A,KF2A+IABS(KFL3)) 
          IF(KFLA.NE.KF1A) KFS=-KFS 
 
C...Splitting of diquark into meson plus new diquark. 
        ELSE 
          KFL1A=MOD(KF1A/1000,10) 
          KFL1B=MOD(KF1A/100,10) 
  120     KFL1D=KFL1A+INT(RLU(0)+0.5)*(KFL1B-KFL1A) 
          KFL1E=KFL1A+KFL1B-KFL1D 
          IF((KFL1D.EQ.3.AND.RLU(0).GT.PARDM).OR.(KFL1E.EQ.3.AND. 
     &    RLU(0).LT.PARDM)) THEN 
            KFL1D=KFL1A+KFL1B-KFL1D 
            KFL1E=KFL1A+KFL1B-KFL1E 
          ENDIF 
          KFL3A=1+INT((2.+PAR2*PAR3M*PARJ(7))*RLU(0)) 
          IF((KFL1E.NE.KFL3A.AND.RLU(0).GT.(1.+PAR4M)/MAX(2.,1.+PAR4M)) 
     &    .OR.(KFL1E.EQ.KFL3A.AND.RLU(0).GT.2./MAX(2.,1.+PAR4M))) 
     &    GOTO 120 
          KFLDS=3 
          IF(KFL1E.NE.KFL3A) KFLDS=2*INT(RLU(0)+1./(1.+PAR4M))+1 
          KFL3=ISIGN(10000+1000*MAX(KFL1E,KFL3A)+100*MIN(KFL1E,KFL3A)+ 
     &    KFLDS,-KFL1) 
          KFLA=MAX(KFL1D,KFL3A) 
          KFLB=MIN(KFL1D,KFL3A) 
          IF(KFLA.NE.KFL1D) KFS=-KFS 
        ENDIF 
 
C...Form meson, with spin and flavour mixing for diagonal states. 
        IF(KFLA.LE.2) KMUL=INT(PARJ(11)+RLU(0)) 
        IF(KFLA.EQ.3) KMUL=INT(PARJ(12)+RLU(0)) 
        IF(KFLA.GE.4) KMUL=INT(PARJ(13)+RLU(0)) 
        IF(KMUL.EQ.0.AND.PARJ(14).GT.0.) THEN 
          IF(RLU(0).LT.PARJ(14)) KMUL=2 
        ELSEIF(KMUL.EQ.1.AND.PARJ(15)+PARJ(16)+PARJ(17).GT.0.) THEN 
          RMUL=RLU(0) 
          IF(RMUL.LT.PARJ(15)) KMUL=3 
          IF(KMUL.EQ.1.AND.RMUL.LT.PARJ(15)+PARJ(16)) KMUL=4 
          IF(KMUL.EQ.1.AND.RMUL.LT.PARJ(15)+PARJ(16)+PARJ(17)) KMUL=5 
        ENDIF 
        KFLS=3 
        IF(KMUL.EQ.0.OR.KMUL.EQ.3) KFLS=1 
        IF(KMUL.EQ.5) KFLS=5 
        IF(KFLA.NE.KFLB) THEN 
          KF=(100*KFLA+10*KFLB+KFLS)*KFS*(-1)**KFLA 
        ELSE 
          RMIX=RLU(0) 
          IMIX=2*KFLA+10*KMUL 
          IF(KFLA.LE.3) KF=110*(1+INT(RMIX+PARF(IMIX-1))+ 
     &    INT(RMIX+PARF(IMIX)))+KFLS 
          IF(KFLA.GE.4) KF=110*KFLA+KFLS 
        ENDIF 
        IF(KMUL.EQ.2.OR.KMUL.EQ.3) KF=KF+ISIGN(10000,KF) 
        IF(KMUL.EQ.4) KF=KF+ISIGN(20000,KF) 
 
C...Optional extra suppression of eta and eta'. 
        IF(KF.EQ.221) THEN 
          IF(RLU(0).GT.PARJ(25)) GOTO 110 
        ELSEIF(KF.EQ.331) THEN 
          IF(RLU(0).GT.PARJ(26)) GOTO 110 
        ENDIF 
 
C...Generate diquark flavour. 
      ELSE 
  130   IF(KF1A.LE.10.AND.KF2A.EQ.0) THEN 
          KFLA=KF1A 
  140     KFLB=1+INT((2.+PAR2*PAR3)*RLU(0)) 
          KFLC=1+INT((2.+PAR2*PAR3)*RLU(0)) 
          KFLDS=1 
          IF(KFLB.GE.KFLC) KFLDS=3 
          IF(KFLDS.EQ.1.AND.PAR4*RLU(0).GT.1.) GOTO 140 
          IF(KFLDS.EQ.3.AND.PAR4.LT.RLU(0)) GOTO 140 
          KFL3=ISIGN(1000*MAX(KFLB,KFLC)+100*MIN(KFLB,KFLC)+KFLDS,KFL1) 
 
C...Take diquark flavour from input. 
        ELSEIF(KF1A.LE.10) THEN 
          KFLA=KF1A 
          KFLB=MOD(KF2A/1000,10) 
          KFLC=MOD(KF2A/100,10) 
          KFLDS=MOD(KF2A,10) 
 
C...Generate (or take from input) quark to go with diquark. 
        ELSE 
          IF(KF2A.EQ.0) KFL3=ISIGN(1+INT((2.+PAR2)*RLU(0)),KFL1) 
          KFLA=KF2A+IABS(KFL3) 
          KFLB=MOD(KF1A/1000,10) 
          KFLC=MOD(KF1A/100,10) 
          KFLDS=MOD(KF1A,10) 
        ENDIF 
 
C...SU(6) factors for formation of baryon. Try again if fails. 
        KBARY=KFLDS 
        IF(KFLDS.EQ.3.AND.KFLB.NE.KFLC) KBARY=5 
        IF(KFLA.NE.KFLB.AND.KFLA.NE.KFLC) KBARY=KBARY+1 
        WT=PARF(60+KBARY)+PARJ(18)*PARF(70+KBARY) 
        IF(MBARY.EQ.1.AND.MSTJ(12).GE.2) THEN 
          WTDQ=PARS0 
          IF(MAX(KFLB,KFLC).EQ.3) WTDQ=PARS1 
          IF(MIN(KFLB,KFLC).EQ.3) WTDQ=PARS2 
          IF(KFLDS.EQ.1) WTDQ=WTDQ/(3.*PAR4M) 
          IF(KFLDS.EQ.1) WT=WT*(1.+WTDQ)/(1.+PARSM/(3.*PAR4M)) 
          IF(KFLDS.EQ.3) WT=WT*(1.+WTDQ)/(1.+PARSM) 
        ENDIF 
        IF(KF2A.EQ.0.AND.WT.LT.RLU(0)) GOTO 130 
 
C...Form baryon. Distinguish Lambda- and Sigmalike baryons. 
        KFLD=MAX(KFLA,KFLB,KFLC) 
        KFLF=MIN(KFLA,KFLB,KFLC) 
        KFLE=KFLA+KFLB+KFLC-KFLD-KFLF 
        KFLS=2 
        IF((PARF(60+KBARY)+PARJ(18)*PARF(70+KBARY))*RLU(0).GT. 
     &  PARF(60+KBARY)) KFLS=4 
        KFLL=0 
        IF(KFLS.EQ.2.AND.KFLD.GT.KFLE.AND.KFLE.GT.KFLF) THEN 
          IF(KFLDS.EQ.1.AND.KFLA.EQ.KFLD) KFLL=1 
          IF(KFLDS.EQ.1.AND.KFLA.NE.KFLD) KFLL=INT(0.25+RLU(0)) 
          IF(KFLDS.EQ.3.AND.KFLA.NE.KFLD) KFLL=INT(0.75+RLU(0)) 
        ENDIF 
        IF(KFLL.EQ.0) KF=ISIGN(1000*KFLD+100*KFLE+10*KFLF+KFLS,KFL1) 
        IF(KFLL.EQ.1) KF=ISIGN(1000*KFLD+100*KFLF+10*KFLE+KFLS,KFL1) 
      ENDIF 
      RETURN 
 
C...Use tabulated probabilities to select new flavour and hadron. 
  150 IF(KTAB2.EQ.0.AND.MSTJ(12).LE.0) THEN 
        KT3L=1 
        KT3U=6 
      ELSEIF(KTAB2.EQ.0.AND.KTAB1.GE.7.AND.MSTJ(12).LE.1) THEN 
        KT3L=1 
        KT3U=6 
      ELSEIF(KTAB2.EQ.0) THEN 
        KT3L=1 
        KT3U=22 
      ELSE 
        KT3L=KTAB2 
        KT3U=KTAB2 
      ENDIF 
      RFL=0. 
      DO 170 KTS=0,2 
      DO 160 KT3=KT3L,KT3U 
      RFL=RFL+PARF(120+80*KTAB1+25*KTS+KT3) 
  160 CONTINUE 
  170 CONTINUE 
      RFL=RLU(0)*RFL 
      DO 190 KTS=0,2 
      KTABS=KTS 
      DO 180 KT3=KT3L,KT3U 
      KTAB3=KT3 
      RFL=RFL-PARF(120+80*KTAB1+25*KTS+KT3) 
      IF(RFL.LE.0.) GOTO 200 
  180 CONTINUE 
  190 CONTINUE 
  200 CONTINUE 
 
C...Reconstruct flavour of produced quark/diquark. 
      IF(KTAB3.LE.6) THEN 
        KFL3A=KTAB3 
        KFL3B=0 
        KFL3=ISIGN(KFL3A,KFL1*(2*KTAB1-13)) 
      ELSE 
        KFL3A=1 
        IF(KTAB3.GE.8) KFL3A=2 
        IF(KTAB3.GE.11) KFL3A=3 
        IF(KTAB3.GE.16) KFL3A=4 
        KFL3B=(KTAB3-6-KFL3A*(KFL3A-2))/2 
        KFL3=1000*KFL3A+100*KFL3B+1 
        IF(KFL3A.EQ.KFL3B.OR.KTAB3.NE.6+KFL3A*(KFL3A-2)+2*KFL3B) KFL3= 
     &  KFL3+2 
        KFL3=ISIGN(KFL3,KFL1*(13-2*KTAB1)) 
      ENDIF 
 
C...Reconstruct meson code. 
      IF(KFL3A.EQ.KFL1A.AND.KFL3B.EQ.KFL1B.AND.(KFL3A.LE.3.OR. 
     &KFL3B.NE.0)) THEN 
        RFL=RLU(0)*(PARF(143+80*KTAB1+25*KTABS)+PARF(144+80*KTAB1+ 
     &  25*KTABS)+PARF(145+80*KTAB1+25*KTABS)) 
        KF=110+2*KTABS+1 
        IF(RFL.GT.PARF(143+80*KTAB1+25*KTABS)) KF=220+2*KTABS+1 
        IF(RFL.GT.PARF(143+80*KTAB1+25*KTABS)+PARF(144+80*KTAB1+ 
     &  25*KTABS)) KF=330+2*KTABS+1 
      ELSEIF(KTAB1.LE.6.AND.KTAB3.LE.6) THEN 
        KFLA=MAX(KTAB1,KTAB3) 
        KFLB=MIN(KTAB1,KTAB3) 
        KFS=ISIGN(1,KFL1) 
        IF(KFLA.NE.KF1A) KFS=-KFS 
        KF=(100*KFLA+10*KFLB+2*KTABS+1)*KFS*(-1)**KFLA 
      ELSEIF(KTAB1.GE.7.AND.KTAB3.GE.7) THEN 
        KFS=ISIGN(1,KFL1) 
        IF(KFL1A.EQ.KFL3A) THEN 
          KFLA=MAX(KFL1B,KFL3B) 
          KFLB=MIN(KFL1B,KFL3B) 
          IF(KFLA.NE.KFL1B) KFS=-KFS 
        ELSEIF(KFL1A.EQ.KFL3B) THEN 
          KFLA=KFL3A 
          KFLB=KFL1B 
          KFS=-KFS 
        ELSEIF(KFL1B.EQ.KFL3A) THEN 
          KFLA=KFL1A 
          KFLB=KFL3B 
        ELSEIF(KFL1B.EQ.KFL3B) THEN 
          KFLA=MAX(KFL1A,KFL3A) 
          KFLB=MIN(KFL1A,KFL3A) 
          IF(KFLA.NE.KFL1A) KFS=-KFS 
        ELSE 
          CALL LUERRM(2,'(LUKFDI:) no matching flavours for qq -> qq') 
          GOTO 100 
        ENDIF 
        KF=(100*KFLA+10*KFLB+2*KTABS+1)*KFS*(-1)**KFLA 
 
C...Reconstruct baryon code. 
      ELSE 
        IF(KTAB1.GE.7) THEN 
          KFLA=KFL3A 
          KFLB=KFL1A 
          KFLC=KFL1B 
        ELSE 
          KFLA=KFL1A 
          KFLB=KFL3A 
          KFLC=KFL3B 
        ENDIF 
        KFLD=MAX(KFLA,KFLB,KFLC) 
        KFLF=MIN(KFLA,KFLB,KFLC) 
        KFLE=KFLA+KFLB+KFLC-KFLD-KFLF 
        IF(KTABS.EQ.0) KF=ISIGN(1000*KFLD+100*KFLF+10*KFLE+2,KFL1) 
        IF(KTABS.GE.1) KF=ISIGN(1000*KFLD+100*KFLE+10*KFLF+2*KTABS,KFL1) 
      ENDIF 
 
C...Check that constructed flavour code is an allowed one. 
      IF(KFL2.NE.0) KFL3=0 
      KC=LUCOMP(KF) 
      IF(KC.EQ.0) THEN 
        CALL LUERRM(2,'(LUKFDI:) user-defined flavour probabilities '// 
     &  'failed') 
        GOTO 100 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUPTDI(KFL,PX,PY) 
 
C...Purpose: to generate transverse momentum according to a Gaussian. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
 
C...Generate p_T and azimuthal angle, gives p_x and p_y. 
      KFLA=IABS(KFL) 
      PT=PARJ(21)*SQRT(-LOG(MAX(1E-10,RLU(0)))) 
      IF(PARJ(23).GT.RLU(0)) PT=PARJ(24)*PT 
      IF(MSTJ(91).EQ.1) PT=PARJ(22)*PT 
      IF(KFLA.EQ.0.AND.MSTJ(13).LE.0) PT=0. 
      PHI=PARU(2)*RLU(0) 
      PX=PT*COS(PHI) 
      PY=PT*SIN(PHI) 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUZDIS(KFL1,KFL2,PR,Z) 
 
C...Purpose: to generate the longitudinal splitting variable z. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT1/,/LUDAT2/ 
 
C...Check if heavy flavour fragmentation. 
      KFLA=IABS(KFL1) 
      KFLB=IABS(KFL2) 
      KFLH=KFLA 
      IF(KFLA.GE.10) KFLH=MOD(KFLA/1000,10) 
 
C...Lund symmetric scaling function: determine parameters of shape. 
      IF(MSTJ(11).EQ.1.OR.(MSTJ(11).EQ.3.AND.KFLH.LE.3).OR. 
     &MSTJ(11).GE.4) THEN 
        FA=PARJ(41) 
        IF(MSTJ(91).EQ.1) FA=PARJ(43) 
        IF(KFLB.GE.10) FA=FA+PARJ(45) 
        FBB=PARJ(42) 
        IF(MSTJ(91).EQ.1) FBB=PARJ(44) 
        FB=FBB*PR 
        FC=1. 
        IF(KFLA.GE.10) FC=FC-PARJ(45) 
        IF(KFLB.GE.10) FC=FC+PARJ(45) 
        IF(MSTJ(11).GE.4.AND.KFLH.GE.4.AND.KFLH.LE.5) THEN 
          FRED=PARJ(46) 
          IF(MSTJ(11).EQ.5.AND.KFLH.EQ.5) FRED=PARJ(47) 
          FC=FC+FRED*FBB*PARF(100+KFLH)**2 
        ELSEIF(MSTJ(11).GE.4.AND.KFLH.GE.6.AND.KFLH.LE.8) THEN 
          FRED=PARJ(46) 
          IF(MSTJ(11).EQ.5) FRED=PARJ(48) 
          FC=FC+FRED*FBB*PMAS(KFLH,1)**2 
        ENDIF 
        MC=1 
        IF(ABS(FC-1.).GT.0.01) MC=2 
 
C...Determine position of maximum. Special cases for a = 0 or a = c. 
        IF(FA.LT.0.02) THEN 
          MA=1 
          ZMAX=1. 
          IF(FC.GT.FB) ZMAX=FB/FC 
        ELSEIF(ABS(FC-FA).LT.0.01) THEN 
          MA=2 
          ZMAX=FB/(FB+FC) 
        ELSE 
          MA=3 
          ZMAX=0.5*(FB+FC-SQRT((FB-FC)**2+4.*FA*FB))/(FC-FA) 
          IF(ZMAX.GT.0.9999.AND.FB.GT.100.) ZMAX=MIN(ZMAX,1.-FA/FB) 
        ENDIF 
 
C...Subdivide z range if distribution very peaked near endpoint. 
        MMAX=2 
        IF(ZMAX.LT.0.1) THEN 
          MMAX=1 
          ZDIV=2.75*ZMAX 
          IF(MC.EQ.1) THEN 
            FINT=1.-LOG(ZDIV) 
          ELSE 
            ZDIVC=ZDIV**(1.-FC) 
            FINT=1.+(1.-1./ZDIVC)/(FC-1.) 
          ENDIF 
        ELSEIF(ZMAX.GT.0.85.AND.FB.GT.1.) THEN 
          MMAX=3 
          FSCB=SQRT(4.+(FC/FB)**2) 
          ZDIV=FSCB-1./ZMAX-(FC/FB)*LOG(ZMAX*0.5*(FSCB+FC/FB)) 
          IF(MA.GE.2) ZDIV=ZDIV+(FA/FB)*LOG(1.-ZMAX) 
          ZDIV=MIN(ZMAX,MAX(0.,ZDIV)) 
          FINT=1.+FB*(1.-ZDIV) 
        ENDIF 
 
C...Choice of z, preweighted for peaks at low or high z. 
  100   Z=RLU(0) 
        FPRE=1. 
        IF(MMAX.EQ.1) THEN 
          IF(FINT*RLU(0).LE.1.) THEN 
            Z=ZDIV*Z 
          ELSEIF(MC.EQ.1) THEN 
            Z=ZDIV**Z 
            FPRE=ZDIV/Z 
          ELSE 
            Z=1./(ZDIVC+Z*(1.-ZDIVC))**(1./(1.-FC)) 
            FPRE=(ZDIV/Z)**FC 
          ENDIF 
        ELSEIF(MMAX.EQ.3) THEN 
          IF(FINT*RLU(0).LE.1.) THEN 
            Z=ZDIV+LOG(Z)/FB 
            FPRE=EXP(FB*(Z-ZDIV)) 
          ELSE 
            Z=ZDIV+Z*(1.-ZDIV) 
          ENDIF 
        ENDIF 
 
C...Weighting according to correct formula. 
        IF(Z.LE.0..OR.Z.GE.1.) GOTO 100 
        FEXP=FC*LOG(ZMAX/Z)+FB*(1./ZMAX-1./Z) 
        IF(MA.GE.2) FEXP=FEXP+FA*LOG((1.-Z)/(1.-ZMAX)) 
        FVAL=EXP(MAX(-50.,MIN(50.,FEXP))) 
        IF(FVAL.LT.RLU(0)*FPRE) GOTO 100 
 
C...Generate z according to Field-Feynman, SLAC, (1-z)**c OR z**c. 
      ELSE 
        FC=PARJ(50+MAX(1,KFLH)) 
        IF(MSTJ(91).EQ.1) FC=PARJ(59) 
  110   Z=RLU(0) 
        IF(FC.GE.0..AND.FC.LE.1.) THEN 
          IF(FC.GT.RLU(0)) Z=1.-Z**(1./3.) 
        ELSEIF(FC.GT.-1.AND.FC.LT.0.) THEN 
          IF(-4.*FC*Z*(1.-Z)**2.LT.RLU(0)*((1.-Z)**2-FC*Z)**2) GOTO 110 
        ELSE 
          IF(FC.GT.0.) Z=1.-Z**(1./FC) 
          IF(FC.LT.0.) Z=Z**(-1./FC) 
        ENDIF 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUSHOW(IP1,IP2,QMAX) 
 
C...Purpose: to generate timelike parton showers from given partons. 
      IMPLICIT DOUBLE PRECISION(D) 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION PMTH(5,50),PS(5),PMA(4),PMSD(4),IEP(4),IPA(4), 
     &KFLA(4),KFLD(4),KFL(4),ITRY(4),ISI(4),ISL(4),DP(4),DPT(5,4), 
     &KSH(0:40),KCII(2),NIIS(2),IIIS(2,2),THEIIS(2,2),PHIIIS(2,2), 
     &ISII(2) 
 
C...Initialization of cutoff masses etc. 
      IF(MSTJ(41).LE.0.OR.(MSTJ(41).EQ.1.AND.QMAX.LE.PARJ(82)).OR. 
     &QMAX.LE.MIN(PARJ(82),PARJ(83))) RETURN 
      DO 100 IFL=0,40 
      KSH(IFL)=0 
  100 CONTINUE 
      KSH(21)=1 
      PMTH(1,21)=ULMASS(21) 
      PMTH(2,21)=SQRT(PMTH(1,21)**2+0.25*PARJ(82)**2) 
      PMTH(3,21)=2.*PMTH(2,21) 
      PMTH(4,21)=PMTH(3,21) 
      PMTH(5,21)=PMTH(3,21) 
      PMTH(1,22)=ULMASS(22) 
      PMTH(2,22)=SQRT(PMTH(1,22)**2+0.25*PARJ(83)**2) 
      PMTH(3,22)=2.*PMTH(2,22) 
      PMTH(4,22)=PMTH(3,22) 
      PMTH(5,22)=PMTH(3,22) 
      PMQTH1=PARJ(82) 
      IF(MSTJ(41).GE.2) PMQTH1=MIN(PARJ(82),PARJ(83)) 
      PMQTH2=PMTH(2,21) 
      IF(MSTJ(41).GE.2) PMQTH2=MIN(PMTH(2,21),PMTH(2,22)) 
      DO 110 IFL=1,8 
      KSH(IFL)=1 
      PMTH(1,IFL)=ULMASS(IFL) 
      PMTH(2,IFL)=SQRT(PMTH(1,IFL)**2+0.25*PMQTH1**2) 
      PMTH(3,IFL)=PMTH(2,IFL)+PMQTH2 
      PMTH(4,IFL)=SQRT(PMTH(1,IFL)**2+0.25*PARJ(82)**2)+PMTH(2,21) 
      PMTH(5,IFL)=SQRT(PMTH(1,IFL)**2+0.25*PARJ(83)**2)+PMTH(2,22) 
  110 CONTINUE 
      DO 120 IFL=11,17,2 
      IF(MSTJ(41).GE.2) KSH(IFL)=1 
      PMTH(1,IFL)=ULMASS(IFL) 
      PMTH(2,IFL)=SQRT(PMTH(1,IFL)**2+0.25*PARJ(83)**2) 
      PMTH(3,IFL)=PMTH(2,IFL)+PMTH(2,22) 
      PMTH(4,IFL)=PMTH(3,IFL) 
      PMTH(5,IFL)=PMTH(3,IFL) 
  120 CONTINUE 
      PT2MIN=MAX(0.5*PARJ(82),1.1*PARJ(81))**2 
      ALAMS=PARJ(81)**2 
      ALFM=LOG(PT2MIN/ALAMS) 
 
C...Store positions of shower initiating partons. 
      IF(IP1.GT.0.AND.IP1.LE.MIN(N,MSTU(4)-MSTU(32)).AND.IP2.EQ.0) THEN 
        NPA=1 
        IPA(1)=IP1 
      ELSEIF(MIN(IP1,IP2).GT.0.AND.MAX(IP1,IP2).LE.MIN(N,MSTU(4)- 
     &MSTU(32))) THEN 
        NPA=2 
        IPA(1)=IP1 
        IPA(2)=IP2 
      ELSEIF(IP1.GT.0.AND.IP1.LE.MIN(N,MSTU(4)-MSTU(32)).AND.IP2.LT.0 
     &.AND.IP2.GE.-3) THEN 
        NPA=IABS(IP2) 
        DO 130 I=1,NPA 
        IPA(I)=IP1+I-1 
  130   CONTINUE 
      ELSE 
        CALL LUERRM(12, 
     &  '(LUSHOW:) failed to reconstruct showering system') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
 
C...Check on phase space available for emission. 
      IREJ=0 
      DO 140 J=1,5 
      PS(J)=0. 
  140 CONTINUE 
      PM=0. 
      DO 160 I=1,NPA 
      KFLA(I)=IABS(K(IPA(I),2)) 
      PMA(I)=P(IPA(I),5) 
C...Special cutoff masses for t, l, h with variable masses.
      IFLA=KFLA(I)
      IF(KFLA(I).GE.6.AND.KFLA(I).LE.8) THEN
        IFLA=37+KFLA(I)+ISIGN(2,K(IPA(I),2))
        PMTH(1,IFLA)=PMA(I)
        PMTH(2,IFLA)=SQRT(PMTH(1,IFLA)**2+0.25*PMQTH1**2) 
        PMTH(3,IFLA)=PMTH(2,IFLA)+PMQTH2 
        PMTH(4,IFLA)=SQRT(PMTH(1,IFLA)**2+0.25*PARJ(82)**2)+PMTH(2,21) 
        PMTH(5,IFLA)=SQRT(PMTH(1,IFLA)**2+0.25*PARJ(83)**2)+PMTH(2,22) 
      ENDIF 
      IF(KFLA(I).LE.40) THEN 
        IF(KSH(KFLA(I)).EQ.1) PMA(I)=PMTH(3,IFLA)
      ENDIF 
      PM=PM+PMA(I) 
      IF(KFLA(I).GT.40) THEN 
        IREJ=IREJ+1 
      ELSE 
        IF(KSH(KFLA(I)).EQ.0.OR.PMA(I).GT.QMAX) IREJ=IREJ+1 
      ENDIF 
      DO 150 J=1,4 
      PS(J)=PS(J)+P(IPA(I),J) 
  150 CONTINUE 
  160 CONTINUE 
      IF(IREJ.EQ.NPA) RETURN 
      PS(5)=SQRT(MAX(0.,PS(4)**2-PS(1)**2-PS(2)**2-PS(3)**2)) 
      IF(NPA.EQ.1) PS(5)=PS(4) 
      IF(PS(5).LE.PM+PMQTH1) RETURN 
 
C...Check if 3-jet matrix elements to be used. 
      M3JC=0 
      IF(NPA.EQ.2.AND.MSTJ(47).GE.1) THEN 
        IF(KFLA(1).GE.1.AND.KFLA(1).LE.8.AND.KFLA(2).GE.1.AND. 
     &  KFLA(2).LE.8) M3JC=1 
        IF((KFLA(1).EQ.11.OR.KFLA(1).EQ.13.OR.KFLA(1).EQ.15.OR. 
     &  KFLA(1).EQ.17).AND.KFLA(2).EQ.KFLA(1)) M3JC=1 
        IF((KFLA(1).EQ.11.OR.KFLA(1).EQ.13.OR.KFLA(1).EQ.15.OR. 
     &  KFLA(1).EQ.17).AND.KFLA(2).EQ.KFLA(1)+1) M3JC=1 
        IF((KFLA(1).EQ.12.OR.KFLA(1).EQ.14.OR.KFLA(1).EQ.16.OR. 
     &  KFLA(1).EQ.18).AND.KFLA(2).EQ.KFLA(1)-1) M3JC=1 
        IF(MSTJ(47).EQ.2.OR.MSTJ(47).EQ.4) M3JC=1 
        M3JCM=0 
        IF(M3JC.EQ.1.AND.MSTJ(47).GE.3.AND.KFLA(1).EQ.KFLA(2)) THEN 
          M3JCM=1 
          QME=(2.*PMTH(1,KFLA(1))/PS(5))**2 
        ENDIF 
      ENDIF 
 
C...Find if interference with initial state partons. 
      MIIS=0 
      IF(MSTJ(50).GE.1.AND.MSTJ(50).LE.3.AND.NPA.EQ.2) MIIS=MSTJ(50) 
      IF(MIIS.NE.0) THEN 
        DO 180 I=1,2 
        KCII(I)=0 
        KCA=LUCOMP(KFLA(I)) 
        IF(KCA.NE.0) KCII(I)=KCHG(KCA,2)*ISIGN(1,K(IPA(I),2)) 
        NIIS(I)=0 
        IF(KCII(I).NE.0) THEN 
          DO 170 J=1,2 
          ICSI=MOD(K(IPA(I),3+J)/MSTU(5),MSTU(5)) 
          IF(ICSI.GT.0.AND.ICSI.NE.IPA(1).AND.ICSI.NE.IPA(2).AND. 
     &    (KCII(I).EQ.(-1)**(J+1).OR.KCII(I).EQ.2)) THEN 
            NIIS(I)=NIIS(I)+1 
            IIIS(I,NIIS(I))=ICSI 
          ENDIF 
  170     CONTINUE 
        ENDIF 
  180   CONTINUE 
        IF(NIIS(1)+NIIS(2).EQ.0) MIIS=0 
      ENDIF 
 
C...Boost interfering initial partons to rest frame 
C...and reconstruct their polar and azimuthal angles. 
      IF(MIIS.NE.0) THEN 
        DO 200 I=1,2 
        DO 190 J=1,5 
        K(N+I,J)=K(IPA(I),J) 
        P(N+I,J)=P(IPA(I),J) 
        V(N+I,J)=0. 
  190   CONTINUE 
  200   CONTINUE 
        DO 220 I=3,2+NIIS(1) 
        DO 210 J=1,5 
        K(N+I,J)=K(IIIS(1,I-2),J) 
        P(N+I,J)=P(IIIS(1,I-2),J) 
        V(N+I,J)=0. 
  210   CONTINUE 
  220   CONTINUE 
        DO 240 I=3+NIIS(1),2+NIIS(1)+NIIS(2) 
        DO 230 J=1,5 
        K(N+I,J)=K(IIIS(2,I-2-NIIS(1)),J) 
        P(N+I,J)=P(IIIS(2,I-2-NIIS(1)),J) 
        V(N+I,J)=0. 
  230   CONTINUE 
  240   CONTINUE 
        CALL LUDBRB(N+1,N+2+NIIS(1)+NIIS(2),0.,0.,-DBLE(PS(1)/PS(4)), 
     &  -DBLE(PS(2)/PS(4)),-DBLE(PS(3)/PS(4))) 
        PHI=ULANGL(P(N+1,1),P(N+1,2)) 
        CALL LUDBRB(N+1,N+2+NIIS(1)+NIIS(2),0.,-PHI,0D0,0D0,0D0) 
        THE=ULANGL(P(N+1,3),P(N+1,1)) 
        CALL LUDBRB(N+1,N+2+NIIS(1)+NIIS(2),-THE,0.,0D0,0D0,0D0) 
        DO 250 I=3,2+NIIS(1) 
        THEIIS(1,I-2)=ULANGL(P(N+I,3),SQRT(P(N+I,1)**2+P(N+I,2)**2)) 
        PHIIIS(1,I-2)=ULANGL(P(N+I,1),P(N+I,2)) 
  250   CONTINUE 
        DO 260 I=3+NIIS(1),2+NIIS(1)+NIIS(2) 
        THEIIS(2,I-2-NIIS(1))=PARU(1)-ULANGL(P(N+I,3), 
     &  SQRT(P(N+I,1)**2+P(N+I,2)**2)) 
        PHIIIS(2,I-2-NIIS(1))=ULANGL(P(N+I,1),P(N+I,2)) 
  260   CONTINUE 
      ENDIF 
 
C...Define imagined single initiator of shower for parton system. 
      NS=N 
      IF(N.GT.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUSHOW:) no more memory left in LUJETS') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      IF(NPA.GE.2) THEN 
        K(N+1,1)=11 
        K(N+1,2)=21 
        K(N+1,3)=0 
        K(N+1,4)=0 
        K(N+1,5)=0 
        P(N+1,1)=0. 
        P(N+1,2)=0. 
        P(N+1,3)=0. 
        P(N+1,4)=PS(5) 
        P(N+1,5)=PS(5) 
        V(N+1,5)=PS(5)**2 
        N=N+1 
      ENDIF 
 
C...Loop over partons that may branch. 
      NEP=NPA 
      IM=NS 
      IF(NPA.EQ.1) IM=NS-1 
  270 IM=IM+1 
      IF(N.GT.NS) THEN 
        IF(IM.GT.N) GOTO 510 
        KFLM=IABS(K(IM,2)) 
        IF(KFLM.GT.40) GOTO 270 
        IF(KSH(KFLM).EQ.0) GOTO 270 
        IFLM=KFLM
        IF(KFLM.GE.6.AND.KFLM.LE.8) IFLM=37+KFLM+ISIGN(2,K(IM,2)) 
        IF(P(IM,5).LT.PMTH(2,IFLM)) GOTO 270 
        IGM=K(IM,3) 
      ELSE 
        IGM=-1 
      ENDIF 
      IF(N+NEP.GT.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUSHOW:) no more memory left in LUJETS') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
 
C...Position of aunt (sister to branching parton). 
C...Origin and flavour of daughters. 
      IAU=0 
      IF(IGM.GT.0) THEN 
        IF(K(IM-1,3).EQ.IGM) IAU=IM-1 
        IF(N.GE.IM+1.AND.K(IM+1,3).EQ.IGM) IAU=IM+1 
      ENDIF 
      IF(IGM.GE.0) THEN 
        K(IM,4)=N+1 
        DO 280 I=1,NEP 
        K(N+I,3)=IM 
  280   CONTINUE 
      ELSE 
        K(N+1,3)=IPA(1) 
      ENDIF 
      IF(IGM.LE.0) THEN 
        DO 290 I=1,NEP 
        K(N+I,2)=K(IPA(I),2) 
  290   CONTINUE 
      ELSEIF(KFLM.NE.21) THEN 
        K(N+1,2)=K(IM,2) 
        K(N+2,2)=K(IM,5) 
      ELSEIF(K(IM,5).EQ.21) THEN 
        K(N+1,2)=21 
        K(N+2,2)=21 
      ELSE 
        K(N+1,2)=K(IM,5) 
        K(N+2,2)=-K(IM,5) 
      ENDIF 
 
C...Reset flags on daughers and tries made. 
      DO 300 IP=1,NEP 
      K(N+IP,1)=3 
      K(N+IP,4)=0 
      K(N+IP,5)=0 
      KFLD(IP)=IABS(K(N+IP,2)) 
      IF(KCHG(LUCOMP(KFLD(IP)),2).EQ.0) K(N+IP,1)=1 
      ITRY(IP)=0 
      ISL(IP)=0 
      ISI(IP)=0 
      IF(KFLD(IP).LE.40) THEN 
        IF(KSH(KFLD(IP)).EQ.1) ISI(IP)=1 
      ENDIF 
  300 CONTINUE 
      ISLM=0 
 
C...Maximum virtuality of daughters. 
      IF(IGM.LE.0) THEN 
        DO 310 I=1,NPA 
        IF(NPA.GE.3) P(N+I,4)=(PS(4)*P(IPA(I),4)-PS(1)*P(IPA(I),1)- 
     &  PS(2)*P(IPA(I),2)-PS(3)*P(IPA(I),3))/PS(5) 
        P(N+I,5)=MIN(QMAX,PS(5)) 
        IF(NPA.GE.3) P(N+I,5)=MIN(P(N+I,5),P(N+I,4)) 
        IF(ISI(I).EQ.0) P(N+I,5)=P(IPA(I),5) 
  310   CONTINUE 
      ELSE 
        IF(MSTJ(43).LE.2) PEM=V(IM,2) 
        IF(MSTJ(43).GE.3) PEM=P(IM,4) 
        P(N+1,5)=MIN(P(IM,5),V(IM,1)*PEM) 
        P(N+2,5)=MIN(P(IM,5),(1.-V(IM,1))*PEM) 
        IF(K(N+2,2).EQ.22) P(N+2,5)=PMTH(1,22) 
      ENDIF 
      DO 320 I=1,NEP 
      PMSD(I)=P(N+I,5) 
      IF(ISI(I).EQ.1) THEN 
        IFLD=KFLD(I)
        IF(KFLD(I).GE.6.AND.KFLD(I).LE.8) IFLD=37+KFLD(I)+
     &  ISIGN(2,K(N+I,2)) 
        IF(P(N+I,5).LE.PMTH(3,IFLD)) P(N+I,5)=PMTH(1,IFLD) 
      ENDIF 
      V(N+I,5)=P(N+I,5)**2 
  320 CONTINUE 
 
C...Choose one of the daughters for evolution. 
  330 INUM=0 
      IF(NEP.EQ.1) INUM=1 
      DO 340 I=1,NEP 
      IF(INUM.EQ.0.AND.ISL(I).EQ.1) INUM=I 
  340 CONTINUE 
      DO 350 I=1,NEP 
      IF(INUM.EQ.0.AND.ITRY(I).EQ.0.AND.ISI(I).EQ.1) THEN 
        IFLD=KFLD(I)
        IF(KFLD(I).GE.6.AND.KFLD(I).LE.8) IFLD=37+KFLD(I)+
     &  ISIGN(2,K(N+I,2)) 
        IF(P(N+I,5).GE.PMTH(2,IFLD)) INUM=I 
      ENDIF 
  350 CONTINUE 
      IF(INUM.EQ.0) THEN 
        RMAX=0. 
        DO 360 I=1,NEP 
        IF(ISI(I).EQ.1.AND.PMSD(I).GE.PMQTH2) THEN 
          RPM=P(N+I,5)/PMSD(I) 
          IFLD=KFLD(I)
          IF(KFLD(I).GE.6.AND.KFLD(I).LE.8) IFLD=37+KFLD(I)+
     &    ISIGN(2,K(N+I,2)) 
          IF(RPM.GT.RMAX.AND.P(N+I,5).GE.PMTH(2,IFLD)) THEN 
            RMAX=RPM 
            INUM=I 
          ENDIF 
        ENDIF 
  360   CONTINUE 
      ENDIF 
 
C...Store information on choice of evolving daughter. 
      INUM=MAX(1,INUM) 
      IEP(1)=N+INUM 
      DO 370 I=2,NEP 
      IEP(I)=IEP(I-1)+1 
      IF(IEP(I).GT.N+NEP) IEP(I)=N+1 
  370 CONTINUE 
      DO 380 I=1,NEP 
      KFL(I)=IABS(K(IEP(I),2)) 
  380 CONTINUE 
      ITRY(INUM)=ITRY(INUM)+1 
      IF(ITRY(INUM).GT.200) THEN 
        CALL LUERRM(14,'(LUSHOW:) caught in infinite loop') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      Z=0.5 
      IF(KFL(1).GT.40) GOTO 430 
      IF(KSH(KFL(1)).EQ.0) GOTO 430 
      IFL=KFL(1)
      IF(KFL(1).GE.6.AND.KFL(1).LE.8) IFL=37+KFL(1)+
     &ISIGN(2,K(IEP(1),2)) 
      IF(P(IEP(1),5).LT.PMTH(2,IFL)) GOTO 430 
 
C...Select side for interference with initial state partons. 
      IF(MIIS.GE.1.AND.IEP(1).LE.NS+3) THEN 
        III=IEP(1)-NS-1 
        ISII(III)=0 
        IF(IABS(KCII(III)).EQ.1.AND.NIIS(III).EQ.1) THEN 
          ISII(III)=1 
        ELSEIF(KCII(III).EQ.2.AND.NIIS(III).EQ.1) THEN 
          IF(RLU(0).GT.0.5) ISII(III)=1 
        ELSEIF(KCII(III).EQ.2.AND.NIIS(III).EQ.2) THEN 
          ISII(III)=1 
          IF(RLU(0).GT.0.5) ISII(III)=2 
        ENDIF 
      ENDIF 
 
C...Calculate allowed z range. 
      IF(NEP.EQ.1) THEN 
        PMED=PS(4) 
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN 
        PMED=P(IM,5) 
      ELSE 
        IF(INUM.EQ.1) PMED=V(IM,1)*PEM 
        IF(INUM.EQ.2) PMED=(1.-V(IM,1))*PEM 
      ENDIF 
      IF(MOD(MSTJ(43),2).EQ.1) THEN 
        ZC=PMTH(2,21)/PMED 
        ZCE=PMTH(2,22)/PMED 
      ELSE 
        ZC=0.5*(1.-SQRT(MAX(0.,1.-(2.*PMTH(2,21)/PMED)**2))) 
        IF(ZC.LT.1E-4) ZC=(PMTH(2,21)/PMED)**2 
        ZCE=0.5*(1.-SQRT(MAX(0.,1.-(2.*PMTH(2,22)/PMED)**2))) 
        IF(ZCE.LT.1E-4) ZCE=(PMTH(2,22)/PMED)**2 
      ENDIF 
      ZC=MIN(ZC,0.491) 
      ZCE=MIN(ZCE,0.491) 
      IF((MSTJ(41).EQ.1.AND.ZC.GT.0.49).OR.(MSTJ(41).GE.2.AND. 
     &MIN(ZC,ZCE).GT.0.49)) THEN 
        P(IEP(1),5)=PMTH(1,IFL) 
        V(IEP(1),5)=P(IEP(1),5)**2 
        GOTO 430 
      ENDIF 
 
C...Integral of Altarelli-Parisi z kernel for QCD. 
      IF(MSTJ(49).EQ.0.AND.KFL(1).EQ.21) THEN 
        FBR=6.*LOG((1.-ZC)/ZC)+MSTJ(45)*(0.5-ZC) 
      ELSEIF(MSTJ(49).EQ.0) THEN 
        FBR=(8./3.)*LOG((1.-ZC)/ZC) 
 
C...Integral of Altarelli-Parisi z kernel for scalar gluon. 
      ELSEIF(MSTJ(49).EQ.1.AND.KFL(1).EQ.21) THEN 
        FBR=(PARJ(87)+MSTJ(45)*PARJ(88))*(1.-2.*ZC) 
      ELSEIF(MSTJ(49).EQ.1) THEN 
        FBR=(1.-2.*ZC)/3. 
        IF(IGM.EQ.0.AND.M3JC.EQ.1) FBR=4.*FBR 
 
C...Integral of Altarelli-Parisi z kernel for Abelian vector gluon. 
      ELSEIF(KFL(1).EQ.21) THEN 
        FBR=6.*MSTJ(45)*(0.5-ZC) 
      ELSE 
        FBR=2.*LOG((1.-ZC)/ZC) 
      ENDIF 
 
C...Reset QCD probability for lepton. 
      IF(KFL(1).GE.11.AND.KFL(1).LE.18) FBR=0. 
 
C...Integral of Altarelli-Parisi kernel for photon emission. 
      IF(MSTJ(41).GE.2.AND.KFL(1).GE.1.AND.KFL(1).LE.18) THEN 
        FBRE=(KCHG(KFL(1),1)/3.)**2*2.*LOG((1.-ZCE)/ZCE) 
        IF(MSTJ(41).EQ.10) FBRE=PARJ(84)*FBRE 
      ENDIF 
 
C...Inner veto algorithm starts. Find maximum mass for evolution. 
  390 PMS=V(IEP(1),5) 
      IF(IGM.GE.0) THEN 
        PM2=0. 
        DO 400 I=2,NEP 
        PM=P(IEP(I),5) 
        IF(KFL(I).LE.40) THEN 
          IFLI=KFL(I)
          IF(KFL(I).GE.6.AND.KFL(I).LE.8) IFLI=37+KFL(I)+
     &    ISIGN(2,K(IEP(I),2)) 
          IF(KSH(KFL(I)).EQ.1) PM=PMTH(2,IFLI) 
        ENDIF 
        PM2=PM2+PM 
  400   CONTINUE 
        PMS=MIN(PMS,(P(IM,5)-PM2)**2) 
      ENDIF 
 
C...Select mass for daughter in QCD evolution. 
      B0=27./6. 
      DO 410 IFF=4,MSTJ(45) 
      IF(PMS.GT.4.*PMTH(2,IFF)**2) B0=(33.-2.*IFF)/6. 
  410 CONTINUE 
      IF(FBR.LT.1E-3) THEN 
        PMSQCD=0. 
      ELSEIF(MSTJ(44).LE.0) THEN 
        PMSQCD=PMS*EXP(MAX(-50.,LOG(RLU(0))*PARU(2)/(PARU(111)*FBR))) 
      ELSEIF(MSTJ(44).EQ.1) THEN 
        PMSQCD=4.*ALAMS*(0.25*PMS/ALAMS)**(RLU(0)**(B0/FBR)) 
      ELSE 
        PMSQCD=PMS*EXP(MAX(-50.,ALFM*B0*LOG(RLU(0))/FBR)) 
      ENDIF 
      IF(ZC.GT.0.49.OR.PMSQCD.LE.PMTH(4,IFL)**2) PMSQCD=PMTH(2,IFL)**2 
      V(IEP(1),5)=PMSQCD 
      MCE=1 
 
C...Select mass for daughter in QED evolution. 
      IF(MSTJ(41).GE.2.AND.KFL(1).GE.1.AND.KFL(1).LE.18) THEN 
        PMSQED=PMS*EXP(MAX(-50.,LOG(RLU(0))*PARU(2)/(PARU(101)*FBRE))) 
        IF(ZCE.GT.0.49.OR.PMSQED.LE.PMTH(5,IFL)**2) PMSQED= 
     &  PMTH(2,IFL)**2 
        IF(PMSQED.GT.PMSQCD) THEN 
          V(IEP(1),5)=PMSQED 
          MCE=2 
        ENDIF 
      ENDIF 
 
C...Check whether daughter mass below cutoff. 
      P(IEP(1),5)=SQRT(V(IEP(1),5)) 
      IF(P(IEP(1),5).LE.PMTH(3,IFL)) THEN 
        P(IEP(1),5)=PMTH(1,IFL) 
        V(IEP(1),5)=P(IEP(1),5)**2 
        GOTO 430 
      ENDIF 
 
C...Select z value of branching: q -> qgamma. 
      IF(MCE.EQ.2) THEN 
        Z=1.-(1.-ZCE)*(ZCE/(1.-ZCE))**RLU(0) 
        IF(1.+Z**2.LT.2.*RLU(0)) GOTO 390 
        K(IEP(1),5)=22 
 
C...Select z value of branching: q -> qg, g -> gg, g -> qqbar. 
      ELSEIF(MSTJ(49).NE.1.AND.KFL(1).NE.21) THEN 
        Z=1.-(1.-ZC)*(ZC/(1.-ZC))**RLU(0) 
        IF(1.+Z**2.LT.2.*RLU(0)) GOTO 390 
        K(IEP(1),5)=21 
      ELSEIF(MSTJ(49).EQ.0.AND.MSTJ(45)*(0.5-ZC).LT.RLU(0)*FBR) THEN 
        Z=(1.-ZC)*(ZC/(1.-ZC))**RLU(0) 
        IF(RLU(0).GT.0.5) Z=1.-Z 
        IF((1.-Z*(1.-Z))**2.LT.RLU(0)) GOTO 390 
        K(IEP(1),5)=21 
      ELSEIF(MSTJ(49).NE.1) THEN 
        Z=ZC+(1.-2.*ZC)*RLU(0) 
        IF(Z**2+(1.-Z)**2.LT.RLU(0)) GOTO 390 
        KFLB=1+INT(MSTJ(45)*RLU(0)) 
        PMQ=4.*PMTH(2,KFLB)**2/V(IEP(1),5) 
        IF(PMQ.GE.1.) GOTO 390 
        PMQ0=4.*PMTH(2,21)**2/V(IEP(1),5) 
        IF(MOD(MSTJ(43),2).EQ.0.AND.(1.+0.5*PMQ)*SQRT(1.-PMQ).LT. 
     &  RLU(0)*(1.+0.5*PMQ0)*SQRT(1.-PMQ0)) GOTO 390 
        K(IEP(1),5)=KFLB 
 
C...Ditto for scalar gluon model. 
      ELSEIF(KFL(1).NE.21) THEN 
        Z=1.-SQRT(ZC**2+RLU(0)*(1.-2.*ZC)) 
        K(IEP(1),5)=21 
      ELSEIF(RLU(0)*(PARJ(87)+MSTJ(45)*PARJ(88)).LE.PARJ(87)) THEN 
        Z=ZC+(1.-2.*ZC)*RLU(0) 
        K(IEP(1),5)=21 
      ELSE 
        Z=ZC+(1.-2.*ZC)*RLU(0) 
        KFLB=1+INT(MSTJ(45)*RLU(0)) 
        PMQ=4.*PMTH(2,KFLB)**2/V(IEP(1),5) 
        IF(PMQ.GE.1.) GOTO 390 
        K(IEP(1),5)=KFLB 
      ENDIF 
      IF(MCE.EQ.1.AND.MSTJ(44).GE.2) THEN 
        IF(Z*(1.-Z)*V(IEP(1),5).LT.PT2MIN) GOTO 390 
        IF(ALFM/LOG(V(IEP(1),5)*Z*(1.-Z)/ALAMS).LT.RLU(0)) GOTO 390 
      ENDIF 
 
C...Check if z consistent with chosen m. 
      IF(KFL(1).EQ.21) THEN 
        KFLGD1=IABS(K(IEP(1),5)) 
        KFLGD2=KFLGD1 
      ELSE 
        KFLGD1=KFL(1) 
        KFLGD2=IABS(K(IEP(1),5)) 
      ENDIF 
      IF(NEP.EQ.1) THEN 
        PED=PS(4) 
      ELSEIF(NEP.GE.3) THEN 
        PED=P(IEP(1),4) 
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN 
        PED=0.5*(V(IM,5)+V(IEP(1),5)-PM2**2)/P(IM,5) 
      ELSE 
        IF(IEP(1).EQ.N+1) PED=V(IM,1)*PEM 
        IF(IEP(1).EQ.N+2) PED=(1.-V(IM,1))*PEM 
      ENDIF 
      IF(MOD(MSTJ(43),2).EQ.1) THEN 
        IFLGD1=KFLGD1
        IF(KFLGD1.GE.6.AND.KFLGD1.LE.8) IFLGD1=IFL
        PMQTH3=0.5*PARJ(82) 
        IF(KFLGD2.EQ.22) PMQTH3=0.5*PARJ(83) 
        PMQ1=(PMTH(1,IFLGD1)**2+PMQTH3**2)/V(IEP(1),5) 
        PMQ2=(PMTH(1,KFLGD2)**2+PMQTH3**2)/V(IEP(1),5) 
        ZD=SQRT(MAX(0.,(1.-V(IEP(1),5)/PED**2)*((1.-PMQ1-PMQ2)**2- 
     &  4.*PMQ1*PMQ2))) 
        ZH=1.+PMQ1-PMQ2 
      ELSE 
        ZD=SQRT(MAX(0.,1.-V(IEP(1),5)/PED**2)) 
        ZH=1. 
      ENDIF 
      ZL=0.5*(ZH-ZD) 
      ZU=0.5*(ZH+ZD) 
      IF(Z.LT.ZL.OR.Z.GT.ZU) GOTO 390 
      IF(KFL(1).EQ.21) V(IEP(1),3)=LOG(ZU*(1.-ZL)/MAX(1E-20,ZL* 
     &(1.-ZU))) 
      IF(KFL(1).NE.21) V(IEP(1),3)=LOG((1.-ZL)/MAX(1E-10,1.-ZU)) 
 
C...Width suppression for q -> q + g.
      IF(MSTJ(40).NE.0.AND.KFL(1).NE.21) THEN
        IF(IGM.EQ.0) THEN
          EGLU=0.5*PS(5)*(1.-Z)*(1.+V(IEP(1),5)/V(NS+1,5))
        ELSE
          EGLU=PMED*(1.-Z)
        ENDIF
        CHI=PARJ(89)**2/(PARJ(89)**2+EGLU**2)
        IF(MSTJ(40).EQ.1) THEN
          IF(CHI.LT.RLU(0)) GOTO 390  
        ELSEIF(MSTJ(40).EQ.2) THEN
          IF(1.-CHI.LT.RLU(0)) GOTO 390
        ENDIF
      ENDIF
 
C...Three-jet matrix element correction. 
      IF(IGM.EQ.0.AND.M3JC.EQ.1) THEN 
        X1=Z*(1.+V(IEP(1),5)/V(NS+1,5)) 
        X2=1.-V(IEP(1),5)/V(NS+1,5) 
        X3=(1.-X1)+(1.-X2) 
        IF(MCE.EQ.2) THEN 
          KI1=K(IPA(INUM),2) 
          KI2=K(IPA(3-INUM),2) 
          QF1=KCHG(IABS(KI1),1)*ISIGN(1,KI1)/3. 
          QF2=KCHG(IABS(KI2),1)*ISIGN(1,KI2)/3. 
          WSHOW=QF1**2*(1.-X1)/X3*(1.+(X1/(2.-X2))**2)+ 
     &    QF2**2*(1.-X2)/X3*(1.+(X2/(2.-X1))**2) 
          WME=(QF1*(1.-X1)/X3-QF2*(1.-X2)/X3)**2*(X1**2+X2**2) 
        ELSEIF(MSTJ(49).NE.1) THEN 
          WSHOW=1.+(1.-X1)/X3*(X1/(2.-X2))**2+ 
     &    (1.-X2)/X3*(X2/(2.-X1))**2 
          WME=X1**2+X2**2 
          IF(M3JCM.EQ.1) WME=WME-QME*X3-0.5*QME**2- 
     &    (0.5*QME+0.25*QME**2)*((1.-X2)/MAX(1E-7,1.-X1)+
     &    (1.-X1)/MAX(1E-7,1.-X2)) 
        ELSE 
          WSHOW=4.*X3*((1.-X1)/(2.-X2)**2+(1.-X2)/(2.-X1)**2) 
          WME=X3**2 
          IF(MSTJ(102).GE.2) WME=X3**2-2.*(1.+X3)*(1.-X1)*(1.-X2)* 
     &    PARJ(171) 
        ENDIF 
        IF(WME.LT.RLU(0)*WSHOW) GOTO 390 
 
C...Impose angular ordering by rejection of nonordered emission. 
      ELSEIF(MCE.EQ.1.AND.IGM.GT.0.AND.MSTJ(42).GE.2) THEN 
        MAOM=1 
        ZM=V(IM,1) 
        IF(IEP(1).EQ.N+2) ZM=1.-V(IM,1) 
        THE2ID=Z*(1.-Z)*(ZM*P(IM,4))**2/V(IEP(1),5) 
        IAOM=IM 
  420   IF(K(IAOM,5).EQ.22) THEN 
          IAOM=K(IAOM,3) 
          IF(K(IAOM,3).LE.NS) MAOM=0 
          IF(MAOM.EQ.1) GOTO 420 
        ENDIF 
        IF(MAOM.EQ.1) THEN 
          THE2IM=V(IAOM,1)*(1.-V(IAOM,1))*P(IAOM,4)**2/V(IAOM,5) 
          IF(THE2ID.LT.THE2IM) GOTO 390 
        ENDIF 
      ENDIF 
 
C...Impose user-defined maximum angle at first branching. 
      IF(MSTJ(48).EQ.1) THEN 
        IF(NEP.EQ.1.AND.IM.EQ.NS) THEN 
          THE2ID=Z*(1.-Z)*PS(4)**2/V(IEP(1),5) 
          IF(THE2ID.LT.1./PARJ(85)**2) GOTO 390 
        ELSEIF(NEP.EQ.2.AND.IEP(1).EQ.NS+2) THEN 
          THE2ID=Z*(1.-Z)*(0.5*P(IM,4))**2/V(IEP(1),5) 
          IF(THE2ID.LT.1./PARJ(85)**2) GOTO 390 
        ELSEIF(NEP.EQ.2.AND.IEP(1).EQ.NS+3) THEN 
          THE2ID=Z*(1.-Z)*(0.5*P(IM,4))**2/V(IEP(1),5) 
          IF(THE2ID.LT.1./PARJ(86)**2) GOTO 390 
        ENDIF 
      ENDIF 
 
C...Impose angular constraint in first branching from interference 
C...with initial state partons. 
      IF(MIIS.GE.2.AND.IEP(1).LE.NS+3) THEN 
        THE2D=MAX((1.-Z)/Z,Z/(1.-Z))*V(IEP(1),5)/(0.5*P(IM,4))**2 
        IF(IEP(1).EQ.NS+2.AND.ISII(1).GE.1) THEN 
          IF(THE2D.GT.THEIIS(1,ISII(1))**2) GOTO 390 
        ELSEIF(IEP(1).EQ.NS+3.AND.ISII(2).GE.1) THEN 
          IF(THE2D.GT.THEIIS(2,ISII(2))**2) GOTO 390 
        ENDIF 
      ENDIF 
 
C...End of inner veto algorithm. Check if only one leg evolved so far. 
  430 V(IEP(1),1)=Z 
      ISL(1)=0 
      ISL(2)=0 
      IF(NEP.EQ.1) GOTO 460 
      IF(NEP.EQ.2.AND.P(IEP(1),5)+P(IEP(2),5).GE.P(IM,5)) GOTO 330 
      DO 440 I=1,NEP 
      IF(ITRY(I).EQ.0.AND.KFLD(I).LE.40) THEN 
        IF(KSH(KFLD(I)).EQ.1) THEN 
          IFLD=KFLD(I)
          IF(KFLD(I).GE.6.AND.KFLD(I).LE.8) IFLD=37+KFLD(I)+
     &    ISIGN(2,K(N+I,2)) 
          IF(P(N+I,5).GE.PMTH(2,IFLD)) GOTO 330 
        ENDIF 
      ENDIF 
  440 CONTINUE 
 
C...Check if chosen multiplet m1,m2,z1,z2 is physical. 
      IF(NEP.EQ.3) THEN 
        PA1S=(P(N+1,4)+P(N+1,5))*(P(N+1,4)-P(N+1,5)) 
        PA2S=(P(N+2,4)+P(N+2,5))*(P(N+2,4)-P(N+2,5)) 
        PA3S=(P(N+3,4)+P(N+3,5))*(P(N+3,4)-P(N+3,5)) 
        PTS=0.25*(2.*PA1S*PA2S+2.*PA1S*PA3S+2.*PA2S*PA3S- 
     &  PA1S**2-PA2S**2-PA3S**2)/PA1S 
        IF(PTS.LE.0.) GOTO 330 
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2.OR.MOD(MSTJ(43),2).EQ.0) THEN 
        DO 450 I1=N+1,N+2 
        KFLDA=IABS(K(I1,2)) 
        IF(KFLDA.GT.40) GOTO 450 
        IF(KSH(KFLDA).EQ.0) GOTO 450 
        IFLDA=KFLDA 
        IF(KFLDA.GE.6.AND.KFLDA.LE.8) IFLDA=37+KFLDA+
     &  ISIGN(2,K(I1,2)) 
        IF(P(I1,5).LT.PMTH(2,IFLDA)) GOTO 450 
        IF(KFLDA.EQ.21) THEN 
          KFLGD1=IABS(K(I1,5)) 
          KFLGD2=KFLGD1 
        ELSE 
          KFLGD1=KFLDA 
          KFLGD2=IABS(K(I1,5)) 
        ENDIF 
        I2=2*N+3-I1 
        IF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN 
          PED=0.5*(V(IM,5)+V(I1,5)-V(I2,5))/P(IM,5) 
        ELSE 
          IF(I1.EQ.N+1) ZM=V(IM,1) 
          IF(I1.EQ.N+2) ZM=1.-V(IM,1) 
          PML=SQRT((V(IM,5)-V(N+1,5)-V(N+2,5))**2- 
     &    4.*V(N+1,5)*V(N+2,5)) 
          PED=PEM*(0.5*(V(IM,5)-PML+V(I1,5)-V(I2,5))+PML*ZM)/V(IM,5) 
        ENDIF 
        IF(MOD(MSTJ(43),2).EQ.1) THEN 
          PMQTH3=0.5*PARJ(82) 
          IF(KFLGD2.EQ.22) PMQTH3=0.5*PARJ(83) 
          IFLGD1=KFLGD1
          IF(KFLGD1.GE.6.AND.KFLGD1.LE.8) IFLGD1=IFLDA
          PMQ1=(PMTH(1,IFLGD1)**2+PMQTH3**2)/V(I1,5) 
          PMQ2=(PMTH(1,KFLGD2)**2+PMQTH3**2)/V(I1,5) 
          ZD=SQRT(MAX(0.,(1.-V(I1,5)/PED**2)*((1.-PMQ1-PMQ2)**2- 
     &    4.*PMQ1*PMQ2))) 
          ZH=1.+PMQ1-PMQ2 
        ELSE 
          ZD=SQRT(MAX(0.,1.-V(I1,5)/PED**2)) 
          ZH=1. 
        ENDIF 
        ZL=0.5*(ZH-ZD) 
        ZU=0.5*(ZH+ZD) 
        IF(I1.EQ.N+1.AND.(V(I1,1).LT.ZL.OR.V(I1,1).GT.ZU)) ISL(1)=1 
        IF(I1.EQ.N+2.AND.(V(I1,1).LT.ZL.OR.V(I1,1).GT.ZU)) ISL(2)=1 
        IF(KFLDA.EQ.21) V(I1,4)=LOG(ZU*(1.-ZL)/MAX(1E-20,ZL*(1.-ZU))) 
        IF(KFLDA.NE.21) V(I1,4)=LOG((1.-ZL)/MAX(1E-10,1.-ZU)) 
  450   CONTINUE 
        IF(ISL(1).EQ.1.AND.ISL(2).EQ.1.AND.ISLM.NE.0) THEN 
          ISL(3-ISLM)=0 
          ISLM=3-ISLM 
        ELSEIF(ISL(1).EQ.1.AND.ISL(2).EQ.1) THEN 
          ZDR1=MAX(0.,V(N+1,3)/MAX(1E-6,V(N+1,4))-1.) 
          ZDR2=MAX(0.,V(N+2,3)/MAX(1E-6,V(N+2,4))-1.) 
          IF(ZDR2.GT.RLU(0)*(ZDR1+ZDR2)) ISL(1)=0 
          IF(ISL(1).EQ.1) ISL(2)=0 
          IF(ISL(1).EQ.0) ISLM=1 
          IF(ISL(2).EQ.0) ISLM=2 
        ENDIF 
        IF(ISL(1).EQ.1.OR.ISL(2).EQ.1) GOTO 330 
      ENDIF 
      IFLD1=KFLD(1)
      IF(KFLD(1).GE.6.AND.KFLD(1).LE.8) IFLD1=37+KFLD(1)+
     &ISIGN(2,K(N+1,2)) 
      IFLD2=KFLD(2)
      IF(KFLD(2).GE.6.AND.KFLD(2).LE.8) IFLD2=37+KFLD(2)+
     &ISIGN(2,K(N+2,2)) 
      IF(IGM.GT.0.AND.MOD(MSTJ(43),2).EQ.1.AND.(P(N+1,5).GE. 
     &PMTH(2,IFLD1).OR.P(N+2,5).GE.PMTH(2,IFLD2))) THEN 
        PMQ1=V(N+1,5)/V(IM,5) 
        PMQ2=V(N+2,5)/V(IM,5) 
        ZD=SQRT(MAX(0.,(1.-V(IM,5)/PEM**2)*((1.-PMQ1-PMQ2)**2- 
     &  4.*PMQ1*PMQ2))) 
        ZH=1.+PMQ1-PMQ2 
        ZL=0.5*(ZH-ZD) 
        ZU=0.5*(ZH+ZD) 
        IF(V(IM,1).LT.ZL.OR.V(IM,1).GT.ZU) GOTO 330 
      ENDIF 
 
C...Accepted branch. Construct four-momentum for initial partons. 
  460 MAZIP=0 
      MAZIC=0 
      IF(NEP.EQ.1) THEN 
        P(N+1,1)=0. 
        P(N+1,2)=0. 
        P(N+1,3)=SQRT(MAX(0.,(P(IPA(1),4)+P(N+1,5))*(P(IPA(1),4)- 
     &  P(N+1,5)))) 
        P(N+1,4)=P(IPA(1),4) 
        V(N+1,2)=P(N+1,4) 
      ELSEIF(IGM.EQ.0.AND.NEP.EQ.2) THEN 
        PED1=0.5*(V(IM,5)+V(N+1,5)-V(N+2,5))/P(IM,5) 
        P(N+1,1)=0. 
        P(N+1,2)=0. 
        P(N+1,3)=SQRT(MAX(0.,(PED1+P(N+1,5))*(PED1-P(N+1,5)))) 
        P(N+1,4)=PED1 
        P(N+2,1)=0. 
        P(N+2,2)=0. 
        P(N+2,3)=-P(N+1,3) 
        P(N+2,4)=P(IM,5)-PED1 
        V(N+1,2)=P(N+1,4) 
        V(N+2,2)=P(N+2,4) 
      ELSEIF(NEP.EQ.3) THEN 
        P(N+1,1)=0. 
        P(N+1,2)=0. 
        P(N+1,3)=SQRT(MAX(0.,PA1S)) 
        P(N+2,1)=SQRT(PTS) 
        P(N+2,2)=0. 
        P(N+2,3)=0.5*(PA3S-PA2S-PA1S)/P(N+1,3) 
        P(N+3,1)=-P(N+2,1) 
        P(N+3,2)=0. 
        P(N+3,3)=-(P(N+1,3)+P(N+2,3)) 
        V(N+1,2)=P(N+1,4) 
        V(N+2,2)=P(N+2,4) 
        V(N+3,2)=P(N+3,4) 
 
C...Construct transverse momentum for ordinary branching in shower. 
      ELSE 
        ZM=V(IM,1) 
        PZM=SQRT(MAX(0.,(PEM+P(IM,5))*(PEM-P(IM,5)))) 
        PMLS=(V(IM,5)-V(N+1,5)-V(N+2,5))**2-4.*V(N+1,5)*V(N+2,5) 
        IF(PZM.LE.0.) THEN 
          PTS=0. 
        ELSEIF(MOD(MSTJ(43),2).EQ.1) THEN 
          PTS=(PEM**2*(ZM*(1.-ZM)*V(IM,5)-(1.-ZM)*V(N+1,5)- 
     &    ZM*V(N+2,5))-0.25*PMLS)/PZM**2 
        ELSE 
          PTS=PMLS*(ZM*(1.-ZM)*PEM**2/V(IM,5)-0.25)/PZM**2 
        ENDIF 
        PT=SQRT(MAX(0.,PTS)) 
 
C...Find coefficient of azimuthal asymmetry due to gluon polarization. 
        HAZIP=0. 
        IF(MSTJ(49).NE.1.AND.MOD(MSTJ(46),2).EQ.1.AND.K(IM,2).EQ.21. 
     &  AND.IAU.NE.0) THEN 
          IF(K(IGM,3).NE.0) MAZIP=1 
          ZAU=V(IGM,1) 
          IF(IAU.EQ.IM+1) ZAU=1.-V(IGM,1) 
          IF(MAZIP.EQ.0) ZAU=0. 
          IF(K(IGM,2).NE.21) THEN 
            HAZIP=2.*ZAU/(1.+ZAU**2) 
          ELSE 
            HAZIP=(ZAU/(1.-ZAU*(1.-ZAU)))**2 
          ENDIF 
          IF(K(N+1,2).NE.21) THEN 
            HAZIP=HAZIP*(-2.*ZM*(1.-ZM))/(1.-2.*ZM*(1.-ZM)) 
          ELSE 
            HAZIP=HAZIP*(ZM*(1.-ZM)/(1.-ZM*(1.-ZM)))**2 
          ENDIF 
        ENDIF 
 
C...Find coefficient of azimuthal asymmetry due to soft gluon 
C...interference. 
        HAZIC=0. 
        IF(MSTJ(49).NE.2.AND.MSTJ(46).GE.2.AND.(K(N+1,2).EQ.21.OR. 
     &  K(N+2,2).EQ.21).AND.IAU.NE.0) THEN 
          IF(K(IGM,3).NE.0) MAZIC=N+1 
          IF(K(IGM,3).NE.0.AND.K(N+1,2).NE.21) MAZIC=N+2 
          IF(K(IGM,3).NE.0.AND.K(N+1,2).EQ.21.AND.K(N+2,2).EQ.21.AND. 
     &    ZM.GT.0.5) MAZIC=N+2 
          IF(K(IAU,2).EQ.22) MAZIC=0 
          ZS=ZM 
          IF(MAZIC.EQ.N+2) ZS=1.-ZM 
          ZGM=V(IGM,1) 
          IF(IAU.EQ.IM-1) ZGM=1.-V(IGM,1) 
          IF(MAZIC.EQ.0) ZGM=1. 
          IF(MAZIC.NE.0) HAZIC=(P(IM,5)/P(IGM,5))*
     &    SQRT((1.-ZS)*(1.-ZGM)/(ZS*ZGM)) 
          HAZIC=MIN(0.95,HAZIC) 
        ENDIF 
      ENDIF 
 
C...Construct kinematics for ordinary branching in shower. 
  470 IF(NEP.EQ.2.AND.IGM.GT.0) THEN 
        IF(MOD(MSTJ(43),2).EQ.1) THEN 
          P(N+1,4)=PEM*V(IM,1) 
        ELSE 
          P(N+1,4)=PEM*(0.5*(V(IM,5)-SQRT(PMLS)+V(N+1,5)-V(N+2,5))+ 
     &    SQRT(PMLS)*ZM)/V(IM,5) 
        ENDIF 
        PHI=PARU(2)*RLU(0) 
        P(N+1,1)=PT*COS(PHI) 
        P(N+1,2)=PT*SIN(PHI) 
        IF(PZM.GT.0.) THEN 
          P(N+1,3)=0.5*(V(N+2,5)-V(N+1,5)-V(IM,5)+2.*PEM*P(N+1,4))/PZM 
        ELSE 
          P(N+1,3)=0. 
        ENDIF 
        P(N+2,1)=-P(N+1,1) 
        P(N+2,2)=-P(N+1,2) 
        P(N+2,3)=PZM-P(N+1,3) 
        P(N+2,4)=PEM-P(N+1,4) 
        IF(MSTJ(43).LE.2) THEN 
          V(N+1,2)=(PEM*P(N+1,4)-PZM*P(N+1,3))/P(IM,5) 
          V(N+2,2)=(PEM*P(N+2,4)-PZM*P(N+2,3))/P(IM,5) 
        ENDIF 
      ENDIF 
 
C...Rotate and boost daughters. 
      IF(IGM.GT.0) THEN 
        IF(MSTJ(43).LE.2) THEN 
          BEX=P(IGM,1)/P(IGM,4) 
          BEY=P(IGM,2)/P(IGM,4) 
          BEZ=P(IGM,3)/P(IGM,4) 
          GA=P(IGM,4)/P(IGM,5) 
          GABEP=GA*(GA*(BEX*P(IM,1)+BEY*P(IM,2)+BEZ*P(IM,3))/(1.+GA)- 
     &    P(IM,4)) 
        ELSE 
          BEX=0. 
          BEY=0. 
          BEZ=0. 
          GA=1. 
          GABEP=0. 
        ENDIF 
        THE=ULANGL(P(IM,3)+GABEP*BEZ,SQRT((P(IM,1)+GABEP*BEX)**2+ 
     &  (P(IM,2)+GABEP*BEY)**2)) 
        PHI=ULANGL(P(IM,1)+GABEP*BEX,P(IM,2)+GABEP*BEY) 
        DO 480 I=N+1,N+2 
        DP(1)=COS(THE)*COS(PHI)*P(I,1)-SIN(PHI)*P(I,2)+ 
     &  SIN(THE)*COS(PHI)*P(I,3) 
        DP(2)=COS(THE)*SIN(PHI)*P(I,1)+COS(PHI)*P(I,2)+ 
     &  SIN(THE)*SIN(PHI)*P(I,3) 
        DP(3)=-SIN(THE)*P(I,1)+COS(THE)*P(I,3) 
        DP(4)=P(I,4) 
        DBP=BEX*DP(1)+BEY*DP(2)+BEZ*DP(3) 
        DGABP=GA*(GA*DBP/(1D0+GA)+DP(4)) 
        P(I,1)=DP(1)+DGABP*BEX 
        P(I,2)=DP(2)+DGABP*BEY 
        P(I,3)=DP(3)+DGABP*BEZ 
        P(I,4)=GA*(DP(4)+DBP) 
  480   CONTINUE 
      ENDIF 
 
C...Weight with azimuthal distribution, if required. 
      IF(MAZIP.NE.0.OR.MAZIC.NE.0) THEN 
        DO 490 J=1,3 
        DPT(1,J)=P(IM,J) 
        DPT(2,J)=P(IAU,J) 
        DPT(3,J)=P(N+1,J) 
  490   CONTINUE 
        DPMA=DPT(1,1)*DPT(2,1)+DPT(1,2)*DPT(2,2)+DPT(1,3)*DPT(2,3) 
        DPMD=DPT(1,1)*DPT(3,1)+DPT(1,2)*DPT(3,2)+DPT(1,3)*DPT(3,3) 
        DPMM=DPT(1,1)**2+DPT(1,2)**2+DPT(1,3)**2 
        DO 500 J=1,3 
        DPT(4,J)=DPT(2,J)-DPMA*DPT(1,J)/DPMM 
        DPT(5,J)=DPT(3,J)-DPMD*DPT(1,J)/DPMM 
  500   CONTINUE 
        DPT(4,4)=SQRT(DPT(4,1)**2+DPT(4,2)**2+DPT(4,3)**2) 
        DPT(5,4)=SQRT(DPT(5,1)**2+DPT(5,2)**2+DPT(5,3)**2) 
        IF(MIN(DPT(4,4),DPT(5,4)).GT.0.1*PARJ(82)) THEN 
          CAD=(DPT(4,1)*DPT(5,1)+DPT(4,2)*DPT(5,2)+ 
     &    DPT(4,3)*DPT(5,3))/(DPT(4,4)*DPT(5,4)) 
          IF(MAZIP.NE.0) THEN 
            IF(1.+HAZIP*(2.*CAD**2-1.).LT.RLU(0)*(1.+ABS(HAZIP))) 
     &      GOTO 470 
          ENDIF 
          IF(MAZIC.NE.0) THEN 
            IF(MAZIC.EQ.N+2) CAD=-CAD 
            IF((1.-HAZIC)*(1.-HAZIC*CAD)/(1.+HAZIC**2-2.*HAZIC*CAD) 
     &      .LT.RLU(0)) GOTO 470 
          ENDIF 
        ENDIF 
      ENDIF 
 
C...Azimuthal anisotropy due to interference with initial state partons. 
      IF(MOD(MIIS,2).EQ.1.AND.IGM.EQ.NS+1.AND.(K(N+1,2).EQ.21.OR. 
     &K(N+2,2).EQ.21)) THEN 
        III=IM-NS-1 
        IF(ISII(III).GE.1) THEN 
          IAZIID=N+1 
          IF(K(N+1,2).NE.21) IAZIID=N+2 
          IF(K(N+1,2).EQ.21.AND.K(N+2,2).EQ.21.AND. 
     &    P(N+1,4).GT.P(N+2,4)) IAZIID=N+2 
          THEIID=ULANGL(P(IAZIID,3),SQRT(P(IAZIID,1)**2+P(IAZIID,2)**2)) 
          IF(III.EQ.2) THEIID=PARU(1)-THEIID 
          PHIIID=ULANGL(P(IAZIID,1),P(IAZIID,2)) 
          HAZII=MIN(0.95,THEIID/THEIIS(III,ISII(III))) 
          CAD=COS(PHIIID-PHIIIS(III,ISII(III))) 
          PHIREL=ABS(PHIIID-PHIIIS(III,ISII(III))) 
          IF(PHIREL.GT.PARU(1)) PHIREL=PARU(2)-PHIREL 
          IF((1.-HAZII)*(1.-HAZII*CAD)/(1.+HAZII**2-2.*HAZII*CAD) 
     &    .LT.RLU(0)) GOTO 470 
        ENDIF 
      ENDIF 
 
C...Continue loop over partons that may branch, until none left. 
      IF(IGM.GE.0) K(IM,1)=14 
      N=N+NEP 
      NEP=2 
      IF(N.GT.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUSHOW:) no more memory left in LUJETS') 
        IF(MSTU(21).GE.1) N=NS 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      GOTO 270 
 
C...Set information on imagined shower initiator. 
  510 IF(NPA.GE.2) THEN 
        K(NS+1,1)=11 
        K(NS+1,2)=94 
        K(NS+1,3)=IP1 
        IF(IP2.GT.0.AND.IP2.LT.IP1) K(NS+1,3)=IP2 
        K(NS+1,4)=NS+2 
        K(NS+1,5)=NS+1+NPA 
        IIM=1 
      ELSE 
        IIM=0 
      ENDIF 
 
C...Reconstruct string drawing information. 
      DO 520 I=NS+1+IIM,N 
      IF(K(I,1).LE.10.AND.K(I,2).EQ.22) THEN 
        K(I,1)=1 
      ELSEIF(K(I,1).LE.10.AND.IABS(K(I,2)).GE.11.AND. 
     &IABS(K(I,2)).LE.18) THEN 
        K(I,1)=1 
      ELSEIF(K(I,1).LE.10) THEN 
        K(I,4)=MSTU(5)*(K(I,4)/MSTU(5)) 
        K(I,5)=MSTU(5)*(K(I,5)/MSTU(5)) 
      ELSEIF(K(MOD(K(I,4),MSTU(5))+1,2).NE.22) THEN 
        ID1=MOD(K(I,4),MSTU(5)) 
        IF(K(I,2).GE.1.AND.K(I,2).LE.8) ID1=MOD(K(I,4),MSTU(5))+1 
        ID2=2*MOD(K(I,4),MSTU(5))+1-ID1 
        K(I,4)=MSTU(5)*(K(I,4)/MSTU(5))+ID1 
        K(I,5)=MSTU(5)*(K(I,5)/MSTU(5))+ID2 
        K(ID1,4)=K(ID1,4)+MSTU(5)*I 
        K(ID1,5)=K(ID1,5)+MSTU(5)*ID2 
        K(ID2,4)=K(ID2,4)+MSTU(5)*ID1 
        K(ID2,5)=K(ID2,5)+MSTU(5)*I 
      ELSE 
        ID1=MOD(K(I,4),MSTU(5)) 
        ID2=ID1+1 
        K(I,4)=MSTU(5)*(K(I,4)/MSTU(5))+ID1 
        K(I,5)=MSTU(5)*(K(I,5)/MSTU(5))+ID1 
        IF(IABS(K(I,2)).LE.10.OR.K(ID1,1).GE.11) THEN 
          K(ID1,4)=K(ID1,4)+MSTU(5)*I 
          K(ID1,5)=K(ID1,5)+MSTU(5)*I 
        ELSE 
          K(ID1,4)=0 
          K(ID1,5)=0 
        ENDIF 
        K(ID2,4)=0 
        K(ID2,5)=0 
      ENDIF 
  520 CONTINUE 
 
C...Transformation from CM frame. 
      IF(NPA.GE.2) THEN 
        BEX=PS(1)/PS(4) 
        BEY=PS(2)/PS(4) 
        BEZ=PS(3)/PS(4) 
        GA=PS(4)/PS(5) 
        GABEP=GA*(GA*(BEX*P(IPA(1),1)+BEY*P(IPA(1),2)+BEZ*P(IPA(1),3)) 
     &  /(1.+GA)-P(IPA(1),4)) 
      ELSE 
        BEX=0. 
        BEY=0. 
        BEZ=0. 
        GABEP=0. 
      ENDIF 
      THE=ULANGL(P(IPA(1),3)+GABEP*BEZ,SQRT((P(IPA(1),1) 
     &+GABEP*BEX)**2+(P(IPA(1),2)+GABEP*BEY)**2)) 
      PHI=ULANGL(P(IPA(1),1)+GABEP*BEX,P(IPA(1),2)+GABEP*BEY) 
      IF(NPA.EQ.3) THEN 
        CHI=ULANGL(COS(THE)*COS(PHI)*(P(IPA(2),1)+GABEP*BEX)+COS(THE)* 
     &  SIN(PHI)*(P(IPA(2),2)+GABEP*BEY)-SIN(THE)*(P(IPA(2),3)+GABEP* 
     &  BEZ),-SIN(PHI)*(P(IPA(2),1)+GABEP*BEX)+COS(PHI)*(P(IPA(2),2)+ 
     &  GABEP*BEY)) 
        MSTU(33)=1 
        CALL LUDBRB(NS+1,N,0.,CHI,0D0,0D0,0D0) 
      ENDIF 
      DBEX=DBLE(BEX) 
      DBEY=DBLE(BEY) 
      DBEZ=DBLE(BEZ) 
      MSTU(33)=1 
      CALL LUDBRB(NS+1,N,THE,PHI,DBEX,DBEY,DBEZ) 
 
C...Decay vertex of shower. 
      DO 540 I=NS+1,N 
      DO 530 J=1,5 
      V(I,J)=V(IP1,J) 
  530 CONTINUE 
  540 CONTINUE 
 
C...Delete trivial shower, else connect initiators. 
      IF(N.EQ.NS+NPA+IIM) THEN 
        N=NS 
      ELSE 
        DO 550 IP=1,NPA 
        K(IPA(IP),1)=14 
        K(IPA(IP),4)=K(IPA(IP),4)+NS+IIM+IP 
        K(IPA(IP),5)=K(IPA(IP),5)+NS+IIM+IP 
        K(NS+IIM+IP,3)=IPA(IP) 
        IF(IIM.EQ.1.AND.MSTU(16).NE.2) K(NS+IIM+IP,3)=NS+1 
        IF(K(NS+IIM+IP,1).NE.1) THEN 
          K(NS+IIM+IP,4)=MSTU(5)*IPA(IP)+K(NS+IIM+IP,4) 
          K(NS+IIM+IP,5)=MSTU(5)*IPA(IP)+K(NS+IIM+IP,5) 
        ENDIF 
  550   CONTINUE 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUBOEI(NSAV) 
 
C...Purpose: to modify event so as to approximately take into account 
C...Bose-Einstein effects according to a simple phenomenological 
C...parametrization. 
      IMPLICIT DOUBLE PRECISION(D) 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUJETS/,/LUDAT1/ 
      DIMENSION DPS(4),KFBE(9),NBE(0:9),BEI(100) 
      DATA KFBE/211,-211,111,321,-321,130,310,221,331/ 
 
C...Boost event to overall CM frame. Calculate CM energy. 
      IF((MSTJ(51).NE.1.AND.MSTJ(51).NE.2).OR.N-NSAV.LE.1) RETURN 
      DO 100 J=1,4 
      DPS(J)=0. 
  100 CONTINUE 
      DO 120 I=1,N 
      KFA=IABS(K(I,2))
      IF(K(I,1).LE.10.AND.((KFA.GT.10.AND.KFA.LE.20).OR.KFA.EQ.22).AND.
     &K(I,3).GT.0) THEN
        KFMA=IABS(K(K(I,3),2))
        IF(KFMA.GT.10.AND.KFMA.LE.80) K(I,1)=-K(I,1)
      ENDIF
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 120 
      DO 110 J=1,4 
      DPS(J)=DPS(J)+P(I,J) 
  110 CONTINUE 
  120 CONTINUE 
      CALL LUDBRB(0,0,0.,0.,-DPS(1)/DPS(4),-DPS(2)/DPS(4), 
     &-DPS(3)/DPS(4)) 
      PECM=0. 
      DO 130 I=1,N 
      IF(K(I,1).GE.1.AND.K(I,1).LE.10) PECM=PECM+P(I,4) 
  130 CONTINUE 
 
C...Reserve copy of particles by species at end of record. 
      NBE(0)=N+MSTU(3) 
      DO 160 IBE=1,MIN(9,MSTJ(52)) 
      NBE(IBE)=NBE(IBE-1) 
      DO 150 I=NSAV+1,N 
      IF(K(I,2).NE.KFBE(IBE)) GOTO 150 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 150 
      IF(NBE(IBE).GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUBOEI:) no more memory left in LUJETS') 
        RETURN 
      ENDIF 
      NBE(IBE)=NBE(IBE)+1 
      K(NBE(IBE),1)=I 
      DO 140 J=1,3 
      P(NBE(IBE),J)=0. 
  140 CONTINUE 
  150 CONTINUE 
  160 CONTINUE 
      IF(NBE(MIN(9,MSTJ(52)))-NBE(0).LE.1) GOTO 280
 
C...Tabulate integral for subsequent momentum shift. 
      DO 220 IBE=1,MIN(9,MSTJ(52)) 
      IF(IBE.NE.1.AND.IBE.NE.4.AND.IBE.LE.7) GOTO 180 
      IF(IBE.EQ.1.AND.MAX(NBE(1)-NBE(0),NBE(2)-NBE(1),NBE(3)-NBE(2)) 
     &.LE.1) GOTO 180 
      IF(IBE.EQ.4.AND.MAX(NBE(4)-NBE(3),NBE(5)-NBE(4),NBE(6)-NBE(5), 
     &NBE(7)-NBE(6)).LE.1) GOTO 180 
      IF(IBE.GE.8.AND.NBE(IBE)-NBE(IBE-1).LE.1) GOTO 180 
      IF(IBE.EQ.1) PMHQ=2.*ULMASS(211) 
      IF(IBE.EQ.4) PMHQ=2.*ULMASS(321) 
      IF(IBE.EQ.8) PMHQ=2.*ULMASS(221) 
      IF(IBE.EQ.9) PMHQ=2.*ULMASS(331) 
      QDEL=0.1*MIN(PMHQ,PARJ(93)) 
      IF(MSTJ(51).EQ.1) THEN 
        NBIN=MIN(100,NINT(9.*PARJ(93)/QDEL)) 
        BEEX=EXP(0.5*QDEL/PARJ(93)) 
        BERT=EXP(-QDEL/PARJ(93)) 
      ELSE 
        NBIN=MIN(100,NINT(3.*PARJ(93)/QDEL)) 
      ENDIF 
      DO 170 IBIN=1,NBIN 
      QBIN=QDEL*(IBIN-0.5) 
      BEI(IBIN)=QDEL*(QBIN**2+QDEL**2/12.)/SQRT(QBIN**2+PMHQ**2) 
      IF(MSTJ(51).EQ.1) THEN 
        BEEX=BEEX*BERT 
        BEI(IBIN)=BEI(IBIN)*BEEX 
      ELSE 
        BEI(IBIN)=BEI(IBIN)*EXP(-(QBIN/PARJ(93))**2) 
      ENDIF 
      IF(IBIN.GE.2) BEI(IBIN)=BEI(IBIN)+BEI(IBIN-1) 
  170 CONTINUE 
 
C...Loop through particle pairs and find old relative momentum. 
  180 DO 210 I1M=NBE(IBE-1)+1,NBE(IBE)-1 
      I1=K(I1M,1) 
      DO 200 I2M=I1M+1,NBE(IBE) 
      I2=K(I2M,1) 
      Q2OLD=MAX(0.,(P(I1,4)+P(I2,4))**2-(P(I1,1)+P(I2,1))**2-(P(I1,2)+ 
     &P(I2,2))**2-(P(I1,3)+P(I2,3))**2-(P(I1,5)+P(I2,5))**2) 
      QOLD=SQRT(Q2OLD) 
 
C...Calculate new relative momentum. 
      IF(QOLD.LT.1E-3*QDEL) THEN 
        GOTO 200 
      ELSEIF(QOLD.LE.QDEL) THEN 
        QMOV=QOLD/3. 
      ELSEIF(QOLD.LT.(NBIN-0.1)*QDEL) THEN 
        RBIN=QOLD/QDEL 
        IBIN=RBIN 
        RINP=(RBIN**3-IBIN**3)/(3*IBIN*(IBIN+1)+1) 
        QMOV=(BEI(IBIN)+RINP*(BEI(IBIN+1)-BEI(IBIN)))* 
     &  SQRT(Q2OLD+PMHQ**2)/Q2OLD 
      ELSE 
        QMOV=BEI(NBIN)*SQRT(Q2OLD+PMHQ**2)/Q2OLD 
      ENDIF 
      Q2NEW=Q2OLD*(QOLD/(QOLD+3.*PARJ(92)*QMOV))**(2./3.) 
 
C...Calculate and save shift to be performed on three-momenta. 
      HC1=(P(I1,4)+P(I2,4))**2-(Q2OLD-Q2NEW) 
      HC2=(Q2OLD-Q2NEW)*(P(I1,4)-P(I2,4))**2 
      HA=0.5*(1.-SQRT(HC1*Q2NEW/(HC1*Q2OLD-HC2))) 
      DO 190 J=1,3 
      PD=HA*(P(I2,J)-P(I1,J)) 
      P(I1M,J)=P(I1M,J)+PD 
      P(I2M,J)=P(I2M,J)-PD 
  190 CONTINUE 
  200 CONTINUE 
  210 CONTINUE 
  220 CONTINUE 
 
C...Shift momenta and recalculate energies. 
      DO 240 IM=NBE(0)+1,NBE(MIN(9,MSTJ(52))) 
      I=K(IM,1) 
      DO 230 J=1,3 
      P(I,J)=P(I,J)+P(IM,J) 
  230 CONTINUE 
      P(I,4)=SQRT(P(I,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2) 
  240 CONTINUE 
 
C...Rescale all momenta for energy conservation. 
      PES=0. 
      PQS=0. 
      DO 250 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 250 
      PES=PES+P(I,4) 
      PQS=PQS+P(I,5)**2/P(I,4) 
  250 CONTINUE 
      FAC=(PECM-PQS)/(PES-PQS) 
      DO 270 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 270 
      DO 260 J=1,3 
      P(I,J)=FAC*P(I,J) 
  260 CONTINUE 
      P(I,4)=SQRT(P(I,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2) 
  270 CONTINUE 
 
C...Boost back to correct reference frame. 
  280 CALL LUDBRB(0,0,0.,0.,DPS(1)/DPS(4),DPS(2)/DPS(4),DPS(3)/DPS(4)) 
      DO 290 I=1,N
      IF(K(I,1).LT.0) K(I,1)=-K(I,1)
  290 CONTINUE
 
      RETURN 
      END 
 
C********************************************************************* 
 
      FUNCTION ULMASS(KF) 
 
C...Purpose: to give the mass of a particle/parton. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT1/,/LUDAT2/ 
 
C...Reset variables. Compressed code. 
      ULMASS=0. 
      KFA=IABS(KF) 
      KC=LUCOMP(KF) 
      IF(KC.EQ.0) RETURN 
      PARF(106)=PMAS(6,1) 
      PARF(107)=PMAS(7,1) 
      PARF(108)=PMAS(8,1) 
 
C...Guarantee use of constituent masses for internal checks. 
      IF((MSTJ(93).EQ.1.OR.MSTJ(93).EQ.2).AND.KFA.LE.10) THEN 
        ULMASS=PARF(100+KFA) 
        IF(MSTJ(93).EQ.2) ULMASS=MAX(0.,ULMASS-PARF(121)) 
 
C...Masses that can be read directly off table. 
      ELSEIF(KFA.LE.100.OR.KC.LE.80.OR.KC.GT.100) THEN 
        ULMASS=PMAS(KC,1) 
 
C...Find constituent partons and their masses. 
      ELSE 
        KFLA=MOD(KFA/1000,10) 
        KFLB=MOD(KFA/100,10) 
        KFLC=MOD(KFA/10,10) 
        KFLS=MOD(KFA,10) 
        KFLR=MOD(KFA/10000,10) 
        PMA=PARF(100+KFLA) 
        PMB=PARF(100+KFLB) 
        PMC=PARF(100+KFLC) 
 
C...Construct masses for various meson, diquark and baryon cases. 
        IF(KFLA.EQ.0.AND.KFLR.EQ.0.AND.KFLS.LE.3) THEN 
          IF(KFLS.EQ.1) PMSPL=-3./(PMB*PMC) 
          IF(KFLS.GE.3) PMSPL=1./(PMB*PMC) 
          ULMASS=PARF(111)+PMB+PMC+PARF(113)*PARF(101)**2*PMSPL 
        ELSEIF(KFLA.EQ.0) THEN 
          KMUL=2 
          IF(KFLS.EQ.1) KMUL=3 
          IF(KFLR.EQ.2) KMUL=4 
          IF(KFLS.EQ.5) KMUL=5 
          ULMASS=PARF(113+KMUL)+PMB+PMC 
        ELSEIF(KFLC.EQ.0) THEN 
          IF(KFLS.EQ.1) PMSPL=-3./(PMA*PMB) 
          IF(KFLS.EQ.3) PMSPL=1./(PMA*PMB) 
          ULMASS=2.*PARF(112)/3.+PMA+PMB+PARF(114)*PARF(101)**2*PMSPL 
          IF(MSTJ(93).EQ.1) ULMASS=PMA+PMB 
          IF(MSTJ(93).EQ.2) ULMASS=MAX(0.,ULMASS-PARF(122)- 
     &    2.*PARF(112)/3.) 
        ELSE 
          IF(KFLS.EQ.2.AND.KFLA.EQ.KFLB) THEN 
            PMSPL=1./(PMA*PMB)-2./(PMA*PMC)-2./(PMB*PMC) 
          ELSEIF(KFLS.EQ.2.AND.KFLB.GE.KFLC) THEN 
            PMSPL=-2./(PMA*PMB)-2./(PMA*PMC)+1./(PMB*PMC) 
          ELSEIF(KFLS.EQ.2) THEN 
            PMSPL=-3./(PMB*PMC) 
          ELSE 
            PMSPL=1./(PMA*PMB)+1./(PMA*PMC)+1./(PMB*PMC) 
          ENDIF 
          ULMASS=PARF(112)+PMA+PMB+PMC+PARF(114)*PARF(101)**2*PMSPL 
        ENDIF 
      ENDIF 
 
C...Optional mass broadening according to truncated Breit-Wigner 
C...(either in m or in m^2). 
      IF(MSTJ(24).GE.1.AND.PMAS(KC,2).GT.1E-4) THEN 
        IF(MSTJ(24).EQ.1.OR.(MSTJ(24).EQ.2.AND.KFA.GT.100)) THEN 
          ULMASS=ULMASS+0.5*PMAS(KC,2)*TAN((2.*RLU(0)-1.)* 
     &    ATAN(2.*PMAS(KC,3)/PMAS(KC,2))) 
        ELSE 
          PM0=ULMASS 
          PMLOW=ATAN((MAX(0.,PM0-PMAS(KC,3))**2-PM0**2)/ 
     &    (PM0*PMAS(KC,2))) 
          PMUPP=ATAN(((PM0+PMAS(KC,3))**2-PM0**2)/(PM0*PMAS(KC,2))) 
          ULMASS=SQRT(MAX(0.,PM0**2+PM0*PMAS(KC,2)*TAN(PMLOW+ 
     &    (PMUPP-PMLOW)*RLU(0)))) 
        ENDIF 
      ENDIF 
      MSTJ(93)=0 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUNAME(KF,CHAU) 
 
C...Purpose: to give the particle/parton name as a character string. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT4/CHAF(500) 
      CHARACTER CHAF*8 
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT4/ 
      CHARACTER CHAU*16 
 
C...Initial values. Charge. Subdivide code. 
      CHAU=' ' 
      KFA=IABS(KF) 
      KC=LUCOMP(KF) 
      IF(KC.EQ.0) RETURN 
      KQ=LUCHGE(KF) 
      KFLA=MOD(KFA/1000,10) 
      KFLB=MOD(KFA/100,10) 
      KFLC=MOD(KFA/10,10) 
      KFLS=MOD(KFA,10) 
      KFLR=MOD(KFA/10000,10) 
 
C...Read out root name and spin for simple particle. 
      IF(KFA.LE.100.OR.(KFA.GT.100.AND.KC.GT.100)) THEN 
        CHAU=CHAF(KC) 
        LEN=0 
        DO 100 LEM=1,8 
        IF(CHAU(LEM:LEM).NE.' ') LEN=LEM 
  100   CONTINUE 
 
C...Construct root name for diquark. Add on spin. 
      ELSEIF(KFLC.EQ.0) THEN 
        CHAU(1:2)=CHAF(KFLA)(1:1)//CHAF(KFLB)(1:1) 
        IF(KFLS.EQ.1) CHAU(3:4)='_0' 
        IF(KFLS.EQ.3) CHAU(3:4)='_1' 
        LEN=4 
 
C...Construct root name for heavy meson. Add on spin and heavy flavour. 
      ELSEIF(KFLA.EQ.0) THEN 
        IF(KFLB.EQ.5) CHAU(1:1)='B' 
        IF(KFLB.EQ.6) CHAU(1:1)='T' 
        IF(KFLB.EQ.7) CHAU(1:1)='L' 
        IF(KFLB.EQ.8) CHAU(1:1)='H' 
        LEN=1 
        IF(KFLR.EQ.0.AND.KFLS.EQ.1) THEN 
        ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.3) THEN 
          CHAU(2:2)='*' 
          LEN=2 
        ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.3) THEN 
          CHAU(2:3)='_1' 
          LEN=3 
        ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.1) THEN 
          CHAU(2:4)='*_0' 
          LEN=4 
        ELSEIF(KFLR.EQ.2) THEN 
          CHAU(2:4)='*_1' 
          LEN=4 
        ELSEIF(KFLS.EQ.5) THEN 
          CHAU(2:4)='*_2' 
          LEN=4 
        ENDIF 
        IF(KFLC.GE.3.AND.KFLR.EQ.0.AND.KFLS.LE.3) THEN 
          CHAU(LEN+1:LEN+2)='_'//CHAF(KFLC)(1:1) 
          LEN=LEN+2 
        ELSEIF(KFLC.GE.3) THEN 
          CHAU(LEN+1:LEN+1)=CHAF(KFLC)(1:1) 
          LEN=LEN+1 
        ENDIF 
 
C...Construct root name and spin for heavy baryon. 
      ELSE 
        IF(KFLB.LE.2.AND.KFLC.LE.2) THEN 
          CHAU='Sigma ' 
          IF(KFLC.GT.KFLB) CHAU='Lambda' 
          IF(KFLS.EQ.4) CHAU='Sigma*' 
          LEN=5 
          IF(CHAU(6:6).NE.' ') LEN=6 
        ELSEIF(KFLB.LE.2.OR.KFLC.LE.2) THEN 
          CHAU='Xi ' 
          IF(KFLA.GT.KFLB.AND.KFLB.GT.KFLC) CHAU='Xi''' 
          IF(KFLS.EQ.4) CHAU='Xi*' 
          LEN=2 
          IF(CHAU(3:3).NE.' ') LEN=3 
        ELSE 
          CHAU='Omega ' 
          IF(KFLA.GT.KFLB.AND.KFLB.GT.KFLC) CHAU='Omega''' 
          IF(KFLS.EQ.4) CHAU='Omega*' 
          LEN=5 
          IF(CHAU(6:6).NE.' ') LEN=6 
        ENDIF 
 
C...Add on heavy flavour content for heavy baryon. 
        CHAU(LEN+1:LEN+2)='_'//CHAF(KFLA)(1:1) 
        LEN=LEN+2 
        IF(KFLB.GE.KFLC.AND.KFLC.GE.4) THEN 
          CHAU(LEN+1:LEN+2)=CHAF(KFLB)(1:1)//CHAF(KFLC)(1:1) 
          LEN=LEN+2 
        ELSEIF(KFLB.GE.KFLC.AND.KFLB.GE.4) THEN 
          CHAU(LEN+1:LEN+1)=CHAF(KFLB)(1:1) 
          LEN=LEN+1 
        ELSEIF(KFLC.GT.KFLB.AND.KFLB.GE.4) THEN 
          CHAU(LEN+1:LEN+2)=CHAF(KFLC)(1:1)//CHAF(KFLB)(1:1) 
          LEN=LEN+2 
        ELSEIF(KFLC.GT.KFLB.AND.KFLC.GE.4) THEN 
          CHAU(LEN+1:LEN+1)=CHAF(KFLC)(1:1) 
          LEN=LEN+1 
        ENDIF 
      ENDIF 
 
C...Add on bar sign for antiparticle (where necessary). 
      IF(KF.GT.0.OR.LEN.EQ.0) THEN 
      ELSEIF(KFA.GT.10.AND.KFA.LE.40.AND.KQ.NE.0.AND.MOD(KQ,3).EQ.0) 
     &THEN 
      ELSEIF(KFA.EQ.89.OR.(KFA.GE.91.AND.KFA.LE.99)) THEN 
      ELSEIF(KFA.GT.100.AND.KFLA.EQ.0.AND.KQ.NE.0) THEN 
      ELSEIF(MSTU(15).LE.1) THEN 
        CHAU(LEN+1:LEN+1)='~' 
        LEN=LEN+1 
      ELSE 
        CHAU(LEN+1:LEN+3)='bar' 
        LEN=LEN+3 
      ENDIF 
 
C...Add on charge where applicable (conventional cases skipped). 
      IF(KQ.EQ.6) CHAU(LEN+1:LEN+2)='++' 
      IF(KQ.EQ.-6) CHAU(LEN+1:LEN+2)='--' 
      IF(KQ.EQ.3) CHAU(LEN+1:LEN+1)='+' 
      IF(KQ.EQ.-3) CHAU(LEN+1:LEN+1)='-' 
      IF(KQ.EQ.0.AND.(KFA.LE.22.OR.LEN.EQ.0)) THEN 
      ELSEIF(KQ.EQ.0.AND.(KFA.GE.81.AND.KFA.LE.100)) THEN 
      ELSEIF(KFA.EQ.28.OR.KFA.EQ.29) THEN 
      ELSEIF(KFA.GT.100.AND.KFLA.EQ.0.AND.KFLB.EQ.KFLC.AND. 
     &KFLB.NE.1) THEN 
      ELSEIF(KQ.EQ.0) THEN 
        CHAU(LEN+1:LEN+1)='0' 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      FUNCTION LUCHGE(KF) 
 
C...Purpose: to give three times the charge for a particle/parton. 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT2/ 
 
C...Initial values. Simple case of direct readout. 
      LUCHGE=0 
      KFA=IABS(KF) 
      KC=LUCOMP(KFA) 
      IF(KC.EQ.0) THEN 
      ELSEIF(KFA.LE.100.OR.KC.LE.80.OR.KC.GT.100) THEN 
        LUCHGE=KCHG(KC,1) 
 
C...Construction from quark content for heavy meson, diquark, baryon. 
      ELSEIF(MOD(KFA/1000,10).EQ.0) THEN 
        LUCHGE=(KCHG(MOD(KFA/100,10),1)-KCHG(MOD(KFA/10,10),1))* 
     &  (-1)**MOD(KFA/100,10) 
      ELSEIF(MOD(KFA/10,10).EQ.0) THEN 
        LUCHGE=KCHG(MOD(KFA/1000,10),1)+KCHG(MOD(KFA/100,10),1) 
      ELSE 
        LUCHGE=KCHG(MOD(KFA/1000,10),1)+KCHG(MOD(KFA/100,10),1)+ 
     &  KCHG(MOD(KFA/10,10),1) 
      ENDIF 
 
C...Add on correct sign. 
      LUCHGE=LUCHGE*ISIGN(1,KF) 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      FUNCTION LUCOMP(KF) 
 
C...Purpose: to compress the standard KF codes for use in mass and decay 
C...arrays; also to check whether a given code actually is defined. 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT2/ 
      DIMENSION KFTAB(25),KCTAB(25) 
      DATA KFTAB/211,111,221,311,321,130,310,213,113,223, 
     &313,323,2112,2212,210,2110,2210,110,220,330,440,30443,30553,0,0/ 
      DATA KCTAB/101,111,112,102,103,221,222,121,131,132, 
     &122,123,332,333,281,282,283,284,285,286,287,231,235,0,0/ 
 
C...Starting values. 
      LUCOMP=0 
      KFA=IABS(KF) 
 
C...Simple cases: direct translation or table. 
      IF(KFA.EQ.0.OR.KFA.GE.100000) THEN 
        RETURN 
      ELSEIF(KFA.LE.100) THEN 
        LUCOMP=KFA 
        IF(KF.LT.0.AND.KCHG(KFA,3).EQ.0) LUCOMP=0 
        RETURN 
      ELSE 
        DO 100 IKF=1,23 
        IF(KFA.EQ.KFTAB(IKF)) THEN 
          LUCOMP=KCTAB(IKF) 
          IF(KF.LT.0.AND.KCHG(LUCOMP,3).EQ.0) LUCOMP=0 
          RETURN 
        ENDIF 
  100   CONTINUE 
      ENDIF 
 
C...Subdivide KF code into constituent pieces. 
      KFLA=MOD(KFA/1000,10) 
      KFLB=MOD(KFA/100,10) 
      KFLC=MOD(KFA/10,10) 
      KFLS=MOD(KFA,10) 
      KFLR=MOD(KFA/10000,10) 
 
C...Mesons. 
      IF(KFA-10000*KFLR.LT.1000) THEN 
        IF(KFLB.EQ.0.OR.KFLB.EQ.9.OR.KFLC.EQ.0.OR.KFLC.EQ.9) THEN 
        ELSEIF(KFLB.LT.KFLC) THEN 
        ELSEIF(KF.LT.0.AND.KFLB.EQ.KFLC) THEN 
        ELSEIF(KFLB.EQ.KFLC) THEN 
          IF(KFLR.EQ.0.AND.KFLS.EQ.1) THEN 
            LUCOMP=110+KFLB 
          ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.3) THEN 
            LUCOMP=130+KFLB 
          ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.3) THEN 
            LUCOMP=150+KFLB 
          ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.1) THEN 
            LUCOMP=170+KFLB 
          ELSEIF(KFLR.EQ.2.AND.KFLS.EQ.3) THEN 
            LUCOMP=190+KFLB 
          ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.5) THEN 
            LUCOMP=210+KFLB 
          ENDIF 
        ELSEIF(KFLB.LE.5) THEN 
          IF(KFLR.EQ.0.AND.KFLS.EQ.1) THEN 
            LUCOMP=100+((KFLB-1)*(KFLB-2))/2+KFLC 
          ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.3) THEN 
            LUCOMP=120+((KFLB-1)*(KFLB-2))/2+KFLC 
          ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.3) THEN 
            LUCOMP=140+((KFLB-1)*(KFLB-2))/2+KFLC 
          ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.1) THEN 
            LUCOMP=160+((KFLB-1)*(KFLB-2))/2+KFLC 
          ELSEIF(KFLR.EQ.2.AND.KFLS.EQ.3) THEN 
            LUCOMP=180+((KFLB-1)*(KFLB-2))/2+KFLC 
          ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.5) THEN 
            LUCOMP=200+((KFLB-1)*(KFLB-2))/2+KFLC 
          ENDIF 
        ELSEIF((KFLS.EQ.1.AND.KFLR.LE.1).OR.(KFLS.EQ.3.AND.KFLR.LE.2) 
     &  .OR.(KFLS.EQ.5.AND.KFLR.EQ.0)) THEN 
          LUCOMP=80+KFLB 
        ENDIF 
 
C...Diquarks. 
      ELSEIF((KFLR.EQ.0.OR.KFLR.EQ.1).AND.KFLC.EQ.0) THEN 
        IF(KFLS.NE.1.AND.KFLS.NE.3) THEN 
        ELSEIF(KFLA.EQ.9.OR.KFLB.EQ.0.OR.KFLB.EQ.9) THEN 
        ELSEIF(KFLA.LT.KFLB) THEN 
        ELSEIF(KFLS.EQ.1.AND.KFLA.EQ.KFLB) THEN 
        ELSE 
          LUCOMP=90 
        ENDIF 
 
C...Spin 1/2 baryons. 
      ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.2) THEN 
        IF(KFLA.EQ.9.OR.KFLB.EQ.0.OR.KFLB.EQ.9.OR.KFLC.EQ.9) THEN 
        ELSEIF(KFLA.LE.KFLC.OR.KFLA.LT.KFLB) THEN 
        ELSEIF(KFLA.GE.6.OR.KFLB.GE.4.OR.KFLC.GE.4) THEN 
          LUCOMP=80+KFLA 
        ELSEIF(KFLB.LT.KFLC) THEN 
          LUCOMP=300+((KFLA+1)*KFLA*(KFLA-1))/6+(KFLC*(KFLC-1))/2+KFLB 
        ELSE 
          LUCOMP=330+((KFLA+1)*KFLA*(KFLA-1))/6+(KFLB*(KFLB-1))/2+KFLC 
        ENDIF 
 
C...Spin 3/2 baryons. 
      ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.4) THEN 
        IF(KFLA.EQ.9.OR.KFLB.EQ.0.OR.KFLB.EQ.9.OR.KFLC.EQ.9) THEN 
        ELSEIF(KFLA.LT.KFLB.OR.KFLB.LT.KFLC) THEN 
        ELSEIF(KFLA.GE.6.OR.KFLB.GE.4) THEN 
          LUCOMP=80+KFLA 
        ELSE 
          LUCOMP=360+((KFLA+1)*KFLA*(KFLA-1))/6+(KFLB*(KFLB-1))/2+KFLC 
        ENDIF 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUERRM(MERR,CHMESS) 
 
C...Purpose: to inform user of errors in program execution. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUJETS/,/LUDAT1/ 
      CHARACTER CHMESS*(*) 
 
C...Write first few warnings, then be silent. 
      IF(MERR.LE.10) THEN 
        MSTU(27)=MSTU(27)+1 
        MSTU(28)=MERR 
        IF(MSTU(25).EQ.1.AND.MSTU(27).LE.MSTU(26)) WRITE(MSTU(11),5000) 
     &  MERR,MSTU(31),CHMESS 
 
C...Write first few errors, then be silent or stop program. 
      ELSEIF(MERR.LE.20) THEN 
        MSTU(23)=MSTU(23)+1 
        MSTU(24)=MERR-10 
        IF(MSTU(21).GE.1.AND.MSTU(23).LE.MSTU(22)) WRITE(MSTU(11),5100) 
     &  MERR-10,MSTU(31),CHMESS 
        IF(MSTU(21).GE.2.AND.MSTU(23).GT.MSTU(22)) THEN 
          WRITE(MSTU(11),5100) MERR-10,MSTU(31),CHMESS 
          WRITE(MSTU(11),5200) 
          IF(MERR.NE.17) CALL LULIST(2) 
          STOP 
        ENDIF 
 
C...Stop program in case of irreparable error. 
      ELSE 
        WRITE(MSTU(11),5300) MERR-20,MSTU(31),CHMESS 
        STOP 
      ENDIF 
 
C...Formats for output. 
 5000 FORMAT(/5X,'Advisory warning type',I2,' given after',I6, 
     &' LUEXEC calls:'/5X,A) 
 5100 FORMAT(/5X,'Error type',I2,' has occured after',I6, 
     &' LUEXEC calls:'/5X,A) 
 5200 FORMAT(5X,'Execution will be stopped after listing of last ', 
     &'event!') 
 5300 FORMAT(/5X,'Fatal error type',I2,' has occured after',I6, 
     &' LUEXEC calls:'/5X,A/5X,'Execution will now be stopped!') 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      FUNCTION ULALEM(Q2) 
 
C...Purpose: to calculate the running alpha_electromagnetic. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
 
C...Calculate real part of photon vacuum polarization. 
C...For leptons simplify by using asymptotic (Q^2 >> m^2) expressions. 
C...For hadrons use parametrization of H. Burkhardt et al. 
C...See R. Kleiss et al, CERN 89-08, vol. 3, pp. 129-131. 
      AEMPI=PARU(101)/(3.*PARU(1)) 
      IF(MSTU(101).LE.0.OR.Q2.LT.2E-6) THEN 
        RPIGG=0. 
      ELSEIF(MSTU(101).EQ.2.AND.Q2.LT.PARU(104)) THEN
        RPIGG=0.
      ELSEIF(MSTU(101).EQ.2) THEN
        RPIGG=1.-PARU(101)/PARU(103) 
      ELSEIF(Q2.LT.0.09) THEN 
        RPIGG=AEMPI*(13.4916+LOG(Q2))+0.00835*LOG(1.+Q2) 
      ELSEIF(Q2.LT.9.) THEN 
        RPIGG=AEMPI*(16.3200+2.*LOG(Q2))+0.00238*LOG(1.+3.927*Q2) 
      ELSEIF(Q2.LT.1E4) THEN 
        RPIGG=AEMPI*(13.4955+3.*LOG(Q2))+0.00165+0.00299*LOG(1.+Q2) 
      ELSE 
        RPIGG=AEMPI*(13.4955+3.*LOG(Q2))+0.00221+0.00293*LOG(1.+Q2) 
      ENDIF 
 
C...Calculate running alpha_em. 
      ULALEM=PARU(101)/(1.-RPIGG) 
      PARU(108)=ULALEM 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      FUNCTION ULALPS(Q2) 
 
C...Purpose: to give the value of alpha_strong. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT1/,/LUDAT2/ 
 
C...Constant alpha_strong trivial. 
      IF(MSTU(111).LE.0) THEN 
        ULALPS=PARU(111) 
        MSTU(118)=MSTU(112) 
        PARU(117)=0. 
        PARU(118)=PARU(111) 
        RETURN 
      ENDIF 
 
C...Find effective Q2, number of flavours and Lambda. 
      Q2EFF=Q2 
      IF(MSTU(115).GE.2) Q2EFF=MAX(Q2,PARU(114)) 
      NF=MSTU(112) 
      ALAM2=PARU(112)**2 
  100 IF(NF.GT.MAX(2,MSTU(113))) THEN 
        Q2THR=PARU(113)*PMAS(NF,1)**2 
        IF(Q2EFF.LT.Q2THR) THEN 
          NF=NF-1 
          ALAM2=ALAM2*(Q2THR/ALAM2)**(2./(33.-2.*NF)) 
          GOTO 100 
        ENDIF 
      ENDIF 
  110 IF(NF.LT.MIN(8,MSTU(114))) THEN 
        Q2THR=PARU(113)*PMAS(NF+1,1)**2 
        IF(Q2EFF.GT.Q2THR) THEN 
          NF=NF+1 
          ALAM2=ALAM2*(ALAM2/Q2THR)**(2./(33.-2.*NF)) 
          GOTO 110 
        ENDIF 
      ENDIF 
      IF(MSTU(115).EQ.1) Q2EFF=Q2EFF+ALAM2 
      PARU(117)=SQRT(ALAM2) 
 
C...Evaluate first or second order alpha_strong. 
      B0=(33.-2.*NF)/6. 
      ALGQ=LOG(MAX(1.0001,Q2EFF/ALAM2)) 
      IF(MSTU(111).EQ.1) THEN 
        ULALPS=MIN(PARU(115),PARU(2)/(B0*ALGQ)) 
      ELSE 
        B1=(153.-19.*NF)/6. 
        ULALPS=MIN(PARU(115),PARU(2)/(B0*ALGQ)*(1.-B1*LOG(ALGQ)/ 
     &  (B0**2*ALGQ))) 
      ENDIF 
      MSTU(118)=NF 
      PARU(118)=ULALPS 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      FUNCTION ULANGL(X,Y) 
 
C...Purpose: to reconstruct an angle from given x and y coordinates. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
 
      ULANGL=0. 
      R=SQRT(X**2+Y**2) 
      IF(R.LT.1E-20) RETURN 
      IF(ABS(X)/R.LT.0.8) THEN 
        ULANGL=SIGN(ACOS(X/R),Y) 
      ELSE 
        ULANGL=ASIN(Y/R) 
        IF(X.LT.0..AND.ULANGL.GE.0.) THEN 
          ULANGL=PARU(1)-ULANGL 
        ELSEIF(X.LT.0.) THEN 
          ULANGL=-PARU(1)-ULANGL 
        ENDIF 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
C********************************************************************* 
 
      SUBROUTINE RLUGET(LFN,MOVE) 
 
C...Purpose: to dump the state of the random number generator on a file 
C...for subsequent startup from this state onwards. 
      COMMON/LUDATR/MRLU(6),RRLU(100) 
      SAVE /LUDATR/ 
      CHARACTER CHERR*8 
 
C...Backspace required number of records (or as many as there are). 
      IF(MOVE.LT.0) THEN 
        NBCK=MIN(MRLU(6),-MOVE) 
        DO 100 IBCK=1,NBCK 
        BACKSPACE(LFN,ERR=110,IOSTAT=IERR) 
  100   CONTINUE 
        MRLU(6)=MRLU(6)-NBCK 
      ENDIF 
 
C...Unformatted write on unit LFN. 
      WRITE(LFN,ERR=110,IOSTAT=IERR) (MRLU(I1),I1=1,5), 
     &(RRLU(I2),I2=1,100) 
      MRLU(6)=MRLU(6)+1 
      RETURN 
 
C...Write error. 
  110 WRITE(CHERR,'(I8)') IERR 
      CALL LUERRM(18,'(RLUGET:) error when accessing file, IOSTAT ='// 
     &CHERR) 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE RLUSET(LFN,MOVE) 
 
C...Purpose: to read a state of the random number generator from a file 
C...for subsequent generation from this state onwards. 
      COMMON/LUDATR/MRLU(6),RRLU(100) 
      SAVE /LUDATR/ 
      CHARACTER CHERR*8 
 
C...Backspace required number of records (or as many as there are). 
      IF(MOVE.LT.0) THEN 
        NBCK=MIN(MRLU(6),-MOVE) 
        DO 100 IBCK=1,NBCK 
        BACKSPACE(LFN,ERR=120,IOSTAT=IERR) 
  100   CONTINUE 
        MRLU(6)=MRLU(6)-NBCK 
      ENDIF 
 
C...Unformatted read from unit LFN. 
      NFOR=1+MAX(0,MOVE) 
      DO 110 IFOR=1,NFOR 
      READ(LFN,ERR=120,IOSTAT=IERR) (MRLU(I1),I1=1,5), 
     &(RRLU(I2),I2=1,100) 
  110 CONTINUE 
      MRLU(6)=MRLU(6)+NFOR 
      RETURN 
 
C...Write error. 
  120 WRITE(CHERR,'(I8)') IERR 
      CALL LUERRM(18,'(RLUSET:) error when accessing file, IOSTAT ='// 
     &CHERR) 
 
      RETURN 
      END 
 
C********************************************************************* 
 
C********************************************************************* 
 
      SUBROUTINE LUEDIT(MEDIT) 
 
C...Purpose: to perform global manipulations on the event record, 
C...in particular to exclude unstable or undetectable partons/particles. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION NS(2),PTS(2),PLS(2) 
 
C...Remove unwanted partons/particles. 
      IF((MEDIT.GE.0.AND.MEDIT.LE.3).OR.MEDIT.EQ.5) THEN 
        IMAX=N 
        IF(MSTU(2).GT.0) IMAX=MSTU(2) 
        I1=MAX(1,MSTU(1))-1 
        DO 110 I=MAX(1,MSTU(1)),IMAX 
        IF(K(I,1).EQ.0.OR.K(I,1).GT.20) GOTO 110 
        IF(MEDIT.EQ.1) THEN 
          IF(K(I,1).GT.10) GOTO 110 
        ELSEIF(MEDIT.EQ.2) THEN 
          IF(K(I,1).GT.10) GOTO 110 
          KC=LUCOMP(K(I,2)) 
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.KC.EQ.18) 
     &    GOTO 110 
        ELSEIF(MEDIT.EQ.3) THEN 
          IF(K(I,1).GT.10) GOTO 110 
          KC=LUCOMP(K(I,2)) 
          IF(KC.EQ.0) GOTO 110 
          IF(KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) GOTO 110 
        ELSEIF(MEDIT.EQ.5) THEN 
          IF(K(I,1).EQ.13.OR.K(I,1).EQ.14) GOTO 110 
          KC=LUCOMP(K(I,2)) 
          IF(KC.EQ.0) GOTO 110 
          IF(K(I,1).GE.11.AND.KCHG(KC,2).EQ.0) GOTO 110 
        ENDIF 
 
C...Pack remaining partons/particles. Origin no longer known. 
        I1=I1+1 
        DO 100 J=1,5 
        K(I1,J)=K(I,J) 
        P(I1,J)=P(I,J) 
        V(I1,J)=V(I,J) 
  100   CONTINUE 
        K(I1,3)=0 
  110   CONTINUE 
        IF(I1.LT.N) MSTU(3)=0 
        IF(I1.LT.N) MSTU(70)=0 
        N=I1 
 
C...Selective removal of class of entries. New position of retained. 
      ELSEIF(MEDIT.GE.11.AND.MEDIT.LE.15) THEN 
        I1=0 
        DO 120 I=1,N 
        K(I,3)=MOD(K(I,3),MSTU(5)) 
        IF(MEDIT.EQ.11.AND.K(I,1).LT.0) GOTO 120 
        IF(MEDIT.EQ.12.AND.K(I,1).EQ.0) GOTO 120 
        IF(MEDIT.EQ.13.AND.(K(I,1).EQ.11.OR.K(I,1).EQ.12.OR. 
     &  K(I,1).EQ.15).AND.K(I,2).NE.94) GOTO 120 
        IF(MEDIT.EQ.14.AND.(K(I,1).EQ.13.OR.K(I,1).EQ.14.OR. 
     &  K(I,2).EQ.94)) GOTO 120 
        IF(MEDIT.EQ.15.AND.K(I,1).GE.21) GOTO 120 
        I1=I1+1 
        K(I,3)=K(I,3)+MSTU(5)*I1 
  120   CONTINUE 
 
C...Find new event history information and replace old. 
        DO 140 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GT.20.OR.K(I,3)/MSTU(5).EQ.0) GOTO 140 
        ID=I 
  130   IM=MOD(K(ID,3),MSTU(5)) 
        IF(MEDIT.EQ.13.AND.IM.GT.0.AND.IM.LE.N) THEN 
          IF((K(IM,1).EQ.11.OR.K(IM,1).EQ.12.OR.K(IM,1).EQ.15).AND. 
     &    K(IM,2).NE.94) THEN 
            ID=IM 
            GOTO 130 
          ENDIF 
        ELSEIF(MEDIT.EQ.14.AND.IM.GT.0.AND.IM.LE.N) THEN 
          IF(K(IM,1).EQ.13.OR.K(IM,1).EQ.14.OR.K(IM,2).EQ.94) THEN 
            ID=IM 
            GOTO 130 
          ENDIF 
        ENDIF 
        K(I,3)=MSTU(5)*(K(I,3)/MSTU(5)) 
        IF(IM.NE.0) K(I,3)=K(I,3)+K(IM,3)/MSTU(5) 
        IF(K(I,1).NE.3.AND.K(I,1).NE.13.AND.K(I,1).NE.14) THEN 
          IF(K(I,4).GT.0.AND.K(I,4).LE.MSTU(4)) K(I,4)= 
     &    K(K(I,4),3)/MSTU(5) 
          IF(K(I,5).GT.0.AND.K(I,5).LE.MSTU(4)) K(I,5)= 
     &    K(K(I,5),3)/MSTU(5) 
        ELSE 
          KCM=MOD(K(I,4)/MSTU(5),MSTU(5)) 
          IF(KCM.GT.0.AND.KCM.LE.MSTU(4)) KCM=K(KCM,3)/MSTU(5) 
          KCD=MOD(K(I,4),MSTU(5)) 
          IF(KCD.GT.0.AND.KCD.LE.MSTU(4)) KCD=K(KCD,3)/MSTU(5) 
          K(I,4)=MSTU(5)**2*(K(I,4)/MSTU(5)**2)+MSTU(5)*KCM+KCD 
          KCM=MOD(K(I,5)/MSTU(5),MSTU(5)) 
          IF(KCM.GT.0.AND.KCM.LE.MSTU(4)) KCM=K(KCM,3)/MSTU(5) 
          KCD=MOD(K(I,5),MSTU(5)) 
          IF(KCD.GT.0.AND.KCD.LE.MSTU(4)) KCD=K(KCD,3)/MSTU(5) 
          K(I,5)=MSTU(5)**2*(K(I,5)/MSTU(5)**2)+MSTU(5)*KCM+KCD 
        ENDIF 
  140   CONTINUE 
 
C...Pack remaining entries. 
        I1=0 
        MSTU90=MSTU(90) 
        MSTU(90)=0 
        DO 170 I=1,N 
        IF(K(I,3)/MSTU(5).EQ.0) GOTO 170 
        I1=I1+1 
        DO 150 J=1,5 
        K(I1,J)=K(I,J) 
        P(I1,J)=P(I,J) 
        V(I1,J)=V(I,J) 
  150   CONTINUE 
        K(I1,3)=MOD(K(I1,3),MSTU(5)) 
        DO 160 IZ=1,MSTU90 
        IF(I.EQ.MSTU(90+IZ)) THEN 
          MSTU(90)=MSTU(90)+1 
          MSTU(90+MSTU(90))=I1 
          PARU(90+MSTU(90))=PARU(90+IZ) 
        ENDIF 
  160   CONTINUE 
  170   CONTINUE 
        IF(I1.LT.N) MSTU(3)=0 
        IF(I1.LT.N) MSTU(70)=0 
        N=I1 
 
C...Fill in some missing daughter pointers (lost in colour flow). 
      ELSEIF(MEDIT.EQ.16) THEN 
        DO 190 I=1,N 
        IF(K(I,1).LE.10.OR.K(I,1).GT.20) GOTO 190 
        IF(K(I,4).NE.0.OR.K(I,5).NE.0) GOTO 190 
C...Find daughters who point to mother.
        DO 180 I1=I+1,N 
        IF(K(I1,3).NE.I) THEN 
        ELSEIF(K(I,4).EQ.0) THEN 
          K(I,4)=I1 
        ELSE 
          K(I,5)=I1 
        ENDIF 
  180   CONTINUE 
        IF(K(I,5).EQ.0) K(I,5)=K(I,4)
        IF(K(I,4).NE.0) GOTO 190
C...Find daughters who point to documentation version of mother.      
        IM=K(I,3)
        IF(IM.LE.0.OR.IM.GE.I) GOTO 190
        IF(K(IM,1).LE.20.OR.K(IM,1).GT.30) GOTO 190  
        IF(K(IM,2).NE.K(I,2).OR.ABS(P(IM,5)-P(I,5)).GT.1E-2) GOTO 190
        DO 182 I1=I+1,N 
        IF(K(I1,3).NE.IM) THEN 
        ELSEIF(K(I,4).EQ.0) THEN 
          K(I,4)=I1 
        ELSE 
          K(I,5)=I1 
        ENDIF 
  182   CONTINUE 
        IF(K(I,5).EQ.0) K(I,5)=K(I,4)
        IF(K(I,4).NE.0) GOTO 190
C...Find daughters who point to documentation daughters who,
C...in their turn, point to documentation mother.
        ID1=IM
        ID2=IM
        DO 184 I1=IM+1,I-1
        IF(K(I1,3).EQ.IM.AND.K(I1,1).GT.20.AND.K(I1,1).LE.30) THEN
          ID2=I1
          IF(ID1.EQ.IM) ID1=I1
        ENDIF
  184   CONTINUE 
        DO 186 I1=I+1,N 
        IF(K(I1,3).NE.ID1.AND.K(I1,3).NE.ID2) THEN 
        ELSEIF(K(I,4).EQ.0) THEN 
          K(I,4)=I1 
        ELSE 
          K(I,5)=I1 
        ENDIF 
  186   CONTINUE 
        IF(K(I,5).EQ.0) K(I,5)=K(I,4)
  190   CONTINUE 
 
C...Save top entries at bottom of LUJETS commonblock. 
      ELSEIF(MEDIT.EQ.21) THEN 
        IF(2*N.GE.MSTU(4)) THEN 
          CALL LUERRM(11,'(LUEDIT:) no more memory left in LUJETS') 
          RETURN 
        ENDIF 
        DO 210 I=1,N 
        DO 200 J=1,5 
        K(MSTU(4)-I,J)=K(I,J) 
        P(MSTU(4)-I,J)=P(I,J) 
        V(MSTU(4)-I,J)=V(I,J) 
  200   CONTINUE 
  210   CONTINUE 
        MSTU(32)=N 
 
C...Restore bottom entries of commonblock LUJETS to top. 
      ELSEIF(MEDIT.EQ.22) THEN 
        DO 230 I=1,MSTU(32) 
        DO 220 J=1,5 
        K(I,J)=K(MSTU(4)-I,J) 
        P(I,J)=P(MSTU(4)-I,J) 
        V(I,J)=V(MSTU(4)-I,J) 
  220   CONTINUE 
  230   CONTINUE 
        N=MSTU(32) 
 
C...Mark primary entries at top of commonblock LUJETS as untreated. 
      ELSEIF(MEDIT.EQ.23) THEN 
        I1=0 
        DO 240 I=1,N 
        KH=K(I,3) 
        IF(KH.GE.1) THEN 
          IF(K(KH,1).GT.20) KH=0 
        ENDIF 
        IF(KH.NE.0) GOTO 250 
        I1=I1+1 
        IF(K(I,1).GT.10.AND.K(I,1).LE.20) K(I,1)=K(I,1)-10 
  240   CONTINUE 
  250   N=I1 
 
C...Place largest axis along z axis and second largest in xy plane. 
      ELSEIF(MEDIT.EQ.31.OR.MEDIT.EQ.32) THEN 
        CALL LUDBRB(1,N+MSTU(3),0.,-ULANGL(P(MSTU(61),1), 
     &  P(MSTU(61),2)),0D0,0D0,0D0) 
        CALL LUDBRB(1,N+MSTU(3),-ULANGL(P(MSTU(61),3), 
     &  P(MSTU(61),1)),0.,0D0,0D0,0D0) 
        CALL LUDBRB(1,N+MSTU(3),0.,-ULANGL(P(MSTU(61)+1,1), 
     &  P(MSTU(61)+1,2)),0D0,0D0,0D0) 
        IF(MEDIT.EQ.31) RETURN 
 
C...Rotate to put slim jet along +z axis. 
        DO 260 IS=1,2 
        NS(IS)=0 
        PTS(IS)=0. 
        PLS(IS)=0. 
  260   CONTINUE 
        DO 270 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 270 
        IF(MSTU(41).GE.2) THEN 
          KC=LUCOMP(K(I,2)) 
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &    KC.EQ.18) GOTO 270 
          IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &    GOTO 270 
        ENDIF 
        IS=2.-SIGN(0.5,P(I,3)) 
        NS(IS)=NS(IS)+1 
        PTS(IS)=PTS(IS)+SQRT(P(I,1)**2+P(I,2)**2) 
  270   CONTINUE 
        IF(NS(1)*PTS(2)**2.LT.NS(2)*PTS(1)**2) 
     &  CALL LUDBRB(1,N+MSTU(3),PARU(1),0.,0D0,0D0,0D0) 
 
C...Rotate to put second largest jet into -z,+x quadrant. 
        DO 280 I=1,N 
        IF(P(I,3).GE.0.) GOTO 280 
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 280 
        IF(MSTU(41).GE.2) THEN 
          KC=LUCOMP(K(I,2)) 
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &    KC.EQ.18) GOTO 280 
          IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &    GOTO 280 
        ENDIF 
        IS=2.-SIGN(0.5,P(I,1)) 
        PLS(IS)=PLS(IS)-P(I,3) 
  280   CONTINUE 
        IF(PLS(2).GT.PLS(1)) CALL LUDBRB(1,N+MSTU(3),0.,PARU(1), 
     &  0D0,0D0,0D0) 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LULIST(MLIST) 
 
C...Purpose: to give program heading, or list an event, or particle 
C...data, or current parameter values. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/ 
      CHARACTER CHAP*16,CHAC*16,CHAN*16,CHAD(5)*16,CHDL(7)*4 
      DIMENSION PS(6) 
      DATA CHDL/'(())',' ','()','!!','<>','==','(==)'/ 
 
C...Initialization printout: version number and date of last change. 
      IF(MLIST.EQ.0.OR.MSTU(12).EQ.1) THEN 
        CALL LULOGO 
        MSTU(12)=0 
        IF(MLIST.EQ.0) RETURN 
      ENDIF 
 
C...List event data, including additional lines after N. 
      IF(MLIST.GE.1.AND.MLIST.LE.3) THEN 
        IF(MLIST.EQ.1) WRITE(MSTU(11),5100) 
        IF(MLIST.EQ.2) WRITE(MSTU(11),5200) 
        IF(MLIST.EQ.3) WRITE(MSTU(11),5300) 
        LMX=12 
        IF(MLIST.GE.2) LMX=16 
        ISTR=0 
        IMAX=N 
        IF(MSTU(2).GT.0) IMAX=MSTU(2) 
        DO 120 I=MAX(1,MSTU(1)),MAX(IMAX,N+MAX(0,MSTU(3))) 
        IF((I.GT.IMAX.AND.I.LE.N).OR.K(I,1).LT.0) GOTO 120 
 
C...Get particle name, pad it and check it is not too long. 
        CALL LUNAME(K(I,2),CHAP) 
        LEN=0 
        DO 100 LEM=1,16 
        IF(CHAP(LEM:LEM).NE.' ') LEN=LEM 
  100   CONTINUE 
        MDL=(K(I,1)+19)/10 
        LDL=0 
        IF(MDL.EQ.2.OR.MDL.GE.8) THEN 
          CHAC=CHAP 
          IF(LEN.GT.LMX) CHAC(LMX:LMX)='?' 
        ELSE 
          LDL=1 
          IF(MDL.EQ.1.OR.MDL.EQ.7) LDL=2 
          IF(LEN.EQ.0) THEN 
            CHAC=CHDL(MDL)(1:2*LDL)//' ' 
          ELSE 
            CHAC=CHDL(MDL)(1:LDL)//CHAP(1:MIN(LEN,LMX-2*LDL))// 
     &      CHDL(MDL)(LDL+1:2*LDL)//' ' 
            IF(LEN+2*LDL.GT.LMX) CHAC(LMX:LMX)='?' 
          ENDIF 
        ENDIF 
 
C...Add information on string connection. 
        IF(K(I,1).EQ.1.OR.K(I,1).EQ.2.OR.K(I,1).EQ.11.OR.K(I,1).EQ.12) 
     &  THEN 
          KC=LUCOMP(K(I,2)) 
          KCC=0 
          IF(KC.NE.0) KCC=KCHG(KC,2) 
          IF(IABS(K(I,2)).EQ.39) THEN 
            IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='X' 
          ELSEIF(KCC.NE.0.AND.ISTR.EQ.0) THEN 
            ISTR=1 
            IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='A' 
          ELSEIF(KCC.NE.0.AND.(K(I,1).EQ.2.OR.K(I,1).EQ.12)) THEN 
            IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='I' 
          ELSEIF(KCC.NE.0) THEN 
            ISTR=0 
            IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='V' 
          ENDIF 
        ENDIF 
 
C...Write data for particle/jet. 
        IF(MLIST.EQ.1.AND.ABS(P(I,4)).LT.9999.) THEN 
          WRITE(MSTU(11),5400) I,CHAC(1:12),(K(I,J1),J1=1,3), 
     &    (P(I,J2),J2=1,5) 
        ELSEIF(MLIST.EQ.1.AND.ABS(P(I,4)).LT.99999.) THEN 
          WRITE(MSTU(11),5500) I,CHAC(1:12),(K(I,J1),J1=1,3), 
     &    (P(I,J2),J2=1,5) 
        ELSEIF(MLIST.EQ.1) THEN 
          WRITE(MSTU(11),5600) I,CHAC(1:12),(K(I,J1),J1=1,3), 
     &    (P(I,J2),J2=1,5) 
        ELSEIF(MSTU(5).EQ.10000.AND.(K(I,1).EQ.3.OR.K(I,1).EQ.13.OR. 
     &  K(I,1).EQ.14)) THEN 
          WRITE(MSTU(11),5700) I,CHAC,(K(I,J1),J1=1,3), 
     &    K(I,4)/100000000,MOD(K(I,4)/10000,10000),MOD(K(I,4),10000), 
     &    K(I,5)/100000000,MOD(K(I,5)/10000,10000),MOD(K(I,5),10000), 
     &    (P(I,J2),J2=1,5) 
        ELSE 
          WRITE(MSTU(11),5800) I,CHAC,(K(I,J1),J1=1,5),(P(I,J2),J2=1,5) 
        ENDIF 
        IF(MLIST.EQ.3) WRITE(MSTU(11),5900) (V(I,J),J=1,5) 
 
C...Insert extra separator lines specified by user. 
        IF(MSTU(70).GE.1) THEN 
          ISEP=0 
          DO 110 J=1,MIN(10,MSTU(70)) 
          IF(I.EQ.MSTU(70+J)) ISEP=1 
  110     CONTINUE 
          IF(ISEP.EQ.1.AND.MLIST.EQ.1) WRITE(MSTU(11),6000) 
          IF(ISEP.EQ.1.AND.MLIST.GE.2) WRITE(MSTU(11),6100) 
        ENDIF 
  120   CONTINUE 
 
C...Sum of charges and momenta. 
        DO 130 J=1,6 
        PS(J)=PLU(0,J) 
  130   CONTINUE 
        IF(MLIST.EQ.1.AND.ABS(PS(4)).LT.9999.) THEN 
          WRITE(MSTU(11),6200) PS(6),(PS(J),J=1,5) 
        ELSEIF(MLIST.EQ.1.AND.ABS(PS(4)).LT.99999.) THEN 
          WRITE(MSTU(11),6300) PS(6),(PS(J),J=1,5) 
        ELSEIF(MLIST.EQ.1) THEN 
          WRITE(MSTU(11),6400) PS(6),(PS(J),J=1,5) 
        ELSE 
          WRITE(MSTU(11),6500) PS(6),(PS(J),J=1,5) 
        ENDIF 
 
C...Give simple list of KF codes defined in program. 
      ELSEIF(MLIST.EQ.11) THEN 
        WRITE(MSTU(11),6600) 
        DO 140 KF=1,40 
        CALL LUNAME(KF,CHAP) 
        CALL LUNAME(-KF,CHAN) 
        IF(CHAP.NE.' '.AND.CHAN.EQ.' ') WRITE(MSTU(11),6700) KF,CHAP 
        IF(CHAN.NE.' ') WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN 
  140   CONTINUE 
        DO 170 KFLS=1,3,2 
        DO 160 KFLA=1,8 
        DO 150 KFLB=1,KFLA-(3-KFLS)/2 
        KF=1000*KFLA+100*KFLB+KFLS 
        CALL LUNAME(KF,CHAP) 
        CALL LUNAME(-KF,CHAN) 
        WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN 
  150   CONTINUE 
  160   CONTINUE 
  170   CONTINUE 
        KF=130 
        CALL LUNAME(KF,CHAP) 
        WRITE(MSTU(11),6700) KF,CHAP 
        KF=310 
        CALL LUNAME(KF,CHAP) 
        WRITE(MSTU(11),6700) KF,CHAP 
        DO 200 KMUL=0,5 
        KFLS=3 
        IF(KMUL.EQ.0.OR.KMUL.EQ.3) KFLS=1 
        IF(KMUL.EQ.5) KFLS=5 
        KFLR=0 
        IF(KMUL.EQ.2.OR.KMUL.EQ.3) KFLR=1 
        IF(KMUL.EQ.4) KFLR=2 
        DO 190 KFLB=1,8 
        DO 180 KFLC=1,KFLB-1 
        KF=10000*KFLR+100*KFLB+10*KFLC+KFLS 
        CALL LUNAME(KF,CHAP) 
        CALL LUNAME(-KF,CHAN) 
        WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN 
  180   CONTINUE 
        KF=10000*KFLR+110*KFLB+KFLS 
        CALL LUNAME(KF,CHAP) 
        WRITE(MSTU(11),6700) KF,CHAP 
  190   CONTINUE 
  200 CONTINUE 
        KF=30443 
        CALL LUNAME(KF,CHAP) 
        WRITE(MSTU(11),6700) KF,CHAP 
        KF=30553 
        CALL LUNAME(KF,CHAP) 
        WRITE(MSTU(11),6700) KF,CHAP 
        DO 240 KFLSP=1,3 
        KFLS=2+2*(KFLSP/3) 
        DO 230 KFLA=1,8 
        DO 220 KFLB=1,KFLA 
        DO 210 KFLC=1,KFLB 
        IF(KFLSP.EQ.1.AND.(KFLA.EQ.KFLB.OR.KFLB.EQ.KFLC)) GOTO 210 
        IF(KFLSP.EQ.2.AND.KFLA.EQ.KFLC) GOTO 210 
        IF(KFLSP.EQ.1) KF=1000*KFLA+100*KFLC+10*KFLB+KFLS 
        IF(KFLSP.GE.2) KF=1000*KFLA+100*KFLB+10*KFLC+KFLS 
        CALL LUNAME(KF,CHAP) 
        CALL LUNAME(-KF,CHAN) 
        WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN 
  210   CONTINUE 
  220   CONTINUE 
  230   CONTINUE 
  240   CONTINUE 
 
C...List parton/particle data table. Check whether to be listed. 
      ELSEIF(MLIST.EQ.12) THEN 
        WRITE(MSTU(11),6800) 
        MSTJ24=MSTJ(24) 
        MSTJ(24)=0 
        KFMAX=30553 
        IF(MSTU(2).NE.0) KFMAX=MSTU(2) 
        DO 270 KF=MAX(1,MSTU(1)),KFMAX 
        KC=LUCOMP(KF) 
        IF(KC.EQ.0) GOTO 270 
        IF(MSTU(14).EQ.0.AND.KF.GT.100.AND.KC.LE.100) GOTO 270 
        IF(MSTU(14).GT.0.AND.KF.GT.100.AND.MAX(MOD(KF/1000,10), 
     &  MOD(KF/100,10)).GT.MSTU(14)) GOTO 270 
        IF(MSTU(14).GT.0.AND.KF.GT.100.AND.KC.EQ.90) GOTO 270 
 
C...Find particle name and mass. Print information. 
        CALL LUNAME(KF,CHAP) 
        IF(KF.LE.100.AND.CHAP.EQ.' '.AND.MDCY(KC,2).EQ.0) GOTO 270 
        CALL LUNAME(-KF,CHAN) 
        PM=ULMASS(KF) 
        WRITE(MSTU(11),6900) KF,KC,CHAP,CHAN,KCHG(KC,1),KCHG(KC,2), 
     &  KCHG(KC,3),PM,PMAS(KC,2),PMAS(KC,3),PMAS(KC,4),MDCY(KC,1) 
 
C...Particle decay: channel number, branching ration, matrix element, 
C...decay products. 
        IF(KF.GT.100.AND.KC.LE.100) GOTO 270 
        DO 260 IDC=MDCY(KC,2),MDCY(KC,2)+MDCY(KC,3)-1 
        DO 250 J=1,5 
        CALL LUNAME(KFDP(IDC,J),CHAD(J)) 
  250   CONTINUE 
        WRITE(MSTU(11),7000) IDC,MDME(IDC,1),MDME(IDC,2),BRAT(IDC), 
     &  (CHAD(J),J=1,5) 
  260   CONTINUE 
  270   CONTINUE 
        MSTJ(24)=MSTJ24 
 
C...List parameter value table. 
      ELSEIF(MLIST.EQ.13) THEN 
        WRITE(MSTU(11),7100) 
        DO 280 I=1,200 
        WRITE(MSTU(11),7200) I,MSTU(I),PARU(I),MSTJ(I),PARJ(I),PARF(I) 
  280   CONTINUE 
      ENDIF 
 
C...Format statements for output on unit MSTU(11) (by default 6). 
 5100 FORMAT(///28X,'Event listing (summary)'//4X,'I  particle/jet KS', 
     &5X,'KF orig    p_x      p_y      p_z       E        m'/) 
 5200 FORMAT(///28X,'Event listing (standard)'//4X,'I  particle/jet', 
     &'  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)', 
     &'       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/) 
 5300 FORMAT(///28X,'Event listing (with vertices)'//4X,'I  particle/j', 
     &'et  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)', 
     &'       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/73X, 
     &'V(I,1)       V(I,2)       V(I,3)       V(I,4)       V(I,5)'/) 
 5400 FORMAT(1X,I4,2X,A12,1X,I2,1X,I6,1X,I4,5F9.3) 
 5500 FORMAT(1X,I4,2X,A12,1X,I2,1X,I6,1X,I4,5F9.2) 
 5600 FORMAT(1X,I4,2X,A12,1X,I2,1X,I6,1X,I4,5F9.1) 
 5700 FORMAT(1X,I4,2X,A16,1X,I3,1X,I8,2X,I4,2(3X,I1,2I4),5F13.5) 
 5800 FORMAT(1X,I4,2X,A16,1X,I3,1X,I8,2X,I4,2(3X,I9),5F13.5) 
 5900 FORMAT(66X,5(1X,F12.3)) 
 6000 FORMAT(1X,78('=')) 
 6100 FORMAT(1X,130('=')) 
 6200 FORMAT(19X,'sum:',F6.2,5X,5F9.3) 
 6300 FORMAT(19X,'sum:',F6.2,5X,5F9.2) 
 6400 FORMAT(19X,'sum:',F6.2,5X,5F9.1) 
 6500 FORMAT(19X,'sum charge:',F6.2,3X,'sum momentum and inv. mass:', 
     &5F13.5) 
 6600 FORMAT(///20X,'List of KF codes in program'/) 
 6700 FORMAT(4X,I6,4X,A16,6X,I6,4X,A16) 
 6800 FORMAT(///30X,'Particle/parton data table'//5X,'KF',5X,'KC',4X, 
     &'particle',8X,'antiparticle',6X,'chg  col  anti',8X,'mass',7X, 
     &'width',7X,'w-cut',5X,'lifetime',1X,'decay'/11X,'IDC',1X,'on/off', 
     &1X,'ME',3X,'Br.rat.',4X,'decay products') 
 6900 FORMAT(/1X,I6,3X,I4,4X,A16,A16,3I5,1X,F12.5,2(1X,F11.5), 
     &2X,F12.5,3X,I2) 
 7000 FORMAT(10X,I4,2X,I3,2X,I3,2X,F8.5,4X,5A16) 
 7100 FORMAT(///20X,'Parameter value table'//4X,'I',3X,'MSTU(I)', 
     &8X,'PARU(I)',3X,'MSTJ(I)',8X,'PARJ(I)',8X,'PARF(I)') 
 7200 FORMAT(1X,I4,1X,I9,1X,F14.5,1X,I9,1X,F14.5,1X,F14.5) 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LULOGO 
 
C...Purpose: to write logo for JETSET and PYTHIA programs. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /LUDAT1/ 
      SAVE /PYPARS/ 
      CHARACTER MONTH(12)*3, LOGO(48)*32, REFER(22)*36, LINE*79, 
     &VERS*1, SUBV*3, DATE*2, YEAR*4 
 
C...Data on months, logo, titles, and references. 
      DATA MONTH/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep', 
     &'Oct','Nov','Dec'/ 
      DATA (LOGO(J),J=1,10)/ 
     &'PPP  Y   Y TTTTT H   H III   A  ', 
     &'P  P  Y Y    T   H   H  I   A A ', 
     &'PPP    Y     T   HHHHH  I  AAAAA', 
     &'P      Y     T   H   H  I  A   A', 
     &'P      Y     T   H   H III A   A', 
     &'JJJJ EEEE TTTTT  SSS  EEEE TTTTT', 
     &'   J E      T   S     E      T  ', 
     &'   J EEE    T    SSS  EEE    T  ', 
     &'J  J E      T       S E      T  ', 
     &' JJ  EEEE   T    SSS  EEEE   T  '/ 
      DATA (LOGO(J),J=11,29)/ 
     &'            *......*            ', 
     &'       *:::!!:::::::::::*       ', 
     &'    *::::::!!::::::::::::::*    ', 
     &'  *::::::::!!::::::::::::::::*  ', 
     &' *:::::::::!!:::::::::::::::::* ', 
     &' *:::::::::!!:::::::::::::::::* ', 
     &'  *::::::::!!::::::::::::::::*! ', 
     &'    *::::::!!::::::::::::::* !! ', 
     &'    !! *:::!!:::::::::::*    !! ', 
     &'    !!     !* -><- *         !! ', 
     &'    !!     !!                !! ', 
     &'    !!     !!                !! ', 
     &'    !!                       !! ', 
     &'    !!        ep             !! ', 
     &'    !!                       !! ', 
     &'    !!                 pp    !! ', 
     &'    !!   e+e-                !! ', 
     &'    !!                       !! ', 
     &'    !!                          '/ 
      DATA (LOGO(J),J=30,48)/ 
     &'Welcome to the Lund Monte Carlo!', 
     &'                                ', 
     &'  This is PYTHIA version x.xxx  ', 
     &'Last date of change: xx xxx 199x', 
     &'                                ', 
     &'  This is JETSET version x.xxx  ', 
     &'Last date of change: xx xxx 199x', 
     &'                                ', 
     &'          Main author:          ', 
     &'       Torbjorn Sjostrand       ', 
     &' Dept. of theoretical physics 2 ', 
     &'       University of Lund       ', 
     &'         Solvegatan 14A         ', 
     &'      S-223 62 Lund, Sweden     ', 
     &'   phone: +46 - 46 - 222 48 16  ', 
     &'   E-mail: torbjorn@thep.lu.se  ', 
     &'                                ', 
     &'  Copyright Torbjorn Sjostrand  ', 
     &'     and CERN, Geneva 1993      '/ 
      DATA (REFER(J),J=1,6)/ 
     &'The latest program versions and docu',
     &'mentation is found on WWW address   ',
     &'http://thep.lu.se/tf2/staff/torbjorn',
     &'/Welcome.html                       ',
     &'                                    ',
     &'                                    '/
      DATA (REFER(J),J=7,22)/ 
     &'When you cite these programs, priori', 
     &'ty should always be given to the    ', 
     &'latest published description. Curren', 
     &'tly this is                         ', 
     &'T. Sjostrand, Computer Physics Commu', 
     &'n. 82 (1994) 74.                    ', 
     &'The most recent long description (un', 
     &'published) is                       ', 
     &'T. Sjostrand, LU TP 95-20 and CERN-T',
     &'H.7112/93 (revised August 1995).    ', 
     &'Also remember that the programs, to ', 
     &'a large extent, represent original  ', 
     &'physics research. Other publications', 
     &' of special relevance to your       ', 
     &'studies may therefore deserve separa', 
     &'te mention.                         '/ 
 
C...Check if PYTHIA linked. 
      IF(MSTP(183)/10.NE.199) THEN 
        LOGO(32)=' Warning: PYTHIA is not loaded! ' 
        LOGO(33)='Did you remember to link PYDATA?' 
      ELSE 
        WRITE(VERS,'(I1)') MSTP(181) 
        LOGO(32)(26:26)=VERS 
        WRITE(SUBV,'(I3)') MSTP(182) 
        LOGO(32)(28:30)=SUBV 
        WRITE(DATE,'(I2)') MSTP(185) 
        LOGO(33)(22:23)=DATE 
        LOGO(33)(25:27)=MONTH(MSTP(184)) 
        WRITE(YEAR,'(I4)') MSTP(183) 
        LOGO(33)(29:32)=YEAR 
      ENDIF 
 
C...Check if JETSET linked. 
      IF(MSTU(183)/10.NE.199) THEN 
        LOGO(35)='  Error: JETSET is not loaded!  ' 
        LOGO(36)='Did you remember to link LUDATA?' 
      ELSE 
        WRITE(VERS,'(I1)') MSTU(181) 
        LOGO(35)(26:26)=VERS 
        WRITE(SUBV,'(I3)') MSTU(182) 
        LOGO(35)(28:30)=SUBV 
        WRITE(DATE,'(I2)') MSTU(185) 
        LOGO(36)(22:23)=DATE 
        LOGO(36)(25:27)=MONTH(MSTU(184)) 
        WRITE(YEAR,'(I4)') MSTU(183) 
        LOGO(36)(29:32)=YEAR 
      ENDIF 
 
C...Loop over lines in header. Define page feed and side borders. 
      DO 100 ILIN=1,48 
      LINE=' ' 
      IF(ILIN.EQ.1) THEN 
        LINE(1:1)='1' 
      ELSE 
        LINE(2:3)='**' 
        LINE(78:79)='**' 
      ENDIF 
 
C...Separator lines and logos. 
      IF(ILIN.EQ.2.OR.ILIN.EQ.3.OR.ILIN.EQ.47.OR.ILIN.EQ.48) THEN 
        LINE(4:77)='***********************************************'// 
     &  '***************************' 
      ELSEIF(ILIN.GE.6.AND.ILIN.LE.10) THEN 
        LINE(6:37)=LOGO(ILIN-5) 
        LINE(44:75)=LOGO(ILIN) 
      ELSEIF(ILIN.GE.13.AND.ILIN.LE.31) THEN 
        LINE(6:37)=LOGO(ILIN-2) 
        LINE(44:75)=LOGO(ILIN+17) 
      ELSEIF(ILIN.GE.34.AND.ILIN.LE.44) THEN 
        LINE(5:40)=REFER(2*ILIN-67) 
        LINE(41:76)=REFER(2*ILIN-66) 
      ENDIF 
 
C...Write lines to appropriate unit. 
      IF(MSTU(183)/10.EQ.199) THEN 
        WRITE(MSTU(11),'(A79)') LINE 
      ELSE 
        WRITE(*,'(A79)') LINE 
      ENDIF 
  100 CONTINUE 
 
C...Check that matching subversions are linked. 
      IF(MSTU(183)/10.EQ.199.AND.MSTP(183)/10.EQ.199) THEN 
        IF(MSTU(182).LT.MSTP(186)) WRITE(MSTU(11), 
     &  '(/'' Warning: JETSET subversion too old for PYTHIA''/)') 
        IF(MSTP(182).LT.MSTU(186)) WRITE(MSTU(11), 
     &  '(/'' Warning: PYTHIA subversion too old for JETSET''/)') 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUUPDA(MUPDA,LFN) 
 
C...Purpose: to facilitate the updating of particle and decay data. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      COMMON/LUDAT4/CHAF(500) 
      CHARACTER CHAF*8 
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/,/LUDAT4/ 
      CHARACTER CHINL*80,CHKC*4,CHVAR(19)*9,CHLIN*72, 
     &CHBLK(20)*72,CHOLD*12,CHTMP*12,CHNEW*12,CHCOM*12 
      DATA CHVAR/ 'KCHG(I,1)','KCHG(I,2)','KCHG(I,3)','PMAS(I,1)', 
     &'PMAS(I,2)','PMAS(I,3)','PMAS(I,4)','MDCY(I,1)','MDCY(I,2)', 
     &'MDCY(I,3)','MDME(I,1)','MDME(I,2)','BRAT(I)  ','KFDP(I,1)', 
     &'KFDP(I,2)','KFDP(I,3)','KFDP(I,4)','KFDP(I,5)','CHAF(I)  '/ 
 
C...Write information on file for editing. 
      IF(MSTU(12).GE.1) CALL LULIST(0) 
      IF(MUPDA.EQ.1) THEN 
        DO 110 KC=1,MSTU(6) 
        WRITE(LFN,5000) KC,CHAF(KC),(KCHG(KC,J1),J1=1,3), 
     &  (PMAS(KC,J2),J2=1,4),MDCY(KC,1) 
        DO 100 IDC=MDCY(KC,2),MDCY(KC,2)+MDCY(KC,3)-1 
        WRITE(LFN,5100) MDME(IDC,1),MDME(IDC,2),BRAT(IDC), 
     &  (KFDP(IDC,J),J=1,5) 
  100   CONTINUE 
  110   CONTINUE 
 
C...Reset variables and read information from edited file. 
      ELSEIF(MUPDA.EQ.2) THEN 
        DO 130 I=1,MSTU(7) 
        MDME(I,1)=1 
        MDME(I,2)=0 
        BRAT(I)=0. 
        DO 120 J=1,5 
        KFDP(I,J)=0 
  120   CONTINUE 
  130   CONTINUE 
        KC=0 
        IDC=0 
        NDC=0 
  140   READ(LFN,5200,END=150) CHINL 
        IF(CHINL(2:5).NE.'    ') THEN 
          CHKC=CHINL(2:5) 
          IF(KC.NE.0) THEN 
            MDCY(KC,2)=0 
            IF(NDC.NE.0) MDCY(KC,2)=IDC+1-NDC 
            MDCY(KC,3)=NDC 
          ENDIF 
          READ(CHKC,5300) KC 
          IF(KC.LE.0.OR.KC.GT.MSTU(6)) CALL LUERRM(27, 
     &    '(LUUPDA:) Read KC code illegal, KC ='//CHKC) 
          READ(CHINL,5000) KCR,CHAF(KC),(KCHG(KC,J1),J1=1,3), 
     &    (PMAS(KC,J2),J2=1,4),MDCY(KC,1) 
          NDC=0 
        ELSE 
          IDC=IDC+1 
          NDC=NDC+1 
          IF(IDC.GE.MSTU(7)) CALL LUERRM(27, 
     &    '(LUUPDA:) Decay data arrays full by KC ='//CHKC) 
          READ(CHINL,5100) MDME(IDC,1),MDME(IDC,2),BRAT(IDC), 
     &    (KFDP(IDC,J),J=1,5) 
        ENDIF 
        GOTO 140 
  150   MDCY(KC,2)=0 
        IF(NDC.NE.0) MDCY(KC,2)=IDC+1-NDC 
        MDCY(KC,3)=NDC 
 
C...Perform possible tests that new information is consistent. 
        MSTJ24=MSTJ(24) 
        MSTJ(24)=0 
        DO 180 KC=1,MSTU(6) 
        WRITE(CHKC,5300) KC 
        IF(MIN(PMAS(KC,1),PMAS(KC,2),PMAS(KC,3),PMAS(KC,1)-PMAS(KC,3), 
     &  PMAS(KC,4)).LT.0..OR.MDCY(KC,3).LT.0) CALL LUERRM(17, 
     &  '(LUUPDA:) Mass/width/life/(# channels) wrong for KC ='//CHKC) 
        BRSUM=0. 
        DO 170 IDC=MDCY(KC,2),MDCY(KC,2)+MDCY(KC,3)-1 
        IF(MDME(IDC,2).GT.80) GOTO 170 
        KQ=KCHG(KC,1) 
        PMS=PMAS(KC,1)-PMAS(KC,3)-PARJ(64) 
        MERR=0 
        DO 160 J=1,5 
        KP=KFDP(IDC,J) 
        IF(KP.EQ.0.OR.KP.EQ.81.OR.IABS(KP).EQ.82) THEN 
        ELSEIF(LUCOMP(KP).EQ.0) THEN 
          MERR=3 
        ELSE 
          KQ=KQ-LUCHGE(KP) 
          PMS=PMS-ULMASS(KP) 
        ENDIF 
  160   CONTINUE 
        IF(KQ.NE.0) MERR=MAX(2,MERR) 
        IF(KFDP(IDC,2).NE.0.AND.(KC.LE.20.OR.KC.GT.40).AND. 
     &  (KC.LE.80.OR.KC.GT.100).AND.MDME(IDC,2).NE.34.AND. 
     &  MDME(IDC,2).NE.61.AND.PMS.LT.0.) MERR=MAX(1,MERR) 
        IF(MERR.EQ.3) CALL LUERRM(17, 
     &  '(LUUPDA:) Unknown particle code in decay of KC ='//CHKC) 
        IF(MERR.EQ.2) CALL LUERRM(17, 
     &  '(LUUPDA:) Charge not conserved in decay of KC ='//CHKC) 
        IF(MERR.EQ.1) CALL LUERRM(7, 
     &  '(LUUPDA:) Kinematically unallowed decay of KC ='//CHKC) 
        BRSUM=BRSUM+BRAT(IDC) 
  170   CONTINUE 
        WRITE(CHTMP,5500) BRSUM 
        IF(ABS(BRSUM).GT.0.0005.AND.ABS(BRSUM-1.).GT.0.0005) CALL 
     &  LUERRM(7,'(LUUPDA:) Sum of branching ratios is '//CHTMP(5:12)// 
     &  ' for KC ='//CHKC) 
  180   CONTINUE 
        MSTJ(24)=MSTJ24 
 
C...Initialize writing of DATA statements for inclusion in program. 
      ELSEIF(MUPDA.EQ.3) THEN 
        DO 250 IVAR=1,19 
        NDIM=MSTU(6) 
        IF(IVAR.GE.11.AND.IVAR.LE.18) NDIM=MSTU(7) 
        NLIN=1 
        CHLIN=' ' 
        CHLIN(7:35)='DATA ('//CHVAR(IVAR)//',I=   1,    )/' 
        LLIN=35 
        CHOLD='START' 
 
C...Loop through variables for conversion to characters. 
        DO 230 IDIM=1,NDIM 
        IF(IVAR.EQ.1) WRITE(CHTMP,5400) KCHG(IDIM,1) 
        IF(IVAR.EQ.2) WRITE(CHTMP,5400) KCHG(IDIM,2) 
        IF(IVAR.EQ.3) WRITE(CHTMP,5400) KCHG(IDIM,3) 
        IF(IVAR.EQ.4) WRITE(CHTMP,5500) PMAS(IDIM,1) 
        IF(IVAR.EQ.5) WRITE(CHTMP,5500) PMAS(IDIM,2) 
        IF(IVAR.EQ.6) WRITE(CHTMP,5500) PMAS(IDIM,3) 
        IF(IVAR.EQ.7) WRITE(CHTMP,5500) PMAS(IDIM,4) 
        IF(IVAR.EQ.8) WRITE(CHTMP,5400) MDCY(IDIM,1) 
        IF(IVAR.EQ.9) WRITE(CHTMP,5400) MDCY(IDIM,2) 
        IF(IVAR.EQ.10) WRITE(CHTMP,5400) MDCY(IDIM,3) 
        IF(IVAR.EQ.11) WRITE(CHTMP,5400) MDME(IDIM,1) 
        IF(IVAR.EQ.12) WRITE(CHTMP,5400) MDME(IDIM,2) 
        IF(IVAR.EQ.13) WRITE(CHTMP,5500) BRAT(IDIM) 
        IF(IVAR.EQ.14) WRITE(CHTMP,5400) KFDP(IDIM,1) 
        IF(IVAR.EQ.15) WRITE(CHTMP,5400) KFDP(IDIM,2) 
        IF(IVAR.EQ.16) WRITE(CHTMP,5400) KFDP(IDIM,3) 
        IF(IVAR.EQ.17) WRITE(CHTMP,5400) KFDP(IDIM,4) 
        IF(IVAR.EQ.18) WRITE(CHTMP,5400) KFDP(IDIM,5) 
        IF(IVAR.EQ.19) CHTMP=CHAF(IDIM) 
 
C...Length of variable, trailing decimal zeros, quotation marks. 
        LLOW=1 
        LHIG=1 
        DO 190 LL=1,12 
        IF(CHTMP(13-LL:13-LL).NE.' ') LLOW=13-LL 
        IF(CHTMP(LL:LL).NE.' ') LHIG=LL 
  190   CONTINUE 
        CHNEW=CHTMP(LLOW:LHIG)//' ' 
        LNEW=1+LHIG-LLOW 
        IF((IVAR.GE.4.AND.IVAR.LE.7).OR.IVAR.EQ.13) THEN 
          LNEW=LNEW+1 
  200     LNEW=LNEW-1 
          IF(CHNEW(LNEW:LNEW).EQ.'0') GOTO 200 
          IF(LNEW.EQ.1) CHNEW(1:2)='0.' 
          IF(LNEW.EQ.1) LNEW=2 
        ELSEIF(IVAR.EQ.19) THEN 
          DO 210 LL=LNEW,1,-1 
          IF(CHNEW(LL:LL).EQ.'''') THEN 
            CHTMP=CHNEW 
            CHNEW=CHTMP(1:LL)//''''//CHTMP(LL+1:11) 
            LNEW=LNEW+1 
          ENDIF 
  210     CONTINUE 
          CHTMP=CHNEW 
          CHNEW(1:LNEW+2)=''''//CHTMP(1:LNEW)//'''' 
          LNEW=LNEW+2 
        ENDIF 
 
C...Form composite character string, often including repetition counter. 
        IF(CHNEW.NE.CHOLD) THEN 
          NRPT=1 
          CHOLD=CHNEW 
          CHCOM=CHNEW 
          LCOM=LNEW 
        ELSE 
          LRPT=LNEW+1 
          IF(NRPT.GE.2) LRPT=LNEW+3 
          IF(NRPT.GE.10) LRPT=LNEW+4 
          IF(NRPT.GE.100) LRPT=LNEW+5 
          IF(NRPT.GE.1000) LRPT=LNEW+6 
          LLIN=LLIN-LRPT 
          NRPT=NRPT+1 
          WRITE(CHTMP,5400) NRPT 
          LRPT=1 
          IF(NRPT.GE.10) LRPT=2 
          IF(NRPT.GE.100) LRPT=3 
          IF(NRPT.GE.1000) LRPT=4 
          CHCOM(1:LRPT+1+LNEW)=CHTMP(13-LRPT:12)//'*'//CHNEW(1:LNEW) 
          LCOM=LRPT+1+LNEW 
        ENDIF 
 
C...Add characters to end of line, to new line (after storing old line), 
C...or to new block of lines (after writing old block). 
        IF(LLIN+LCOM.LE.70) THEN 
          CHLIN(LLIN+1:LLIN+LCOM+1)=CHCOM(1:LCOM)//',' 
          LLIN=LLIN+LCOM+1 
        ELSEIF(NLIN.LE.19) THEN 
          CHLIN(LLIN+1:72)=' ' 
          CHBLK(NLIN)=CHLIN 
          NLIN=NLIN+1 
          CHLIN(6:6+LCOM+1)='&'//CHCOM(1:LCOM)//',' 
          LLIN=6+LCOM+1 
        ELSE 
          CHLIN(LLIN:72)='/'//' ' 
          CHBLK(NLIN)=CHLIN 
          WRITE(CHTMP,5400) IDIM-NRPT 
          CHBLK(1)(30:33)=CHTMP(9:12) 
          DO 220 ILIN=1,NLIN 
          WRITE(LFN,5600) CHBLK(ILIN) 
  220     CONTINUE 
          NLIN=1 
          CHLIN=' ' 
          CHLIN(7:35+LCOM+1)='DATA ('//CHVAR(IVAR)//',I=    ,    )/'// 
     &    CHCOM(1:LCOM)//',' 
          WRITE(CHTMP,5400) IDIM-NRPT+1 
          CHLIN(25:28)=CHTMP(9:12) 
          LLIN=35+LCOM+1 
        ENDIF 
  230   CONTINUE 
 
C...Write final block of lines. 
        CHLIN(LLIN:72)='/'//' ' 
        CHBLK(NLIN)=CHLIN 
        WRITE(CHTMP,5400) NDIM 
        CHBLK(1)(30:33)=CHTMP(9:12) 
        DO 240 ILIN=1,NLIN 
        WRITE(LFN,5600) CHBLK(ILIN) 
  240   CONTINUE 
  250   CONTINUE 
      ENDIF 
 
C...Formats for reading and writing particle data. 
 5000 FORMAT(1X,I4,2X,A8,3I3,3F12.5,2X,F12.5,I3) 
 5100 FORMAT(5X,2I5,F12.5,5I8) 
 5200 FORMAT(A80) 
 5300 FORMAT(I4) 
 5400 FORMAT(I12) 
 5500 FORMAT(F12.5) 
 5600 FORMAT(A72) 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      FUNCTION KLU(I,J) 
 
C...Purpose: to provide various integer-valued event related data. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Default value. For I=0 number of entries, number of stable entries 
C...or 3 times total charge. 
      KLU=0 
      IF(I.LT.0.OR.I.GT.MSTU(4).OR.J.LE.0) THEN 
      ELSEIF(I.EQ.0.AND.J.EQ.1) THEN 
        KLU=N 
      ELSEIF(I.EQ.0.AND.(J.EQ.2.OR.J.EQ.6)) THEN 
        DO 100 I1=1,N 
        IF(J.EQ.2.AND.K(I1,1).GE.1.AND.K(I1,1).LE.10) KLU=KLU+1 
        IF(J.EQ.6.AND.K(I1,1).GE.1.AND.K(I1,1).LE.10) KLU=KLU+ 
     &  LUCHGE(K(I1,2)) 
  100   CONTINUE 
      ELSEIF(I.EQ.0) THEN 
 
C...For I > 0 direct readout of K matrix or charge. 
      ELSEIF(J.LE.5) THEN 
        KLU=K(I,J) 
      ELSEIF(J.EQ.6) THEN 
        KLU=LUCHGE(K(I,2)) 
 
C...Status (existing/fragmented/decayed), parton/hadron separation. 
      ELSEIF(J.LE.8) THEN 
        IF(K(I,1).GE.1.AND.K(I,1).LE.10) KLU=1 
        IF(J.EQ.8) KLU=KLU*K(I,2) 
      ELSEIF(J.LE.12) THEN 
        KFA=IABS(K(I,2)) 
        KC=LUCOMP(KFA) 
        KQ=0 
        IF(KC.NE.0) KQ=KCHG(KC,2) 
        IF(J.EQ.9.AND.KC.NE.0.AND.KQ.NE.0) KLU=K(I,2) 
        IF(J.EQ.10.AND.KC.NE.0.AND.KQ.EQ.0) KLU=K(I,2) 
        IF(J.EQ.11) KLU=KC 
        IF(J.EQ.12) KLU=KQ*ISIGN(1,K(I,2)) 
 
C...Heaviest flavour in hadron/diquark. 
      ELSEIF(J.EQ.13) THEN 
        KFA=IABS(K(I,2)) 
        KLU=MOD(KFA/100,10)*(-1)**MOD(KFA/100,10) 
        IF(KFA.LT.10) KLU=KFA 
        IF(MOD(KFA/1000,10).NE.0) KLU=MOD(KFA/1000,10) 
        KLU=KLU*ISIGN(1,K(I,2)) 
 
C...Particle history: generation, ancestor, rank. 
      ELSEIF(J.LE.15) THEN 
        I2=I 
        I1=I 
  110   KLU=KLU+1 
        I2=I1 
        I1=K(I1,3) 
        IF(I1.GT.0.AND.K(I1,1).GT.0.AND.K(I1,1).LE.20) GOTO 110 
        IF(J.EQ.15) KLU=I2 
      ELSEIF(J.EQ.16) THEN 
        KFA=IABS(K(I,2))
        IF(K(I,1).LE.20.AND.((KFA.GE.11.AND.KFA.LE.20).OR.KFA.EQ.22.OR.        
     &  (KFA.GT.100.AND.MOD(KFA/10,10).NE.0))) THEN  
          I1=I
  120     I2=I1 
          I1=K(I1,3)
          IF(I1.GT.0) THEN
            KFAM=IABS(K(I1,2))
            ILP=1
            IF(KFAM.NE.0.AND.KFAM.LE.10) ILP=0
            IF(KFAM.EQ.21.OR.KFAM.EQ.91.OR.KFAM.EQ.92.OR.KFAM.EQ.93) 
     &      ILP=0
            IF(KFAM.GT.100.AND.MOD(KFAM/10,10).EQ.0) ILP=0
            IF(ILP.EQ.1) GOTO 120
          ENDIF
          IF(K(I1,1).EQ.12) THEN
            DO 130 I3=I1+1,I2 
            IF(K(I3,3).EQ.K(I2,3).AND.K(I3,2).NE.91.AND.K(I3,2).NE.92
     &      .AND.K(I3,2).NE.93) KLU=KLU+1
  130       CONTINUE
          ELSE
            I3=I2
  140       KLU=KLU+1
            I3=I3+1
            IF(I3.LT.N.AND.K(I3,3).EQ.K(I2,3)) GOTO 140           
          ENDIF 
        ENDIF 
 
C...Particle coming from collapsing jet system or not. 
      ELSEIF(J.EQ.17) THEN 
        I1=I 
  150   KLU=KLU+1 
        I3=I1 
        I1=K(I1,3) 
        I0=MAX(1,I1) 
        KC=LUCOMP(K(I0,2)) 
        IF(I1.EQ.0.OR.K(I0,1).LE.0.OR.K(I0,1).GT.20.OR.KC.EQ.0) THEN 
          IF(KLU.EQ.1) KLU=-1 
          IF(KLU.GT.1) KLU=0 
          RETURN 
        ENDIF 
        IF(KCHG(KC,2).EQ.0) GOTO 150 
        IF(K(I1,1).NE.12) KLU=0 
        IF(K(I1,1).NE.12) RETURN 
        I2=I1 
  160   I2=I2+1 
        IF(I2.LT.N.AND.K(I2,1).NE.11) GOTO 160 
        K3M=K(I3-1,3) 
        IF(K3M.GE.I1.AND.K3M.LE.I2) KLU=0 
        K3P=K(I3+1,3) 
        IF(I3.LT.N.AND.K3P.GE.I1.AND.K3P.LE.I2) KLU=0 
 
C...Number of decay products. Colour flow. 
      ELSEIF(J.EQ.18) THEN 
        IF(K(I,1).EQ.11.OR.K(I,1).EQ.12) KLU=MAX(0,K(I,5)-K(I,4)+1) 
        IF(K(I,4).EQ.0.OR.K(I,5).EQ.0) KLU=0 
      ELSEIF(J.LE.22) THEN 
        IF(K(I,1).NE.3.AND.K(I,1).NE.13.AND.K(I,1).NE.14) RETURN 
        IF(J.EQ.19) KLU=MOD(K(I,4)/MSTU(5),MSTU(5)) 
        IF(J.EQ.20) KLU=MOD(K(I,5)/MSTU(5),MSTU(5)) 
        IF(J.EQ.21) KLU=MOD(K(I,4),MSTU(5)) 
        IF(J.EQ.22) KLU=MOD(K(I,5),MSTU(5)) 
      ELSE 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      FUNCTION PLU(I,J) 
 
C...Purpose: to provide various real-valued event related data. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION PSUM(4) 
 
C...Set default value. For I = 0 sum of momenta or charges, 
C...or invariant mass of system. 
      PLU=0. 
      IF(I.LT.0.OR.I.GT.MSTU(4).OR.J.LE.0) THEN 
      ELSEIF(I.EQ.0.AND.J.LE.4) THEN 
        DO 100 I1=1,N 
        IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PLU=PLU+P(I1,J) 
  100   CONTINUE 
      ELSEIF(I.EQ.0.AND.J.EQ.5) THEN 
        DO 120 J1=1,4 
        PSUM(J1)=0. 
        DO 110 I1=1,N 
        IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PSUM(J1)=PSUM(J1)+P(I1,J1) 
  110   CONTINUE 
  120 CONTINUE 
        PLU=SQRT(MAX(0.,PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2)) 
      ELSEIF(I.EQ.0.AND.J.EQ.6) THEN 
        DO 130 I1=1,N 
        IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PLU=PLU+LUCHGE(K(I1,2))/3. 
  130   CONTINUE 
      ELSEIF(I.EQ.0) THEN 
 
C...Direct readout of P matrix. 
      ELSEIF(J.LE.5) THEN 
        PLU=P(I,J) 
 
C...Charge, total momentum, transverse momentum, transverse mass. 
      ELSEIF(J.LE.12) THEN 
        IF(J.EQ.6) PLU=LUCHGE(K(I,2))/3. 
        IF(J.EQ.7.OR.J.EQ.8) PLU=P(I,1)**2+P(I,2)**2+P(I,3)**2 
        IF(J.EQ.9.OR.J.EQ.10) PLU=P(I,1)**2+P(I,2)**2 
        IF(J.EQ.11.OR.J.EQ.12) PLU=P(I,5)**2+P(I,1)**2+P(I,2)**2 
        IF(J.EQ.8.OR.J.EQ.10.OR.J.EQ.12) PLU=SQRT(PLU) 
 
C...Theta and phi angle in radians or degrees. 
      ELSEIF(J.LE.16) THEN 
        IF(J.LE.14) PLU=ULANGL(P(I,3),SQRT(P(I,1)**2+P(I,2)**2)) 
        IF(J.GE.15) PLU=ULANGL(P(I,1),P(I,2)) 
        IF(J.EQ.14.OR.J.EQ.16) PLU=PLU*180./PARU(1) 
 
C...True rapidity, rapidity with pion mass, pseudorapidity. 
      ELSEIF(J.LE.19) THEN 
        PMR=0. 
        IF(J.EQ.17) PMR=P(I,5) 
        IF(J.EQ.18) PMR=ULMASS(211) 
        PR=MAX(1E-20,PMR**2+P(I,1)**2+P(I,2)**2) 
        PLU=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/SQRT(PR), 
     &  1E20)),P(I,3)) 
 
C...Energy and momentum fractions (only to be used in CM frame). 
      ELSEIF(J.LE.25) THEN 
        IF(J.EQ.20) PLU=2.*SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)/PARU(21) 
        IF(J.EQ.21) PLU=2.*P(I,3)/PARU(21) 
        IF(J.EQ.22) PLU=2.*SQRT(P(I,1)**2+P(I,2)**2)/PARU(21) 
        IF(J.EQ.23) PLU=2.*P(I,4)/PARU(21) 
        IF(J.EQ.24) PLU=(P(I,4)+P(I,3))/PARU(21) 
        IF(J.EQ.25) PLU=(P(I,4)-P(I,3))/PARU(21) 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUSPHE(SPH,APL) 
 
C...Purpose: to perform sphericity tensor analysis to give sphericity, 
C...aplanarity and the related event axes. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION SM(3,3),SV(3,3) 
 
C...Calculate matrix to be diagonalized. 
      NP=0 
      DO 110 J1=1,3 
      DO 100 J2=J1,3 
      SM(J1,J2)=0. 
  100 CONTINUE 
  110 CONTINUE 
      PS=0. 
      DO 140 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 140 
      IF(MSTU(41).GE.2) THEN 
        KC=LUCOMP(K(I,2)) 
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &  KC.EQ.18) GOTO 140 
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &  GOTO 140 
      ENDIF 
      NP=NP+1 
      PA=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2) 
      PWT=1. 
      IF(ABS(PARU(41)-2.).GT.0.001) PWT=MAX(1E-10,PA)**(PARU(41)-2.) 
      DO 130 J1=1,3 
      DO 120 J2=J1,3 
      SM(J1,J2)=SM(J1,J2)+PWT*P(I,J1)*P(I,J2) 
  120 CONTINUE 
  130 CONTINUE 
      PS=PS+PWT*PA**2 
  140 CONTINUE 
 
C...Very low multiplicities (0 or 1) not considered. 
      IF(NP.LE.1) THEN 
        CALL LUERRM(8,'(LUSPHE:) too few particles for analysis') 
        SPH=-1. 
        APL=-1. 
        RETURN 
      ENDIF 
      DO 160 J1=1,3 
      DO 150 J2=J1,3 
      SM(J1,J2)=SM(J1,J2)/PS 
  150 CONTINUE 
  160 CONTINUE 
 
C...Find eigenvalues to matrix (third degree equation). 
      SQ=(SM(1,1)*SM(2,2)+SM(1,1)*SM(3,3)+SM(2,2)*SM(3,3)-SM(1,2)**2- 
     &SM(1,3)**2-SM(2,3)**2)/3.-1./9. 
      SR=-0.5*(SQ+1./9.+SM(1,1)*SM(2,3)**2+SM(2,2)*SM(1,3)**2+SM(3,3)* 
     &SM(1,2)**2-SM(1,1)*SM(2,2)*SM(3,3))+SM(1,2)*SM(1,3)*SM(2,3)+1./27. 
      SP=COS(ACOS(MAX(MIN(SR/SQRT(-SQ**3),1.),-1.))/3.) 
      P(N+1,4)=1./3.+SQRT(-SQ)*MAX(2.*SP,SQRT(3.*(1.-SP**2))-SP) 
      P(N+3,4)=1./3.+SQRT(-SQ)*MIN(2.*SP,-SQRT(3.*(1.-SP**2))-SP) 
      P(N+2,4)=1.-P(N+1,4)-P(N+3,4) 
      IF(P(N+2,4).LT.1E-5) THEN 
        CALL LUERRM(8,'(LUSPHE:) all particles back-to-back') 
        SPH=-1. 
        APL=-1. 
        RETURN 
      ENDIF 
 
C...Find first and last eigenvector by solving equation system. 
      DO 240 I=1,3,2 
      DO 180 J1=1,3 
      SV(J1,J1)=SM(J1,J1)-P(N+I,4) 
      DO 170 J2=J1+1,3 
      SV(J1,J2)=SM(J1,J2) 
      SV(J2,J1)=SM(J1,J2) 
  170 CONTINUE 
  180 CONTINUE 
      SMAX=0. 
      DO 200 J1=1,3 
      DO 190 J2=1,3 
      IF(ABS(SV(J1,J2)).LE.SMAX) GOTO 190 
      JA=J1 
      JB=J2 
      SMAX=ABS(SV(J1,J2)) 
  190 CONTINUE 
  200 CONTINUE 
      SMAX=0. 
      DO 220 J3=JA+1,JA+2 
      J1=J3-3*((J3-1)/3) 
      RL=SV(J1,JB)/SV(JA,JB) 
      DO 210 J2=1,3 
      SV(J1,J2)=SV(J1,J2)-RL*SV(JA,J2) 
      IF(ABS(SV(J1,J2)).LE.SMAX) GOTO 210 
      JC=J1 
      SMAX=ABS(SV(J1,J2)) 
  210 CONTINUE 
  220 CONTINUE 
      JB1=JB+1-3*(JB/3) 
      JB2=JB+2-3*((JB+1)/3) 
      P(N+I,JB1)=-SV(JC,JB2) 
      P(N+I,JB2)=SV(JC,JB1) 
      P(N+I,JB)=-(SV(JA,JB1)*P(N+I,JB1)+SV(JA,JB2)*P(N+I,JB2))/ 
     &SV(JA,JB) 
      PA=SQRT(P(N+I,1)**2+P(N+I,2)**2+P(N+I,3)**2) 
      SGN=(-1.)**INT(RLU(0)+0.5) 
      DO 230 J=1,3 
      P(N+I,J)=SGN*P(N+I,J)/PA 
  230 CONTINUE 
  240 CONTINUE 
 
C...Middle axis orthogonal to other two. Fill other codes. 
      SGN=(-1.)**INT(RLU(0)+0.5) 
      P(N+2,1)=SGN*(P(N+1,2)*P(N+3,3)-P(N+1,3)*P(N+3,2)) 
      P(N+2,2)=SGN*(P(N+1,3)*P(N+3,1)-P(N+1,1)*P(N+3,3)) 
      P(N+2,3)=SGN*(P(N+1,1)*P(N+3,2)-P(N+1,2)*P(N+3,1)) 
      DO 260 I=1,3 
      K(N+I,1)=31 
      K(N+I,2)=95 
      K(N+I,3)=I 
      K(N+I,4)=0 
      K(N+I,5)=0 
      P(N+I,5)=0. 
      DO 250 J=1,5 
      V(I,J)=0. 
  250 CONTINUE 
  260 CONTINUE 
 
C...Calculate sphericity and aplanarity. Select storing option. 
      SPH=1.5*(P(N+2,4)+P(N+3,4)) 
      APL=1.5*P(N+3,4) 
      MSTU(61)=N+1 
      MSTU(62)=NP 
      IF(MSTU(43).LE.1) MSTU(3)=3 
      IF(MSTU(43).GE.2) N=N+3 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUTHRU(THR,OBL) 
 
C...Purpose: to perform thrust analysis to give thrust, oblateness 
C...and the related event axes. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION TDI(3),TPR(3) 
 
C...Take copy of particles that are to be considered in thrust analysis. 
      NP=0 
      PS=0. 
      DO 100 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 100 
      IF(MSTU(41).GE.2) THEN 
        KC=LUCOMP(K(I,2)) 
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &  KC.EQ.18) GOTO 100 
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &  GOTO 100 
      ENDIF 
      IF(N+NP+MSTU(44)+15.GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUTHRU:) no more memory left in LUJETS') 
        THR=-2. 
        OBL=-2. 
        RETURN 
      ENDIF 
      NP=NP+1 
      K(N+NP,1)=23 
      P(N+NP,1)=P(I,1) 
      P(N+NP,2)=P(I,2) 
      P(N+NP,3)=P(I,3) 
      P(N+NP,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2) 
      P(N+NP,5)=1. 
      IF(ABS(PARU(42)-1.).GT.0.001) P(N+NP,5)=P(N+NP,4)**(PARU(42)-1.) 
      PS=PS+P(N+NP,4)*P(N+NP,5) 
  100 CONTINUE 
 
C...Very low multiplicities (0 or 1) not considered. 
      IF(NP.LE.1) THEN 
        CALL LUERRM(8,'(LUTHRU:) too few particles for analysis') 
        THR=-1. 
        OBL=-1. 
        RETURN 
      ENDIF 
 
C...Loop over thrust and major. T axis along z direction in latter case. 
      DO 320 ILD=1,2 
      IF(ILD.EQ.2) THEN 
        K(N+NP+1,1)=31 
        PHI=ULANGL(P(N+NP+1,1),P(N+NP+1,2)) 
        MSTU(33)=1 
        CALL LUDBRB(N+1,N+NP+1,0.,-PHI,0D0,0D0,0D0) 
        THE=ULANGL(P(N+NP+1,3),P(N+NP+1,1)) 
        CALL LUDBRB(N+1,N+NP+1,-THE,0.,0D0,0D0,0D0) 
      ENDIF 
 
C...Find and order particles with highest p (pT for major). 
      DO 110 ILF=N+NP+4,N+NP+MSTU(44)+4 
      P(ILF,4)=0. 
  110 CONTINUE 
      DO 160 I=N+1,N+NP 
      IF(ILD.EQ.2) P(I,4)=SQRT(P(I,1)**2+P(I,2)**2) 
      DO 130 ILF=N+NP+MSTU(44)+3,N+NP+4,-1 
      IF(P(I,4).LE.P(ILF,4)) GOTO 140 
      DO 120 J=1,5 
      P(ILF+1,J)=P(ILF,J) 
  120 CONTINUE 
  130 CONTINUE 
      ILF=N+NP+3 
  140 DO 150 J=1,5 
      P(ILF+1,J)=P(I,J) 
  150 CONTINUE 
  160 CONTINUE 
 
C...Find and order initial axes with highest thrust (major). 
      DO 170 ILG=N+NP+MSTU(44)+5,N+NP+MSTU(44)+15 
      P(ILG,4)=0. 
  170 CONTINUE 
      NC=2**(MIN(MSTU(44),NP)-1) 
      DO 250 ILC=1,NC 
      DO 180 J=1,3 
      TDI(J)=0. 
  180 CONTINUE 
      DO 200 ILF=1,MIN(MSTU(44),NP) 
      SGN=P(N+NP+ILF+3,5) 
      IF(2**ILF*((ILC+2**(ILF-1)-1)/2**ILF).GE.ILC) SGN=-SGN 
      DO 190 J=1,4-ILD 
      TDI(J)=TDI(J)+SGN*P(N+NP+ILF+3,J) 
  190 CONTINUE 
  200 CONTINUE 
      TDS=TDI(1)**2+TDI(2)**2+TDI(3)**2 
      DO 220 ILG=N+NP+MSTU(44)+MIN(ILC,10)+4,N+NP+MSTU(44)+5,-1 
      IF(TDS.LE.P(ILG,4)) GOTO 230 
      DO 210 J=1,4 
      P(ILG+1,J)=P(ILG,J) 
  210 CONTINUE 
  220 CONTINUE 
      ILG=N+NP+MSTU(44)+4 
  230 DO 240 J=1,3 
      P(ILG+1,J)=TDI(J) 
  240 CONTINUE 
      P(ILG+1,4)=TDS 
  250 CONTINUE 
 
C...Iterate direction of axis until stable maximum. 
      P(N+NP+ILD,4)=0. 
      ILG=0 
  260 ILG=ILG+1 
      THP=0. 
  270 THPS=THP 
      DO 280 J=1,3 
      IF(THP.LE.1E-10) TDI(J)=P(N+NP+MSTU(44)+4+ILG,J) 
      IF(THP.GT.1E-10) TDI(J)=TPR(J) 
      TPR(J)=0. 
  280 CONTINUE 
      DO 300 I=N+1,N+NP 
      SGN=SIGN(P(I,5),TDI(1)*P(I,1)+TDI(2)*P(I,2)+TDI(3)*P(I,3)) 
      DO 290 J=1,4-ILD 
      TPR(J)=TPR(J)+SGN*P(I,J) 
  290 CONTINUE 
  300 CONTINUE 
      THP=SQRT(TPR(1)**2+TPR(2)**2+TPR(3)**2)/PS 
      IF(THP.GE.THPS+PARU(48)) GOTO 270 
 
C...Save good axis. Try new initial axis until a number of tries agree. 
      IF(THP.LT.P(N+NP+ILD,4)-PARU(48).AND.ILG.LT.MIN(10,NC)) GOTO 260 
      IF(THP.GT.P(N+NP+ILD,4)+PARU(48)) THEN 
        IAGR=0 
        SGN=(-1.)**INT(RLU(0)+0.5) 
        DO 310 J=1,3 
        P(N+NP+ILD,J)=SGN*TPR(J)/(PS*THP) 
  310   CONTINUE 
        P(N+NP+ILD,4)=THP 
        P(N+NP+ILD,5)=0. 
      ENDIF 
      IAGR=IAGR+1 
      IF(IAGR.LT.MSTU(45).AND.ILG.LT.MIN(10,NC)) GOTO 260 
  320 CONTINUE 
 
C...Find minor axis and value by orthogonality. 
      SGN=(-1.)**INT(RLU(0)+0.5) 
      P(N+NP+3,1)=-SGN*P(N+NP+2,2) 
      P(N+NP+3,2)=SGN*P(N+NP+2,1) 
      P(N+NP+3,3)=0. 
      THP=0. 
      DO 330 I=N+1,N+NP 
      THP=THP+P(I,5)*ABS(P(N+NP+3,1)*P(I,1)+P(N+NP+3,2)*P(I,2)) 
  330 CONTINUE 
      P(N+NP+3,4)=THP/PS 
      P(N+NP+3,5)=0. 
 
C...Fill axis information. Rotate back to original coordinate system. 
      DO 350 ILD=1,3 
      K(N+ILD,1)=31 
      K(N+ILD,2)=96 
      K(N+ILD,3)=ILD 
      K(N+ILD,4)=0 
      K(N+ILD,5)=0 
      DO 340 J=1,5 
      P(N+ILD,J)=P(N+NP+ILD,J) 
      V(N+ILD,J)=0. 
  340 CONTINUE 
  350 CONTINUE 
      CALL LUDBRB(N+1,N+3,THE,PHI,0D0,0D0,0D0) 
 
C...Calculate thrust and oblateness. Select storing option. 
      THR=P(N+1,4) 
      OBL=P(N+2,4)-P(N+3,4) 
      MSTU(61)=N+1 
      MSTU(62)=NP 
      IF(MSTU(43).LE.1) MSTU(3)=3 
      IF(MSTU(43).GE.2) N=N+3 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUCLUS(NJET) 
 
C...Purpose: to subdivide the particle content of an event into 
C...jets/clusters. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION PS(5) 
      SAVE NSAV,NP,PS,PSS,RINIT,NPRE,NREM 
 
C...Functions: distance measure in pT or (pseudo)mass. 
      R2T(I1,I2)=(P(I1,5)*P(I2,5)-P(I1,1)*P(I2,1)-P(I1,2)*P(I2,2)- 
     &P(I1,3)*P(I2,3))*2.*P(I1,5)*P(I2,5)/(0.0001+P(I1,5)+P(I2,5))**2 
      R2M(I1,I2)=2.*P(I1,4)*P(I2,4)*(1.-(P(I1,1)*P(I2,1)+P(I1,2)* 
     &P(I2,2)+P(I1,3)*P(I2,3))/(P(I1,5)*P(I2,5))) 
 
C...If first time, reset. If reentering, skip preliminaries. 
      IF(MSTU(48).LE.0) THEN 
        NP=0 
        DO 100 J=1,5 
        PS(J)=0. 
  100   CONTINUE 
        PSS=0. 
      ELSE 
        NJET=NSAV 
        IF(MSTU(43).GE.2) N=N-NJET 
        DO 110 I=N+1,N+NJET 
        P(I,5)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2) 
  110   CONTINUE 
        IF(MSTU(46).LE.3) R2ACC=PARU(44)**2 
        IF(MSTU(46).GE.4) R2ACC=PARU(45)*PS(5)**2 
        NLOOP=0 
        GOTO 300 
      ENDIF 
 
C...Find which particles are to be considered in cluster search. 
      DO 140 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 140 
      IF(MSTU(41).GE.2) THEN 
        KC=LUCOMP(K(I,2)) 
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &  KC.EQ.18) GOTO 140 
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &  GOTO 140 
      ENDIF 
      IF(N+2*NP.GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUCLUS:) no more memory left in LUJETS') 
        NJET=-1 
        RETURN 
      ENDIF 
 
C...Take copy of these particles, with space left for jets later on. 
      NP=NP+1 
      K(N+NP,3)=I 
      DO 120 J=1,5 
      P(N+NP,J)=P(I,J) 
  120 CONTINUE 
      IF(MSTU(42).EQ.0) P(N+NP,5)=0. 
      IF(MSTU(42).EQ.1.AND.K(I,2).NE.22) P(N+NP,5)=PMAS(101,1) 
      P(N+NP,4)=SQRT(P(N+NP,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2) 
      P(N+NP,5)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2) 
      DO 130 J=1,4 
      PS(J)=PS(J)+P(N+NP,J) 
  130 CONTINUE 
      PSS=PSS+P(N+NP,5) 
  140 CONTINUE 
      DO 160 I=N+1,N+NP 
      K(I+NP,3)=K(I,3) 
      DO 150 J=1,5 
      P(I+NP,J)=P(I,J) 
  150 CONTINUE 
  160 CONTINUE 
      PS(5)=SQRT(MAX(0.,PS(4)**2-PS(1)**2-PS(2)**2-PS(3)**2)) 
 
C...Very low multiplicities not considered. 
      IF(NP.LT.MSTU(47)) THEN 
        CALL LUERRM(8,'(LUCLUS:) too few particles for analysis') 
        NJET=-1 
        RETURN 
      ENDIF 
 
C...Find precluster configuration. If too few jets, make harder cuts. 
      NLOOP=0 
      IF(MSTU(46).LE.3) R2ACC=PARU(44)**2 
      IF(MSTU(46).GE.4) R2ACC=PARU(45)*PS(5)**2 
      RINIT=1.25*PARU(43) 
      IF(NP.LE.MSTU(47)+2) RINIT=0. 
  170 RINIT=0.8*RINIT 
      NPRE=0 
      NREM=NP 
      DO 180 I=N+NP+1,N+2*NP 
      K(I,4)=0 
  180 CONTINUE 
 
C...Sum up small momentum region. Jet if enough absolute momentum. 
      IF(MSTU(46).LE.2) THEN 
        DO 190 J=1,4 
        P(N+1,J)=0. 
  190   CONTINUE 
        DO 210 I=N+NP+1,N+2*NP 
        IF(P(I,5).GT.2.*RINIT) GOTO 210 
        NREM=NREM-1 
        K(I,4)=1 
        DO 200 J=1,4 
        P(N+1,J)=P(N+1,J)+P(I,J) 
  200   CONTINUE 
  210   CONTINUE 
        P(N+1,5)=SQRT(P(N+1,1)**2+P(N+1,2)**2+P(N+1,3)**2) 
        IF(P(N+1,5).GT.2.*RINIT) NPRE=1 
        IF(RINIT.GE.0.2*PARU(43).AND.NPRE+NREM.LT.MSTU(47)) GOTO 170 
        IF(NREM.EQ.0) GOTO 170 
      ENDIF 
 
C...Find fastest remaining particle. 
  220 NPRE=NPRE+1 
      PMAX=0. 
      DO 230 I=N+NP+1,N+2*NP 
      IF(K(I,4).NE.0.OR.P(I,5).LE.PMAX) GOTO 230 
      IMAX=I 
      PMAX=P(I,5) 
  230 CONTINUE 
      DO 240 J=1,5 
      P(N+NPRE,J)=P(IMAX,J) 
  240 CONTINUE 
      NREM=NREM-1 
      K(IMAX,4)=NPRE 
 
C...Sum up precluster around it according to pT separation. 
      IF(MSTU(46).LE.2) THEN 
        DO 260 I=N+NP+1,N+2*NP 
        IF(K(I,4).NE.0) GOTO 260 
        R2=R2T(I,IMAX) 
        IF(R2.GT.RINIT**2) GOTO 260 
        NREM=NREM-1 
        K(I,4)=NPRE 
        DO 250 J=1,4 
        P(N+NPRE,J)=P(N+NPRE,J)+P(I,J) 
  250   CONTINUE 
  260   CONTINUE 
        P(N+NPRE,5)=SQRT(P(N+NPRE,1)**2+P(N+NPRE,2)**2+P(N+NPRE,3)**2) 
 
C...Sum up precluster around it according to mass separation. 
      ELSE 
  270   IMIN=0 
        R2MIN=RINIT**2 
        DO 280 I=N+NP+1,N+2*NP 
        IF(K(I,4).NE.0) GOTO 280 
        R2=R2M(I,N+NPRE) 
        IF(R2.GE.R2MIN) GOTO 280 
        IMIN=I 
        R2MIN=R2 
  280   CONTINUE 
        IF(IMIN.NE.0) THEN 
          DO 290 J=1,4 
          P(N+NPRE,J)=P(N+NPRE,J)+P(IMIN,J) 
  290     CONTINUE 
          P(N+NPRE,5)=SQRT(P(N+NPRE,1)**2+P(N+NPRE,2)**2+P(N+NPRE,3)**2) 
          NREM=NREM-1 
          K(IMIN,4)=NPRE 
          GOTO 270 
        ENDIF 
      ENDIF 
 
C...Check if more preclusters to be found. Start over if too few. 
      IF(RINIT.GE.0.2*PARU(43).AND.NPRE+NREM.LT.MSTU(47)) GOTO 170 
      IF(NREM.GT.0) GOTO 220 
      NJET=NPRE 
 
C...Reassign all particles to nearest jet. Sum up new jet momenta. 
  300 TSAV=0. 
      PSJT=0. 
  310 IF(MSTU(46).LE.1) THEN 
        DO 330 I=N+1,N+NJET 
        DO 320 J=1,4 
        V(I,J)=0. 
  320   CONTINUE 
  330 CONTINUE 
        DO 360 I=N+NP+1,N+2*NP 
        R2MIN=PSS**2 
        DO 340 IJET=N+1,N+NJET 
        IF(P(IJET,5).LT.RINIT) GOTO 340 
        R2=R2T(I,IJET) 
        IF(R2.GE.R2MIN) GOTO 340 
        IMIN=IJET 
        R2MIN=R2 
  340   CONTINUE 
        K(I,4)=IMIN-N 
        DO 350 J=1,4 
        V(IMIN,J)=V(IMIN,J)+P(I,J) 
  350   CONTINUE 
  360   CONTINUE 
        PSJT=0. 
        DO 380 I=N+1,N+NJET 
        DO 370 J=1,4 
        P(I,J)=V(I,J) 
  370   CONTINUE 
        P(I,5)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2) 
        PSJT=PSJT+P(I,5) 
  380   CONTINUE 
      ENDIF 
 
C...Find two closest jets. 
      R2MIN=2.*MAX(R2ACC,PS(5)**2) 
      DO 400 ITRY1=N+1,N+NJET-1 
      DO 390 ITRY2=ITRY1+1,N+NJET 
      IF(MSTU(46).LE.2) R2=R2T(ITRY1,ITRY2) 
      IF(MSTU(46).GE.3) R2=R2M(ITRY1,ITRY2) 
      IF(R2.GE.R2MIN) GOTO 390 
      IMIN1=ITRY1 
      IMIN2=ITRY2 
      R2MIN=R2 
  390 CONTINUE 
  400 CONTINUE 
 
C...If allowed, join two closest jets and start over. 
      IF(NJET.GT.MSTU(47).AND.R2MIN.LT.R2ACC) THEN 
        IREC=MIN(IMIN1,IMIN2) 
        IDEL=MAX(IMIN1,IMIN2) 
        DO 410 J=1,4 
        P(IREC,J)=P(IMIN1,J)+P(IMIN2,J) 
  410   CONTINUE 
        P(IREC,5)=SQRT(P(IREC,1)**2+P(IREC,2)**2+P(IREC,3)**2) 
        DO 430 I=IDEL+1,N+NJET 
        DO 420 J=1,5 
        P(I-1,J)=P(I,J) 
  420   CONTINUE 
  430 CONTINUE 
        IF(MSTU(46).GE.2) THEN 
          DO 440 I=N+NP+1,N+2*NP 
          IORI=N+K(I,4) 
          IF(IORI.EQ.IDEL) K(I,4)=IREC-N 
          IF(IORI.GT.IDEL) K(I,4)=K(I,4)-1 
  440     CONTINUE 
        ENDIF 
        NJET=NJET-1 
        GOTO 300 
 
C...Divide up broad jet if empty cluster in list of final ones. 
      ELSEIF(NJET.EQ.MSTU(47).AND.MSTU(46).LE.1.AND.NLOOP.LE.2) THEN 
        DO 450 I=N+1,N+NJET 
        K(I,5)=0 
  450   CONTINUE 
        DO 460 I=N+NP+1,N+2*NP 
        K(N+K(I,4),5)=K(N+K(I,4),5)+1 
  460   CONTINUE 
        IEMP=0 
        DO 470 I=N+1,N+NJET 
        IF(K(I,5).EQ.0) IEMP=I 
  470   CONTINUE 
        IF(IEMP.NE.0) THEN 
          NLOOP=NLOOP+1 
          ISPL=0 
          R2MAX=0. 
          DO 480 I=N+NP+1,N+2*NP 
          IF(K(N+K(I,4),5).LE.1.OR.P(I,5).LT.RINIT) GOTO 480 
          IJET=N+K(I,4) 
          R2=R2T(I,IJET) 
          IF(R2.LE.R2MAX) GOTO 480 
          ISPL=I 
          R2MAX=R2 
  480     CONTINUE 
          IF(ISPL.NE.0) THEN 
            IJET=N+K(ISPL,4) 
            DO 490 J=1,4 
            P(IEMP,J)=P(ISPL,J) 
            P(IJET,J)=P(IJET,J)-P(ISPL,J) 
  490       CONTINUE 
            P(IEMP,5)=P(ISPL,5) 
            P(IJET,5)=SQRT(P(IJET,1)**2+P(IJET,2)**2+P(IJET,3)**2) 
            IF(NLOOP.LE.2) GOTO 300 
          ENDIF 
        ENDIF 
      ENDIF 
 
C...If generalized thrust has not yet converged, continue iteration. 
      IF(MSTU(46).LE.1.AND.NLOOP.LE.2.AND.PSJT/PSS.GT.TSAV+PARU(48)) 
     &THEN 
        TSAV=PSJT/PSS 
        GOTO 310 
      ENDIF 
 
C...Reorder jets according to energy. 
      DO 510 I=N+1,N+NJET 
      DO 500 J=1,5 
      V(I,J)=P(I,J) 
  500 CONTINUE 
  510 CONTINUE 
      DO 540 INEW=N+1,N+NJET 
      PEMAX=0. 
      DO 520 ITRY=N+1,N+NJET 
      IF(V(ITRY,4).LE.PEMAX) GOTO 520 
      IMAX=ITRY 
      PEMAX=V(ITRY,4) 
  520 CONTINUE 
      K(INEW,1)=31 
      K(INEW,2)=97 
      K(INEW,3)=INEW-N 
      K(INEW,4)=0 
      DO 530 J=1,5 
      P(INEW,J)=V(IMAX,J) 
  530 CONTINUE 
      V(IMAX,4)=-1. 
      K(IMAX,5)=INEW 
  540 CONTINUE 
 
C...Clean up particle-jet assignments and jet information. 
      DO 550 I=N+NP+1,N+2*NP 
      IORI=K(N+K(I,4),5) 
      K(I,4)=IORI-N 
      IF(K(K(I,3),1).NE.3) K(K(I,3),4)=IORI-N 
      K(IORI,4)=K(IORI,4)+1 
  550 CONTINUE 
      IEMP=0 
      PSJT=0. 
      DO 570 I=N+1,N+NJET 
      K(I,5)=0 
      PSJT=PSJT+P(I,5) 
      P(I,5)=SQRT(MAX(P(I,4)**2-P(I,5)**2,0.)) 
      DO 560 J=1,5 
      V(I,J)=0. 
  560 CONTINUE 
      IF(K(I,4).EQ.0) IEMP=I 
  570 CONTINUE 
 
C...Select storing option. Output variables. Check for failure. 
      MSTU(61)=N+1 
      MSTU(62)=NP 
      MSTU(63)=NPRE 
      PARU(61)=PS(5) 
      PARU(62)=PSJT/PSS 
      PARU(63)=SQRT(R2MIN) 
      IF(NJET.LE.1) PARU(63)=0. 
      IF(IEMP.NE.0) THEN 
        CALL LUERRM(8,'(LUCLUS:) failed to reconstruct as requested') 
        NJET=-1 
      ENDIF 
      IF(MSTU(43).LE.1) MSTU(3)=NJET 
      IF(MSTU(43).GE.2) N=N+NJET 
      NSAV=NJET 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUCELL(NJET) 
 
C...Purpose: to provide a simple way of jet finding in an eta-phi-ET 
C...coordinate frame, as used for calorimeters at hadron colliders. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Loop over all particles. Find cell that was hit by given particle. 
      PTLRAT=1./SINH(PARU(51))**2 
      NP=0 
      NC=N 
      DO 110 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 110 
      IF(P(I,1)**2+P(I,2)**2.LE.PTLRAT*P(I,3)**2) GOTO 110 
      IF(MSTU(41).GE.2) THEN 
        KC=LUCOMP(K(I,2)) 
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &  KC.EQ.18) GOTO 110 
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &  GOTO 110 
      ENDIF 
      NP=NP+1 
      PT=SQRT(P(I,1)**2+P(I,2)**2) 
      ETA=SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
      IETA=MAX(1,MIN(MSTU(51),1+INT(MSTU(51)*0.5*(ETA/PARU(51)+1.)))) 
      PHI=ULANGL(P(I,1),P(I,2)) 
      IPHI=MAX(1,MIN(MSTU(52),1+INT(MSTU(52)*0.5*(PHI/PARU(1)+1.)))) 
      IETPH=MSTU(52)*IETA+IPHI 
 
C...Add to cell already hit, or book new cell. 
      DO 100 IC=N+1,NC 
      IF(IETPH.EQ.K(IC,3)) THEN 
        K(IC,4)=K(IC,4)+1 
        P(IC,5)=P(IC,5)+PT 
        GOTO 110 
      ENDIF 
  100 CONTINUE 
      IF(NC.GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUCELL:) no more memory left in LUJETS') 
        NJET=-2 
        RETURN 
      ENDIF 
      NC=NC+1 
      K(NC,3)=IETPH 
      K(NC,4)=1 
      K(NC,5)=2 
      P(NC,1)=(PARU(51)/MSTU(51))*(2*IETA-1-MSTU(51)) 
      P(NC,2)=(PARU(1)/MSTU(52))*(2*IPHI-1-MSTU(52)) 
      P(NC,5)=PT 
  110 CONTINUE 
 
C...Smear true bin content by calorimeter resolution. 
      IF(MSTU(53).GE.1) THEN 
        DO 130 IC=N+1,NC 
        PEI=P(IC,5) 
        IF(MSTU(53).EQ.2) PEI=P(IC,5)*COSH(P(IC,1)) 
  120   PEF=PEI+PARU(55)*SQRT(-2.*LOG(MAX(1E-10,RLU(0)))*PEI)* 
     &  COS(PARU(2)*RLU(0)) 
        IF(PEF.LT.0..OR.PEF.GT.PARU(56)*PEI) GOTO 120 
        P(IC,5)=PEF 
        IF(MSTU(53).EQ.2) P(IC,5)=PEF/COSH(P(IC,1)) 
  130   CONTINUE 
      ENDIF 
 
C...Remove cells below threshold. 
      IF(PARU(58).GT.0.) THEN 
        NCC=NC 
        NC=N 
        DO 140 IC=N+1,NCC 
        IF(P(IC,5).GT.PARU(58)) THEN 
          NC=NC+1 
          K(NC,3)=K(IC,3) 
          K(NC,4)=K(IC,4) 
          K(NC,5)=K(IC,5) 
          P(NC,1)=P(IC,1) 
          P(NC,2)=P(IC,2) 
          P(NC,5)=P(IC,5) 
        ENDIF 
  140   CONTINUE 
      ENDIF 
 
C...Find initiator cell: the one with highest pT of not yet used ones. 
      NJ=NC 
  150 ETMAX=0. 
      DO 160 IC=N+1,NC 
      IF(K(IC,5).NE.2) GOTO 160 
      IF(P(IC,5).LE.ETMAX) GOTO 160 
      ICMAX=IC 
      ETA=P(IC,1) 
      PHI=P(IC,2) 
      ETMAX=P(IC,5) 
  160 CONTINUE 
      IF(ETMAX.LT.PARU(52)) GOTO 220 
      IF(NJ.GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUCELL:) no more memory left in LUJETS') 
        NJET=-2 
        RETURN 
      ENDIF 
      K(ICMAX,5)=1 
      NJ=NJ+1 
      K(NJ,4)=0 
      K(NJ,5)=1 
      P(NJ,1)=ETA 
      P(NJ,2)=PHI 
      P(NJ,3)=0. 
      P(NJ,4)=0. 
      P(NJ,5)=0. 
 
C...Sum up unused cells within required distance of initiator. 
      DO 170 IC=N+1,NC 
      IF(K(IC,5).EQ.0) GOTO 170 
      IF(ABS(P(IC,1)-ETA).GT.PARU(54)) GOTO 170 
      DPHIA=ABS(P(IC,2)-PHI) 
      IF(DPHIA.GT.PARU(54).AND.DPHIA.LT.PARU(2)-PARU(54)) GOTO 170 
      PHIC=P(IC,2) 
      IF(DPHIA.GT.PARU(1)) PHIC=PHIC+SIGN(PARU(2),PHI) 
      IF((P(IC,1)-ETA)**2+(PHIC-PHI)**2.GT.PARU(54)**2) GOTO 170 
      K(IC,5)=-K(IC,5) 
      K(NJ,4)=K(NJ,4)+K(IC,4) 
      P(NJ,3)=P(NJ,3)+P(IC,5)*P(IC,1) 
      P(NJ,4)=P(NJ,4)+P(IC,5)*PHIC 
      P(NJ,5)=P(NJ,5)+P(IC,5) 
  170 CONTINUE 
 
C...Reject cluster below minimum ET, else accept. 
      IF(P(NJ,5).LT.PARU(53)) THEN 
        NJ=NJ-1 
        DO 180 IC=N+1,NC 
        IF(K(IC,5).LT.0) K(IC,5)=-K(IC,5) 
  180   CONTINUE 
      ELSEIF(MSTU(54).LE.2) THEN 
        P(NJ,3)=P(NJ,3)/P(NJ,5) 
        P(NJ,4)=P(NJ,4)/P(NJ,5) 
        IF(ABS(P(NJ,4)).GT.PARU(1)) P(NJ,4)=P(NJ,4)-SIGN(PARU(2), 
     &  P(NJ,4)) 
        DO 190 IC=N+1,NC 
        IF(K(IC,5).LT.0) K(IC,5)=0 
  190   CONTINUE 
      ELSE 
        DO 200 J=1,4 
        P(NJ,J)=0. 
  200   CONTINUE 
        DO 210 IC=N+1,NC 
        IF(K(IC,5).GE.0) GOTO 210 
        P(NJ,1)=P(NJ,1)+P(IC,5)*COS(P(IC,2)) 
        P(NJ,2)=P(NJ,2)+P(IC,5)*SIN(P(IC,2)) 
        P(NJ,3)=P(NJ,3)+P(IC,5)*SINH(P(IC,1)) 
        P(NJ,4)=P(NJ,4)+P(IC,5)*COSH(P(IC,1)) 
        K(IC,5)=0 
  210   CONTINUE 
      ENDIF 
      GOTO 150 
 
C...Arrange clusters in falling ET sequence. 
  220 DO 250 I=1,NJ-NC 
      ETMAX=0. 
      DO 230 IJ=NC+1,NJ 
      IF(K(IJ,5).EQ.0) GOTO 230 
      IF(P(IJ,5).LT.ETMAX) GOTO 230 
      IJMAX=IJ 
      ETMAX=P(IJ,5) 
  230 CONTINUE 
      K(IJMAX,5)=0 
      K(N+I,1)=31 
      K(N+I,2)=98 
      K(N+I,3)=I 
      K(N+I,4)=K(IJMAX,4) 
      K(N+I,5)=0 
      DO 240 J=1,5 
      P(N+I,J)=P(IJMAX,J) 
      V(N+I,J)=0. 
  240 CONTINUE 
  250 CONTINUE 
      NJET=NJ-NC 
 
C...Convert to massless or massive four-vectors. 
      IF(MSTU(54).EQ.2) THEN 
        DO 260 I=N+1,N+NJET 
        ETA=P(I,3) 
        P(I,1)=P(I,5)*COS(P(I,4)) 
        P(I,2)=P(I,5)*SIN(P(I,4)) 
        P(I,3)=P(I,5)*SINH(ETA) 
        P(I,4)=P(I,5)*COSH(ETA) 
        P(I,5)=0. 
  260   CONTINUE 
      ELSEIF(MSTU(54).GE.3) THEN 
        DO 270 I=N+1,N+NJET 
        P(I,5)=SQRT(MAX(0.,P(I,4)**2-P(I,1)**2-P(I,2)**2-P(I,3)**2)) 
  270   CONTINUE 
      ENDIF 
 
C...Information about storage. 
      MSTU(61)=N+1 
      MSTU(62)=NP 
      MSTU(63)=NC-N 
      IF(MSTU(43).LE.1) MSTU(3)=NJET 
      IF(MSTU(43).GE.2) N=N+NJET 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUJMAS(PMH,PML) 
 
C...Purpose: to determine, approximately, the two jet masses that 
C...minimize the sum m_H^2 + m_L^2, a la Clavelli and Wyler. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION SM(3,3),SAX(3),PS(3,5) 
 
C...Reset. 
      NP=0 
      DO 120 J1=1,3 
      DO 100 J2=J1,3 
      SM(J1,J2)=0. 
  100 CONTINUE 
      DO 110 J2=1,4 
      PS(J1,J2)=0. 
  110 CONTINUE 
  120 CONTINUE 
      PSS=0. 
 
C...Take copy of particles that are to be considered in mass analysis. 
      DO 170 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 170 
      IF(MSTU(41).GE.2) THEN 
        KC=LUCOMP(K(I,2)) 
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &  KC.EQ.18) GOTO 170 
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &  GOTO 170 
      ENDIF 
      IF(N+NP+1.GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUJMAS:) no more memory left in LUJETS') 
        PMH=-2. 
        PML=-2. 
        RETURN 
      ENDIF 
      NP=NP+1 
      DO 130 J=1,5 
      P(N+NP,J)=P(I,J) 
  130 CONTINUE 
      IF(MSTU(42).EQ.0) P(N+NP,5)=0. 
      IF(MSTU(42).EQ.1.AND.K(I,2).NE.22) P(N+NP,5)=PMAS(101,1) 
      P(N+NP,4)=SQRT(P(N+NP,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2) 
 
C...Fill information in sphericity tensor and total momentum vector. 
      DO 150 J1=1,3 
      DO 140 J2=J1,3 
      SM(J1,J2)=SM(J1,J2)+P(I,J1)*P(I,J2) 
  140 CONTINUE 
  150 CONTINUE 
      PSS=PSS+(P(I,1)**2+P(I,2)**2+P(I,3)**2) 
      DO 160 J=1,4 
      PS(3,J)=PS(3,J)+P(N+NP,J) 
  160 CONTINUE 
  170 CONTINUE 
 
C...Very low multiplicities (0 or 1) not considered. 
      IF(NP.LE.1) THEN 
        CALL LUERRM(8,'(LUJMAS:) too few particles for analysis') 
        PMH=-1. 
        PML=-1. 
        RETURN 
      ENDIF 
      PARU(61)=SQRT(MAX(0.,PS(3,4)**2-PS(3,1)**2-PS(3,2)**2-PS(3,3)**2)) 
 
C...Find largest eigenvalue to matrix (third degree equation). 
      DO 190 J1=1,3 
      DO 180 J2=J1,3 
      SM(J1,J2)=SM(J1,J2)/PSS 
  180 CONTINUE 
  190 CONTINUE 
      SQ=(SM(1,1)*SM(2,2)+SM(1,1)*SM(3,3)+SM(2,2)*SM(3,3)-SM(1,2)**2- 
     &SM(1,3)**2-SM(2,3)**2)/3.-1./9. 
      SR=-0.5*(SQ+1./9.+SM(1,1)*SM(2,3)**2+SM(2,2)*SM(1,3)**2+SM(3,3)* 
     &SM(1,2)**2-SM(1,1)*SM(2,2)*SM(3,3))+SM(1,2)*SM(1,3)*SM(2,3)+1./27. 
      SP=COS(ACOS(MAX(MIN(SR/SQRT(-SQ**3),1.),-1.))/3.) 
      SMA=1./3.+SQRT(-SQ)*MAX(2.*SP,SQRT(3.*(1.-SP**2))-SP) 
 
C...Find largest eigenvector by solving equation system. 
      DO 210 J1=1,3 
      SM(J1,J1)=SM(J1,J1)-SMA 
      DO 200 J2=J1+1,3 
      SM(J2,J1)=SM(J1,J2) 
  200 CONTINUE 
  210 CONTINUE 
      SMAX=0. 
      DO 230 J1=1,3 
      DO 220 J2=1,3 
      IF(ABS(SM(J1,J2)).LE.SMAX) GOTO 220 
      JA=J1 
      JB=J2 
      SMAX=ABS(SM(J1,J2)) 
  220 CONTINUE 
  230 CONTINUE 
      SMAX=0. 
      DO 250 J3=JA+1,JA+2 
      J1=J3-3*((J3-1)/3) 
      RL=SM(J1,JB)/SM(JA,JB) 
      DO 240 J2=1,3 
      SM(J1,J2)=SM(J1,J2)-RL*SM(JA,J2) 
      IF(ABS(SM(J1,J2)).LE.SMAX) GOTO 240 
      JC=J1 
      SMAX=ABS(SM(J1,J2)) 
  240 CONTINUE 
  250 CONTINUE 
      JB1=JB+1-3*(JB/3) 
      JB2=JB+2-3*((JB+1)/3) 
      SAX(JB1)=-SM(JC,JB2) 
      SAX(JB2)=SM(JC,JB1) 
      SAX(JB)=-(SM(JA,JB1)*SAX(JB1)+SM(JA,JB2)*SAX(JB2))/SM(JA,JB) 
 
C...Divide particles into two initial clusters by hemisphere. 
      DO 270 I=N+1,N+NP 
      PSAX=P(I,1)*SAX(1)+P(I,2)*SAX(2)+P(I,3)*SAX(3) 
      IS=1 
      IF(PSAX.LT.0.) IS=2 
      K(I,3)=IS 
      DO 260 J=1,4 
      PS(IS,J)=PS(IS,J)+P(I,J) 
  260 CONTINUE 
  270 CONTINUE 
      PMS=MAX(1E-10,PS(1,4)**2-PS(1,1)**2-PS(1,2)**2-PS(1,3)**2)+ 
     &MAX(1E-10,PS(2,4)**2-PS(2,1)**2-PS(2,2)**2-PS(2,3)**2) 
 
C...Reassign one particle at a time; find maximum decrease of m^2 sum. 
  280 PMD=0. 
      IM=0 
      DO 290 J=1,4 
      PS(3,J)=PS(1,J)-PS(2,J) 
  290 CONTINUE 
      DO 300 I=N+1,N+NP 
      PPS=P(I,4)*PS(3,4)-P(I,1)*PS(3,1)-P(I,2)*PS(3,2)-P(I,3)*PS(3,3) 
      IF(K(I,3).EQ.1) PMDI=2.*(P(I,5)**2-PPS) 
      IF(K(I,3).EQ.2) PMDI=2.*(P(I,5)**2+PPS) 
      IF(PMDI.LT.PMD) THEN 
        PMD=PMDI 
        IM=I 
      ENDIF 
  300 CONTINUE 
 
C...Loop back if significant reduction in sum of m^2. 
      IF(PMD.LT.-PARU(48)*PMS) THEN 
        PMS=PMS+PMD 
        IS=K(IM,3) 
        DO 310 J=1,4 
        PS(IS,J)=PS(IS,J)-P(IM,J) 
        PS(3-IS,J)=PS(3-IS,J)+P(IM,J) 
  310   CONTINUE 
        K(IM,3)=3-IS 
        GOTO 280 
      ENDIF 
 
C...Final masses and output. 
      MSTU(61)=N+1 
      MSTU(62)=NP 
      PS(1,5)=SQRT(MAX(0.,PS(1,4)**2-PS(1,1)**2-PS(1,2)**2-PS(1,3)**2)) 
      PS(2,5)=SQRT(MAX(0.,PS(2,4)**2-PS(2,1)**2-PS(2,2)**2-PS(2,3)**2)) 
      PMH=MAX(PS(1,5),PS(2,5)) 
      PML=MIN(PS(1,5),PS(2,5)) 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUFOWO(H10,H20,H30,H40) 
 
C...Purpose: to calculate the first few Fox-Wolfram moments. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Copy momenta for particles and calculate H0. 
      NP=0 
      H0=0. 
      HD=0. 
      DO 110 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 110 
      IF(MSTU(41).GE.2) THEN 
        KC=LUCOMP(K(I,2)) 
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &  KC.EQ.18) GOTO 110 
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &  GOTO 110 
      ENDIF 
      IF(N+NP.GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUFOWO:) no more memory left in LUJETS') 
        H10=-1. 
        H20=-1. 
        H30=-1. 
        H40=-1. 
        RETURN 
      ENDIF 
      NP=NP+1 
      DO 100 J=1,3 
      P(N+NP,J)=P(I,J) 
  100 CONTINUE 
      P(N+NP,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2) 
      H0=H0+P(N+NP,4) 
      HD=HD+P(N+NP,4)**2 
  110 CONTINUE 
      H0=H0**2 
 
C...Very low multiplicities (0 or 1) not considered. 
      IF(NP.LE.1) THEN 
        CALL LUERRM(8,'(LUFOWO:) too few particles for analysis') 
        H10=-1. 
        H20=-1. 
        H30=-1. 
        H40=-1. 
        RETURN 
      ENDIF 
 
C...Calculate H1 - H4. 
      H10=0. 
      H20=0. 
      H30=0. 
      H40=0. 
      DO 130 I1=N+1,N+NP 
      DO 120 I2=I1+1,N+NP 
      CTHE=(P(I1,1)*P(I2,1)+P(I1,2)*P(I2,2)+P(I1,3)*P(I2,3))/ 
     &(P(I1,4)*P(I2,4)) 
      H10=H10+P(I1,4)*P(I2,4)*CTHE 
      H20=H20+P(I1,4)*P(I2,4)*(1.5*CTHE**2-0.5) 
      H30=H30+P(I1,4)*P(I2,4)*(2.5*CTHE**3-1.5*CTHE) 
      H40=H40+P(I1,4)*P(I2,4)*(4.375*CTHE**4-3.75*CTHE**2+0.375) 
  120 CONTINUE 
  130 CONTINUE 
 
C...Calculate H1/H0 - H4/H0. Output. 
      MSTU(61)=N+1 
      MSTU(62)=NP 
      H10=(HD+2.*H10)/H0 
      H20=(HD+2.*H20)/H0 
      H30=(HD+2.*H30)/H0 
      H40=(HD+2.*H40)/H0 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUTABU(MTABU) 
 
C...Purpose: to evaluate various properties of an event, with 
C...statistics accumulated during the course of the run and 
C...printed at the end. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/ 
      DIMENSION KFIS(100,2),NPIS(100,0:10),KFFS(400),NPFS(400,4), 
     &FEVFM(10,4),FM1FM(3,10,4),FM2FM(3,10,4),FMOMA(4),FMOMS(4), 
     &FEVEE(50),FE1EC(50),FE2EC(50),FE1EA(25),FE2EA(25), 
     &KFDM(8),KFDC(200,0:8),NPDC(200) 
      SAVE NEVIS,NKFIS,KFIS,NPIS,NEVFS,NPRFS,NFIFS,NCHFS,NKFFS, 
     &KFFS,NPFS,NEVFM,NMUFM,FM1FM,FM2FM,NEVEE,FE1EC,FE2EC,FE1EA, 
     &FE2EA,NEVDC,NKFDC,NREDC,KFDC,NPDC 
      CHARACTER CHAU*16,CHIS(2)*12,CHDC(8)*12 
      DATA NEVIS/0/,NKFIS/0/,NEVFS/0/,NPRFS/0/,NFIFS/0/,NCHFS/0/, 
     &NKFFS/0/,NEVFM/0/,NMUFM/0/,FM1FM/120*0./,FM2FM/120*0./, 
     &NEVEE/0/,FE1EC/50*0./,FE2EC/50*0./,FE1EA/25*0./,FE2EA/25*0./, 
     &NEVDC/0/,NKFDC/0/,NREDC/0/ 
 
C...Reset statistics on initial parton state. 
      IF(MTABU.EQ.10) THEN 
        NEVIS=0 
        NKFIS=0 
 
C...Identify and order flavour content of initial state. 
      ELSEIF(MTABU.EQ.11) THEN 
        NEVIS=NEVIS+1 
        KFM1=2*IABS(MSTU(161)) 
        IF(MSTU(161).GT.0) KFM1=KFM1-1 
        KFM2=2*IABS(MSTU(162)) 
        IF(MSTU(162).GT.0) KFM2=KFM2-1 
        KFMN=MIN(KFM1,KFM2) 
        KFMX=MAX(KFM1,KFM2) 
        DO 100 I=1,NKFIS 
        IF(KFMN.EQ.KFIS(I,1).AND.KFMX.EQ.KFIS(I,2)) THEN 
          IKFIS=-I 
          GOTO 110 
        ELSEIF(KFMN.LT.KFIS(I,1).OR.(KFMN.EQ.KFIS(I,1).AND. 
     &  KFMX.LT.KFIS(I,2))) THEN 
          IKFIS=I 
          GOTO 110 
        ENDIF 
  100   CONTINUE 
        IKFIS=NKFIS+1 
  110   IF(IKFIS.LT.0) THEN 
          IKFIS=-IKFIS 
        ELSE 
          IF(NKFIS.GE.100) RETURN 
          DO 130 I=NKFIS,IKFIS,-1 
          KFIS(I+1,1)=KFIS(I,1) 
          KFIS(I+1,2)=KFIS(I,2) 
          DO 120 J=0,10 
          NPIS(I+1,J)=NPIS(I,J) 
  120     CONTINUE 
  130   CONTINUE 
          NKFIS=NKFIS+1 
          KFIS(IKFIS,1)=KFMN 
          KFIS(IKFIS,2)=KFMX 
          DO 140 J=0,10 
          NPIS(IKFIS,J)=0 
  140     CONTINUE 
        ENDIF 
        NPIS(IKFIS,0)=NPIS(IKFIS,0)+1 
 
C...Count number of partons in initial state. 
        NP=0 
        DO 160 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GT.12) THEN 
        ELSEIF(IABS(K(I,2)).GT.80.AND.IABS(K(I,2)).LE.100) THEN 
        ELSEIF(IABS(K(I,2)).GT.100.AND.MOD(IABS(K(I,2))/10,10).NE.0) 
     &  THEN 
        ELSE 
          IM=I 
  150     IM=K(IM,3) 
          IF(IM.LE.0.OR.IM.GT.N) THEN 
            NP=NP+1 
          ELSEIF(K(IM,1).LE.0.OR.K(IM,1).GT.20) THEN 
            NP=NP+1 
          ELSEIF(IABS(K(IM,2)).GT.80.AND.IABS(K(IM,2)).LE.100) THEN 
          ELSEIF(IABS(K(IM,2)).GT.100.AND.MOD(IABS(K(IM,2))/10,10).NE.0) 
     &    THEN 
          ELSE 
            GOTO 150 
          ENDIF 
        ENDIF 
  160   CONTINUE 
        NPCO=MAX(NP,1) 
        IF(NP.GE.6) NPCO=6 
        IF(NP.GE.8) NPCO=7 
        IF(NP.GE.11) NPCO=8 
        IF(NP.GE.16) NPCO=9 
        IF(NP.GE.26) NPCO=10 
        NPIS(IKFIS,NPCO)=NPIS(IKFIS,NPCO)+1 
        MSTU(62)=NP 
 
C...Write statistics on initial parton state. 
      ELSEIF(MTABU.EQ.12) THEN 
        FAC=1./MAX(1,NEVIS) 
        WRITE(MSTU(11),5000) NEVIS 
        DO 170 I=1,NKFIS 
        KFMN=KFIS(I,1) 
        IF(KFMN.EQ.0) KFMN=KFIS(I,2) 
        KFM1=(KFMN+1)/2 
        IF(2*KFM1.EQ.KFMN) KFM1=-KFM1 
        CALL LUNAME(KFM1,CHAU) 
        CHIS(1)=CHAU(1:12) 
        IF(CHAU(13:13).NE.' ') CHIS(1)(12:12)='?' 
        KFMX=KFIS(I,2) 
        IF(KFIS(I,1).EQ.0) KFMX=0 
        KFM2=(KFMX+1)/2 
        IF(2*KFM2.EQ.KFMX) KFM2=-KFM2 
        CALL LUNAME(KFM2,CHAU) 
        CHIS(2)=CHAU(1:12) 
        IF(CHAU(13:13).NE.' ') CHIS(2)(12:12)='?' 
        WRITE(MSTU(11),5100) CHIS(1),CHIS(2),FAC*NPIS(I,0), 
     &  (NPIS(I,J)/FLOAT(NPIS(I,0)),J=1,10) 
  170   CONTINUE 
 
C...Copy statistics on initial parton state into /LUJETS/. 
      ELSEIF(MTABU.EQ.13) THEN 
        FAC=1./MAX(1,NEVIS) 
        DO 190 I=1,NKFIS 
        KFMN=KFIS(I,1) 
        IF(KFMN.EQ.0) KFMN=KFIS(I,2) 
        KFM1=(KFMN+1)/2 
        IF(2*KFM1.EQ.KFMN) KFM1=-KFM1 
        KFMX=KFIS(I,2) 
        IF(KFIS(I,1).EQ.0) KFMX=0 
        KFM2=(KFMX+1)/2 
        IF(2*KFM2.EQ.KFMX) KFM2=-KFM2 
        K(I,1)=32 
        K(I,2)=99 
        K(I,3)=KFM1 
        K(I,4)=KFM2 
        K(I,5)=NPIS(I,0) 
        DO 180 J=1,5 
        P(I,J)=FAC*NPIS(I,J) 
        V(I,J)=FAC*NPIS(I,J+5) 
  180   CONTINUE 
  190   CONTINUE 
        N=NKFIS 
        DO 200 J=1,5 
        K(N+1,J)=0 
        P(N+1,J)=0. 
        V(N+1,J)=0. 
  200   CONTINUE 
        K(N+1,1)=32 
        K(N+1,2)=99 
        K(N+1,5)=NEVIS 
        MSTU(3)=1 
 
C...Reset statistics on number of particles/partons. 
      ELSEIF(MTABU.EQ.20) THEN 
        NEVFS=0 
        NPRFS=0 
        NFIFS=0 
        NCHFS=0 
        NKFFS=0 
 
C...Identify whether particle/parton is primary or not. 
      ELSEIF(MTABU.EQ.21) THEN 
        NEVFS=NEVFS+1 
        MSTU(62)=0 
        DO 260 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GT.20.OR.K(I,1).EQ.13) GOTO 260 
        MSTU(62)=MSTU(62)+1 
        KC=LUCOMP(K(I,2)) 
        MPRI=0 
        IF(K(I,3).LE.0.OR.K(I,3).GT.N) THEN 
          MPRI=1 
        ELSEIF(K(K(I,3),1).LE.0.OR.K(K(I,3),1).GT.20) THEN 
          MPRI=1 
        ELSEIF(K(K(I,3),2).GE.91.AND.K(K(I,3),2).LE.93) THEN 
          MPRI=1 
        ELSEIF(KC.EQ.0) THEN 
        ELSEIF(K(K(I,3),1).EQ.13) THEN 
          IM=K(K(I,3),3) 
          IF(IM.LE.0.OR.IM.GT.N) THEN 
            MPRI=1 
          ELSEIF(K(IM,1).LE.0.OR.K(IM,1).GT.20) THEN 
            MPRI=1 
          ENDIF 
        ELSEIF(KCHG(KC,2).EQ.0) THEN 
          KCM=LUCOMP(K(K(I,3),2)) 
          IF(KCM.NE.0) THEN 
            IF(KCHG(KCM,2).NE.0) MPRI=1 
          ENDIF 
        ENDIF 
        IF(KC.NE.0.AND.MPRI.EQ.1) THEN 
          IF(KCHG(KC,2).EQ.0) NPRFS=NPRFS+1 
        ENDIF 
        IF(K(I,1).LE.10) THEN 
          NFIFS=NFIFS+1 
          IF(LUCHGE(K(I,2)).NE.0) NCHFS=NCHFS+1 
        ENDIF 
 
C...Fill statistics on number of particles/partons in event. 
        KFA=IABS(K(I,2)) 
        KFS=3-ISIGN(1,K(I,2))-MPRI 
        DO 210 IP=1,NKFFS 
        IF(KFA.EQ.KFFS(IP)) THEN 
          IKFFS=-IP 
          GOTO 220 
        ELSEIF(KFA.LT.KFFS(IP)) THEN 
          IKFFS=IP 
          GOTO 220 
        ENDIF 
  210   CONTINUE 
        IKFFS=NKFFS+1 
  220   IF(IKFFS.LT.0) THEN 
          IKFFS=-IKFFS 
        ELSE 
          IF(NKFFS.GE.400) RETURN 
          DO 240 IP=NKFFS,IKFFS,-1 
          KFFS(IP+1)=KFFS(IP) 
          DO 230 J=1,4 
          NPFS(IP+1,J)=NPFS(IP,J) 
  230     CONTINUE 
  240   CONTINUE 
          NKFFS=NKFFS+1 
          KFFS(IKFFS)=KFA 
          DO 250 J=1,4 
          NPFS(IKFFS,J)=0 
  250     CONTINUE 
        ENDIF 
        NPFS(IKFFS,KFS)=NPFS(IKFFS,KFS)+1 
  260   CONTINUE 
 
C...Write statistics on particle/parton composition of events. 
      ELSEIF(MTABU.EQ.22) THEN 
        FAC=1./MAX(1,NEVFS) 
        WRITE(MSTU(11),5200) NEVFS,FAC*NPRFS,FAC*NFIFS,FAC*NCHFS 
        DO 270 I=1,NKFFS 
        CALL LUNAME(KFFS(I),CHAU) 
        KC=LUCOMP(KFFS(I)) 
        MDCYF=0 
        IF(KC.NE.0) MDCYF=MDCY(KC,1) 
        WRITE(MSTU(11),5300) KFFS(I),CHAU,MDCYF,(FAC*NPFS(I,J),J=1,4), 
     &  FAC*(NPFS(I,1)+NPFS(I,2)+NPFS(I,3)+NPFS(I,4)) 
  270   CONTINUE 
 
C...Copy particle/parton composition information into /LUJETS/. 
      ELSEIF(MTABU.EQ.23) THEN 
        FAC=1./MAX(1,NEVFS) 
        DO 290 I=1,NKFFS 
        K(I,1)=32 
        K(I,2)=99 
        K(I,3)=KFFS(I) 
        K(I,4)=0 
        K(I,5)=NPFS(I,1)+NPFS(I,2)+NPFS(I,3)+NPFS(I,4) 
        DO 280 J=1,4 
        P(I,J)=FAC*NPFS(I,J) 
        V(I,J)=0. 
  280   CONTINUE 
        P(I,5)=FAC*K(I,5) 
        V(I,5)=0. 
  290   CONTINUE 
        N=NKFFS 
        DO 300 J=1,5 
        K(N+1,J)=0 
        P(N+1,J)=0. 
        V(N+1,J)=0. 
  300   CONTINUE 
        K(N+1,1)=32 
        K(N+1,2)=99 
        K(N+1,5)=NEVFS 
        P(N+1,1)=FAC*NPRFS 
        P(N+1,2)=FAC*NFIFS 
        P(N+1,3)=FAC*NCHFS 
        MSTU(3)=1 
 
C...Reset factorial moments statistics. 
      ELSEIF(MTABU.EQ.30) THEN 
        NEVFM=0 
        NMUFM=0 
        DO 330 IM=1,3 
        DO 320 IB=1,10 
        DO 310 IP=1,4 
        FM1FM(IM,IB,IP)=0. 
        FM2FM(IM,IB,IP)=0. 
  310   CONTINUE 
  320   CONTINUE 
  330   CONTINUE 
 
C...Find particles to include, with (pion,pseudo)rapidity and azimuth. 
      ELSEIF(MTABU.EQ.31) THEN 
        NEVFM=NEVFM+1 
        NLOW=N+MSTU(3) 
        NUPP=NLOW 
        DO 410 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 410 
        IF(MSTU(41).GE.2) THEN 
          KC=LUCOMP(K(I,2)) 
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &    KC.EQ.18) GOTO 410 
          IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &    GOTO 410 
        ENDIF 
        PMR=0. 
        IF(MSTU(42).EQ.1.AND.K(I,2).NE.22) PMR=ULMASS(211) 
        IF(MSTU(42).GE.2) PMR=P(I,5) 
        PR=MAX(1E-20,PMR**2+P(I,1)**2+P(I,2)**2) 
        YETA=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/SQRT(PR), 
     &  1E20)),P(I,3)) 
        IF(ABS(YETA).GT.PARU(57)) GOTO 410 
        PHI=ULANGL(P(I,1),P(I,2)) 
        IYETA=512.*(YETA+PARU(57))/(2.*PARU(57)) 
        IYETA=MAX(0,MIN(511,IYETA)) 
        IPHI=512.*(PHI+PARU(1))/PARU(2) 
        IPHI=MAX(0,MIN(511,IPHI)) 
        IYEP=0 
        DO 340 IB=0,9 
        IYEP=IYEP+4**IB*(2*MOD(IYETA/2**IB,2)+MOD(IPHI/2**IB,2)) 
  340   CONTINUE 
 
C...Order particles in (pseudo)rapidity and/or azimuth. 
        IF(NUPP.GT.MSTU(4)-5-MSTU(32)) THEN 
          CALL LUERRM(11,'(LUTABU:) no more memory left in LUJETS') 
          RETURN 
        ENDIF 
        NUPP=NUPP+1 
        IF(NUPP.EQ.NLOW+1) THEN 
          K(NUPP,1)=IYETA 
          K(NUPP,2)=IPHI 
          K(NUPP,3)=IYEP 
        ELSE 
          DO 350 I1=NUPP-1,NLOW+1,-1 
          IF(IYETA.GE.K(I1,1)) GOTO 360 
          K(I1+1,1)=K(I1,1) 
  350     CONTINUE 
  360     K(I1+1,1)=IYETA 
          DO 370 I1=NUPP-1,NLOW+1,-1 
          IF(IPHI.GE.K(I1,2)) GOTO 380 
          K(I1+1,2)=K(I1,2) 
  370     CONTINUE 
  380     K(I1+1,2)=IPHI 
          DO 390 I1=NUPP-1,NLOW+1,-1 
          IF(IYEP.GE.K(I1,3)) GOTO 400 
          K(I1+1,3)=K(I1,3) 
  390     CONTINUE 
  400     K(I1+1,3)=IYEP 
        ENDIF 
  410   CONTINUE 
        K(NUPP+1,1)=2**10 
        K(NUPP+1,2)=2**10 
        K(NUPP+1,3)=4**10 
 
C...Calculate sum of factorial moments in event. 
        DO 480 IM=1,3 
        DO 430 IB=1,10 
        DO 420 IP=1,4 
        FEVFM(IB,IP)=0. 
  420   CONTINUE 
  430   CONTINUE 
        DO 450 IB=1,10 
        IF(IM.LE.2) IBIN=2**(10-IB) 
        IF(IM.EQ.3) IBIN=4**(10-IB) 
        IAGR=K(NLOW+1,IM)/IBIN 
        NAGR=1 
        DO 440 I=NLOW+2,NUPP+1 
        ICUT=K(I,IM)/IBIN 
        IF(ICUT.EQ.IAGR) THEN 
          NAGR=NAGR+1 
        ELSE 
          IF(NAGR.EQ.1) THEN 
          ELSEIF(NAGR.EQ.2) THEN 
            FEVFM(IB,1)=FEVFM(IB,1)+2. 
          ELSEIF(NAGR.EQ.3) THEN 
            FEVFM(IB,1)=FEVFM(IB,1)+6. 
            FEVFM(IB,2)=FEVFM(IB,2)+6. 
          ELSEIF(NAGR.EQ.4) THEN 
            FEVFM(IB,1)=FEVFM(IB,1)+12. 
            FEVFM(IB,2)=FEVFM(IB,2)+24. 
            FEVFM(IB,3)=FEVFM(IB,3)+24. 
          ELSE 
            FEVFM(IB,1)=FEVFM(IB,1)+NAGR*(NAGR-1.) 
            FEVFM(IB,2)=FEVFM(IB,2)+NAGR*(NAGR-1.)*(NAGR-2.) 
            FEVFM(IB,3)=FEVFM(IB,3)+NAGR*(NAGR-1.)*(NAGR-2.)*(NAGR-3.) 
            FEVFM(IB,4)=FEVFM(IB,4)+NAGR*(NAGR-1.)*(NAGR-2.)*(NAGR-3.)* 
     &      (NAGR-4.) 
          ENDIF 
          IAGR=ICUT 
          NAGR=1 
        ENDIF 
  440   CONTINUE 
  450   CONTINUE 
 
C...Add results to total statistics. 
        DO 470 IB=10,1,-1 
        DO 460 IP=1,4 
        IF(FEVFM(1,IP).LT.0.5) THEN 
          FEVFM(IB,IP)=0. 
        ELSEIF(IM.LE.2) THEN 
          FEVFM(IB,IP)=2.**((IB-1)*IP)*FEVFM(IB,IP)/FEVFM(1,IP) 
        ELSE 
          FEVFM(IB,IP)=4.**((IB-1)*IP)*FEVFM(IB,IP)/FEVFM(1,IP) 
        ENDIF 
        FM1FM(IM,IB,IP)=FM1FM(IM,IB,IP)+FEVFM(IB,IP) 
        FM2FM(IM,IB,IP)=FM2FM(IM,IB,IP)+FEVFM(IB,IP)**2 
  460   CONTINUE 
  470   CONTINUE 
  480   CONTINUE 
        NMUFM=NMUFM+(NUPP-NLOW) 
        MSTU(62)=NUPP-NLOW 
 
C...Write accumulated statistics on factorial moments. 
      ELSEIF(MTABU.EQ.32) THEN 
        FAC=1./MAX(1,NEVFM) 
        IF(MSTU(42).LE.0) WRITE(MSTU(11),5400) NEVFM,'eta' 
        IF(MSTU(42).EQ.1) WRITE(MSTU(11),5400) NEVFM,'ypi' 
        IF(MSTU(42).GE.2) WRITE(MSTU(11),5400) NEVFM,'y  ' 
        DO 510 IM=1,3 
        WRITE(MSTU(11),5500) 
        DO 500 IB=1,10 
        BYETA=2.*PARU(57) 
        IF(IM.NE.2) BYETA=BYETA/2**(IB-1) 
        BPHI=PARU(2) 
        IF(IM.NE.1) BPHI=BPHI/2**(IB-1) 
        IF(IM.LE.2) BNAVE=FAC*NMUFM/FLOAT(2**(IB-1)) 
        IF(IM.EQ.3) BNAVE=FAC*NMUFM/FLOAT(4**(IB-1)) 
        DO 490 IP=1,4 
        FMOMA(IP)=FAC*FM1FM(IM,IB,IP) 
        FMOMS(IP)=SQRT(MAX(0.,FAC*(FAC*FM2FM(IM,IB,IP)-FMOMA(IP)**2))) 
  490   CONTINUE 
        WRITE(MSTU(11),5600) BYETA,BPHI,BNAVE,(FMOMA(IP),FMOMS(IP), 
     &  IP=1,4) 
  500   CONTINUE 
  510   CONTINUE 
 
C...Copy statistics on factorial moments into /LUJETS/. 
      ELSEIF(MTABU.EQ.33) THEN 
        FAC=1./MAX(1,NEVFM) 
        DO 540 IM=1,3 
        DO 530 IB=1,10 
        I=10*(IM-1)+IB 
        K(I,1)=32 
        K(I,2)=99 
        K(I,3)=1 
        IF(IM.NE.2) K(I,3)=2**(IB-1) 
        K(I,4)=1 
        IF(IM.NE.1) K(I,4)=2**(IB-1) 
        K(I,5)=0 
        P(I,1)=2.*PARU(57)/K(I,3) 
        V(I,1)=PARU(2)/K(I,4) 
        DO 520 IP=1,4 
        P(I,IP+1)=FAC*FM1FM(IM,IB,IP) 
        V(I,IP+1)=SQRT(MAX(0.,FAC*(FAC*FM2FM(IM,IB,IP)-P(I,IP+1)**2))) 
  520   CONTINUE 
  530   CONTINUE 
  540   CONTINUE 
        N=30 
        DO 550 J=1,5 
        K(N+1,J)=0 
        P(N+1,J)=0. 
        V(N+1,J)=0. 
  550   CONTINUE 
        K(N+1,1)=32 
        K(N+1,2)=99 
        K(N+1,5)=NEVFM 
        MSTU(3)=1 
 
C...Reset statistics on Energy-Energy Correlation. 
      ELSEIF(MTABU.EQ.40) THEN 
        NEVEE=0 
        DO 560 J=1,25 
        FE1EC(J)=0. 
        FE2EC(J)=0. 
        FE1EC(51-J)=0. 
        FE2EC(51-J)=0. 
        FE1EA(J)=0. 
        FE2EA(J)=0. 
  560   CONTINUE 
 
C...Find particles to include, with proper assumed mass. 
      ELSEIF(MTABU.EQ.41) THEN 
        NEVEE=NEVEE+1 
        NLOW=N+MSTU(3) 
        NUPP=NLOW 
        ECM=0. 
        DO 570 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 570 
        IF(MSTU(41).GE.2) THEN 
          KC=LUCOMP(K(I,2)) 
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &    KC.EQ.18) GOTO 570 
          IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &    GOTO 570 
        ENDIF 
        PMR=0. 
        IF(MSTU(42).EQ.1.AND.K(I,2).NE.22) PMR=ULMASS(211) 
        IF(MSTU(42).GE.2) PMR=P(I,5) 
        IF(NUPP.GT.MSTU(4)-5-MSTU(32)) THEN 
          CALL LUERRM(11,'(LUTABU:) no more memory left in LUJETS') 
          RETURN 
        ENDIF 
        NUPP=NUPP+1 
        P(NUPP,1)=P(I,1) 
        P(NUPP,2)=P(I,2) 
        P(NUPP,3)=P(I,3) 
        P(NUPP,4)=SQRT(PMR**2+P(I,1)**2+P(I,2)**2+P(I,3)**2) 
        P(NUPP,5)=MAX(1E-10,SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)) 
        ECM=ECM+P(NUPP,4) 
  570   CONTINUE 
        IF(NUPP.EQ.NLOW) RETURN 
 
C...Analyze Energy-Energy Correlation in event. 
        FAC=(2./ECM**2)*50./PARU(1) 
        DO 580 J=1,50 
        FEVEE(J)=0. 
  580   CONTINUE 
        DO 600 I1=NLOW+2,NUPP 
        DO 590 I2=NLOW+1,I1-1 
        CTHE=(P(I1,1)*P(I2,1)+P(I1,2)*P(I2,2)+P(I1,3)*P(I2,3))/ 
     &  (P(I1,5)*P(I2,5)) 
        THE=ACOS(MAX(-1.,MIN(1.,CTHE))) 
        ITHE=MAX(1,MIN(50,1+INT(50.*THE/PARU(1)))) 
        FEVEE(ITHE)=FEVEE(ITHE)+FAC*P(I1,4)*P(I2,4) 
  590   CONTINUE 
  600   CONTINUE 
        DO 610 J=1,25 
        FE1EC(J)=FE1EC(J)+FEVEE(J) 
        FE2EC(J)=FE2EC(J)+FEVEE(J)**2 
        FE1EC(51-J)=FE1EC(51-J)+FEVEE(51-J) 
        FE2EC(51-J)=FE2EC(51-J)+FEVEE(51-J)**2 
        FE1EA(J)=FE1EA(J)+(FEVEE(51-J)-FEVEE(J)) 
        FE2EA(J)=FE2EA(J)+(FEVEE(51-J)-FEVEE(J))**2 
  610   CONTINUE 
        MSTU(62)=NUPP-NLOW 
 
C...Write statistics on Energy-Energy Correlation. 
      ELSEIF(MTABU.EQ.42) THEN 
        FAC=1./MAX(1,NEVEE) 
        WRITE(MSTU(11),5700) NEVEE 
        DO 620 J=1,25 
        FEEC1=FAC*FE1EC(J) 
        FEES1=SQRT(MAX(0.,FAC*(FAC*FE2EC(J)-FEEC1**2))) 
        FEEC2=FAC*FE1EC(51-J) 
        FEES2=SQRT(MAX(0.,FAC*(FAC*FE2EC(51-J)-FEEC2**2))) 
        FEECA=FAC*FE1EA(J) 
        FEESA=SQRT(MAX(0.,FAC*(FAC*FE2EA(J)-FEECA**2))) 
        WRITE(MSTU(11),5800) 3.6*(J-1),3.6*J,FEEC1,FEES1,FEEC2,FEES2, 
     &  FEECA,FEESA 
  620   CONTINUE 
 
C...Copy statistics on Energy-Energy Correlation into /LUJETS/. 
      ELSEIF(MTABU.EQ.43) THEN 
        FAC=1./MAX(1,NEVEE) 
        DO 630 I=1,25 
        K(I,1)=32 
        K(I,2)=99 
        K(I,3)=0 
        K(I,4)=0 
        K(I,5)=0 
        P(I,1)=FAC*FE1EC(I) 
        V(I,1)=SQRT(MAX(0.,FAC*(FAC*FE2EC(I)-P(I,1)**2))) 
        P(I,2)=FAC*FE1EC(51-I) 
        V(I,2)=SQRT(MAX(0.,FAC*(FAC*FE2EC(51-I)-P(I,2)**2))) 
        P(I,3)=FAC*FE1EA(I) 
        V(I,3)=SQRT(MAX(0.,FAC*(FAC*FE2EA(I)-P(I,3)**2))) 
        P(I,4)=PARU(1)*(I-1)/50. 
        P(I,5)=PARU(1)*I/50. 
        V(I,4)=3.6*(I-1) 
        V(I,5)=3.6*I 
  630   CONTINUE 
        N=25 
        DO 640 J=1,5 
        K(N+1,J)=0 
        P(N+1,J)=0. 
        V(N+1,J)=0. 
  640   CONTINUE 
        K(N+1,1)=32 
        K(N+1,2)=99 
        K(N+1,5)=NEVEE 
        MSTU(3)=1 
 
C...Reset statistics on decay channels. 
      ELSEIF(MTABU.EQ.50) THEN 
        NEVDC=0 
        NKFDC=0 
        NREDC=0 
 
C...Identify and order flavour content of final state. 
      ELSEIF(MTABU.EQ.51) THEN 
        NEVDC=NEVDC+1 
        NDS=0 
        DO 670 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GE.6) GOTO 670 
        NDS=NDS+1 
        IF(NDS.GT.8) THEN 
          NREDC=NREDC+1 
          RETURN 
        ENDIF 
        KFM=2*IABS(K(I,2)) 
        IF(K(I,2).LT.0) KFM=KFM-1 
        DO 650 IDS=NDS-1,1,-1 
        IIN=IDS+1 
        IF(KFM.LT.KFDM(IDS)) GOTO 660 
        KFDM(IDS+1)=KFDM(IDS) 
  650   CONTINUE 
        IIN=1 
  660   KFDM(IIN)=KFM 
  670   CONTINUE 
 
C...Find whether old or new final state. 
        DO 690 IDC=1,NKFDC 
        IF(NDS.LT.KFDC(IDC,0)) THEN 
          IKFDC=IDC 
          GOTO 700 
        ELSEIF(NDS.EQ.KFDC(IDC,0)) THEN 
          DO 680 I=1,NDS 
          IF(KFDM(I).LT.KFDC(IDC,I)) THEN 
            IKFDC=IDC 
            GOTO 700 
          ELSEIF(KFDM(I).GT.KFDC(IDC,I)) THEN 
            GOTO 690 
          ENDIF 
  680     CONTINUE 
          IKFDC=-IDC 
          GOTO 700 
        ENDIF 
  690   CONTINUE 
        IKFDC=NKFDC+1 
  700   IF(IKFDC.LT.0) THEN 
          IKFDC=-IKFDC 
        ELSEIF(NKFDC.GE.200) THEN 
          NREDC=NREDC+1 
          RETURN 
        ELSE 
          DO 720 IDC=NKFDC,IKFDC,-1 
          NPDC(IDC+1)=NPDC(IDC) 
          DO 710 I=0,8 
          KFDC(IDC+1,I)=KFDC(IDC,I) 
  710     CONTINUE 
  720     CONTINUE 
          NKFDC=NKFDC+1 
          KFDC(IKFDC,0)=NDS 
          DO 730 I=1,NDS 
          KFDC(IKFDC,I)=KFDM(I) 
  730     CONTINUE 
          NPDC(IKFDC)=0 
        ENDIF 
        NPDC(IKFDC)=NPDC(IKFDC)+1 
 
C...Write statistics on decay channels. 
      ELSEIF(MTABU.EQ.52) THEN 
        FAC=1./MAX(1,NEVDC) 
        WRITE(MSTU(11),5900) NEVDC 
        DO 750 IDC=1,NKFDC 
        DO 740 I=1,KFDC(IDC,0) 
        KFM=KFDC(IDC,I) 
        KF=(KFM+1)/2 
        IF(2*KF.NE.KFM) KF=-KF 
        CALL LUNAME(KF,CHAU) 
        CHDC(I)=CHAU(1:12) 
        IF(CHAU(13:13).NE.' ') CHDC(I)(12:12)='?' 
  740   CONTINUE 
        WRITE(MSTU(11),6000) FAC*NPDC(IDC),(CHDC(I),I=1,KFDC(IDC,0)) 
  750   CONTINUE 
        IF(NREDC.NE.0) WRITE(MSTU(11),6100) FAC*NREDC 
 
C...Copy statistics on decay channels into /LUJETS/. 
      ELSEIF(MTABU.EQ.53) THEN 
        FAC=1./MAX(1,NEVDC) 
        DO 780 IDC=1,NKFDC 
        K(IDC,1)=32 
        K(IDC,2)=99 
        K(IDC,3)=0 
        K(IDC,4)=0 
        K(IDC,5)=KFDC(IDC,0) 
        DO 760 J=1,5 
        P(IDC,J)=0. 
        V(IDC,J)=0. 
  760   CONTINUE 
        DO 770 I=1,KFDC(IDC,0) 
        KFM=KFDC(IDC,I) 
        KF=(KFM+1)/2 
        IF(2*KF.NE.KFM) KF=-KF 
        IF(I.LE.5) P(IDC,I)=KF 
        IF(I.GE.6) V(IDC,I-5)=KF 
  770   CONTINUE 
        V(IDC,5)=FAC*NPDC(IDC) 
  780   CONTINUE 
        N=NKFDC 
        DO 790 J=1,5 
        K(N+1,J)=0 
        P(N+1,J)=0. 
        V(N+1,J)=0. 
  790   CONTINUE 
        K(N+1,1)=32 
        K(N+1,2)=99 
        K(N+1,5)=NEVDC 
        V(N+1,5)=FAC*NREDC 
        MSTU(3)=1 
      ENDIF 
 
C...Format statements for output on unit MSTU(11) (default 6). 
 5000 FORMAT(///20X,'Event statistics - initial state'/ 
     &20X,'based on an analysis of ',I6,' events'// 
     &3X,'Main flavours after',8X,'Fraction',4X,'Subfractions ', 
     &'according to fragmenting system multiplicity'/ 
     &4X,'hard interaction',24X,'1',7X,'2',7X,'3',7X,'4',7X,'5', 
     &6X,'6-7',5X,'8-10',3X,'11-15',3X,'16-25',4X,'>25'/) 
 5100 FORMAT(3X,A12,1X,A12,F10.5,1X,10F8.4) 
 5200 FORMAT(///20X,'Event statistics - final state'/ 
     &20X,'based on an analysis of ',I7,' events'// 
     &5X,'Mean primary multiplicity =',F10.4/ 
     &5X,'Mean final   multiplicity =',F10.4/ 
     &5X,'Mean charged multiplicity =',F10.4// 
     &5X,'Number of particles produced per event (directly and via ', 
     &'decays/branchings)'/ 
     &5X,'KF    Particle/jet  MDCY',10X,'Particles',13X,'Antiparticles', 
     &8X,'Total'/35X,'prim        seco        prim        seco'/) 
 5300 FORMAT(1X,I6,4X,A16,I2,5(1X,F11.6)) 
 5400 FORMAT(///20X,'Factorial moments analysis of multiplicity'/ 
     &20X,'based on an analysis of ',I6,' events'// 
     &3X,'delta-',A3,' delta-phi     <n>/bin',10X,'<F2>',18X,'<F3>', 
     &18X,'<F4>',18X,'<F5>'/35X,4('     value     error  ')) 
 5500 FORMAT(10X) 
 5600 FORMAT(2X,2F10.4,F12.4,4(F12.4,F10.4)) 
 5700 FORMAT(///20X,'Energy-Energy Correlation and Asymmetry'/ 
     &20X,'based on an analysis of ',I6,' events'// 
     &2X,'theta range',8X,'EEC(theta)',8X,'EEC(180-theta)',7X, 
     &'EECA(theta)'/2X,'in degrees ',3('      value    error')/) 
 5800 FORMAT(2X,F4.1,' - ',F4.1,3(F11.4,F9.4)) 
 5900 FORMAT(///20X,'Decay channel analysis - final state'/ 
     &20X,'based on an analysis of ',I6,' events'// 
     &2X,'Probability',10X,'Complete final state'/) 
 6000 FORMAT(2X,F9.5,5X,8(A12,1X)) 
 6100 FORMAT(2X,F9.5,5X,'into other channels (more than 8 particles ', 
     &'or table overflow)') 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUEEVT(KFL,ECM) 
 
C...Purpose: to handle the generation of an e+e- annihilation jet event. 
      IMPLICIT DOUBLE PRECISION(D) 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Check input parameters. 
      IF(MSTU(12).GE.1) CALL LULIST(0) 
      IF(KFL.LT.0.OR.KFL.GT.8) THEN 
        CALL LUERRM(16,'(LUEEVT:) called with unknown flavour code') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      IF(KFL.LE.5) ECMMIN=PARJ(127)+2.02*PARF(100+MAX(1,KFL)) 
      IF(KFL.GE.6) ECMMIN=PARJ(127)+2.02*PMAS(KFL,1) 
      IF(ECM.LT.ECMMIN) THEN 
        CALL LUERRM(16,'(LUEEVT:) called with too small CM energy') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
 
C...Check consistency of MSTJ options set. 
      IF(MSTJ(109).EQ.2.AND.MSTJ(110).NE.1) THEN 
        CALL LUERRM(6, 
     &  '(LUEEVT:) MSTJ(109) value requires MSTJ(110) = 1') 
        MSTJ(110)=1 
      ENDIF 
      IF(MSTJ(109).EQ.2.AND.MSTJ(111).NE.0) THEN 
        CALL LUERRM(6, 
     &  '(LUEEVT:) MSTJ(109) value requires MSTJ(111) = 0') 
        MSTJ(111)=0 
      ENDIF 
 
C...Initialize alpha_strong and total cross-section. 
      MSTU(111)=MSTJ(108) 
      IF(MSTJ(108).EQ.2.AND.(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.1)) 
     &MSTU(111)=1 
      PARU(112)=PARJ(121) 
      IF(MSTU(111).EQ.2) PARU(112)=PARJ(122) 
      IF(MSTJ(116).GT.0.AND.(MSTJ(116).GE.2.OR.ABS(ECM-PARJ(151)).GE. 
     &PARJ(139).OR.10*MSTJ(102)+KFL.NE.MSTJ(119))) CALL LUXTOT(KFL,ECM, 
     &XTOT) 
      IF(MSTJ(116).GE.3) MSTJ(116)=1 
      PARJ(171)=0. 
 
C...Add initial e+e- to event record (documentation only). 
      NTRY=0 
  100 NTRY=NTRY+1 
      IF(NTRY.GT.100) THEN 
        CALL LUERRM(14,'(LUEEVT:) caught in an infinite loop') 
        RETURN 
      ENDIF 
      MSTU(24)=0 
      NC=0 
      IF(MSTJ(115).GE.2) THEN 
        NC=NC+2 
        CALL LU1ENT(NC-1,11,0.5*ECM,0.,0.) 
        K(NC-1,1)=21 
        CALL LU1ENT(NC,-11,0.5*ECM,PARU(1),0.) 
        K(NC,1)=21 
      ENDIF 
 
C...Radiative photon (in initial state). 
      MK=0 
      ECMC=ECM 
      IF(MSTJ(107).GE.1.AND.MSTJ(116).GE.1) CALL LURADK(ECM,MK,PAK, 
     &THEK,PHIK,ALPK) 
      IF(MK.EQ.1) ECMC=SQRT(ECM*(ECM-2.*PAK)) 
      IF(MSTJ(115).GE.1.AND.MK.EQ.1) THEN 
        NC=NC+1 
        CALL LU1ENT(NC,22,PAK,THEK,PHIK) 
        K(NC,3)=MIN(MSTJ(115)/2,1) 
      ENDIF 
 
C...Virtual exchange boson (gamma or Z0). 
      IF(MSTJ(115).GE.3) THEN 
        NC=NC+1 
        KF=22 
        IF(MSTJ(102).EQ.2) KF=23 
        MSTU10=MSTU(10) 
        MSTU(10)=1 
        P(NC,5)=ECMC 
        CALL LU1ENT(NC,KF,ECMC,0.,0.) 
        K(NC,1)=21 
        K(NC,3)=1 
        MSTU(10)=MSTU10 
      ENDIF 
 
C...Choice of flavour and jet configuration. 
      CALL LUXKFL(KFL,ECM,ECMC,KFLC) 
      IF(KFLC.EQ.0) GOTO 100 
      CALL LUXJET(ECMC,NJET,CUT) 
      KFLN=21 
      IF(NJET.EQ.4) CALL LUX4JT(NJET,CUT,KFLC,ECMC,KFLN,X1,X2,X4, 
     &X12,X14) 
      IF(NJET.EQ.3) CALL LUX3JT(NJET,CUT,KFLC,ECMC,X1,X3) 
      IF(NJET.EQ.2) MSTJ(120)=1 
 
C...Fill jet configuration and origin. 
      IF(NJET.EQ.2.AND.MSTJ(101).NE.5) CALL LU2ENT(NC+1,KFLC,-KFLC,ECMC) 
      IF(NJET.EQ.2.AND.MSTJ(101).EQ.5) CALL LU2ENT(-(NC+1),KFLC,-KFLC, 
     &ECMC) 
      IF(NJET.EQ.3) CALL LU3ENT(NC+1,KFLC,21,-KFLC,ECMC,X1,X3) 
      IF(NJET.EQ.4.AND.KFLN.EQ.21) CALL LU4ENT(NC+1,KFLC,KFLN,KFLN, 
     &-KFLC,ECMC,X1,X2,X4,X12,X14) 
      IF(NJET.EQ.4.AND.KFLN.NE.21) CALL LU4ENT(NC+1,KFLC,-KFLN,KFLN, 
     &-KFLC,ECMC,X1,X2,X4,X12,X14) 
      IF(MSTU(24).NE.0) GOTO 100 
      DO 110 IP=NC+1,N 
      K(IP,3)=K(IP,3)+MIN(MSTJ(115)/2,1)+(MSTJ(115)/3)*(NC-1) 
  110 CONTINUE 
 
C...Angular orientation according to matrix element. 
      IF(MSTJ(106).EQ.1) THEN 
        CALL LUXDIF(NC,NJET,KFLC,ECMC,CHI,THE,PHI) 
        CALL LUDBRB(NC+1,N,0.,CHI,0D0,0D0,0D0) 
        CALL LUDBRB(NC+1,N,THE,PHI,0D0,0D0,0D0) 
      ENDIF 
 
C...Rotation and boost from radiative photon. 
      IF(MK.EQ.1) THEN 
        DBEK=-PAK/(ECM-PAK) 
        NMIN=NC+1-MSTJ(115)/3 
        CALL LUDBRB(NMIN,N,0.,-PHIK,0D0,0D0,0D0) 
        CALL LUDBRB(NMIN,N,ALPK,0.,DBEK*SIN(THEK),0D0,DBEK*COS(THEK)) 
        CALL LUDBRB(NMIN,N,0.,PHIK,0D0,0D0,0D0) 
      ENDIF 
 
C...Generate parton shower. Rearrange along strings and check. 
      IF(MSTJ(101).EQ.5) THEN 
        CALL LUSHOW(N-1,N,ECMC) 
        MSTJ14=MSTJ(14) 
        IF(MSTJ(105).EQ.-1) MSTJ(14)=-1 
        IF(MSTJ(105).GE.0) MSTU(28)=0 
        CALL LUPREP(0) 
        MSTJ(14)=MSTJ14 
        IF(MSTJ(105).GE.0.AND.MSTU(28).NE.0) GOTO 100 
      ENDIF 
 
C...Fragmentation/decay generation. Information for LUTABU. 
      IF(MSTJ(105).EQ.1) CALL LUEXEC 
      MSTU(161)=KFLC 
      MSTU(162)=-KFLC 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUXTOT(KFL,ECM,XTOT) 
 
C...Purpose: to calculate total cross-section, including initial 
C...state radiation effects. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT1/,/LUDAT2/ 
 
C...Status, (optimized) Q^2 scale, alpha_strong. 
      PARJ(151)=ECM 
      MSTJ(119)=10*MSTJ(102)+KFL 
      IF(MSTJ(111).EQ.0) THEN 
        Q2R=ECM**2 
      ELSEIF(MSTU(111).EQ.0) THEN 
        PARJ(168)=MIN(1.,MAX(PARJ(128),EXP(-12.*PARU(1)/ 
     &  ((33.-2.*MSTU(112))*PARU(111))))) 
        Q2R=PARJ(168)*ECM**2 
      ELSE 
        PARJ(168)=MIN(1.,MAX(PARJ(128),PARU(112)/ECM, 
     &  (2.*PARU(112)/ECM)**2)) 
        Q2R=PARJ(168)*ECM**2 
      ENDIF 
      ALSPI=ULALPS(Q2R)/PARU(1) 
 
C...QCD corrections factor in R. 
      IF(MSTJ(101).EQ.0.OR.MSTJ(109).EQ.1) THEN 
        RQCD=1. 
      ELSEIF(IABS(MSTJ(101)).EQ.1.AND.MSTJ(109).EQ.0) THEN 
        RQCD=1.+ALSPI 
      ELSEIF(MSTJ(109).EQ.0) THEN 
        RQCD=1.+ALSPI+(1.986-0.115*MSTU(118))*ALSPI**2 
        IF(MSTJ(111).EQ.1) RQCD=MAX(1.,RQCD+(33.-2.*MSTU(112))/12.* 
     &  LOG(PARJ(168))*ALSPI**2) 
      ELSEIF(IABS(MSTJ(101)).EQ.1) THEN 
        RQCD=1.+(3./4.)*ALSPI 
      ELSE 
        RQCD=1.+(3./4.)*ALSPI-(3./32.+0.519*MSTU(118))*ALSPI**2 
      ENDIF 
 
C...Calculate Z0 width if default value not acceptable. 
      IF(MSTJ(102).GE.3) THEN 
        RVA=3.*(3.+(4.*PARU(102)-1.)**2)+6.*RQCD*(2.+(1.-8.*PARU(102)/ 
     &  3.)**2+(4.*PARU(102)/3.-1.)**2) 
        DO 100 KFLC=5,6 
        VQ=1. 
        IF(MOD(MSTJ(103),2).EQ.1) VQ=SQRT(MAX(0.,1.-(2.*ULMASS(KFLC)/ 
     &  ECM)**2)) 
        IF(KFLC.EQ.5) VF=4.*PARU(102)/3.-1. 
        IF(KFLC.EQ.6) VF=1.-8.*PARU(102)/3. 
        RVA=RVA+3.*RQCD*(0.5*VQ*(3.-VQ**2)*VF**2+VQ**3) 
  100   CONTINUE 
        PARJ(124)=PARU(101)*PARJ(123)*RVA/(48.*PARU(102)*(1.-PARU(102))) 
      ENDIF 
 
C...Calculate propagator and related constants for QFD case. 
      POLL=1.-PARJ(131)*PARJ(132) 
      IF(MSTJ(102).GE.2) THEN 
        SFF=1./(16.*PARU(102)*(1.-PARU(102))) 
        SFW=ECM**4/((ECM**2-PARJ(123)**2)**2+(PARJ(123)*PARJ(124))**2) 
        SFI=SFW*(1.-(PARJ(123)/ECM)**2) 
        VE=4.*PARU(102)-1. 
        SF1I=SFF*(VE*POLL+PARJ(132)-PARJ(131)) 
        SF1W=SFF**2*((VE**2+1.)*POLL+2.*VE*(PARJ(132)-PARJ(131))) 
        HF1I=SFI*SF1I 
        HF1W=SFW*SF1W 
      ENDIF 
 
C...Loop over different flavours: charge, velocity. 
      RTOT=0. 
      RQQ=0. 
      RQV=0. 
      RVA=0. 
      DO 110 KFLC=1,MAX(MSTJ(104),KFL) 
      IF(KFL.GT.0.AND.KFLC.NE.KFL) GOTO 110 
      MSTJ(93)=1 
      PMQ=ULMASS(KFLC) 
      IF(ECM.LT.2.*PMQ+PARJ(127)) GOTO 110 
      QF=KCHG(KFLC,1)/3. 
      VQ=1. 
      IF(MOD(MSTJ(103),2).EQ.1) VQ=SQRT(1.-(2.*PMQ/ECM)**2) 
 
C...Calculate R and sum of charges for QED or QFD case. 
      RQQ=RQQ+3.*QF**2*POLL 
      IF(MSTJ(102).LE.1) THEN 
        RTOT=RTOT+3.*0.5*VQ*(3.-VQ**2)*QF**2*POLL 
      ELSE 
        VF=SIGN(1.,QF)-4.*QF*PARU(102) 
        RQV=RQV-6.*QF*VF*SF1I 
        RVA=RVA+3.*(VF**2+1.)*SF1W 
        RTOT=RTOT+3.*(0.5*VQ*(3.-VQ**2)*(QF**2*POLL-2.*QF*VF*HF1I+ 
     &  VF**2*HF1W)+VQ**3*HF1W) 
      ENDIF 
  110 CONTINUE 
      RSUM=RQQ 
      IF(MSTJ(102).GE.2) RSUM=RQQ+SFI*RQV+SFW*RVA 
 
C...Calculate cross-section, including QCD corrections. 
      PARJ(141)=RQQ 
      PARJ(142)=RTOT 
      PARJ(143)=RTOT*RQCD 
      PARJ(144)=PARJ(143) 
      PARJ(145)=PARJ(141)*86.8/ECM**2 
      PARJ(146)=PARJ(142)*86.8/ECM**2 
      PARJ(147)=PARJ(143)*86.8/ECM**2 
      PARJ(148)=PARJ(147) 
      PARJ(157)=RSUM*RQCD 
      PARJ(158)=0. 
      PARJ(159)=0. 
      XTOT=PARJ(147) 
      IF(MSTJ(107).LE.0) RETURN 
 
C...Virtual cross-section. 
      XKL=PARJ(135) 
      XKU=MIN(PARJ(136),1.-(2.*PARJ(127)/ECM)**2) 
      ALE=2.*LOG(ECM/ULMASS(11))-1. 
      SIGV=ALE/3.+2.*LOG(ECM**2/(ULMASS(13)*ULMASS(15)))/3.-4./3.+ 
     &1.526*LOG(ECM**2/0.932) 
 
C...Soft and hard radiative cross-section in QED case. 
      IF(MSTJ(102).LE.1) THEN 
        SIGV=1.5*ALE-0.5+PARU(1)**2/3.+2.*SIGV 
        SIGS=ALE*(2.*LOG(XKL)-LOG(1.-XKL)-XKL) 
        SIGH=ALE*(2.*LOG(XKU/XKL)-LOG((1.-XKU)/(1.-XKL))-(XKU-XKL)) 
 
C...Soft and hard radiative cross-section in QFD case. 
      ELSE 
        SZM=1.-(PARJ(123)/ECM)**2 
        SZW=PARJ(123)*PARJ(124)/ECM**2 
        PARJ(161)=-RQQ/RSUM 
        PARJ(162)=-(RQQ+RQV+RVA)/RSUM 
        PARJ(163)=(RQV*(1.-0.5*SZM-SFI)+RVA*(1.5-SZM-SFW))/RSUM 
        PARJ(164)=(RQV*SZW**2*(1.-2.*SFW)+RVA*(2.*SFI+SZW**2-4.+3.*SZM- 
     &  SZM**2))/(SZW*RSUM) 
        SIGV=1.5*ALE-0.5+PARU(1)**2/3.+((2.*RQQ+SFI*RQV)/RSUM)*SIGV+ 
     &  (SZW*SFW*RQV/RSUM)*PARU(1)*20./9. 
        SIGS=ALE*(2.*LOG(XKL)+PARJ(161)*LOG(1.-XKL)+PARJ(162)*XKL+ 
     &  PARJ(163)*LOG(((XKL-SZM)**2+SZW**2)/(SZM**2+SZW**2))+ 
     &  PARJ(164)*(ATAN((XKL-SZM)/SZW)-ATAN(-SZM/SZW))) 
        SIGH=ALE*(2.*LOG(XKU/XKL)+PARJ(161)*LOG((1.-XKU)/(1.-XKL))+ 
     &  PARJ(162)*(XKU-XKL)+PARJ(163)*LOG(((XKU-SZM)**2+SZW**2)/ 
     &  ((XKL-SZM)**2+SZW**2))+PARJ(164)*(ATAN((XKU-SZM)/SZW)- 
     &  ATAN((XKL-SZM)/SZW))) 
      ENDIF 
 
C...Total cross-section and fraction of hard photon events. 
      PARJ(160)=SIGH/(PARU(1)/PARU(101)+SIGV+SIGS+SIGH) 
      PARJ(157)=RSUM*(1.+(PARU(101)/PARU(1))*(SIGV+SIGS+SIGH))*RQCD 
      PARJ(144)=PARJ(157) 
      PARJ(148)=PARJ(144)*86.8/ECM**2 
      XTOT=PARJ(148) 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LURADK(ECM,MK,PAK,THEK,PHIK,ALPK) 
 
C...Purpose: to generate initial state photon radiation. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
 
C...Function: cumulative hard photon spectrum in QFD case. 
      FXK(XX)=2.*LOG(XX)+PARJ(161)*LOG(1.-XX)+PARJ(162)*XX+ 
     &PARJ(163)*LOG((XX-SZM)**2+SZW**2)+PARJ(164)*ATAN((XX-SZM)/SZW) 
 
C...Determine whether radiative photon or not. 
      MK=0 
      PAK=0. 
      IF(PARJ(160).LT.RLU(0)) RETURN 
      MK=1 
 
C...Photon energy range. Find photon momentum in QED case. 
      XKL=PARJ(135) 
      XKU=MIN(PARJ(136),1.-(2.*PARJ(127)/ECM)**2) 
      IF(MSTJ(102).LE.1) THEN 
  100   XK=1./(1.+(1./XKL-1.)*((1./XKU-1.)/(1./XKL-1.))**RLU(0)) 
        IF(1.+(1.-XK)**2.LT.2.*RLU(0)) GOTO 100 
 
C...Ditto in QFD case, by numerical inversion of integrated spectrum. 
      ELSE 
        SZM=1.-(PARJ(123)/ECM)**2 
        SZW=PARJ(123)*PARJ(124)/ECM**2 
        FXKL=FXK(XKL) 
        FXKU=FXK(XKU) 
        FXKD=1E-4*(FXKU-FXKL) 
        FXKR=FXKL+RLU(0)*(FXKU-FXKL) 
        NXK=0 
  110   NXK=NXK+1 
        XK=0.5*(XKL+XKU) 
        FXKV=FXK(XK) 
        IF(FXKV.GT.FXKR) THEN 
          XKU=XK 
          FXKU=FXKV 
        ELSE 
          XKL=XK 
          FXKL=FXKV 
        ENDIF 
        IF(NXK.LT.15.AND.FXKU-FXKL.GT.FXKD) GOTO 110 
        XK=XKL+(XKU-XKL)*(FXKR-FXKL)/(FXKU-FXKL) 
      ENDIF 
      PAK=0.5*ECM*XK 
 
C...Photon polar and azimuthal angle. 
      PME=2.*(ULMASS(11)/ECM)**2 
  120 CTHM=PME*(2./PME)**RLU(0) 
      IF(1.-(XK**2*CTHM*(1.-0.5*CTHM)+2.*(1.-XK)*PME/MAX(PME, 
     &CTHM*(1.-0.5*CTHM)))/(1.+(1.-XK)**2).LT.RLU(0)) GOTO 120 
      CTHE=1.-CTHM 
      IF(RLU(0).GT.0.5) CTHE=-CTHE 
      STHE=SQRT(MAX(0.,(CTHM-PME)*(2.-CTHM))) 
      THEK=ULANGL(CTHE,STHE) 
      PHIK=PARU(2)*RLU(0) 
 
C...Rotation angle for hadronic system. 
      SGN=1. 
      IF(0.5*(2.-XK*(1.-CTHE))**2/((2.-XK)**2+(XK*CTHE)**2).GT. 
     &RLU(0)) SGN=-1. 
      ALPK=ASIN(SGN*STHE*(XK-SGN*(2.*SQRT(1.-XK)-2.+XK)*CTHE)/ 
     &(2.-XK*(1.-SGN*CTHE))) 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUXKFL(KFL,ECM,ECMC,KFLC) 
 
C...Purpose: to select flavour for produced qqbar pair. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT1/,/LUDAT2/ 
 
C...Calculate maximum weight in QED or QFD case. 
      IF(MSTJ(102).LE.1) THEN 
        RFMAX=4./9. 
      ELSE 
        POLL=1.-PARJ(131)*PARJ(132) 
        SFF=1./(16.*PARU(102)*(1.-PARU(102))) 
        SFW=ECMC**4/((ECMC**2-PARJ(123)**2)**2+(PARJ(123)*PARJ(124))**2) 
        SFI=SFW*(1.-(PARJ(123)/ECMC)**2) 
        VE=4.*PARU(102)-1. 
        HF1I=SFI*SFF*(VE*POLL+PARJ(132)-PARJ(131)) 
        HF1W=SFW*SFF**2*((VE**2+1.)*POLL+2.*VE*(PARJ(132)-PARJ(131))) 
        RFMAX=MAX(4./9.*POLL-4./3.*(1.-8.*PARU(102)/3.)*HF1I+ 
     &  ((1.-8.*PARU(102)/3.)**2+1.)*HF1W,1./9.*POLL+2./3.* 
     &  (-1.+4.*PARU(102)/3.)*HF1I+((-1.+4.*PARU(102)/3.)**2+1.)*HF1W) 
      ENDIF 
 
C...Choose flavour. Gives charge and velocity. 
      NTRY=0 
  100 NTRY=NTRY+1 
      IF(NTRY.GT.100) THEN 
        CALL LUERRM(14,'(LUXKFL:) caught in an infinite loop') 
        KFLC=0 
        RETURN 
      ENDIF 
      KFLC=KFL 
      IF(KFL.LE.0) KFLC=1+INT(MSTJ(104)*RLU(0)) 
      MSTJ(93)=1 
      PMQ=ULMASS(KFLC) 
      IF(ECM.LT.2.*PMQ+PARJ(127)) GOTO 100 
      QF=KCHG(KFLC,1)/3. 
      VQ=1. 
      IF(MOD(MSTJ(103),2).EQ.1) VQ=SQRT(MAX(0.,1.-(2.*PMQ/ECMC)**2)) 
 
C...Calculate weight in QED or QFD case. 
      IF(MSTJ(102).LE.1) THEN 
        RF=QF**2 
        RFV=0.5*VQ*(3.-VQ**2)*QF**2 
      ELSE 
        VF=SIGN(1.,QF)-4.*QF*PARU(102) 
        RF=QF**2*POLL-2.*QF*VF*HF1I+(VF**2+1.)*HF1W 
        RFV=0.5*VQ*(3.-VQ**2)*(QF**2*POLL-2.*QF*VF*HF1I+VF**2*HF1W)+ 
     &  VQ**3*HF1W 
        IF(RFV.GT.0.) PARJ(171)=MIN(1.,VQ**3*HF1W/RFV) 
      ENDIF 
 
C...Weighting or new event (radiative photon). Cross-section update. 
      IF(KFL.LE.0.AND.RF.LT.RLU(0)*RFMAX) GOTO 100 
      PARJ(158)=PARJ(158)+1. 
      IF(ECMC.LT.2.*PMQ+PARJ(127).OR.RFV.LT.RLU(0)*RF) KFLC=0 
      IF(MSTJ(107).LE.0.AND.KFLC.EQ.0) GOTO 100 
      IF(KFLC.NE.0) PARJ(159)=PARJ(159)+1. 
      PARJ(144)=PARJ(157)*PARJ(159)/PARJ(158) 
      PARJ(148)=PARJ(144)*86.8/ECM**2 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUXJET(ECM,NJET,CUT) 
 
C...Purpose: to select number of jets in matrix element approach. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      DIMENSION ZHUT(5) 
 
C...Relative three-jet rate in Zhu second order parametrization. 
      DATA ZHUT/3.0922, 6.2291, 7.4782, 7.8440, 8.2560/ 
 
C...Trivial result for two-jets only, including parton shower. 
      IF(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.5) THEN 
        CUT=0. 
 
C...QCD and Abelian vector gluon theory: Q^2 for jet rate and R. 
      ELSEIF(MSTJ(109).EQ.0.OR.MSTJ(109).EQ.2) THEN 
        CF=4./3. 
        IF(MSTJ(109).EQ.2) CF=1. 
        IF(MSTJ(111).EQ.0) THEN 
          Q2=ECM**2 
          Q2R=ECM**2 
        ELSEIF(MSTU(111).EQ.0) THEN 
          PARJ(169)=MIN(1.,PARJ(129)) 
          Q2=PARJ(169)*ECM**2 
          PARJ(168)=MIN(1.,MAX(PARJ(128),EXP(-12.*PARU(1)/ 
     &    ((33.-2.*MSTU(112))*PARU(111))))) 
          Q2R=PARJ(168)*ECM**2 
        ELSE 
          PARJ(169)=MIN(1.,MAX(PARJ(129),(2.*PARU(112)/ECM)**2)) 
          Q2=PARJ(169)*ECM**2 
          PARJ(168)=MIN(1.,MAX(PARJ(128),PARU(112)/ECM, 
     &    (2.*PARU(112)/ECM)**2)) 
          Q2R=PARJ(168)*ECM**2 
        ENDIF 
 
C...alpha_strong for R and R itself. 
        ALSPI=(3./4.)*CF*ULALPS(Q2R)/PARU(1) 
        IF(IABS(MSTJ(101)).EQ.1) THEN 
          RQCD=1.+ALSPI 
        ELSEIF(MSTJ(109).EQ.0) THEN 
          RQCD=1.+ALSPI+(1.986-0.115*MSTU(118))*ALSPI**2 
          IF(MSTJ(111).EQ.1) RQCD=MAX(1.,RQCD+(33.-2.*MSTU(112))/12.* 
     &    LOG(PARJ(168))*ALSPI**2) 
        ELSE 
          RQCD=1.+ALSPI-(3./32.+0.519*MSTU(118))*(4.*ALSPI/3.)**2 
        ENDIF 
 
C...alpha_strong for jet rate. Initial value for y cut. 
        ALSPI=(3./4.)*CF*ULALPS(Q2)/PARU(1) 
        CUT=MAX(0.001,PARJ(125),(PARJ(126)/ECM)**2) 
        IF(IABS(MSTJ(101)).LE.1.OR.(MSTJ(109).EQ.0.AND.MSTJ(111).EQ.0)) 
     &  CUT=MAX(CUT,EXP(-SQRT(0.75/ALSPI))/2.) 
        IF(MSTJ(110).EQ.2) CUT=MAX(0.01,MIN(0.05,CUT)) 
 
C...Parametrization of first order three-jet cross-section. 
  100   IF(MSTJ(101).EQ.0.OR.CUT.GE.0.25) THEN 
          PARJ(152)=0. 
        ELSE 
          PARJ(152)=(2.*ALSPI/3.)*((3.-6.*CUT+2.*LOG(CUT))* 
     &    LOG(CUT/(1.-2.*CUT))+(2.5+1.5*CUT-6.571)*(1.-3.*CUT)+ 
     &    5.833*(1.-3.*CUT)**2-3.894*(1.-3.*CUT)**3+ 
     &    1.342*(1.-3.*CUT)**4)/RQCD 
          IF(MSTJ(109).EQ.2.AND.(MSTJ(101).EQ.2.OR.MSTJ(101).LE.-2)) 
     &    PARJ(152)=0. 
        ENDIF 
 
C...Parametrization of second order three-jet cross-section. 
        IF(IABS(MSTJ(101)).LE.1.OR.MSTJ(101).EQ.3.OR.MSTJ(109).EQ.2.OR. 
     &  CUT.GE.0.25) THEN 
          PARJ(153)=0. 
        ELSEIF(MSTJ(110).LE.1) THEN 
          CT=LOG(1./CUT-2.) 
          PARJ(153)=ALSPI**2*CT**2*(2.419+0.5989*CT+0.6782*CT**2- 
     &    0.2661*CT**3+0.01159*CT**4)/RQCD 
 
C...Interpolation in second/first order ratio for Zhu parametrization. 
        ELSEIF(MSTJ(110).EQ.2) THEN 
          IZA=0 
          DO 110 IY=1,5 
          IF(ABS(CUT-0.01*IY).LT.0.0001) IZA=IY 
  110     CONTINUE 
          IF(IZA.NE.0) THEN 
            ZHURAT=ZHUT(IZA) 
          ELSE 
            IZ=100.*CUT 
            ZHURAT=ZHUT(IZ)+(100.*CUT-IZ)*(ZHUT(IZ+1)-ZHUT(IZ)) 
          ENDIF 
          PARJ(153)=ALSPI*PARJ(152)*ZHURAT 
        ENDIF 
 
C...Shift in second order three-jet cross-section with optimized Q^2. 
        IF(MSTJ(111).EQ.1.AND.IABS(MSTJ(101)).GE.2.AND.MSTJ(101).NE.3. 
     &  AND.CUT.LT.0.25) PARJ(153)=PARJ(153)+(33.-2.*MSTU(112))/12.* 
     &  LOG(PARJ(169))*ALSPI*PARJ(152) 
 
C...Parametrization of second order four-jet cross-section. 
        IF(IABS(MSTJ(101)).LE.1.OR.CUT.GE.0.125) THEN 
          PARJ(154)=0. 
        ELSE 
          CT=LOG(1./CUT-5.) 
          IF(CUT.LE.0.018) THEN 
            XQQGG=6.349-4.330*CT+0.8304*CT**2 
            IF(MSTJ(109).EQ.2) XQQGG=(4./3.)**2*(3.035-2.091*CT+ 
     &      0.4059*CT**2) 
            XQQQQ=1.25*(-0.1080+0.01486*CT+0.009364*CT**2) 
            IF(MSTJ(109).EQ.2) XQQQQ=8.*XQQQQ 
          ELSE 
            XQQGG=-0.09773+0.2959*CT-0.2764*CT**2+0.08832*CT**3 
            IF(MSTJ(109).EQ.2) XQQGG=(4./3.)**2*(-0.04079+0.1340*CT- 
     &      0.1326*CT**2+0.04365*CT**3) 
            XQQQQ=1.25*(0.003661-0.004888*CT-0.001081*CT**2+0.002093* 
     &      CT**3) 
            IF(MSTJ(109).EQ.2) XQQQQ=8.*XQQQQ 
          ENDIF 
          PARJ(154)=ALSPI**2*CT**2*(XQQGG+XQQQQ)/RQCD 
          PARJ(155)=XQQQQ/(XQQGG+XQQQQ) 
        ENDIF 
 
C...If negative three-jet rate, change y' optimization parameter. 
        IF(MSTJ(111).EQ.1.AND.PARJ(152)+PARJ(153).LT.0..AND. 
     &  PARJ(169).LT.0.99) THEN 
          PARJ(169)=MIN(1.,1.2*PARJ(169)) 
          Q2=PARJ(169)*ECM**2 
          ALSPI=(3./4.)*CF*ULALPS(Q2)/PARU(1) 
          GOTO 100 
        ENDIF 
 
C...If too high cross-section, use harder cuts, or fail. 
        IF(PARJ(152)+PARJ(153)+PARJ(154).GE.1) THEN 
          IF(MSTJ(110).EQ.2.AND.CUT.GT.0.0499.AND.MSTJ(111).EQ.1.AND. 
     &    PARJ(169).LT.0.99) THEN 
            PARJ(169)=MIN(1.,1.2*PARJ(169)) 
            Q2=PARJ(169)*ECM**2 
            ALSPI=(3./4.)*CF*ULALPS(Q2)/PARU(1) 
            GOTO 100 
          ELSEIF(MSTJ(110).EQ.2.AND.CUT.GT.0.0499) THEN 
            CALL LUERRM(26, 
     &      '(LUXJET:) no allowed y cut value for Zhu parametrization') 
          ENDIF 
          CUT=0.26*(4.*CUT)**(PARJ(152)+PARJ(153)+PARJ(154))**(-1./3.) 
          IF(MSTJ(110).EQ.2) CUT=MAX(0.01,MIN(0.05,CUT)) 
          GOTO 100 
        ENDIF 
 
C...Scalar gluon (first order only). 
      ELSE 
        ALSPI=ULALPS(ECM**2)/PARU(1) 
        CUT=MAX(0.001,PARJ(125),(PARJ(126)/ECM)**2,EXP(-3./ALSPI)) 
        PARJ(152)=0. 
        IF(CUT.LT.0.25) PARJ(152)=(ALSPI/3.)*((1.-2.*CUT)* 
     &  LOG((1.-2.*CUT)/CUT)+0.5*(9.*CUT**2-1.)) 
        PARJ(153)=0. 
        PARJ(154)=0. 
      ENDIF 
 
C...Select number of jets. 
      PARJ(150)=CUT 
      IF(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.5) THEN 
        NJET=2 
      ELSEIF(MSTJ(101).LE.0) THEN 
        NJET=MIN(4,2-MSTJ(101)) 
      ELSE 
        RNJ=RLU(0) 
        NJET=2 
        IF(PARJ(152)+PARJ(153)+PARJ(154).GT.RNJ) NJET=3 
        IF(PARJ(154).GT.RNJ) NJET=4 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUX3JT(NJET,CUT,KFL,ECM,X1,X2) 
 
C...Purpose: to select the kinematical variables of three-jet events. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      DIMENSION ZHUP(5,12) 
 
C...Coefficients of Zhu second order parametrization. 
      DATA ((ZHUP(IC1,IC2),IC2=1,12),IC1=1,5)/ 
     &    18.29,    89.56,    4.541,   -52.09,   -109.8,    24.90, 
     &    11.63,    3.683,    17.50, 0.002440,   -1.362,  -0.3537, 
     &    11.42,    6.299,   -22.55,   -8.915,    59.25,   -5.855, 
     &   -32.85,   -1.054,   -16.90, 0.006489,  -0.8156,  0.01095, 
     &    7.847,   -3.964,   -35.83,    1.178,    29.39,   0.2806, 
     &    47.82,   -12.36,   -56.72,  0.04054,  -0.4365,   0.6062, 
     &    5.441,   -56.89,   -50.27,    15.13,    114.3,   -18.19, 
     &    97.05,   -1.890,   -139.9,  0.08153,  -0.4984,   0.9439, 
     &   -17.65,    51.44,   -58.32,    70.95,   -255.7,   -78.99, 
     &    476.9,    29.65,   -239.3,   0.4745,   -1.174,    6.081/ 
 
C...Dilogarithm of x for x<0.5 (x>0.5 obtained by analytic trick). 
      DILOG(X)=X+X**2/4.+X**3/9.+X**4/16.+X**5/25.+X**6/36.+X**7/49. 
 
C...Event type. Mass effect factors and other common constants. 
      MSTJ(120)=2 
      MSTJ(121)=0 
      PMQ=ULMASS(KFL) 
      QME=(2.*PMQ/ECM)**2 
      IF(MSTJ(109).NE.1) THEN 
        CUTL=LOG(CUT) 
        CUTD=LOG(1./CUT-2.) 
        IF(MSTJ(109).EQ.0) THEN 
          CF=4./3. 
          CN=3. 
          TR=2. 
          WTMX=MIN(20.,37.-6.*CUTD) 
          IF(MSTJ(110).EQ.2) WTMX=2.*(7.5+80.*CUT) 
        ELSE 
          CF=1. 
          CN=0. 
          TR=12. 
          WTMX=0. 
        ENDIF 
 
C...Alpha_strong and effects of optimized Q^2 scale. Maximum weight. 
        ALS2PI=PARU(118)/PARU(2) 
        WTOPT=0. 
        IF(MSTJ(111).EQ.1) WTOPT=(33.-2.*MSTU(112))/6.*LOG(PARJ(169))* 
     &  ALS2PI 
        WTMAX=MAX(0.,1.+WTOPT+ALS2PI*WTMX) 
 
C...Choose three-jet events in allowed region. 
  100   NJET=3 
  110   Y13L=CUTL+CUTD*RLU(0) 
        Y23L=CUTL+CUTD*RLU(0) 
        Y13=EXP(Y13L) 
        Y23=EXP(Y23L) 
        Y12=1.-Y13-Y23 
        IF(Y12.LE.CUT) GOTO 110 
        IF(Y13**2+Y23**2+2.*Y12.LE.2.*RLU(0)) GOTO 110 
 
C...Second order corrections. 
        IF(MSTJ(101).EQ.2.AND.MSTJ(110).LE.1) THEN 
          Y12L=LOG(Y12) 
          Y13M=LOG(1.-Y13) 
          Y23M=LOG(1.-Y23) 
          Y12M=LOG(1.-Y12) 
          IF(Y13.LE.0.5) Y13I=DILOG(Y13) 
          IF(Y13.GE.0.5) Y13I=1.644934-Y13L*Y13M-DILOG(1.-Y13) 
          IF(Y23.LE.0.5) Y23I=DILOG(Y23) 
          IF(Y23.GE.0.5) Y23I=1.644934-Y23L*Y23M-DILOG(1.-Y23) 
          IF(Y12.LE.0.5) Y12I=DILOG(Y12) 
          IF(Y12.GE.0.5) Y12I=1.644934-Y12L*Y12M-DILOG(1.-Y12) 
          WT1=(Y13**2+Y23**2+2.*Y12)/(Y13*Y23) 
          WT2=CF*(-2.*(CUTL-Y12L)**2-3.*CUTL-1.+3.289868+ 
     &    2.*(2.*CUTL-Y12L)*CUT/Y12)+ 
     &    CN*((CUTL-Y12L)**2-(CUTL-Y13L)**2-(CUTL-Y23L)**2-11.*CUTL/6.+ 
     &    67./18.+1.644934-(2.*CUTL-Y12L)*CUT/Y12+(2.*CUTL-Y13L)* 
     &    CUT/Y13+(2.*CUTL-Y23L)*CUT/Y23)+ 
     &    TR*(2.*CUTL/3.-10./9.)+ 
     &    CF*(Y12/(Y12+Y13)+Y12/(Y12+Y23)+(Y12+Y23)/Y13+(Y12+Y13)/Y23+ 
     &    Y13L*(4.*Y12**2+2.*Y12*Y13+4.*Y12*Y23+Y13*Y23)/(Y12+Y23)**2+ 
     &    Y23L*(4.*Y12**2+2.*Y12*Y23+4.*Y12*Y13+Y13*Y23)/(Y12+Y13)**2)/ 
     &    WT1+ 
     &    CN*(Y13L*Y13/(Y12+Y23)+Y23L*Y23/(Y12+Y13))/WT1+ 
     &    (CN-2.*CF)*((Y12**2+(Y12+Y13)**2)*(Y12L*Y23L-Y12L*Y12M-Y23L* 
     &    Y23M+1.644934-Y12I-Y23I)/(Y13*Y23)+(Y12**2+(Y12+Y23)**2)* 
     &    (Y12L*Y13L-Y12L*Y12M-Y13L*Y13M+1.644934-Y12I-Y13I)/ 
     &    (Y13*Y23)+(Y13**2+Y23**2)/(Y13*Y23*(Y13+Y23))- 
     &    2.*Y12L*Y12**2/(Y13+Y23)**2-4.*Y12L*Y12/(Y13+Y23))/WT1- 
     &    CN*(Y13L*Y23L-Y13L*Y13M-Y23L*Y23M+1.644934-Y13I-Y23I) 
          IF(1.+WTOPT+ALS2PI*WT2.LE.0.) MSTJ(121)=1 
          IF(1.+WTOPT+ALS2PI*WT2.LE.WTMAX*RLU(0)) GOTO 110 
          PARJ(156)=(WTOPT+ALS2PI*WT2)/(1.+WTOPT+ALS2PI*WT2) 
 
        ELSEIF(MSTJ(101).EQ.2.AND.MSTJ(110).EQ.2) THEN 
C...Second order corrections; Zhu parametrization of ERT. 
          ZX=(Y23-Y13)**2 
          ZY=1.-Y12 
          IZA=0 
          DO 120 IY=1,5 
          IF(ABS(CUT-0.01*IY).LT.0.0001) IZA=IY 
  120     CONTINUE 
          IF(IZA.NE.0) THEN 
            IZ=IZA 
            WT2=ZHUP(IZ,1)+ZHUP(IZ,2)*ZX+ZHUP(IZ,3)*ZX**2+(ZHUP(IZ,4)+ 
     &      ZHUP(IZ,5)*ZX)*ZY+(ZHUP(IZ,6)+ZHUP(IZ,7)*ZX)*ZY**2+ 
     &      (ZHUP(IZ,8)+ZHUP(IZ,9)*ZX)*ZY**3+ZHUP(IZ,10)/(ZX-ZY**2)+ 
     &      ZHUP(IZ,11)/(1.-ZY)+ZHUP(IZ,12)/ZY 
          ELSE 
            IZ=100.*CUT 
            WTL=ZHUP(IZ,1)+ZHUP(IZ,2)*ZX+ZHUP(IZ,3)*ZX**2+(ZHUP(IZ,4)+ 
     &      ZHUP(IZ,5)*ZX)*ZY+(ZHUP(IZ,6)+ZHUP(IZ,7)*ZX)*ZY**2+ 
     &      (ZHUP(IZ,8)+ZHUP(IZ,9)*ZX)*ZY**3+ZHUP(IZ,10)/(ZX-ZY**2)+ 
     &      ZHUP(IZ,11)/(1.-ZY)+ZHUP(IZ,12)/ZY 
            IZ=IZ+1 
            WTU=ZHUP(IZ,1)+ZHUP(IZ,2)*ZX+ZHUP(IZ,3)*ZX**2+(ZHUP(IZ,4)+ 
     &      ZHUP(IZ,5)*ZX)*ZY+(ZHUP(IZ,6)+ZHUP(IZ,7)*ZX)*ZY**2+ 
     &      (ZHUP(IZ,8)+ZHUP(IZ,9)*ZX)*ZY**3+ZHUP(IZ,10)/(ZX-ZY**2)+ 
     &      ZHUP(IZ,11)/(1.-ZY)+ZHUP(IZ,12)/ZY 
            WT2=WTL+(WTU-WTL)*(100.*CUT+1.-IZ) 
          ENDIF 
          IF(1.+WTOPT+2.*ALS2PI*WT2.LE.0.) MSTJ(121)=1 
          IF(1.+WTOPT+2.*ALS2PI*WT2.LE.WTMAX*RLU(0)) GOTO 110 
          PARJ(156)=(WTOPT+2.*ALS2PI*WT2)/(1.+WTOPT+2.*ALS2PI*WT2) 
        ENDIF 
 
C...Impose mass cuts (gives two jets). For fixed jet number new try. 
        X1=1.-Y23 
        X2=1.-Y13 
        X3=1.-Y12 
        IF(4.*Y23*Y13*Y12/X3**2.LE.QME) NJET=2 
        IF(MOD(MSTJ(103),4).GE.2.AND.IABS(MSTJ(101)).LE.1.AND.QME*X3+ 
     &  0.5*QME**2+(0.5*QME+0.25*QME**2)*((1.-X2)/(1.-X1)+ 
     &  (1.-X1)/(1.-X2)).GT.(X1**2+X2**2)*RLU(0)) NJET=2 
        IF(MSTJ(101).EQ.-1.AND.NJET.EQ.2) GOTO 100 
 
C...Scalar gluon model (first order only, no mass effects). 
      ELSE 
  130   NJET=3 
  140   X3=SQRT(4.*CUT**2+RLU(0)*((1.-CUT)**2-4.*CUT**2)) 
        IF(LOG((X3-CUT)/CUT).LE.RLU(0)*LOG((1.-2.*CUT)/CUT)) GOTO 140 
        YD=SIGN(2.*CUT*((X3-CUT)/CUT)**RLU(0)-X3,RLU(0)-0.5) 
        X1=1.-0.5*(X3+YD) 
        X2=1.-0.5*(X3-YD) 
        IF(4.*(1.-X1)*(1.-X2)*(1.-X3)/X3**2.LE.QME) NJET=2 
        IF(MSTJ(102).GE.2) THEN 
          IF(X3**2-2.*(1.+X3)*(1.-X1)*(1.-X2)*PARJ(171).LT. 
     &    X3**2*RLU(0)) NJET=2 
        ENDIF 
        IF(MSTJ(101).EQ.-1.AND.NJET.EQ.2) GOTO 130 
      ENDIF 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUX4JT(NJET,CUT,KFL,ECM,KFLN,X1,X2,X4,X12,X14) 
 
C...Purpose: to select the kinematical variables of four-jet events. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      DIMENSION WTA(4),WTB(4),WTC(4),WTD(4),WTE(4) 
 
C...Common constants. Colour factors for QCD and Abelian gluon theory. 
      PMQ=ULMASS(KFL) 
      QME=(2.*PMQ/ECM)**2 
      CT=LOG(1./CUT-5.) 
      IF(MSTJ(109).EQ.0) THEN 
        CF=4./3. 
        CN=3. 
        TR=2.5 
      ELSE 
        CF=1. 
        CN=0. 
        TR=15. 
      ENDIF 
 
C...Choice of process (qqbargg or qqbarqqbar). 
  100 NJET=4 
      IT=1 
      IF(PARJ(155).GT.RLU(0)) IT=2 
      IF(MSTJ(101).LE.-3) IT=-MSTJ(101)-2 
      IF(IT.EQ.1) WTMX=0.7/CUT**2 
      IF(IT.EQ.1.AND.MSTJ(109).EQ.2) WTMX=0.6/CUT**2 
      IF(IT.EQ.2) WTMX=0.1125*CF*TR/CUT**2 
      ID=1 
 
C...Sample the five kinematical variables (for qqgg preweighted in y34). 
  110 Y134=3.*CUT+(1.-6.*CUT)*RLU(0) 
      Y234=3.*CUT+(1.-6.*CUT)*RLU(0) 
      IF(IT.EQ.1) Y34=(1.-5.*CUT)*EXP(-CT*RLU(0)) 
      IF(IT.EQ.2) Y34=CUT+(1.-6.*CUT)*RLU(0) 
      IF(Y34.LE.Y134+Y234-1..OR.Y34.GE.Y134*Y234) GOTO 110 
      VT=RLU(0) 
      CP=COS(PARU(1)*RLU(0)) 
      Y14=(Y134-Y34)*VT 
      Y13=Y134-Y14-Y34 
      VB=Y34*(1.-Y134-Y234+Y34)/((Y134-Y34)*(Y234-Y34)) 
      Y24=0.5*(Y234-Y34)*(1.-4.*SQRT(MAX(0.,VT*(1.-VT)*VB*(1.-VB)))* 
     &CP-(1.-2.*VT)*(1.-2.*VB)) 
      Y23=Y234-Y34-Y24 
      Y12=1.-Y134-Y23-Y24 
      IF(MIN(Y12,Y13,Y14,Y23,Y24).LE.CUT) GOTO 110 
      Y123=Y12+Y13+Y23 
      Y124=Y12+Y14+Y24 
 
C...Calculate matrix elements for qqgg or qqqq process. 
      IC=0 
      WTTOT=0. 
  120 IC=IC+1 
      IF(IT.EQ.1) THEN 
        WTA(IC)=(Y12*Y34**2-Y13*Y24*Y34+Y14*Y23*Y34+3.*Y12*Y23*Y34+ 
     &  3.*Y12*Y14*Y34+4.*Y12**2*Y34-Y13*Y23*Y24+2.*Y12*Y23*Y24- 
     &  Y13*Y14*Y24-2.*Y12*Y13*Y24+2.*Y12**2*Y24+Y14*Y23**2+2.*Y12* 
     &  Y23**2+Y14**2*Y23+4.*Y12*Y14*Y23+4.*Y12**2*Y23+2.*Y12*Y14**2+ 
     &  2.*Y12*Y13*Y14+4.*Y12**2*Y14+2.*Y12**2*Y13+2.*Y12**3)/(2.*Y13* 
     &  Y134*Y234*Y24)+(Y24*Y34+Y12*Y34+Y13*Y24-Y14*Y23+Y12*Y13)/(Y13* 
     &  Y134**2)+2.*Y23*(1.-Y13)/(Y13*Y134*Y24)+Y34/(2.*Y13*Y24) 
        WTB(IC)=(Y12*Y24*Y34+Y12*Y14*Y34-Y13*Y24**2+Y13*Y14*Y24+2.*Y12* 
     &  Y14*Y24)/(Y13*Y134*Y23*Y14)+Y12*(1.+Y34)*Y124/(Y134*Y234*Y14* 
     &  Y24)-(2.*Y13*Y24+Y14**2+Y13*Y23+2.*Y12*Y13)/(Y13*Y134*Y14)+ 
     &  Y12*Y123*Y124/(2.*Y13*Y14*Y23*Y24) 
        WTC(IC)=-(5.*Y12*Y34**2+2.*Y12*Y24*Y34+2.*Y12*Y23*Y34+2.*Y12* 
     &  Y14*Y34+2.*Y12*Y13*Y34+4.*Y12**2*Y34-Y13*Y24**2+Y14*Y23*Y24+ 
     &  Y13*Y23*Y24+Y13*Y14*Y24-Y12*Y14*Y24-Y13**2*Y24-3.*Y12*Y13*Y24- 
     &  Y14*Y23**2-Y14**2*Y23+Y13*Y14*Y23-3.*Y12*Y14*Y23-Y12*Y13*Y23)/ 
     &  (4.*Y134*Y234*Y34**2)+(3.*Y12*Y34**2-3.*Y13*Y24*Y34+3.*Y12*Y24* 
     &  Y34+3.*Y14*Y23*Y34-Y13*Y24**2-Y12*Y23*Y34+6.*Y12*Y14*Y34+2.*Y12* 
     &  Y13*Y34-2.*Y12**2*Y34+Y14*Y23*Y24-3.*Y13*Y23*Y24-2.*Y13*Y14* 
     &  Y24+4.*Y12*Y14*Y24+2.*Y12*Y13*Y24+3.*Y14*Y23**2+2.*Y14**2*Y23+ 
     &  2.*Y14**2*Y12+2.*Y12**2*Y14+6.*Y12*Y14*Y23-2.*Y12*Y13**2- 
     &  2.*Y12**2*Y13)/(4.*Y13*Y134*Y234*Y34) 
        WTC(IC)=WTC(IC)+(2.*Y12*Y34**2-2.*Y13*Y24*Y34+Y12*Y24*Y34+ 
     &  4.*Y13*Y23*Y34+4.*Y12*Y14*Y34+2.*Y12*Y13*Y34+2.*Y12**2*Y34- 
     &  Y13*Y24**2+3.*Y14*Y23*Y24+4.*Y13*Y23*Y24-2.*Y13*Y14*Y24+ 
     &  4.*Y12*Y14*Y24+2.*Y12*Y13*Y24+2.*Y14*Y23**2+4.*Y13*Y23**2+ 
     &  2.*Y13*Y14*Y23+2.*Y12*Y14*Y23+4.*Y12*Y13*Y23+2.*Y12*Y14**2+4.* 
     &  Y12**2*Y13+4.*Y12*Y13*Y14+2.*Y12**2*Y14)/(4.*Y13*Y134*Y24*Y34)- 
     &  (Y12*Y34**2-2.*Y14*Y24*Y34-2.*Y13*Y24*Y34-Y14*Y23*Y34+Y13*Y23* 
     &  Y34+Y12*Y14*Y34+2.*Y12*Y13*Y34-2.*Y14**2*Y24-4.*Y13*Y14*Y24- 
     &  4.*Y13**2*Y24-Y14**2*Y23-Y13**2*Y23+Y12*Y13*Y14-Y12*Y13**2)/ 
     &  (2.*Y13*Y34*Y134**2)+(Y12*Y34**2-4.*Y14*Y24*Y34-2.*Y13*Y24*Y34- 
     &  2.*Y14*Y23*Y34-4.*Y13*Y23*Y34-4.*Y12*Y14*Y34-4.*Y12*Y13*Y34- 
     &  2.*Y13*Y14*Y24+2.*Y13**2*Y24+2.*Y14**2*Y23-2.*Y13*Y14*Y23- 
     &  Y12*Y14**2-6.*Y12*Y13*Y14-Y12*Y13**2)/(4.*Y34**2*Y134**2) 
        WTTOT=WTTOT+Y34*CF*(CF*WTA(IC)+(CF-0.5*CN)*WTB(IC)+CN*WTC(IC))/ 
     &  8. 
      ELSE 
        WTD(IC)=(Y13*Y23*Y34+Y12*Y23*Y34-Y12**2*Y34+Y13*Y23*Y24+2.*Y12* 
     &  Y23*Y24-Y14*Y23**2+Y12*Y13*Y24+Y12*Y14*Y23+Y12*Y13*Y14)/(Y13**2* 
     &  Y123**2)-(Y12*Y34**2-Y13*Y24*Y34+Y12*Y24*Y34-Y14*Y23*Y34-Y12* 
     &  Y23*Y34-Y13*Y24**2+Y14*Y23*Y24-Y13*Y23*Y24-Y13**2*Y24+Y14* 
     &  Y23**2)/(Y13**2*Y123*Y134)+(Y13*Y14*Y12+Y34*Y14*Y12-Y34**2*Y12+ 
     &  Y13*Y14*Y24+2.*Y34*Y14*Y24-Y23*Y14**2+Y34*Y13*Y24+Y34*Y23*Y14+ 
     &  Y34*Y13*Y23)/(Y13**2*Y134**2)-(Y34*Y12**2-Y13*Y24*Y12+Y34*Y24* 
     &  Y12-Y23*Y14*Y12-Y34*Y14*Y12-Y13*Y24**2+Y23*Y14*Y24-Y13*Y14*Y24- 
     &  Y13**2*Y24+Y23*Y14**2)/(Y13**2*Y134*Y123) 
        WTE(IC)=(Y12*Y34*(Y23-Y24+Y14+Y13)+Y13*Y24**2-Y14*Y23*Y24+Y13* 
     &  Y23*Y24+Y13*Y14*Y24+Y13**2*Y24-Y14*Y23*(Y14+Y23+Y13))/(Y13*Y23* 
     &  Y123*Y134)-Y12*(Y12*Y34-Y23*Y24-Y13*Y24-Y14*Y23-Y14*Y13)/(Y13* 
     &  Y23*Y123**2)-(Y14+Y13)*(Y24+Y23)*Y34/(Y13*Y23*Y134*Y234)+ 
     &  (Y12*Y34*(Y14-Y24+Y23+Y13)+Y13*Y24**2-Y23*Y14*Y24+Y13*Y14*Y24+ 
     &  Y13*Y23*Y24+Y13**2*Y24-Y23*Y14*(Y14+Y23+Y13))/(Y13*Y14*Y134* 
     &  Y123)-Y34*(Y34*Y12-Y14*Y24-Y13*Y24-Y23*Y14-Y23*Y13)/(Y13*Y14* 
     &  Y134**2)-(Y23+Y13)*(Y24+Y14)*Y12/(Y13*Y14*Y123*Y124) 
        WTTOT=WTTOT+CF*(TR*WTD(IC)+(CF-0.5*CN)*WTE(IC))/16. 
      ENDIF 
 
C...Permutations of momenta in matrix element. Weighting. 
  130 IF(IC.EQ.1.OR.IC.EQ.3.OR.ID.EQ.2.OR.ID.EQ.3) THEN 
        YSAV=Y13 
        Y13=Y14 
        Y14=YSAV 
        YSAV=Y23 
        Y23=Y24 
        Y24=YSAV 
        YSAV=Y123 
        Y123=Y124 
        Y124=YSAV 
      ENDIF 
      IF(IC.EQ.2.OR.IC.EQ.4.OR.ID.EQ.3.OR.ID.EQ.4) THEN 
        YSAV=Y13 
        Y13=Y23 
        Y23=YSAV 
        YSAV=Y14 
        Y14=Y24 
        Y24=YSAV 
        YSAV=Y134 
        Y134=Y234 
        Y234=YSAV 
      ENDIF 
      IF(IC.LE.3) GOTO 120 
      IF(ID.EQ.1.AND.WTTOT.LT.RLU(0)*WTMX) GOTO 110 
      IC=5 
 
C...qqgg events: string configuration and event type. 
      IF(IT.EQ.1) THEN 
        IF(MSTJ(109).EQ.0.AND.ID.EQ.1) THEN 
          PARJ(156)=Y34*(2.*(WTA(1)+WTA(2)+WTA(3)+WTA(4))+4.*(WTC(1)+ 
     &    WTC(2)+WTC(3)+WTC(4)))/(9.*WTTOT) 
          IF(WTA(2)+WTA(4)+2.*(WTC(2)+WTC(4)).GT.RLU(0)*(WTA(1)+WTA(2)+ 
     &    WTA(3)+WTA(4)+2.*(WTC(1)+WTC(2)+WTC(3)+WTC(4)))) ID=2 
          IF(ID.EQ.2) GOTO 130 
        ELSEIF(MSTJ(109).EQ.2.AND.ID.EQ.1) THEN 
          PARJ(156)=Y34*(WTA(1)+WTA(2)+WTA(3)+WTA(4))/(8.*WTTOT) 
          IF(WTA(2)+WTA(4).GT.RLU(0)*(WTA(1)+WTA(2)+WTA(3)+WTA(4))) ID=2 
          IF(ID.EQ.2) GOTO 130 
        ENDIF 
        MSTJ(120)=3 
        IF(MSTJ(109).EQ.0.AND.0.5*Y34*(WTC(1)+WTC(2)+WTC(3)+WTC(4)).GT. 
     &  RLU(0)*WTTOT) MSTJ(120)=4 
        KFLN=21 
 
C...Mass cuts. Kinematical variables out. 
        IF(Y12.LE.CUT+QME) NJET=2 
        IF(NJET.EQ.2) GOTO 150 
        Q12=0.5*(1.-SQRT(1.-QME/Y12)) 
        X1=1.-(1.-Q12)*Y234-Q12*Y134 
        X4=1.-(1.-Q12)*Y134-Q12*Y234 
        X2=1.-Y124 
        X12=(1.-Q12)*Y13+Q12*Y23 
        X14=Y12-0.5*QME 
        IF(Y134*Y234/((1.-X1)*(1.-X4)).LE.RLU(0)) NJET=2 
 
C...qqbarqqbar events: string configuration, choose new flavour. 
      ELSE 
        IF(ID.EQ.1) THEN 
          WTR=RLU(0)*(WTD(1)+WTD(2)+WTD(3)+WTD(4)) 
          IF(WTR.LT.WTD(2)+WTD(3)+WTD(4)) ID=2 
          IF(WTR.LT.WTD(3)+WTD(4)) ID=3 
          IF(WTR.LT.WTD(4)) ID=4 
          IF(ID.GE.2) GOTO 130 
        ENDIF 
        MSTJ(120)=5 
        PARJ(156)=CF*TR*(WTD(1)+WTD(2)+WTD(3)+WTD(4))/(16.*WTTOT) 
  140   KFLN=1+INT(5.*RLU(0)) 
        IF(KFLN.NE.KFL.AND.0.2*PARJ(156).LE.RLU(0)) GOTO 140 
        IF(KFLN.EQ.KFL.AND.1.-0.8*PARJ(156).LE.RLU(0)) GOTO 140 
        IF(KFLN.GT.MSTJ(104)) NJET=2 
        PMQN=ULMASS(KFLN) 
        QMEN=(2.*PMQN/ECM)**2 
 
C...Mass cuts. Kinematical variables out. 
        IF(Y24.LE.CUT+QME.OR.Y13.LE.1.1*QMEN) NJET=2 
        IF(NJET.EQ.2) GOTO 150 
        Q24=0.5*(1.-SQRT(1.-QME/Y24)) 
        Q13=0.5*(1.-SQRT(1.-QMEN/Y13)) 
        X1=1.-(1.-Q24)*Y123-Q24*Y134 
        X4=1.-(1.-Q24)*Y134-Q24*Y123 
        X2=1.-(1.-Q13)*Y234-Q13*Y124 
        X12=(1.-Q24)*((1.-Q13)*Y14+Q13*Y34)+Q24*((1.-Q13)*Y12+Q13*Y23) 
        X14=Y24-0.5*QME 
        X34=(1.-Q24)*((1.-Q13)*Y23+Q13*Y12)+Q24*((1.-Q13)*Y34+Q13*Y14) 
        IF(PMQ**2+PMQN**2+MIN(X12,X34)*ECM**2.LE. 
     &  (PARJ(127)+PMQ+PMQN)**2) NJET=2 
        IF(Y123*Y134/((1.-X1)*(1.-X4)).LE.RLU(0)) NJET=2 
      ENDIF 
  150 IF(MSTJ(101).LE.-2.AND.NJET.EQ.2) GOTO 100 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUXDIF(NC,NJET,KFL,ECM,CHI,THE,PHI) 
 
C...Purpose: to give the angular orientation of events. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Charge. Factors depending on polarization for QED case. 
      QF=KCHG(KFL,1)/3. 
      POLL=1.-PARJ(131)*PARJ(132) 
      POLD=PARJ(132)-PARJ(131) 
      IF(MSTJ(102).LE.1.OR.MSTJ(109).EQ.1) THEN 
        HF1=POLL 
        HF2=0. 
        HF3=PARJ(133)**2 
        HF4=0. 
 
C...Factors depending on flavour, energy and polarization for QFD case. 
      ELSE 
        SFF=1./(16.*PARU(102)*(1.-PARU(102))) 
        SFW=ECM**4/((ECM**2-PARJ(123)**2)**2+(PARJ(123)*PARJ(124))**2) 
        SFI=SFW*(1.-(PARJ(123)/ECM)**2) 
        AE=-1. 
        VE=4.*PARU(102)-1. 
        AF=SIGN(1.,QF) 
        VF=AF-4.*QF*PARU(102) 
        HF1=QF**2*POLL-2.*QF*VF*SFI*SFF*(VE*POLL-AE*POLD)+ 
     &  (VF**2+AF**2)*SFW*SFF**2*((VE**2+AE**2)*POLL-2.*VE*AE*POLD) 
        HF2=-2.*QF*AF*SFI*SFF*(AE*POLL-VE*POLD)+2.*VF*AF*SFW*SFF**2* 
     &  (2.*VE*AE*POLL-(VE**2+AE**2)*POLD) 
        HF3=PARJ(133)**2*(QF**2-2.*QF*VF*SFI*SFF*VE+(VF**2+AF**2)* 
     &  SFW*SFF**2*(VE**2-AE**2)) 
        HF4=-PARJ(133)**2*2.*QF*VF*SFW*(PARJ(123)*PARJ(124)/ECM**2)* 
     &  SFF*AE 
      ENDIF 
 
C...Mass factor. Differential cross-sections for two-jet events. 
      SQ2=SQRT(2.) 
      QME=0. 
      IF(MSTJ(103).GE.4.AND.IABS(MSTJ(101)).LE.1.AND.MSTJ(102).LE.1.AND. 
     &MSTJ(109).NE.1) QME=(2.*ULMASS(KFL)/ECM)**2 
      IF(NJET.EQ.2) THEN 
        SIGU=4.*SQRT(1.-QME) 
        SIGL=2.*QME*SQRT(1.-QME) 
        SIGT=0. 
        SIGI=0. 
        SIGA=0. 
        SIGP=4. 
 
C...Kinematical variables. Reduce four-jet event to three-jet one. 
      ELSE 
        IF(NJET.EQ.3) THEN 
          X1=2.*P(NC+1,4)/ECM 
          X2=2.*P(NC+3,4)/ECM 
        ELSE 
          ECMR=P(NC+1,4)+P(NC+4,4)+SQRT((P(NC+2,1)+P(NC+3,1))**2+ 
     &    (P(NC+2,2)+P(NC+3,2))**2+(P(NC+2,3)+P(NC+3,3))**2) 
          X1=2.*P(NC+1,4)/ECMR 
          X2=2.*P(NC+4,4)/ECMR 
        ENDIF 
 
C...Differential cross-sections for three-jet (or reduced four-jet). 
        XQ=(1.-X1)/(1.-X2) 
        CT12=(X1*X2-2.*X1-2.*X2+2.+QME)/SQRT((X1**2-QME)*(X2**2-QME)) 
        ST12=SQRT(1.-CT12**2) 
        IF(MSTJ(109).NE.1) THEN 
          SIGU=2.*X1**2+X2**2*(1.+CT12**2)-QME*(3.+CT12**2-X1-X2)- 
     &    QME*X1/XQ+0.5*QME*((X2**2-QME)*ST12**2-2.*X2)*XQ 
          SIGL=(X2*ST12)**2-QME*(3.-CT12**2-2.5*(X1+X2)+X1*X2+QME)+ 
     &    0.5*QME*(X1**2-X1-QME)/XQ+0.5*QME*((X2**2-QME)*CT12**2-X2)*XQ 
          SIGT=0.5*(X2**2-QME-0.5*QME*(X2**2-QME)/XQ)*ST12**2 
          SIGI=((1.-0.5*QME*XQ)*(X2**2-QME)*ST12*CT12+QME*(1.-X1-X2+ 
     &    0.5*X1*X2+0.5*QME)*ST12/CT12)/SQ2 
          SIGA=X2**2*ST12/SQ2 
          SIGP=2.*(X1**2-X2**2*CT12) 
 
C...Differential cross-sect for scalar gluons (no mass effects). 
        ELSE 
          X3=2.-X1-X2 
          XT=X2*ST12 
          CT13=SQRT(MAX(0.,1.-(XT/X3)**2)) 
          SIGU=(1.-PARJ(171))*(X3**2-0.5*XT**2)+ 
     &    PARJ(171)*(X3**2-0.5*XT**2-4.*(1.-X1)*(1.-X2)**2/X1) 
          SIGL=(1.-PARJ(171))*0.5*XT**2+ 
     &    PARJ(171)*0.5*(1.-X1)**2*XT**2 
          SIGT=(1.-PARJ(171))*0.25*XT**2+ 
     &    PARJ(171)*0.25*XT**2*(1.-2.*X1) 
          SIGI=-(0.5/SQ2)*((1.-PARJ(171))*XT*X3*CT13+ 
     &    PARJ(171)*XT*((1.-2.*X1)*X3*CT13-X1*(X1-X2))) 
          SIGA=(0.25/SQ2)*XT*(2.*(1.-X1)-X1*X3) 
          SIGP=X3**2-2.*(1.-X1)*(1.-X2)/X1 
        ENDIF 
      ENDIF 
 
C...Upper bounds for differential cross-section. 
      HF1A=ABS(HF1) 
      HF2A=ABS(HF2) 
      HF3A=ABS(HF3) 
      HF4A=ABS(HF4) 
      SIGMAX=(2.*HF1A+HF3A+HF4A)*ABS(SIGU)+2.*(HF1A+HF3A+HF4A)* 
     &ABS(SIGL)+2.*(HF1A+2.*HF3A+2.*HF4A)*ABS(SIGT)+2.*SQ2* 
     &(HF1A+2.*HF3A+2.*HF4A)*ABS(SIGI)+4.*SQ2*HF2A*ABS(SIGA)+ 
     &2.*HF2A*ABS(SIGP) 
 
C...Generate angular orientation according to differential cross-sect. 
  100 CHI=PARU(2)*RLU(0) 
      CTHE=2.*RLU(0)-1. 
      PHI=PARU(2)*RLU(0) 
      CCHI=COS(CHI) 
      SCHI=SIN(CHI) 
      C2CHI=COS(2.*CHI) 
      S2CHI=SIN(2.*CHI) 
      THE=ACOS(CTHE) 
      STHE=SIN(THE) 
      C2PHI=COS(2.*(PHI-PARJ(134))) 
      S2PHI=SIN(2.*(PHI-PARJ(134))) 
      SIG=((1.+CTHE**2)*HF1+STHE**2*(C2PHI*HF3-S2PHI*HF4))*SIGU+ 
     &2.*(STHE**2*HF1-STHE**2*(C2PHI*HF3-S2PHI*HF4))*SIGL+ 
     &2.*(STHE**2*C2CHI*HF1+((1.+CTHE**2)*C2CHI*C2PHI-2.*CTHE*S2CHI* 
     &S2PHI)*HF3-((1.+CTHE**2)*C2CHI*S2PHI+2.*CTHE*S2CHI*C2PHI)*HF4)* 
     &SIGT-2.*SQ2*(2.*STHE*CTHE*CCHI*HF1-2.*STHE*(CTHE*CCHI*C2PHI- 
     &SCHI*S2PHI)*HF3+2.*STHE*(CTHE*CCHI*S2PHI+SCHI*C2PHI)*HF4)*SIGI+ 
     &4.*SQ2*STHE*CCHI*HF2*SIGA+2.*CTHE*HF2*SIGP 
      IF(SIG.LT.SIGMAX*RLU(0)) GOTO 100 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUONIA(KFL,ECM) 
 
C...Purpose: to generate Upsilon and toponium decays into three 
C...gluons or two gluons and a photon. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Printout. Check input parameters. 
      IF(MSTU(12).GE.1) CALL LULIST(0) 
      IF(KFL.LT.0.OR.KFL.GT.8) THEN 
        CALL LUERRM(16,'(LUONIA:) called with unknown flavour code') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
      IF(ECM.LT.PARJ(127)+2.02*PARF(101)) THEN 
        CALL LUERRM(16,'(LUONIA:) called with too small CM energy') 
        IF(MSTU(21).GE.1) RETURN 
      ENDIF 
 
C...Initial e+e- and onium state (optional). 
      NC=0 
      IF(MSTJ(115).GE.2) THEN 
        NC=NC+2 
        CALL LU1ENT(NC-1,11,0.5*ECM,0.,0.) 
        K(NC-1,1)=21 
        CALL LU1ENT(NC,-11,0.5*ECM,PARU(1),0.) 
        K(NC,1)=21 
      ENDIF 
      KFLC=IABS(KFL) 
      IF(MSTJ(115).GE.3.AND.KFLC.GE.5) THEN 
        NC=NC+1 
        KF=110*KFLC+3 
        MSTU10=MSTU(10) 
        MSTU(10)=1 
        P(NC,5)=ECM 
        CALL LU1ENT(NC,KF,ECM,0.,0.) 
        K(NC,1)=21 
        K(NC,3)=1 
        MSTU(10)=MSTU10 
      ENDIF 
 
C...Choose x1 and x2 according to matrix element. 
      NTRY=0 
  100 X1=RLU(0) 
      X2=RLU(0) 
      X3=2.-X1-X2 
      IF(X3.GE.1..OR.((1.-X1)/(X2*X3))**2+((1.-X2)/(X1*X3))**2+ 
     &((1.-X3)/(X1*X2))**2.LE.2.*RLU(0)) GOTO 100 
      NTRY=NTRY+1 
      NJET=3 
      IF(MSTJ(101).LE.4) CALL LU3ENT(NC+1,21,21,21,ECM,X1,X3) 
      IF(MSTJ(101).GE.5) CALL LU3ENT(-(NC+1),21,21,21,ECM,X1,X3) 
 
C...Photon-gluon-gluon events. Small system modifications. Jet origin. 
      MSTU(111)=MSTJ(108) 
      IF(MSTJ(108).EQ.2.AND.(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.1)) 
     &MSTU(111)=1 
      PARU(112)=PARJ(121) 
      IF(MSTU(111).EQ.2) PARU(112)=PARJ(122) 
      QF=0. 
      IF(KFLC.NE.0) QF=KCHG(KFLC,1)/3. 
      RGAM=7.2*QF**2*PARU(101)/ULALPS(ECM**2) 
      MK=0 
      ECMC=ECM 
      IF(RLU(0).GT.RGAM/(1.+RGAM)) THEN 
        IF(1.-MAX(X1,X2,X3).LE.MAX((PARJ(126)/ECM)**2,PARJ(125))) 
     &  NJET=2 
        IF(NJET.EQ.2.AND.MSTJ(101).LE.4) CALL LU2ENT(NC+1,21,21,ECM) 
        IF(NJET.EQ.2.AND.MSTJ(101).GE.5) CALL LU2ENT(-(NC+1),21,21,ECM) 
      ELSE 
        MK=1 
        ECMC=SQRT(1.-X1)*ECM 
        IF(ECMC.LT.2.*PARJ(127)) GOTO 100 
        K(NC+1,1)=1 
        K(NC+1,2)=22 
        K(NC+1,4)=0 
        K(NC+1,5)=0 
        IF(MSTJ(101).GE.5) K(NC+2,4)=MSTU(5)*(NC+3) 
        IF(MSTJ(101).GE.5) K(NC+2,5)=MSTU(5)*(NC+3) 
        IF(MSTJ(101).GE.5) K(NC+3,4)=MSTU(5)*(NC+2) 
        IF(MSTJ(101).GE.5) K(NC+3,5)=MSTU(5)*(NC+2) 
        NJET=2 
        IF(ECMC.LT.4.*PARJ(127)) THEN 
          MSTU10=MSTU(10) 
          MSTU(10)=1 
          P(NC+2,5)=ECMC 
          CALL LU1ENT(NC+2,83,0.5*(X2+X3)*ECM,PARU(1),0.) 
          MSTU(10)=MSTU10 
          NJET=0 
        ENDIF 
      ENDIF 
      DO 110 IP=NC+1,N 
      K(IP,3)=K(IP,3)+(MSTJ(115)/2)+(KFLC/5)*(MSTJ(115)/3)*(NC-1) 
  110 CONTINUE 
 
C...Differential cross-sections. Upper limit for cross-section. 
      IF(MSTJ(106).EQ.1) THEN 
        SQ2=SQRT(2.) 
        HF1=1.-PARJ(131)*PARJ(132) 
        HF3=PARJ(133)**2 
        CT13=(X1*X3-2.*X1-2.*X3+2.)/(X1*X3) 
        ST13=SQRT(1.-CT13**2) 
        SIGL=0.5*X3**2*((1.-X2)**2+(1.-X3)**2)*ST13**2 
        SIGU=(X1*(1.-X1))**2+(X2*(1.-X2))**2+(X3*(1.-X3))**2-SIGL 
        SIGT=0.5*SIGL 
        SIGI=(SIGL*CT13/ST13+0.5*X1*X3*(1.-X2)**2*ST13)/SQ2 
        SIGMAX=(2.*HF1+HF3)*ABS(SIGU)+2.*(HF1+HF3)*ABS(SIGL)+2.*(HF1+ 
     &  2.*HF3)*ABS(SIGT)+2.*SQ2*(HF1+2.*HF3)*ABS(SIGI) 
 
C...Angular orientation of event. 
  120   CHI=PARU(2)*RLU(0) 
        CTHE=2.*RLU(0)-1. 
        PHI=PARU(2)*RLU(0) 
        CCHI=COS(CHI) 
        SCHI=SIN(CHI) 
        C2CHI=COS(2.*CHI) 
        S2CHI=SIN(2.*CHI) 
        THE=ACOS(CTHE) 
        STHE=SIN(THE) 
        C2PHI=COS(2.*(PHI-PARJ(134))) 
        S2PHI=SIN(2.*(PHI-PARJ(134))) 
        SIG=((1.+CTHE**2)*HF1+STHE**2*C2PHI*HF3)*SIGU+2.*(STHE**2*HF1- 
     &  STHE**2*C2PHI*HF3)*SIGL+2.*(STHE**2*C2CHI*HF1+((1.+CTHE**2)* 
     &  C2CHI*C2PHI-2.*CTHE*S2CHI*S2PHI)*HF3)*SIGT-2.*SQ2*(2.*STHE*CTHE* 
     &  CCHI*HF1-2.*STHE*(CTHE*CCHI*C2PHI-SCHI*S2PHI)*HF3)*SIGI 
        IF(SIG.LT.SIGMAX*RLU(0)) GOTO 120 
        CALL LUDBRB(NC+1,N,0.,CHI,0D0,0D0,0D0) 
        CALL LUDBRB(NC+1,N,THE,PHI,0D0,0D0,0D0) 
      ENDIF 
 
C...Generate parton shower. Rearrange along strings and check. 
      IF(MSTJ(101).GE.5.AND.NJET.GE.2) THEN 
        CALL LUSHOW(NC+MK+1,-NJET,ECMC) 
        MSTJ14=MSTJ(14) 
        IF(MSTJ(105).EQ.-1) MSTJ(14)=-1 
        IF(MSTJ(105).GE.0) MSTU(28)=0 
        CALL LUPREP(0) 
        MSTJ(14)=MSTJ14 
        IF(MSTJ(105).GE.0.AND.MSTU(28).NE.0) GOTO 100 
      ENDIF 
 
C...Generate fragmentation. Information for LUTABU: 
      IF(MSTJ(105).EQ.1) CALL LUEXEC 
      MSTU(161)=110*KFLC+3 
      MSTU(162)=0 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUHEPC(MCONV) 
 
C...Purpose: to convert JETSET event record contents to or from 
C...the standard event record commonblock. 
C...Note that HEPEVT is in double precision according to LEP 2 standard.
      PARAMETER (NMXHEP=2000) 
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP), 
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP) 
      DOUBLE PRECISION PHEP,VHEP
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /HEPEVT/ 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Conversion from JETSET to standard, the easy part. 
      IF(MCONV.EQ.1) THEN 
        NEVHEP=0 
        IF(N.GT.NMXHEP) CALL LUERRM(8, 
     &  '(LUHEPC:) no more space in /HEPEVT/') 
        NHEP=MIN(N,NMXHEP) 
        DO 140 I=1,NHEP 
        ISTHEP(I)=0 
        IF(K(I,1).GE.1.AND.K(I,1).LE.10) ISTHEP(I)=1 
        IF(K(I,1).GE.11.AND.K(I,1).LE.20) ISTHEP(I)=2 
        IF(K(I,1).GE.21.AND.K(I,1).LE.30) ISTHEP(I)=3 
        IF(K(I,1).GE.31.AND.K(I,1).LE.100) ISTHEP(I)=K(I,1) 
        IDHEP(I)=K(I,2) 
        JMOHEP(1,I)=K(I,3) 
        JMOHEP(2,I)=0 
        IF(K(I,1).NE.3.AND.K(I,1).NE.13.AND.K(I,1).NE.14) THEN 
          JDAHEP(1,I)=K(I,4) 
          JDAHEP(2,I)=K(I,5) 
        ELSE 
          JDAHEP(1,I)=0 
          JDAHEP(2,I)=0 
        ENDIF 
        DO 100 J=1,5 
        PHEP(J,I)=P(I,J) 
  100   CONTINUE 
        DO 110 J=1,4 
        VHEP(J,I)=V(I,J) 
  110   CONTINUE 
 
C...Check if new event (from pileup). 
        IF(I.EQ.1) THEN 
          INEW=1 
        ELSE 
          IF(K(I,1).EQ.21.AND.K(I-1,1).NE.21) INEW=I 
        ENDIF 
 
C...Fill in missing mother information. 
        IF(I.GE.INEW+2.AND.K(I,1).EQ.21.AND.K(I,3).EQ.0) THEN 
          IMO1=I-2 
          IF(I.GE.INEW+3.AND.K(I-1,1).EQ.21.AND.K(I-1,3).EQ.0) 
     &    IMO1=IMO1-1 
          JMOHEP(1,I)=IMO1 
          JMOHEP(2,I)=IMO1+1 
        ELSEIF(K(I,2).GE.91.AND.K(I,2).LE.93) THEN 
          I1=K(I,3)-1 
  120     I1=I1+1 
          IF(I1.GE.I) CALL LUERRM(8, 
     &    '(LUHEPC:) translation of inconsistent event history') 
          IF(I1.LT.I.AND.K(I1,1).NE.1.AND.K(I1,1).NE.11) GOTO 120 
          KC=LUCOMP(K(I1,2)) 
          IF(I1.LT.I.AND.KC.EQ.0) GOTO 120 
          IF(I1.LT.I.AND.KCHG(KC,2).EQ.0) GOTO 120 
          JMOHEP(2,I)=I1 
        ELSEIF(K(I,2).EQ.94) THEN 
          NJET=2 
          IF(NHEP.GE.I+3.AND.K(I+3,3).LE.I) NJET=3 
          IF(NHEP.GE.I+4.AND.K(I+4,3).LE.I) NJET=4 
          JMOHEP(2,I)=MOD(K(I+NJET,4)/MSTU(5),MSTU(5)) 
          IF(JMOHEP(2,I).EQ.JMOHEP(1,I)) JMOHEP(2,I)= 
     &    MOD(K(I+1,4)/MSTU(5),MSTU(5)) 
        ENDIF 
 
C...Fill in missing daughter information. 
        IF(K(I,2).EQ.94.AND.MSTU(16).NE.2) THEN 
          DO 130 I1=JDAHEP(1,I),JDAHEP(2,I) 
          I2=MOD(K(I1,4)/MSTU(5),MSTU(5)) 
          JDAHEP(1,I2)=I 
  130     CONTINUE 
        ENDIF 
        IF(K(I,2).GE.91.AND.K(I,2).LE.94) GOTO 140 
        I1=JMOHEP(1,I) 
        IF(I1.LE.0.OR.I1.GT.NHEP) GOTO 140 
        IF(K(I1,1).NE.13.AND.K(I1,1).NE.14) GOTO 140 
        IF(JDAHEP(1,I1).EQ.0) THEN 
          JDAHEP(1,I1)=I 
        ELSE 
          JDAHEP(2,I1)=I 
        ENDIF 
  140   CONTINUE 
        DO 150 I=1,NHEP 
        IF(K(I,1).NE.13.AND.K(I,1).NE.14) GOTO 150 
        IF(JDAHEP(2,I).EQ.0) JDAHEP(2,I)=JDAHEP(1,I) 
  150   CONTINUE 
 
C...Conversion from standard to JETSET, the easy part. 
      ELSE 
        IF(NHEP.GT.MSTU(4)) CALL LUERRM(8, 
     &  '(LUHEPC:) no more space in /LUJETS/') 
        N=MIN(NHEP,MSTU(4)) 
        NKQ=0 
        KQSUM=0 
        DO 180 I=1,N 
        K(I,1)=0 
        IF(ISTHEP(I).EQ.1) K(I,1)=1 
        IF(ISTHEP(I).EQ.2) K(I,1)=11 
        IF(ISTHEP(I).EQ.3) K(I,1)=21 
        K(I,2)=IDHEP(I) 
        K(I,3)=JMOHEP(1,I) 
        K(I,4)=JDAHEP(1,I) 
        K(I,5)=JDAHEP(2,I) 
        DO 160 J=1,5 
        P(I,J)=PHEP(J,I) 
  160   CONTINUE 
        DO 170 J=1,4 
        V(I,J)=VHEP(J,I) 
  170   CONTINUE 
        V(I,5)=0. 
        IF(ISTHEP(I).EQ.2.AND.PHEP(4,I).GT.PHEP(5,I)) THEN 
          I1=JDAHEP(1,I) 
          IF(I1.GT.0.AND.I1.LE.NHEP) V(I,5)=(VHEP(4,I1)-VHEP(4,I))* 
     &    PHEP(5,I)/PHEP(4,I) 
        ENDIF 
 
C...Fill in missing information on colour connection in jet systems. 
        IF(ISTHEP(I).EQ.1) THEN 
          KC=LUCOMP(K(I,2)) 
          KQ=0 
          IF(KC.NE.0) KQ=KCHG(KC,2)*ISIGN(1,K(I,2)) 
          IF(KQ.NE.0) NKQ=NKQ+1 
          IF(KQ.NE.2) KQSUM=KQSUM+KQ 
          IF(KQ.NE.0.AND.KQSUM.NE.0) THEN 
            K(I,1)=2 
          ELSEIF(KQ.EQ.2.AND.I.LT.N) THEN 
            IF(K(I+1,2).EQ.21) K(I,1)=2 
          ENDIF 
        ENDIF 
  180   CONTINUE 
        IF(NKQ.EQ.1.OR.KQSUM.NE.0) CALL LUERRM(8, 
     &  '(LUHEPC:) input parton configuration not colour singlet') 
      ENDIF 
 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUTEST(MTEST) 
 
C...Purpose: to provide a simple program (disguised as subroutine) to 
C...run at installation as a check that the program works as intended. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUJETS/,/LUDAT1/ 
      DIMENSION PSUM(5),PINI(6),PFIN(6) 
 
C...Loop over events to be generated. 
      IF(MTEST.GE.1) CALL LUTABU(20) 
      NERR=0 
      DO 180 IEV=1,600 
 
C...Reset parameter values. Switch on some nonstandard features. 
      MSTJ(1)=1 
      MSTJ(3)=0 
      MSTJ(11)=1 
      MSTJ(42)=2 
      MSTJ(43)=4 
      MSTJ(44)=2 
      PARJ(17)=0.1 
      PARJ(22)=1.5 
      PARJ(43)=1. 
      PARJ(54)=-0.05 
      MSTJ(101)=5 
      MSTJ(104)=5 
      MSTJ(105)=0 
      MSTJ(107)=1 
      IF(IEV.EQ.301.OR.IEV.EQ.351.OR.IEV.EQ.401) MSTJ(116)=3 
 
C...Ten events each for some single jets configurations. 
      IF(IEV.LE.50) THEN 
        ITY=(IEV+9)/10 
        MSTJ(3)=-1 
        IF(ITY.EQ.3.OR.ITY.EQ.4) MSTJ(11)=2 
        IF(ITY.EQ.1) CALL LU1ENT(1,1,15.,0.,0.) 
        IF(ITY.EQ.2) CALL LU1ENT(1,3101,15.,0.,0.) 
        IF(ITY.EQ.3) CALL LU1ENT(1,-2203,15.,0.,0.) 
        IF(ITY.EQ.4) CALL LU1ENT(1,-4,30.,0.,0.) 
        IF(ITY.EQ.5) CALL LU1ENT(1,21,15.,0.,0.) 
 
C...Ten events each for some simple jet systems; string fragmentation. 
      ELSEIF(IEV.LE.130) THEN 
        ITY=(IEV-41)/10 
        IF(ITY.EQ.1) CALL LU2ENT(1,1,-1,40.) 
        IF(ITY.EQ.2) CALL LU2ENT(1,4,-4,30.) 
        IF(ITY.EQ.3) CALL LU2ENT(1,2,2103,100.) 
        IF(ITY.EQ.4) CALL LU2ENT(1,21,21,40.) 
        IF(ITY.EQ.5) CALL LU3ENT(1,2101,21,-3203,30.,0.6,0.8) 
        IF(ITY.EQ.6) CALL LU3ENT(1,5,21,-5,40.,0.9,0.8) 
        IF(ITY.EQ.7) CALL LU3ENT(1,21,21,21,60.,0.7,0.5) 
        IF(ITY.EQ.8) CALL LU4ENT(1,2,21,21,-2,40.,0.4,0.64,0.6,0.12,0.2) 
 
C...Seventy events with independent fragmentation and momentum cons. 
      ELSEIF(IEV.LE.200) THEN 
        ITY=1+(IEV-131)/16 
        MSTJ(2)=1+MOD(IEV-131,4) 
        MSTJ(3)=1+MOD((IEV-131)/4,4) 
        IF(ITY.EQ.1) CALL LU2ENT(1,4,-5,40.) 
        IF(ITY.EQ.2) CALL LU3ENT(1,3,21,-3,40.,0.9,0.4) 
        IF(ITY.EQ.3) CALL LU4ENT(1,2,21,21,-2,40.,0.4,0.64,0.6,0.12,0.2) 
        IF(ITY.GE.4) CALL LU4ENT(1,2,-3,3,-2,40.,0.4,0.64,0.6,0.12,0.2) 
 
C...A hundred events with random jets (check invariant mass). 
      ELSEIF(IEV.LE.300) THEN 
  100   DO 110 J=1,5 
        PSUM(J)=0. 
  110   CONTINUE 
        NJET=2.+6.*RLU(0) 
        DO 130 I=1,NJET 
        KFL=21 
        IF(I.EQ.1) KFL=INT(1.+4.*RLU(0)) 
        IF(I.EQ.NJET) KFL=-INT(1.+4.*RLU(0)) 
        EJET=5.+20.*RLU(0) 
        THETA=ACOS(2.*RLU(0)-1.) 
        PHI=6.2832*RLU(0) 
        IF(I.LT.NJET) CALL LU1ENT(-I,KFL,EJET,THETA,PHI) 
        IF(I.EQ.NJET) CALL LU1ENT(I,KFL,EJET,THETA,PHI) 
        IF(I.EQ.1.OR.I.EQ.NJET) MSTJ(93)=1 
        IF(I.EQ.1.OR.I.EQ.NJET) PSUM(5)=PSUM(5)+ULMASS(KFL) 
        DO 120 J=1,4 
        PSUM(J)=PSUM(J)+P(I,J) 
  120   CONTINUE 
  130   CONTINUE 
        IF(PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2.LT. 
     &  (PSUM(5)+PARJ(32))**2) GOTO 100 
 
C...Fifty e+e- continuum events with matrix elements. 
      ELSEIF(IEV.LE.350) THEN 
        MSTJ(101)=2 
        CALL LUEEVT(0,40.) 
 
C...Fifty e+e- continuum event with varying shower options. 
      ELSEIF(IEV.LE.400) THEN 
        MSTJ(42)=1+MOD(IEV,2) 
        MSTJ(43)=1+MOD(IEV/2,4) 
        MSTJ(44)=MOD(IEV/8,3) 
        CALL LUEEVT(0,90.) 
 
C...Fifty e+e- continuum events with coherent shower, including top. 
      ELSEIF(IEV.LE.450) THEN 
        MSTJ(104)=6 
        CALL LUEEVT(0,500.) 
 
C...Fifty Upsilon decays to ggg or gammagg with coherent shower. 
      ELSEIF(IEV.LE.500) THEN 
        CALL LUONIA(5,9.46) 
 
C...One decay each for some heavy mesons. 
      ELSEIF(IEV.LE.560) THEN 
        ITY=IEV-501 
        KFLS=2*(ITY/20)+1 
        KFLB=8-MOD(ITY/5,4) 
        KFLC=KFLB-MOD(ITY,5) 
        CALL LU1ENT(1,100*KFLB+10*KFLC+KFLS,0.,0.,0.) 
 
C...One decay each for some heavy baryons. 
      ELSEIF(IEV.LE.600) THEN 
        ITY=IEV-561 
        KFLS=2*(ITY/20)+2 
        KFLA=8-MOD(ITY/5,4) 
        KFLB=KFLA-MOD(ITY,5) 
        KFLC=MAX(1,KFLB-1) 
        CALL LU1ENT(1,1000*KFLA+100*KFLB+10*KFLC+KFLS,0.,0.,0.) 
      ENDIF 
 
C...Generate event. Find total momentum, energy and charge. 
      DO 140 J=1,4 
      PINI(J)=PLU(0,J) 
  140 CONTINUE 
      PINI(6)=PLU(0,6) 
      CALL LUEXEC 
      DO 150 J=1,4 
      PFIN(J)=PLU(0,J) 
  150 CONTINUE 
      PFIN(6)=PLU(0,6) 
 
C...Check conservation of energy, momentum and charge; 
C...usually exact, but only approximate for single jets. 
      MERR=0 
      IF(IEV.LE.50) THEN 
        IF((PFIN(1)-PINI(1))**2+(PFIN(2)-PINI(2))**2.GE.4.) MERR=MERR+1 
        EPZREM=PINI(4)+PINI(3)-PFIN(4)-PFIN(3) 
        IF(EPZREM.LT.0..OR.EPZREM.GT.2.*PARJ(31)) MERR=MERR+1 
        IF(ABS(PFIN(6)-PINI(6)).GT.2.1) MERR=MERR+1 
      ELSE 
        DO 160 J=1,4 
        IF(ABS(PFIN(J)-PINI(J)).GT.0.0001*PINI(4)) MERR=MERR+1 
  160   CONTINUE 
        IF(ABS(PFIN(6)-PINI(6)).GT.0.1) MERR=MERR+1 
      ENDIF 
      IF(MERR.NE.0) WRITE(MSTU(11),5000) (PINI(J),J=1,4),PINI(6), 
     &(PFIN(J),J=1,4),PFIN(6) 
 
C...Check that all KF codes are known ones, and that partons/particles 
C...satisfy energy-momentum-mass relation. Store particle statistics. 
      DO 170 I=1,N 
      IF(K(I,1).GT.20) GOTO 170 
      IF(LUCOMP(K(I,2)).EQ.0) THEN 
        WRITE(MSTU(11),5100) I 
        MERR=MERR+1 
      ENDIF 
      PD=P(I,4)**2-P(I,1)**2-P(I,2)**2-P(I,3)**2-P(I,5)**2 
      IF(ABS(PD).GT.MAX(0.1,0.001*P(I,4)**2).OR.P(I,4).LT.0.) THEN 
        WRITE(MSTU(11),5200) I 
        MERR=MERR+1 
      ENDIF 
  170 CONTINUE 
      IF(MTEST.GE.1) CALL LUTABU(21) 
 
C...List all erroneous events and some normal ones. 
      IF(MERR.NE.0.OR.MSTU(24).NE.0.OR.MSTU(28).NE.0) THEN 
        CALL LULIST(2) 
      ELSEIF(MTEST.GE.1.AND.MOD(IEV-5,100).EQ.0) THEN 
        CALL LULIST(1) 
      ENDIF 
 
C...Stop execution if too many errors. 
      IF(MERR.NE.0) NERR=NERR+1 
      IF(NERR.GE.10) THEN 
        WRITE(MSTU(11),5300) IEV 
        STOP 
      ENDIF 
  180 CONTINUE 
 
C...Summarize result of run. 
      IF(MTEST.GE.1) CALL LUTABU(22) 
      IF(NERR.EQ.0) WRITE(MSTU(11),5400) 
      IF(NERR.GT.0) WRITE(MSTU(11),5500) NERR 
 
C...Reset commonblock variables changed during run. 
      MSTJ(2)=3 
      PARJ(17)=0. 
      PARJ(22)=1. 
      PARJ(43)=0.5 
      PARJ(54)=0. 
      MSTJ(105)=1 
      MSTJ(107)=0 
 
C...Format statements for output. 
 5000 FORMAT(/' Momentum, energy and/or charge were not conserved ', 
     &'in following event'/' sum of',9X,'px',11X,'py',11X,'pz',11X, 
     &'E',8X,'charge'/' before',2X,4(1X,F12.5),1X,F8.2/' after',3X, 
     &4(1X,F12.5),1X,F8.2) 
 5100 FORMAT(/5X,'Entry no.',I4,' in following event not known code') 
 5200 FORMAT(/5X,'Entry no.',I4,' in following event has faulty ', 
     &'kinematics') 
 5300 FORMAT(/5X,'Ten errors experienced by event ',I3/ 
     &5X,'Something is seriously wrong! Execution stopped now!') 
 5400 FORMAT(//5X,'End result of LUTEST: no errors detected.') 
 5500 FORMAT(//5X,'End result of LUTEST:',I2,' errors detected.'/ 
     &5X,'This should not have happened!') 
 
      RETURN 
      END 
 
C********************************************************************* 
 
      BLOCK DATA LUDATA 
 
C...Purpose: to give default values to parameters and particle and 
C...decay data. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      COMMON/LUDAT4/CHAF(500) 
      CHARACTER CHAF*8 
      COMMON/LUDATR/MRLU(6),RRLU(100) 
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/,/LUDAT4/,/LUDATR/ 
 
C...LUDAT1, containing status codes and most parameters. 
      DATA MSTU/ 
     &    0,    0,    0, 4000,10000,  500, 2000,    0,    0,    2, 
     1    6,    1,    1,    0,    1,    1,    0,    0,    0,    0, 
     2    2,   10,    0,    0,    1,   10,    0,    0,    0,    0, 
     3    0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     4    2,    2,    1,    4,    2,    1,    1,    0,    0,    0, 
     5   25,   24,    0,    1,    0,    0,    0,    0,    0,    0, 
     6    0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     7  30*0, 
     &    1,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     1    1,    5,    3,    5,    0,    0,    0,    0,    0,    0, 
     2  60*0, 
     8    7,  408, 1995,   08,   23,  700,    0,    0,    0,    0, 
     9    0,    0,    0,    0,    0,    0,    0,    0,    0,    0/ 
      DATA PARU/ 
     & 3.1415927, 6.2831854, 0.1973, 5.068, 0.3894, 2.568,   4*0., 
     1 0.001, 0.09, 0.01,  0.,   0.,   0.,   0.,   0.,   0.,   0., 
     2   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     3   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     4  2.0,  1.0, 0.25,  2.5, 0.05,   0.,   0., 0.0001, 0.,   0., 
     5  2.5,  1.5,  7.0,  1.0,  0.5,  2.0,  3.2,   0.,   0.,   0., 
     6  40*0., 
     & 0.00729735, 0.232, 0.007764, 1.0, 1.16639E-5, 0., 0., 0., 
     &   0.,   0., 
     1 0.20, 0.25,  1.0,  4.0,  10.,   0.,   0.,   0.,   0.,   0., 
     2 -0.693, -1.0, 0.387, 1.0, -0.08, -1.0, 1.0, 1.0, 1.0,   0., 
     3  1.0, -1.0,  1.0, -1.0,  1.0,   0.,   0.,   0.,   0.,   0., 
     4  5.0,  1.0,  1.0,   0.,  1.0,  1.0,   0.,   0.,   0.,   0., 
     5  1.0,   0.,   0.,   0., 1000., 1.0,  1.0,  1.0,  1.0,   0., 
     6  1.0,  1.0,  1.0,  1.0,  1.0,   0.,   0.,   0.,   0.,   0., 
     7  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,   0.,   0.,   0., 
     8  1.0,  1.0,  1.0,  0.0,  0.0,  1.0,  1.0,  0.0,  0.0,   0., 
     9   0.,   0.,   0.,   0.,  1.0,   0.,   0.,   0.,   0.,   0./ 
      DATA MSTJ/ 
     &    1,    3,    0,    0,    0,    0,    0,    0,    0,    0, 
     1    4,    2,    0,    1,    0,    0,    0,    0,    0,    0, 
     2    2,    1,    1,    2,    1,    2,    2,    0,    0,    0, 
     3    0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 
     4    2,    2,    4,    2,    5,    3,    3,    0,    0,    0, 
     5    0,    3,    0,    0,    0,    0,    0,    0,    0,    0, 
     6  40*0, 
     &    5,    2,    7,    5,    1,    1,    0,    2,    0,    2, 
     1    0,    0,    0,    0,    1,    1,    0,    0,    0,    0, 
     2  80*0/ 
      DATA PARJ/ 
     & 0.10, 0.30, 0.40, 0.05, 0.50, 0.50, 0.50,   0.,   0.,   0., 
     1 0.50, 0.60, 0.75,   0.,   0.,   0.,   0.,  1.0,  1.0,   0., 
     2 0.36,  1.0, 0.01,  2.0,  1.0,  0.4,   0.,   0.,   0.,   0., 
     3 0.10,  1.0,  0.8,  1.5,   0.,  2.0,  0.2,  2.5,  0.6,   0., 
     4  0.3, 0.58,  0.5,  0.9,  0.5,  1.0,  1.0,  1.0,   0.,   0., 
     5 0.77,0.77,0.77,-0.05,-0.005,-0.00001,-0.00001,-0.00001,1.0,0., 
     6  4.5,  0.7,  0., 0.003,  0.5,  0.5,   0.,   0.,   0.,   0., 
     7  10., 1000., 100., 1000., 0.,  0.7,  10.,   0.,   0.,   0., 
     8 0.29,  1.0,  1.0,   0.,  10.,  10.,   0.,   0.,   0.,   0., 
     9 0.02,  1.0,  0.2,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     &   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     1   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     2  1.0, 0.25,91.187,2.489, 0.01, 2.0,  1.0, 0.25,0.002,   0., 
     3   0.,   0.,   0.,   0., 0.01, 0.99,   0.,   0.,  0.2,   0., 
     4  60*0./ 
 
C...LUDAT2, with particle data and flavour treatment parameters. 
      DATA (KCHG(I,1),I=   1, 500)/-1,2,-1,2,-1,2,-1,2,2*0,-3,0,-3,0, 
     &-3,0,-3,6*0,3,9*0,3,2*0,3,0,-1,44*0,2,-1,2,-1,2,3,11*0,3,0,2*3,0, 
     &3,0,3,0,3,10*0,3,0,2*3,0,3,0,3,0,3,10*0,3,0,2*3,0,3,0,3,0,3,10*0, 
     &3,0,2*3,0,3,0,3,0,3,10*0,3,0,2*3,0,3,0,3,0,3,10*0,3,0,2*3,0,3,0, 
     &3,0,3,70*0,3,0,3,28*0,3,2*0,3,8*0,-3,8*0,3,0,-3,0,3,-3,3*0,3,6,0, 
     &3,5*0,-3,0,3,-3,0,-3,4*0,-3,0,3,6,-3,0,3,-3,0,-3,0,3,6,0,3,5*0, 
     &-3,0,3,-3,0,-3,114*0/ 
      DATA (KCHG(I,2),I=   1, 500)/8*1,12*0,2,16*0,2,1,50*0,-1,410*0/ 
      DATA (KCHG(I,3),I=   1, 500)/8*1,2*0,8*1,5*0,1,9*0,1,2*0,1,0,2*1, 
     &41*0,1,0,7*1,10*0,10*1,10*0,10*1,10*0,10*1,10*0,10*1,10*0,10*1, 
     &10*0,10*1,70*0,3*1,22*0,1,5*0,1,0,2*1,6*0,1,0,2*1,6*0,2*1,0,5*1, 
     &0,6*1,4*0,6*1,4*0,16*1,4*0,6*1,114*0/ 
      DATA (PMAS(I,1),I=   1, 500)/0.0099,0.0056,0.199,1.35,5.,160., 
     &2*250.,2*0.,0.00051,0.,0.1057,0.,1.777,0.,250.,5*0.,91.187,80.25, 
     &80.,6*0.,500.,900.,500.,3*300.,350.,200.,5000.,60*0.,0.1396, 
     &0.4977,0.4936,1.8693,1.8645,1.9688,5.2787,5.2786,5.47972,6.594, 
     &0.135,0.5475,0.9578,2.9788,9.4,320.,2*500.,2*0.,0.7669,0.8961, 
     &0.8916,2.0101,2.0071,2.11,2*5.325,5.5068,6.602,0.7683,0.782, 
     &1.0194,3.0969,9.4603,320.,2*500.,2*0.,1.232,2*1.29,2*2.424,2.536, 
     &2*5.73,5.97,7.3,1.232,1.17,1.4,3.46,9.875,320.,2*500.,2*0.,0.983, 
     &2*1.429,2*2.272,2.5,2*5.68,5.92,7.25,0.9827,1.,1.4,3.4151,9.8598, 
     &320.,2*500.,2*0.,1.26,2*1.402,2*2.372,2.56,2*5.78,6.02,7.3,1.26, 
     &1.282,1.42,3.5106,9.8919,320.,2*500.,2*0.,1.318,1.432,1.425, 
     &2*2.46,2.61,2*5.83,6.07,7.35,1.318,1.275,1.525,3.5562,9.9132, 
     &320.,2*500.,2*0.,2*0.4977,8*0.,3.686,3*0.,10.0233,70*0.,1.1156, 
     &5*0.,2.2849,0.,2.473,2.466,6*0.,5.641,0.,2*5.84,6*0.,0.9396, 
     &0.9383,0.,1.1974,1.1926,1.1894,1.3213,1.3149,0.,2.4525,2.4529, 
     &2.4527,2*2.55,2.73,4*0.,3*5.8,2*5.96,6.12,4*0.,1.234,1.233,1.232, 
     &1.231,1.3872,1.3837,1.3828,1.535,1.5318,1.6724,3*2.5,2*2.63,2.8, 
     &4*0.,3*5.81,2*5.97,6.13,114*0./ 
      DATA (PMAS(I,2),I=   1, 500)/22*0.,2.489,2.066,88*0.,0.0002, 
     &0.001,6*0.,0.149,0.0505,0.0498,7*0.,0.151,0.00843,0.0044,7*0., 
     &0.155,2*0.09,2*0.02,0.,4*0.05,0.155,0.36,0.08,2*0.01,5*0.,0.057, 
     &2*0.287,7*0.05,0.057,0.,0.25,0.014,6*0.,0.4,2*0.174,7*0.05,0.4, 
     &0.024,0.06,0.0009,6*0.,0.11,0.109,0.098,2*0.019,5*0.02,0.11, 
     &0.185,0.076,0.002,146*0.,4*0.12,0.0394,0.036,0.0358,0.0099, 
     &0.0091,131*0./ 
      DATA (PMAS(I,3),I=   1, 500)/22*0.,2*20.,88*0.,0.002,0.005,6*0., 
     &0.4,2*0.2,7*0.,0.4,0.1,0.015,7*0.,0.25,0.005,0.01,2*0.08,0., 
     &4*0.1,0.25,0.2,0.001,2*0.02,5*0.,0.05,2*0.4,6*0.1,2*0.05,0.,0.35, 
     &0.05,6*0.,3*0.3,2*0.1,0.03,4*0.1,0.3,0.05,0.02,0.001,6*0.,0.25, 
     &4*0.12,5*0.05,0.25,0.17,0.2,0.01,146*0.,4*0.14,0.04,2*0.035, 
     &2*0.05,131*0./ 
      DATA (PMAS(I,4),I=   1, 500)/12*0.,658650.,0.,0.0914,68*0.,0.1, 
     &0.387,15*0.,7804.,0.,3709.,0.32,0.1259,0.135,3*0.387,0.15,110*0., 
     &15500.,26.75,83*0.,78.88,5*0.,0.057,0.,0.025,0.09,6*0.,0.387,0., 
     &2*0.387,9*0.,44.3,0.,23.95,49.1,86.9,6*0.,0.13,9*0.,0.387,13*0., 
     &24.60001,130*0./ 
      DATA PARF/ 
     &  0.5, 0.25,  0.5, 0.25,   1.,  0.5,   0.,   0.,   0.,   0., 
     1  0.5,   0.,  0.5,   0.,   1.,   1.,   0.,   0.,   0.,   0., 
     2  0.5,   0.,  0.5,   0.,   1.,   1.,   0.,   0.,   0.,   0., 
     3  0.5,   0.,  0.5,   0.,   1.,   1.,   0.,   0.,   0.,   0., 
     4  0.5,   0.,  0.5,   0.,   1.,   1.,   0.,   0.,   0.,   0., 
     5  0.5,   0.,  0.5,   0.,   1.,   1.,   0.,   0.,   0.,   0., 
     6 0.75,  0.5,   0., 0.1667, 0.0833, 0.1667, 0., 0., 0.,   0., 
     7   0.,   0.,   1., 0.3333, 0.6667, 0.3333, 0., 0., 0.,   0., 
     8   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     9   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     & 0.325, 0.325, 0.5, 1.6,  5.0,   0.,   0.,   0.,   0.,   0., 
     1   0., 0.11, 0.16, 0.048, 0.50, 0.45, 0.55, 0.60,  0.,   0., 
     2  0.2,  0.1,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 
     3  1870*0./ 
      DATA ((VCKM(I,J),J=1,4),I=1,4)/ 
     1  0.95113,  0.04884,  0.00003,  0.00000, 
     2  0.04884,  0.94940,  0.00176,  0.00000, 
     3  0.00003,  0.00176,  0.99821,  0.00000, 
     4  0.00000,  0.00000,  0.00000,  1.00000/ 
 
C...LUDAT3, with particle decay parameters and data. 
      DATA (MDCY(I,1),I=   1, 500)/5*0,3*1,6*0,1,0,1,5*0,3*1,6*0,1,0,1, 
     &2*0,4*1,42*0,7*1,12*0,1,0,15*1,2*0,18*1,2*0,18*1,2*0,18*1,2*0, 
     &18*1,2*0,18*1,3*0,1,8*0,1,3*0,1,70*0,1,5*0,1,0,2*1,6*0,1,0,2*1, 
     &9*0,5*1,0,6*1,4*0,6*1,4*0,16*1,4*0,6*1,114*0/ 
      DATA (MDCY(I,2),I=   1, 500)/1,9,17,25,33,41,50,60,2*0,70,74,76, 
     &81,83,124,126,132,2*0,135,144,156,172,192,6*0,209,0,231,254,274, 
     &292,301,304,305,42*0,314,315,319,328,331,336,338,11*0,358,359, 
     &361,367,430,491,524,560,596,635,666,668,675,681,682,683,684,685, 
     &2*0,686,688,691,694,697,699,700,701,702,703,704,708,713,721,724, 
     &733,734,735,2*0,736,737,742,747,749,751,753,755,757,759,761,762, 
     &765,769,770,771,772,773,2*0,774,775,777,779,781,783,785,787,789, 
     &791,793,794,799,804,806,808,809,810,2*0,811,813,815,817,819,821, 
     &823,825,827,829,831,833,846,850,852,854,855,856,2*0,857,863,873, 
     &884,892,900,904,912,920,924,928,936,945,951,953,955,956,957,2*0, 
     &958,966,8*0,968,3*0,979,70*0,993,5*0,997,0,1073,1074,6*0,1075,0, 
     &1092,1093,9*0,1094,1096,1097,1100,1101,0,1103,1104,1105,1106, 
     &1107,1108,4*0,1109,1110,1111,1112,1113,1114,4*0,1115,1116,1119, 
     &1122,1123,1126,1129,1132,1134,1136,1140,1141,1142,1143,1145,1147, 
     &4*0,1148,1149,1150,1151,1152,1153,114*0/ 
      DATA (MDCY(I,3),I=   1, 500)/5*8,9,2*10,2*0,4,2,5,2,41,2,6,3,2*0, 
     &9,12,16,20,17,6*0,22,0,23,20,18,9,3,1,9,42*0,1,4,9,3,5,2,20,11*0, 
     &1,2,6,63,61,33,2*36,39,31,2,7,6,5*1,2*0,2,3*3,2,5*1,4,5,8,3,9, 
     &3*1,2*0,1,2*5,7*2,1,3,4,5*1,2*0,1,9*2,1,2*5,2*2,3*1,2*0,11*2,13, 
     &4,2*2,3*1,2*0,6,10,11,2*8,4,2*8,2*4,8,9,6,2*2,3*1,2*0,8,2,8*0,11, 
     &3*0,14,70*0,4,5*0,76,0,2*1,6*0,17,0,2*1,9*0,2,1,3,1,2,0,6*1,4*0, 
     &6*1,4*0,1,2*3,1,3*3,2*2,4,3*1,2*2,1,4*0,6*1,114*0/ 
      DATA (MDME(I,1),I=   1,2000)/6*1,-1,7*1,-1,7*1,-1,7*1,-1,7*1,-1, 
     &7*1,-1,1,-1,8*1,2*-1,8*1,2*-1,61*1,-1,2*1,-1,6*1,2*-1,7*1,2*-1, 
     &3*1,-1,6*1,2*-1,6*1,2*-1,3*1,-1,3*1,-1,3*1,5*-1,3*1,-1,6*1,2*-1, 
     &3*1,-1,11*1,2*-1,6*1,8*-1,3*1,-1,3*1,-1,3*1,5*-1,3*1,4*-1,6*1, 
     &2*-1,3*1,-1,5*1,-1,8*1,2*-1,3*1,-1,9*1,-1,3*1,-1,9*1,2*-1,2*1,-1, 
     &16*1,-1,2*1,3*-1,1665*1/ 
      DATA (MDME(I,2),I=   1,2000)/75*102,42,6*102,2*42,2*0,7*41,2*0, 
     &24*41,6*102,45,29*102,8*32,8*0,16*32,4*0,8*32,4*0,32,4*0,8*32, 
     &14*0,16*32,7*0,8*32,4*0,32,7*0,8*32,4*0,32,5*0,4*32,5*0,3*32,0, 
     &6*32,3*0,12,2*42,2*11,9*42,2*45,31,2*45,2*33,31,2*45,20*46,7*0, 
     &24*42,41*0,16*42,46*0,10*42,20*0,2*13,14*42,16*0,48,3*13,16*42, 
     &16*0,48,3*13,16*42,19*0,48,3*13,2*42,0,2*11,28*42,0,2,4*0,2,8*0, 
     &12,32,86,87,88,3,0,2*3,0,2*3,0,2*3,0,3,6*0,3,3*0,1,0,3,2*0,2*3, 
     &3*0,1,4*0,12,3*0,4*32,2*4,86,87,88,33*0,12,32,86,87,88,31*0,12,0, 
     &32,86,87,88,40*0,12,0,32,86,87,88,95*0,12,0,32,86,87,88,2*0,4*42, 
     &6*0,12,11*0,4*32,2*4,9*0,14*42,52*0,10*13,2*84,3*42,8*0,48,3*13, 
     &2*42,2*85,14*0,84,5*0,85,886*0/ 
      DATA (BRAT(I)  ,I=   1, 439)/75*0.,1.,6*0.,0.179,0.178,0.116, 
     &0.235,0.005,0.056,0.018,0.023,0.011,2*0.004,0.0067,0.014,2*0.002, 
     &2*0.001,0.0022,0.054,0.002,0.016,0.005,0.011,0.0101,5*0.006, 
     &0.002,2*0.001,5*0.002,6*0.,1.,29*0.,0.15394,0.11936,0.15394, 
     &0.11926,0.15254,3*0.,0.03368,0.06664,0.03368,0.06664,0.03368, 
     &0.06664,2*0.,0.3214,0.0165,2*0.,0.0165,0.3207,2*0.,0.00001, 
     &0.00059,6*0.,3*0.1081,3*0.,0.0003,0.048,0.8705,4*0.,0.0002, 
     &0.0603,0.,0.0199,0.0008,3*0.,0.143,0.111,0.143,0.111,0.143,0.085, 
     &2*0.,0.03,0.058,0.03,0.058,0.03,0.058,8*0.,0.25,0.01,2*0.,0.01, 
     &0.25,4*0.,0.24,5*0.,3*0.08,6*0.,0.01,0.08,0.82,5*0.,0.09,11*0., 
     &0.01,0.08,0.82,5*0.,0.09,9*0.,1.,6*0.,0.01,0.98,0.01,1.,4*0.215, 
     &2*0.,2*0.07,0.,1.,2*0.08,0.76,0.08,2*0.105,0.04,0.5,0.08,0.14, 
     &0.01,0.015,0.005,1.,3*0.,1.,4*0.,1.,0.25,0.01,2*0.,0.01,0.25, 
     &4*0.,0.24,5*0.,3*0.08,0.,1.,2*0.5,0.635,0.212,0.056,0.017,0.048, 
     &0.032,0.07,0.065,2*0.005,2*0.011,5*0.001,0.07,0.065,2*0.005, 
     &2*0.011,5*0.001,0.026,0.019,0.066,0.041,0.045,0.076,0.0073, 
     &2*0.0047,0.026,0.001,0.0006,0.0066,0.005,2*0.003,2*0.0006, 
     &2*0.001,0.006,0.005,0.012,0.0057,0.067,0.008,0.0022,0.027,0.004, 
     &0.019,0.012,0.002,0.009,0.0218,0.001,0.022,0.087,0.001,0.0019, 
     &0.0015,0.0028,0.034,0.027,2*0.002,2*0.004,2*0.002,0.034,0.027/ 
      DATA (BRAT(I)  ,I= 440, 655)/2*0.002,2*0.004,2*0.002,0.0365, 
     &0.045,0.073,0.062,3*0.021,0.0061,0.015,0.025,0.0088,0.074,0.0109, 
     &0.0041,0.002,0.0035,0.0011,0.001,0.0027,2*0.0016,0.0018,0.011, 
     &0.0063,0.0052,0.018,0.016,0.0034,0.0036,0.0009,0.0006,0.015, 
     &0.0923,0.018,0.022,0.0077,0.009,0.0075,0.024,0.0085,0.067,0.0511, 
     &0.017,0.0004,0.0028,0.01,2*0.02,0.03,2*0.005,2*0.02,0.03,2*0.005, 
     &0.015,0.037,0.028,0.079,0.095,0.052,0.0078,4*0.001,0.028,0.033, 
     &0.026,0.05,0.01,4*0.005,0.25,0.0952,0.02,0.055,2*0.005,0.008, 
     &0.012,0.02,0.055,2*0.005,0.008,0.012,0.01,0.03,0.0035,0.011, 
     &0.0055,0.0042,0.009,0.018,0.015,0.0185,0.0135,0.025,0.0004, 
     &0.0007,0.0008,0.0014,0.0019,0.0025,0.4291,0.08,0.07,0.02,0.015, 
     &0.005,0.02,0.055,2*0.005,0.008,0.012,0.02,0.055,2*0.005,0.008, 
     &0.012,0.01,0.03,0.0035,0.011,0.0055,0.0042,0.009,0.018,0.015, 
     &0.0185,0.0135,0.025,0.0004,0.0007,0.0008,0.0014,0.0019,0.0025, 
     &0.4291,0.08,0.07,0.02,0.015,0.005,0.02,0.055,2*0.005,0.008,0.012, 
     &0.02,0.055,2*0.005,0.008,0.012,0.01,0.03,0.0035,0.011,0.0055, 
     &0.0042,0.009,0.018,0.015,0.0185,0.0135,0.025,2*0.0002,0.0007, 
     &2*0.0004,0.0014,0.001,0.0009,0.0025,0.4291,0.08,0.07,0.02,0.015, 
     &0.005,0.047,0.122,0.006,0.012,0.035,0.012,0.035,0.003,0.007,0.15, 
     &0.037,0.008,0.002,0.05,0.015,0.003,0.001,0.014,0.042,0.014,0.042/ 
      DATA (BRAT(I)  ,I= 656, 931)/0.24,0.065,0.012,0.003,0.001,0.002, 
     &0.001,0.002,0.014,0.003,0.988,0.012,0.389,0.319,0.2367,0.049, 
     &0.005,0.001,0.0003,0.441,0.206,0.3,0.03,0.022,0.001,5*1.,0.99955, 
     &0.00045,0.665,0.333,0.002,0.666,0.333,0.001,0.65,0.3,0.05,0.56, 
     &0.44,5*1.,0.99912,0.00079,0.00005,0.00004,0.888,0.085,0.021, 
     &2*0.003,0.49,0.344,3*0.043,0.023,0.013,0.001,0.0627,0.0597, 
     &0.8776,3*0.027,0.015,0.045,0.015,0.045,0.77,0.029,4*1.,0.28,0.14, 
     &0.313,0.157,0.11,0.28,0.14,0.313,0.157,0.11,0.667,0.333,0.667, 
     &0.333,2*0.5,0.667,0.333,0.667,0.333,4*0.5,1.,0.333,0.334,0.333, 
     &4*0.25,6*1.,0.667,0.333,0.667,0.333,0.667,0.333,0.667,0.333, 
     &2*0.5,0.667,0.333,0.667,0.333,4*0.5,1.,0.52,0.26,0.11,2*0.055, 
     &0.62,0.31,0.035,2*0.0175,0.007,0.993,0.02,0.98,3*1.,2*0.5,0.667, 
     &0.333,0.667,0.333,0.667,0.333,0.667,0.333,2*0.5,0.667,0.333, 
     &0.667,0.333,6*0.5,3*0.12,0.097,0.043,4*0.095,4*0.03,4*0.25,0.273, 
     &0.727,0.35,0.65,3*1.,2*0.35,0.144,0.105,0.048,0.003,0.333,0.166, 
     &0.168,0.084,0.087,0.043,0.059,2*0.029,0.002,0.332,0.166,0.168, 
     &0.084,0.086,0.043,0.059,2*0.029,2*0.002,0.3,0.15,0.16,0.08,0.13, 
     &0.06,0.08,0.04,0.3,0.15,0.16,0.08,0.13,0.06,0.08,0.04,2*0.3, 
     &2*0.2,0.3,0.15,0.16,0.08,0.13,0.06,0.08,0.04,0.3,0.15,0.16,0.08, 
     &0.13,0.06,0.08,0.04,2*0.3,2*0.2,2*0.3,2*0.2,2*0.35,0.144,0.105/ 
      DATA (BRAT(I)  ,I= 932,2000)/0.024,2*0.012,0.003,0.566,0.283, 
     &0.069,0.028,0.023,2*0.0115,0.005,0.003,0.356,2*0.178,0.28, 
     &2*0.004,0.135,0.865,0.22,0.78,3*1.,0.217,0.124,2*0.193,2*0.135, 
     &0.002,0.001,0.686,0.314,2*0.0083,0.1866,0.324,0.184,0.027,0.001, 
     &0.093,0.087,0.078,0.0028,3*0.014,0.008,0.024,0.008,0.024,0.425, 
     &0.02,0.185,0.088,0.043,0.067,0.066,0.641,0.357,2*0.001,0.018, 
     &2*0.005,0.003,0.002,2*0.006,0.018,2*0.005,0.003,0.002,2*0.006, 
     &0.0066,0.025,0.016,0.0088,2*0.005,0.0058,0.005,0.0055,4*0.004, 
     &2*0.002,2*0.004,0.003,0.002,2*0.003,3*0.002,2*0.001,0.002, 
     &2*0.001,2*0.002,0.0013,0.0018,5*0.001,4*0.003,2*0.005,2*0.002, 
     &2*0.001,2*0.002,2*0.001,0.2432,0.057,2*0.035,0.15,2*0.075,0.03, 
     &2*0.015,2*1.,2*0.105,0.04,0.0077,0.02,0.0235,0.0285,0.0435, 
     &0.0011,0.0022,0.0044,0.4291,0.08,0.07,0.02,0.015,0.005,2*1., 
     &0.999,0.001,1.,0.516,0.483,0.001,1.,0.995,0.005,13*1.,0.331, 
     &0.663,0.006,0.663,0.331,0.006,1.,0.88,2*0.06,0.88,2*0.06,0.88, 
     &2*0.06,0.667,2*0.333,0.667,0.676,0.234,0.085,0.005,3*1.,4*0.5, 
     &7*1.,847*0./ 
      DATA (KFDP(I,1),I=   1, 507)/21,22,23,4*-24,25,21,22,23,4*24,25, 
     &21,22,23,4*-24,25,21,22,23,4*24,25,21,22,23,4*-24,25,21,22,23, 
     &4*24,25,37,21,22,23,4*-24,25,2*-37,21,22,23,4*24,25,2*37,22,23, 
     &-24,25,23,24,-12,22,23,-24,25,23,24,-12,-14,35*16,22,23,-24,25, 
     &23,24,-89,22,23,-24,25,-37,23,24,37,1,2,3,4,5,6,7,8,21,1,2,3,4,5, 
     &6,7,8,11,13,15,17,1,2,3,4,5,6,7,8,11,12,13,14,15,16,17,18,4*-1, 
     &4*-3,4*-5,4*-7,-11,-13,-15,-17,1,2,3,4,5,6,7,8,11,13,15,17,21, 
     &2*22,23,24,1,2,3,4,5,6,7,8,11,12,13,14,15,16,17,18,24,37,2*23,25, 
     &35,4*-1,4*-3,4*-5,4*-7,-11,-13,-15,-17,3*24,1,2,3,4,5,6,7,8,11, 
     &13,15,17,21,2*22,23,24,23,25,36,1,2,3,4,5,6,7,8,11,13,15,17,21, 
     &2*22,23,24,23,-1,-3,-5,-7,-11,-13,-15,-17,24,5,6,21,2,1,2,3,4,5, 
     &6,11,13,15,82,-11,-13,2*2,-12,-14,-16,2*-2,2*-4,-2,-4,2*89,37, 
     &2*-89,2*5,-37,2*89,4*-1,4*-3,4*-5,4*-7,-11,-13,-15,-17,-13,130, 
     &310,-13,3*211,12,14,11*-11,11*-13,-311,-313,-311,-313,-20313, 
     &2*-311,-313,-311,-313,2*111,2*221,2*331,2*113,2*223,2*333,-311, 
     &-313,2*-321,211,-311,-321,333,-311,-313,-321,211,2*-321,2*-311, 
     &-321,211,113,8*-11,8*-13,-321,-323,-321,-323,-311,2*-313,-311, 
     &-313,2*-311,-321,-10323,-321,-323,-321,-311,2*-313,211,111,333, 
     &3*-321,-311,-313,-321,-313,310,333,211,2*-321,-311,-313,-311,211, 
     &-321,3*-311,211,113,321,-15,5*-11,5*-13,221,331,333,221,331,333/ 
      DATA (KFDP(I,1),I= 508, 924)/10221,211,213,211,213,321,323,321, 
     &323,2212,221,331,333,221,2*2,6*12,6*14,2*16,3*-411,3*-413,2*-411, 
     &2*-413,2*441,2*443,2*20443,2*2,2*4,2,4,6*12,6*14,2*16,3*-421, 
     &3*-423,2*-421,2*-423,2*441,2*443,2*20443,2*2,2*4,2,4,6*12,6*14, 
     &2*16,3*-431,3*-433,2*-431,2*-433,3*441,3*443,3*20443,2*2,2*4,2,4, 
     &16,2*4,2*12,2*14,2*16,4*2,4*4,2*-11,2*-13,2*-1,2*-3,2*-11,2*-13, 
     &2*-1,3*22,111,211,2*22,211,22,211,111,3*22,111,82,21,3*0,2*211, 
     &321,3*311,2*321,421,2*411,2*421,431,511,521,531,541,211,111,13, 
     &11,211,22,211,2*111,321,130,-213,113,213,211,22,111,11,13,82,11, 
     &13,15,1,2,3,4,21,22,3*0,223,321,311,323,313,2*311,321,313,323, 
     &321,423,2*413,2*423,413,523,2*513,2*523,2*513,523,223,213,113, 
     &-213,313,-313,323,-323,82,21,3*0,221,321,2*311,321,421,2*411,421, 
     &411,421,521,2*511,2*521,2*511,521,221,211,111,321,130,310,211, 
     &111,321,130,310,443,82,553,21,3*0,113,213,323,2*313,323,423, 
     &2*413,2*423,413,523,2*513,2*523,2*513,523,213,-213,10211,10111, 
     &-10211,2*221,213,2*113,-213,2*321,2*311,313,-313,323,-323,443,82, 
     &553,21,3*0,213,113,221,223,321,211,321,311,323,313,323,313,321, 
     &4*311,321,313,323,313,323,311,4*321,421,411,423,413,423,413,421, 
     &2*411,421,413,423,413,423,411,2*421,411,423,413,521,511,523,513, 
     &523,513,521,2*511,521,513,523,513,523,511,2*521,511,523,513,511/ 
      DATA (KFDP(I,1),I= 925,2000)/521,513,523,213,-213,221,223,321, 
     &130,310,111,211,111,2*211,321,130,310,221,111,321,130,310,221, 
     &211,111,443,82,553,21,3*0,111,211,-12,12,-14,14,211,111,211,111, 
     &11,13,82,4*443,10441,20443,445,441,11,13,15,1,2,3,4,21,22,2*553, 
     &10551,20553,555,2212,2*2112,-12,7*-11,7*-13,2*2224,2*2212,2*2214, 
     &2*3122,2*3212,2*3214,5*3222,4*3224,2*3322,3324,2*2224,7*2212, 
     &5*2214,2*2112,2*2114,2*3122,2*3212,2*3214,2*3222,2*3224,4*2,3, 
     &2*2,1,2*2,2*0,-12,-14,-16,5*4122,441,443,20443,2*-2,2*-4,-2,-4, 
     &2*0,2112,-12,3122,2212,2112,2212,3*3122,3*4122,4132,4232,0, 
     &3*5122,5132,5232,0,2112,2212,2*2112,2212,2112,2*2212,3122,3212, 
     &3112,3122,3222,3112,3122,3222,3212,3322,3312,3322,3312,3122,3322, 
     &3312,-12,3*4122,2*4132,2*4232,4332,3*5122,5132,5232,5332,847*0/ 
      DATA (KFDP(I,2),I=   1, 476)/3*1,2,4,6,8,1,3*2,1,3,5,7,2,3*3,2,4, 
     &6,8,3,3*4,1,3,5,7,4,3*5,2,4,6,8,5,3*6,1,3,5,7,6,5,3*7,2,4,6,8,7, 
     &4,6,3*8,1,3,5,7,8,5,7,2*11,12,11,12,2*11,2*13,14,13,14,13,11,13, 
     &-211,-213,-211,-213,-211,-213,3*-211,-321,-323,-321,-323,3*-321, 
     &4*-211,-213,-211,-213,-211,-213,-211,-213,-211,-213,6*-211,2*15, 
     &16,15,16,15,18,2*17,18,17,2*18,2*17,-1,-2,-3,-4,-5,-6,-7,-8,21, 
     &-1,-2,-3,-4,-5,-6,-7,-8,-11,-13,-15,-17,-1,-2,-3,-4,-5,-6,-7,-8, 
     &-11,-12,-13,-14,-15,-16,-17,-18,2,4,6,8,2,4,6,8,2,4,6,8,2,4,6,8, 
     &12,14,16,18,-1,-2,-3,-4,-5,-6,-7,-8,-11,-13,-15,-17,21,22,2*23, 
     &-24,-1,-2,-3,-4,-5,-6,-7,-8,-11,-12,-13,-14,-15,-16,-17,-18,-24, 
     &-37,22,25,2*36,2,4,6,8,2,4,6,8,2,4,6,8,2,4,6,8,12,14,16,18,23,22, 
     &25,-1,-2,-3,-4,-5,-6,-7,-8,-11,-13,-15,-17,21,22,2*23,-24,2*25, 
     &36,-1,-2,-3,-4,-5,-6,-7,-8,-11,-13,-15,-17,21,22,2*23,-24,25,2,4, 
     &6,8,12,14,16,18,25,-5,-6,21,11,-3,-4,-5,-6,-7,-8,-13,-15,-17,-82, 
     &12,14,-1,-3,11,13,15,1,4,3,4,1,3,5,3,5,6,4,21,22,4,7,5,2,4,6,8,2, 
     &4,6,8,2,4,6,8,2,4,6,8,12,14,16,18,14,2*0,14,111,211,111,-11,-13, 
     &11*12,11*14,2*211,2*213,211,20213,2*321,2*323,211,213,211,213, 
     &211,213,211,213,211,213,211,213,3*211,213,211,2*321,8*211,2*113, 
     &2*211,8*12,8*14,2*211,2*213,2*111,221,2*113,223,333,20213,211, 
     &2*321,323,2*311,313,-211,111,113,2*211,321,2*211,311,321,310,211/ 
      DATA (KFDP(I,2),I= 477, 857)/-211,4*211,321,4*211,113,2*211,-321, 
     &16,5*12,5*14,3*211,3*213,211,2*111,2*113,2*-311,2*-313,-2112, 
     &3*321,323,2*-1,6*-11,6*-13,2*-15,211,213,20213,211,213,20213,431, 
     &433,431,433,311,313,311,313,311,313,-1,-4,-3,-4,-1,-3,6*-11, 
     &6*-13,2*-15,211,213,20213,211,213,20213,431,433,431,433,321,323, 
     &321,323,321,323,-1,-4,-3,-4,-1,-3,6*-11,6*-13,2*-15,211,213, 
     &20213,211,213,20213,431,433,431,433,221,331,333,221,331,333,221, 
     &331,333,-1,-4,-3,-4,-1,-3,-15,-3,-1,2*-11,2*-13,2*-15,-1,-4,-3, 
     &-4,-3,-4,-1,-4,2*12,2*14,2,3,2,3,2*12,2*14,2,1,22,11,22,111,-211, 
     &211,11,-211,13,-211,111,113,223,22,111,-82,21,3*0,111,22,-211, 
     &111,22,211,111,22,211,111,22,111,6*22,-211,22,-13,-11,-211,111, 
     &-211,2*111,-321,310,211,111,2*-211,221,22,-11,-13,-82,-11,-13, 
     &-15,-1,-2,-3,-4,2*21,3*0,211,-213,113,-211,111,223,213,113,211, 
     &111,223,211,111,-211,111,321,311,-211,111,211,111,-321,-311,411, 
     &421,111,-211,111,211,-311,311,-321,321,-82,21,3*0,211,-211,111, 
     &211,111,211,111,-211,111,311,321,-211,111,211,111,-321,-311,411, 
     &421,111,-211,111,-321,130,310,-211,111,-321,130,310,22,-82,22,21, 
     &3*0,211,111,-211,111,211,111,211,111,-211,111,321,311,-211,111, 
     &211,111,-321,-311,411,421,-211,211,-211,111,2*211,111,-211,211, 
     &111,211,-321,2*-311,-321,-311,311,-321,321,22,-82,22,21,3*0,111/ 
      DATA (KFDP(I,2),I= 858,2000)/3*211,-311,22,-211,111,-211,111, 
     &-211,211,-213,113,223,221,211,111,211,111,2*211,213,113,223,221, 
     &22,211,111,211,111,4*211,-211,111,-211,111,-211,211,-211,211,321, 
     &311,321,311,-211,111,-211,111,-211,211,-211,2*211,111,211,111, 
     &4*211,-321,-311,-321,-311,411,421,411,421,-211,211,111,211,-321, 
     &130,310,22,-211,111,2*-211,-321,130,310,221,111,-321,130,310,221, 
     &-211,111,22,-82,22,21,3*0,111,-211,11,-11,13,-13,-211,111,-211, 
     &111,-11,-13,-82,211,111,221,111,4*22,-11,-13,-15,-1,-2,-3,-4, 
     &2*21,211,111,3*22,-211,111,22,11,7*12,7*14,-321,-323,-311,-313, 
     &-311,-313,211,213,211,213,211,213,111,221,331,113,223,111,221, 
     &113,223,321,323,321,-211,-213,111,221,331,113,223,333,10221,111, 
     &221,331,113,223,211,213,211,213,321,323,321,323,321,323,311,313, 
     &311,313,2*-1,-3,-1,2203,3201,3203,2203,2101,2103,2*0,11,13,15, 
     &-211,-213,-20213,-431,-433,3*3122,1,4,3,4,1,3,2*0,-211,11,22,111, 
     &211,22,-211,111,22,-211,111,211,2*22,0,-211,111,211,2*22,0, 
     &2*-211,111,22,111,211,22,211,2*-211,2*111,-211,2*211,111,211, 
     &-211,2*111,211,-321,-211,111,11,-211,111,211,111,22,111,2*22, 
     &-211,111,211,3*22,847*0/ 
      DATA (KFDP(I,3),I=   1, 944)/75*0,14,6*0,2*16,2*0,5*111,310,130, 
     &2*0,2*111,310,130,321,113,211,223,221,2*113,2*211,2*223,2*221, 
     &2*113,221,113,2*213,-213,195*0,4*3,4*4,1,4,3,2*2,10*81,25*0,-211, 
     &3*111,-311,-313,-311,-321,-313,-323,111,221,331,113,223,-311, 
     &-313,-311,-321,-313,-323,111,221,331,113,223,22*0,111,113,2*211, 
     &-211,-311,211,111,3*211,-211,7*211,-321,-323,-311,-321,-313,-323, 
     &-211,-213,-321,-323,-311,-321,-313,-323,-211,-213,22*0,111,113, 
     &-311,2*-211,211,-211,310,-211,2*111,211,2*-211,-321,-211,2*211, 
     &-211,111,-211,2*211,0,221,331,333,321,311,221,331,333,321,311, 
     &20*0,3,0,-411,-413,-10413,-10411,-20413,-415,-411,-413,-10413, 
     &-10411,-20413,-415,-411,-413,16*0,-4,-1,-4,-3,2*-2,-421,-423, 
     &-10423,-10421,-20423,-425,-421,-423,-10423,-10421,-20423,-425, 
     &-421,-423,16*0,-4,-1,-4,-3,2*-2,-431,-433,-10433,-10431,-20433, 
     &-435,-431,-433,-10433,-10431,-20433,-435,-431,-433,19*0,-4,-1,-4, 
     &-3,2*-2,3*0,441,443,441,443,441,443,-4,-1,-4,-3,-4,-3,-4,-1,531, 
     &533,531,533,3,2,3,2,511,513,511,513,1,2,0,-11,0,2*111,-211,-11, 
     &11,-13,2*221,3*0,111,27*0,111,2*0,22,111,5*0,111,12*0,2*21,103*0, 
     &-211,2*111,-211,3*111,-211,111,211,14*0,111,6*0,111,-211,8*0,111, 
     &-211,9*0,111,-211,111,-211,4*0,111,-211,111,-211,8*0,111,-211, 
     &111,-211,4*0,111,-211,111,-211,11*0,-211,6*0,111,211,4*0,111/ 
      DATA (KFDP(I,3),I= 945,2000)/13*0,2*111,211,-211,211,-211,7*0, 
     &-211,111,13*0,2*21,-211,111,6*0,2212,3122,3212,3214,2112,2114, 
     &2212,2112,3122,3212,3214,2112,2114,2212,2112,52*0,3*3,1,8*0, 
     &3*4122,8*0,4,1,4,3,2*2,3*0,2112,43*0,3322,861*0/ 
      DATA (KFDP(I,4),I=   1,2000)/88*0,3*111,8*0,-211,0,-211,3*0,111, 
     &2*-211,0,111,0,2*111,113,221,111,-213,-211,211,195*0,13*81,41*0, 
     &111,211,111,211,7*0,111,211,111,211,35*0,2*-211,2*111,211,111, 
     &-211,2*211,2*-211,2*0,-211,111,-211,111,4*0,-211,111,-211,111, 
     &34*0,111,-211,3*111,3*-211,2*111,3*-211,4*0,-321,-311,3*0,-321, 
     &-311,20*0,-3,31*0,6*1,30*0,6*2,33*0,6*3,9*0,8*4,4*0,4*-5,4*0, 
     &2*-5,7*0,-11,264*0,111,-211,4*0,111,57*0,-211,111,5*0,-211,111, 
     &52*0,2101,2103,2*2101,19*0,6*2101,909*0/ 
      DATA (KFDP(I,5),I=   1,2000)/90*0,111,16*0,111,7*0,111,0,2*111, 
     &303*0,-211,2*111,-211,111,-211,111,54*0,111,-211,3*111,-211,111, 
     &1510*0/ 
 
C...LUDAT4, with character strings. 
      DATA (CHAF(I)  ,I=   1, 281)/'d','u','s','c','b','t','l','h', 
     &2*' ','e','nu_e','mu','nu_mu','tau','nu_tau','chi','nu_chi', 
     &2*' ','g','gamma','Z','W','H',2*' ','reggeon','pomeron',2*' ', 
     &'Z''','Z"','W''','H''','A','H','eta_tech','LQ_ue','R',40*' ', 
     &'specflav','rndmflav','phasespa','c-hadron','b-hadron', 
     &'t-hadron','l-hadron','h-hadron','Wvirt','diquark','cluster', 
     &'string','indep.','CMshower','SPHEaxis','THRUaxis','CLUSjet', 
     &'CELLjet','table',' ','pi',2*'K',2*'D','D_s',2*'B','B_s','B_c', 
     &'pi','eta','eta''','eta_c','eta_b','eta_t','eta_l','eta_h',2*' ', 
     &'rho',2*'K*',2*'D*','D*_s',2*'B*','B*_s','B*_c','rho','omega', 
     &'phi','J/psi','Upsilon','Theta','Theta_l','Theta_h',2*' ','b_1', 
     &2*'K_1',2*'D_1','D_1s',2*'B_1','B_1s','B_1c','b_1','h_1','h''_1', 
     &'h_1c','h_1b','h_1t','h_1l','h_1h',2*' ','a_0',2*'K*_0',2*'D*_0', 
     &'D*_0s',2*'B*_0','B*_0s','B*_0c','a_0','f_0','f''_0','chi_0c', 
     &'chi_0b','chi_0t','chi_0l','chi_0h',2*' ','a_1',2*'K*_1', 
     &2*'D*_1','D*_1s',2*'B*_1','B*_1s','B*_1c','a_1','f_1','f''_1', 
     &'chi_1c','chi_1b','chi_1t','chi_1l','chi_1h',2*' ','a_2', 
     &2*'K*_2',2*'D*_2','D*_2s',2*'B*_2','B*_2s','B*_2c','a_2','f_2', 
     &'f''_2','chi_2c','chi_2b','chi_2t','chi_2l','chi_2h',2*' ','K_L', 
     &'K_S',8*' ','psi''',3*' ','Upsilon''',45*' ','pi_diffr'/ 
      DATA (CHAF(I)  ,I= 282, 500)/'n_diffr','p_diffr','rho_diff', 
     &'omega_di','phi_diff','J/psi_di',18*' ','Lambda',5*' ', 
     &'Lambda_c',' ',2*'Xi_c',6*' ','Lambda_b',' ',2*'Xi_b',6*' ','n', 
     &'p',' ',3*'Sigma',2*'Xi',' ',3*'Sigma_c',2*'Xi''_c','Omega_c', 
     &4*' ',3*'Sigma_b',2*'Xi''_b','Omega_b',4*' ',4*'Delta', 
     &3*'Sigma*',2*'Xi*','Omega',3*'Sigma*_c',2*'Xi*_c','Omega*_c', 
     &4*' ',3*'Sigma*_b',2*'Xi*_b','Omega*_b',114*' '/ 
 
C...LUDATR, with initial values for the random number generator. 
      DATA MRLU/19780503,0,0,97,33,0/ 
 
      END 
 
C********************************************************************* 
 
      SUBROUTINE LUTAUD(ITAU,IORIG,KFORIG,NDECAY) 
 
C...Dummy routine, to be replaced by user, to handle the decay of a 
C...polarized tau lepton. 
C...Input: 
C...ITAU is the position where the decaying tau is stored in /LUJETS/. 
C...IORIG is the position where the mother of the tau is stored; 
C...     is 0 when the mother is not stored. 
C...KFORIG is the flavour of the mother of the tau; 
C...     is 0 when the mother is not known. 
C...Note that IORIG=0 does not necessarily imply KFORIG=0; 
C...     e.g. in B hadron semileptonic decays the W  propagator 
C...     is not explicitly stored but the W code is still unambiguous. 
C...Output: 
C...NDECAY is the number of decay products in the current tau decay. 
C...These decay products should be added to the /LUJETS/ common block, 
C...in positions N+1 through N+NDECAY. For each product I you must 
C...give the flavour codes K(I,2) and the five-momenta P(I,1), P(I,2), 
C...P(I,3), P(I,4) and P(I,5). The rest will be stored automatically. 
 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUJETS/,/LUDAT1/ 
 
C...Stop program if this routine is ever called. 
C...You should not copy these lines to your own routine. 
      NDECAY=ITAU+IORIG+KFORIG      
      WRITE(MSTU(11),5000) 
      IF(RLU(0).LT.10.) STOP 
 
C...Format for error printout. 
 5000 FORMAT(1X,'Error: you did not link your LUTAUD routine ', 
     &'correctly.'/1X,'Dummy routine in JETSET file called instead.'/ 
     &1X,'Execution stopped!') 
 
 
      RETURN 
      END 
