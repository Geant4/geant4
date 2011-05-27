c      ###################################################################
c      03/2002 version INCL4.1 prete a l'usage
c      ###################################################################
c
C
C Version avec traitement de la surface Wood-Saxon et Consistent
C   Dynamical Pauli Principle (CDPP).                          (avec J. Cugnon)   
C        
C         version avec W.S., traitement Cugnon R(q),.............. A.B. 4/2000
C
C On installe le calcul pour les participants seulement...      (avec J. Cugnon) 
C                                                                  A.B. 12/2000
C
C Version avec deuton incident possible (C; Volant, J. Cugnon, A. Boudard) 4/2001
C              densite du deuton dans l'espace des impulsions selon pot de Paris.
C
C NOSURF=-2 for a standard INCL4 calculation (it means with W.S. nuclear density 
C                                and stopping time as 70*(A/208)**0.16 )    
c  This version has a unique random generator. One event can be reproduced
C  if ial and IYV(1...19) are printed during the test and then given
C  as input for initialisation of the common hazard.
C
C The pion emission is reduced: ND->NN (*3),  dynamical width of the Delta,
C        but the sigma(piN) IS NOT multiplied by 3.                 June 2001
C
C Corrections for composite projectiles (geometry and time matrix) june/2001
C
c   delta production: correction of the angular distribution 02/09/02 JC/CV
c
c   Modif du temps d'arret pour pi incident (et surface) et R0 des cibles 6<A<18
C								AB,CV,JC  mars/2003
C  11/2003:    Numerics (proposed by J. Hendricks) (TIME: TA.NE.0, energy/momentum
C              
C              Additional input F(9)->F(15)
C              Kveux: swipping of calculations for the avatar ntuple
C              Angular momentum carried by the remnant (L) from the distance
C			between remnant and actual target c.m., and recoil
C			momentum computed from outside (Pin - Sig(Pout))
    
        subroutine pnu(IBERT,F,NOPART,KIND,EP,ALPHA,BETA,GAM,IZREM,
     -           IAREM,ESREM,ERECREM,ALREM,BEREM,GAREM,BIMPACT,L)
C                                                                       P-N00010
C$$$$$$$$$$$$$$$$$-----PROTON-NUCLEUS-----$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ P-N00020
c
c      LIEGE INC-model as a subroutine 
c      The Liege INC model has been applied to several systems and
c      in several conditions. Refinements are still in progress
c
c      PLEASE refer to this version as INCL4.1 in order to avoid 
c      confusion when comparing to the results of your colleagues.
c
c      DIFFERENT from INCL2.0 in the sense that the cascade is stopped when
c      the excitation energy vanishes if this occurs before the predetermined
c      stopping time (herein denoted as temfin)
c      Special are taken to avoid emission of slow particles due to the 
c      imperfect Pauli blocking
c
c      PLEASE notice: There are basically only two parameters in the 
c      model: the average potential depth, denoted as V0, and the time 
c      at which the cascade is stopped denoted as temfin. In this program,
c      the "standard" values (those introduced in the ref NPA620(1997)475)
c      are V0=40MeV and temfin=1.25*some function of 
c      incident energy and impact parameter. You may, of course, change 
c      these parameters V0 and the numerical coefficient in temfin, within 
c      reasonable limits (i.e. V0 cannot be lower than 38MeV; V0=45MeV is 
c      recommended for heavy nuclei). If you do, PLEASE indicate your choice,
c      once again, for the same reason as above.
c
c      The description of the cascade model(incl2.0) can be found in:
c        J.C., C.VOLANT & S.VUILLIER, NPA620(1997)475
c      It is basically the same model as described in J.C. NPA462(1987)751
c      (version 7 (in our jargon), sketched below)
c      + a refinement of the parametrization of the cross-sections,
c      based on J.C., D. L'HOTE, J.VANDERMEULEN, NIM B111(1996)215 and
c      J.C., S.LERAY, E.MARTINEZ, Y.PATIN & S.VUILLIER PRC56(1998)2431
c
c      technical notes: 
c      1.for the parametrizations of cross  sections, see
c       notes of 4/10/96, 9/10/96, 31/12/97 and 13/10/98
c      2.temfin=1.25*... 18/6/98
c      3.sepa in concordance with v0-tf 2/7/98
c      4.special care for stopping the cascade before t=temfin 27/04/99
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                                       P-N00030
C     VERSION 7:  2/2/93                                                P-N00040
C
C++++++++++++ DESCRIPTION OF INPUT AND OUTPUT+++++++++++++++++++++++++++++
C
C **INPUT DATA**
C
C  IBERT=O IN THE FIRST CALL
C        1 IN THE SUBSEQUENT CALLS
C
C  F= (REAL) ARRAY OF DIMENSION 8
C
C    F(1)= A (TARGET)
C    F(2)= Z (TARGET)
C    F(3)= KINETIC ENERGY (IN MEV) OF THE INCIDENT PARTICLE
C    F(4)= SUPPOSED TO BE THE MINIMUM  PROTON ENERGY REQUIRED TO LEAVE
C          THE TARGET. IN THIS CASCADE MODEL, IT IS ZERO
C    F(5)= Nuclear potential V0 (standard value=45 MeV for heavy nuclei)
C    F(6)= Rescale the cascade duration (the standard value t0 is MULTIPLIED
c 	   by this value. F(6)=1. is the standard)
C    F(7)= TYPE OF INCIDENT PARTICLE
C           1.00 FOR PROTON
C           2.00 FOR NEUTRON
C           3.00 FOR PI+
C           4.00 FOR PI0
C           5.00 FOR PI-
C           6.00 FOR DEUTERON
C           7.00 FOR TRITON
C           8.00 FOR HE3
C           9.00 FOR HE4
C    F(8)= SUPPOSED TO BE THE MINIMUM NEUTRON ENERGY REQUIRED TO LEAVE
C          THE TARGET. IN THIS CASCADE MODEL, IT IS ZERO
C
C                                 NOSURF=1 Sharp density (hard sphere), 
C                                 NOSURF=0 Wood-Saxon density, stopping time "70" 
C                                                without B (impact) dependence.
C                                 NOSURF=-1 Wood-Saxon density, stopping time "70" 
C                                                with B (impact) dependence
C                                 (on peut toujours nenormaliser ces fonctions 
C                                  de temps avec le facteur F(6): t=t0*F(6) ) 
C                                XFOISA      Rmaxws = R0 + XFOISA*A
C					     Bmax = Rmaxws for pions and nucleons
C					     Bmax = Rmaxws + rms1t (data) for composits
C         Pauli strict (1) or statistic (0) or without pauli (2):     NPAULSTR
C
C    F(9)= imat, target material identifier for the right choose of Sax.-Wood density
c
c common/hazard/ial,IY1,IY2,IY3,IY4,IY5,IY6,IY7,IY8,IY9,IY10,
c     s               IY11,IY12,IY13,IY14,IY15,IY16,IY17,IY18,IY19
c ......20 numbers in the main routine to initialize the random numbers
c
C **OUTPUT DATA**
C
C  NOPART=-1 PSEUDO REACTION (VOID EVENT)
C          0 ABSORPTION
C         >0 TRUE EVENT, = NUMBER OF PARTICLES EMITTED (EXCLUDING THE REMNANT)
C
C  FOR N=1,NOPART:
C  KIND(N)= TYPE OF PARTICLES (SAME CONVENTION AS FOR F(7), BUT IN INTEGERS)
C           in this version extended to large composite projectiles, KIND
C           is negative for projectile spectators (not entering the
C           potential).
C
C  EP(N)=  KINETIC ENERGY
C
C  ALPHA(N),BETA(N),GAM(N)= DIRECTION COSINES
C
C  IZREM= Z (REMNANT)
C
C  IAREM= A (REMNANT)
C
C  ESREM= EXCITATION ENERGY OF THE REMNANT
C
C  ERECREM= RECOIL ENERGY OF THE REMNANT
C
C  ALREM,BEREM,GAREM=DIRECTION COSINES OF THE REMNANT
C
C  BIMPACT impact parameter
C
C  L Intrinsic momentum of the remnant in units of h/2pi=hbar=197.328
C+++++++++ DESCRIPTION OF THE INC MODEL ++++++++++++++++++++++++++++++++++++
C
C     MODEL DESCRIBED IN J.CUGNON (NP A462(1987)751)                    P-N00050
C     =MODEL (DR) OF J.CUGNON,D.KINET,J.VANDERMEULEN(NP A379(1982)567)  P-N00060
C                                                                       P-N00110
C        +REFLECTION OR TRANSMISSION ON THE POTENTIAL WALL              P-N00120
C              (THE POTENTIAL DEPTH IS THE SAME FOR NUCLEONS & DELTA'S) P-N00130
C              (CONTAINS A COULOMB BARRIER)                             P-N00140
C                                                                       P-N00150
C        +ABSORPTION OF THE PION ABOVE THE (3,3) RESONANCE (NOT IN      P-N00160
C                VERSION 2)                                             P-N00170
C                                                                       P-N00180
C        +POSSIBLE PAULI BLOCKING OF TWO BODY COLLISIONS                P-N00190
C        +POSSIBLE PAULI BLOCKING OF DELTA DECAY                        P-N00200
C                 THE PAULI BLOCKING IS APPLIED TO THE NUCLEONS         P-N00210
C                 ONLY.THE PAULI BLOCKING FACTORS ARE EVALUATED         P-N00220
C                 BY COUNTING THE NUCLEONS INSIDE A VOLUME IN           P-N00230
C                 PHASE SPACE.THE EXTENSION OF THIS VOLUME IS           P-N00240
C                 OF THE ORDER OF H**3                                  P-N00250
C                                                                       P-N00260
C    ADDITIONAL FEATURES:                                               P-N00270
C                                                                       P-N00280
C        +ISOSPIN (WITH NEW PN BACKWARD-FORWARD ASYMMETRY)              P-N00290
C                                                                       P-N00300
C        +"LONGITUDINAL GROWTH" OF THE BARYONS (NOT ACTIVATED HERE)
C
C        + PARTICLE #1 IS ALWAYS THE FASTEST PARTICLE IN THE Z-DIRECTIONP-N00100
C                                 (NOT ACTIVATED HERE)
C        +SIMPLIFIED NEUTRON EVAPORATION AT THE END OF THE CASCADE      P-N00310
C                                 (NOT PRESENT HERE)
C                 
C        +POSSIBLE CONSERVATION OF ANGULAR MOMENTUM (NOT ACTIVATED
C                    HERE, COPIED FROM P_NUCJ)
C        P-NU7=SAME AS P-NU6 + EVAPORATION                              P-N00330
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++P-N00340
 
      DIMENSION f(15),kind(300),ep(300),alpha(300),beta(300),gam(300)
      DIMENSION bmass(300)
      COMMON/hazard/ial,IY1,IY2,IY3,IY4,IY5,IY6,IY7,IY8,IY9,IY10,
     s               IY11,IY12,IY13,IY14,IY15,IY16,IY17,IY18,IY19
      COMMON/kind/kindf7
      DIMENSION IND(20000),JND(20000)                                   P-N00350
      DIMENSION INDIC(3000)                                             P-N00360
      DIMENSION NPAR(625,15),NIMP(600,15),NNCO(15),NIMPP(600)           P-N00370
      DIMENSION NENTR(10,8,20,10),NOUT1(15),NOUT2(15)                   P-N00380
 
c     DIMENSION TEM(15),NSR(40),NSP(40),NSR1(40),NSP1(40)               P-N00400
      DIMENSION TEM(15),NSR1(40),NSP1(40)                               P-N00400
      DIMENSION T(200),LINE(132),Q1(200),Q2(200),Q3(200),Q4(200),NC(300)P-N00410
      DIMENSION Y1(200),Y2(200),Y3(200),YM(200),IPI(200)                P-N00420
      DIMENSION NRNN(15),NRND(15),NRDD(15),NRDN(15),NRDP(15),NRPD(15),NCP-N00430
     -DD(15),NPAUL1(15),NPAUL2(15)                                      P-N00440
      DIMENSION NPDIR(600)                                              P-N00450
      DIMENSION NEJ(6,15),NRES(6,15),NPIA(6,15),NCHPRO(15),NDEL(15)     P-N00460
      DIMENSION EDEP1(15),EDEP2(15),EDEP3(15),EDEP4(15),NINT(15)        P-N00470
     -,EPAR1(15),EPAR2(15),EPAR3(15),EPAR4(15),ENPI(15),E1(15),EZ3(15)  P-N00480
      DIMENSION IHF1(50),IHF2(50),IHF3(50),IHF4(50),IHP(2,100),IHC(50), P-N00490
     -IHE(2,100),IHF5(100),IHREM(100,100)
     
      DIMENSION JPARTICIP(300),eps_c(12),p3_c(12)
     
C For clusters:
      REAL*4 pout(3), lout(3), pout_clst(3,100), lout_clst(3,100)
      DIMENSION id_nuc_clst(4),id_clst(4,100),k_clst(100),e_clst(100),
     s          t_clst(100)

C Dialogue with INCL: function R(q/pf) for each nucleus
      COMMON/SAXW/ XX(30,500),YY(30,500),SS(30,500),NBPINTER,IMAT
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX      

C RMS espace R, espace P, Fermi momentum and energy for light gauss nuc.      
      COMMON/light_gaus_nuc/rms1t(9),pf1t(9),pfln(9),tfln(9),vnuc(9)
      
C Fermi 2 param from A=19 to 28, modified harm oscil A=6 to 18
C (H. De Vries et al. At. Data and Nuc. Data Tab. 36 (1987) 495)
      COMMON/light_nuc/R_light_nuc(30),a_light_nuc(30)      

C Pour version Be (Parametres du projectile)
      COMMON/V_BE/ia_be,iz_be,rms_be,pms_be,bind_be

C Common for study of AVATARS through an NTUPLE optionally produced
      REAL*4 Bavat,TIMEavat,ENERGYavat
      INTEGER Bloc_Paul,Bloc_CDPP,GO_OUT,avm,DEL1avat,DEL2avat
      PARAMETER (avm=1000)
      COMMON/VAR_AVAT/Kveux,Bavat,NOPARTavat,NCOLavat,
     s R1_in(3),R1_first_avat(3),
     s EPSd(250),EPS2(250),EPS4(250),EPS6(250),EPSf(250),
     s NB_AVAT,
     s TIMEavat(avm),L1avat(avm),L2avat(avm),JPARTL1(avm),JPARTL2(avm),
     s DEL1avat(avm),DEL2avat(avm),ENERGYavat(avm),
     s Bloc_Paul(avm),Bloc_CDPP(avm),GO_OUT(avm)

      DIMENSION npproj(300)
      COMMON/SPL2/ XSP(100),YSP(100),AD(100),BD(100),CD(100),NDEUT
c deutons      
c     DIMENSION LTT(15)                                                 P-N00520
      COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TAP-N00530
      COMMON/BL2/CROIS(19900),K,IND,JND                                 P-N00540
      COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N00550
      COMMON/BL4/TMAX5                                                  P-N00560
      COMMON/BL5/TLG(300),NESC(300)                                     P-N00570
      COMMON/BL6/XX10,ISA                                               P-N00580
      COMMON/BL8/RATHR,RAMASS
      COMMON/bl9/hel(300),l1,l2
      COMMON/bl10/ri4,rs4,r2i,r2s,PF
      COMMON/PAUL/CT0,CT1,CT2,CT3,CT4,CT5,CT6,PR,PR2,XRR,XRR2,
     s            CP0,CP1,CP2,CP3,CP4,CP5,CP6
     
      DIMENSION ia1t(9),iz1t(9),fmpinct(9)      

      DATA LINE/132*1H*/,HC,FMP,FMD,FMPI/197.328,938.2796,1232.,138.00/ P-N00590
      DATA ia1t,iz1t,fmpinct/1,1,0,0,0,2,3,3,4,1,0,1,0,-1,1,1,2,2,
     -938.2796,938.2796,138.0,138.0,138.0,1874.35,2806.8,2806.8,3727./
C      DATA rms1t,pf1t/0.,0.,0.,0.,0.,2.10,1.80,1.80,1.63,
C     -0.,0.,0.,0.,0.,77.,110.,110.,153./


c deutons
      DATA NENTR/16000*0/                                               P-N00620
C  Parameters of the WS density for the target nucleus
C      DATA R0,ADIF,RMAXWS,DRWS/6.624,0.549,13.5,0.5/       !208 Pb
C      DATA R0,ADIF,RMAXWS,DRWS/4.11,0.558,10.,0.25/        !56 Fe  
C      DATA R0,ADIF,RMAXWS,DRWS/6.38,0.535,11.,0.25/        !197 Au  
      
      EXTERNAL TIME                
      external mglms

C                                                                       P-N00630
CCC   RELATIVISTIC ENERGY AND MASS FUNCTIONS                            P-N00640
C                                                                       P-N00650
      W(A,B,C,D)=SQRT(A*A+B*B+C*C+D*D)                                  P-N00660
      AM(A,B,C,D)=SQRT(D*D-A*A-B*B-C*C)                                 P-N00670
C                                                                       P-N00680
CCC   FORMATS                                                           P-N00690
C                                                                       P-N00700
   10 FORMAT (A60)                                                      P-N00710
   11 FORMAT (3I3,2F10.3)                                               P-N00720
   12 FORMAT (6F10.3)                                                   P-N00730
   13 FORMAT (6I3)                                                      P-N00740
   14 FORMAT (10I4)                                                     P-N00750
  100 FORMAT (50H B LARGER THAN R1+R2 , NO COLLISION               )    P-N00760
  101 FORMAT (38H INITIAL VALUE OF THE RANDOM VARIABLE=,I11/)           P-N00770
  102 FORMAT (19H  K1 K2 K3 K4 K5 K6/1H ,6I3///)                        P-N00780
  103 FORMAT (1H , 'COMPOSITION OF THE EJECTILES'/                      P-N00790
     -'    N      P     D-    D0    D+    D++'/6F7.2/)                  P-N00800
  104 FORMAT (1H ,15H INCIDENT MASS=,I3,15H INCID. CHARGE=,I3,/18H TARGEP-N00810
     -T MASS=     ,I3,15H TARGET CHARGE=,I3,                 /19H   INCI
     -DENT ENERGY=,F8.2,6H MEV  ,13H   BETA(INC)=,F6.4,14H   GAMMA(INC)=P-N00820
     -,F8.4//)                                                          P-N00830
  105 FORMAT (18H IMPACT PARAMETER=,F6.2,' tirage b',/)                 P-N00840
  106 FORMAT (9H BETA(1)=,F6.4,10H  BETA(2)=,F6.4,11H  GAMMA(1)=,F8.4,1 P-N00850
     -1H  GAMMA(2)=,F8.4//)                                             P-N00860
  107 FORMAT (1H ,132A1)                                                P-N00870
  108 FORMAT (50H  PARAMETERS OF THE GENERAL FRAME                 /)   P-N00880
  109 FORMAT (54H NO MORE COLLISIONS , SYSTEM EXPANDING FOR EVER       ,P-N00890
     -10X,22H NUMBER OF COLLISIONS=,I4,10X, 6H TIME=,F8.4)              P-N00900
  110 FORMAT (22H NUMBER OF COLLISIONS=,I5,17H NUMBER OF PIONS=,I5,21H EP-N00910
     -NERGY OF THE PIONS=,E11.4,4H MEV,28H ENERGY OF THE DIRECT PIONS=,EP-N00920
     -11.4,4H MEV)                                                      P-N00930
  111 FORMAT (19H INCIDENT MOMENTUM=,F8.2,6H MEV/C/)                    P-N00940
  112 FORMAT (55H FINAL LEADING PARTICLE KIN. ENERGY DISTRIBUTION, STEP=P-N00950
     -,F7.2,4H MEV/)                                                    P-N00960
  113 FORMAT (' FINAL LEADING  PARTICLE MOMENTUM DISTRIBUTION, STEP IN  P-N00970
     - X,Y,Z=',3(F6.2,2X),6H MEV/C/'  DISTRIBUTIONS IN X AND Y ARE      P-N00980
     - CENTERED AT 0',/'  DISTRIBUTION IN Z STARTS AT 0',/)             P-N00990
  114 FORMAT (1H1)                                                      P-N01000
  115 FORMAT (22H NUMBER OF COLLISIONS=,I4,/47H NUMBER OF COLLISIONS UNDP-N01010
     -ERGONE BY EACH NUCLEON)                                           P-N01020
  116 FORMAT (40I3)                                                     P-N01030
  117 FORMAT (16H   (2) EJECTILES///)                                   P-N01040
  118 FORMAT (54H DISTRIBUTION OF THE # OF PRIMARY COLLISIONS(0,1,2,.../P-N01050
     -20I4)                                                             P-N01060
  119 FORMAT (39H ENERGY SPECTRUM OF THE REMNANTS ,STEP=,F8.2,4H MEV/)  P-N01070
  120 FORMAT (26H VALUES OF THE MESH POINTS//6H XMIN=,F8.2,8H   XMAX=,F8P-N01080
     -.2,6H   DX=,F8.2/F10.2)                                           P-N01090
  121 FORMAT (1H ,F10.2,1X,40I3)                                        P-N01100
  122 FORMAT (/71H THE ACTUAL DENSITY IS OBTAINED BY DIVIDING BY NRUN ANP-N01110
     -D BY THE VOLUME V)                                                P-N01120
  123 FORMAT (81H PLOT OF THE NUCLEON MOMENTUM DISTRIBUTION IN THE RAPIDP-N01130
     -ITY-P(PERPENDICULAR) PLANE)                                       P-N01140
  124 FORMAT (26H VALUES OF THE MESH POINTS//6H  PP1=,F8.2,6H  PP2=,F8.2P-N01150
     -,6H  DPE=,F8.2/F10.2)                                             P-N01160
  125 FORMAT (//17H NUMBER OF DELTA=,F6.2/22H NUMBER OF COLLISIONS=,F6.2P-N01170
     -/17H NUMBER OF PIONS=,F6.2/15H NUMBER OF NES=,F6.2/16H NUMBER OF PP-N01180
     -ROD=,F6.2/18H NUMBER OF N-D ES=,F6.2/15H NUMBER OF REC=,F6.2/15H NP-N01190
     -UMBER OF DEC=,F6.2/15H NUMBER OF ABS=,F6.2/15H NUMBER OF DES=,F6.2P-N01200
     -/29H # OF PAULI BLOCKED BB COLL.=,F6.2/34H # OF PAULI BLOCKED DELTP-N01210
     -A DECAYS= ,F6.2/)                                                 P-N01220
  126 FORMAT (81H PLOT OF THE    PION MOMENTUM DISTRIBUTION IN THE RAPIDP-N01230
     -ITY-P(PERPENDICULAR) PLANE)                                       P-N01240
  127 FORMAT (1H ,I5,28H PARTICLES OUTSIDE THE FRAME/)                  P-N01250
  128 FORMAT (1H1,55X,6H TIME=,F8.2,5H FM/C///)                         P-N01260
  129 FORMAT (1H )                                                      P-N01270
  130 FORMAT (18H RAPIDITY SPECTRUM/)                                   P-N01280
  131 FORMAT (1H ,20I4)                                                 P-N01290
  132 FORMAT (/17H P(PERP) SPECTRUM/)                                   P-N01300
  133 FORMAT (11H ASSYMETRY=,F8.4)                                      P-N01310
  134 FORMAT(17H NON DIRECT PIONS)                                      P-N01320
  135 FORMAT(13H DIRECT PIONS)                                          P-N01330
  136 FORMAT (' PROB.OF HAVING THE LEADING PARTICLE (LARGEST ENERGY) DIFP-N01340
     -F.FROM THE INCOMING "PROTON"(PART.#1)=',F8.2/)                    P-N01350
  137 FORMAT (34H AVERAGE # OF PRIMARY COLLISIONS= ,F8.2/61H PRIMARY COLP-N01360
     -LISIONS = COLLISIONS MADE BY THE LEADING PARTICLE/)               P-N01370
  138 FORMAT (42H FINAL LEADING PARTICLE MOMENTUM(X,Y,Z,T)=,4(2X,E11.4))P-N01380
  139 FORMAT (48H FINAL "PROTON" ENERGY(PART.#1 , LEADING PART.)=,2(E11.P-N01390
     -4,2X))                                                            P-N01400
  140 FORMAT (21H MASS OF THE REMNANT=,F8.4/34H MOMENTUM OF THE REMNANT P-N01410
     - (X,Y,Z)=,3(E11.4,2X),6H MEV/C/39H ENERGY (KINETIC+MASS) OF THE REP-N01420
     -MNANT= ,E11.4,4H MEV//7H   E*= ,E11.4,4H MEV,4X,E11.4,4H MEV/)    P-N01430
  141 FORMAT (30X,19H  FINAL  QUANTITIES////                            P-N01440
     -21H  (1) INCIDENT PROTON//)                                       P-N01450
  142 FORMAT (27H TOTAL ENERGY OF THE PIONS=,F8.2,4H MEV/)              P-N01460
  143 FORMAT (23H # OF EJECTED NUCLEONS=,F8.2,                          P-N01470
     -43H ENERGY-MOMENTUM OF THE EJECTILES(X,Y,Z,E)=,4(E11.4,3X),/      P-N01480
     -16H KINETIC ENERGY=,E11.4/)                                       P-N01490
  144 FORMAT (1H ,F10.2,1X,20I4)                                        P-N01500
  145 FORMAT (39H TOTAL ENERGY OF THE "PROTON"(PART.#1)=,F10.2//        P-N01510
     -18H # OF PARTICIPANTS,F8.2/)                                      P-N01520
  146 FORMAT (1H ,F10.3,15H NUCLEONS    S=,F10.3/)                      P-N01530
  147 FORMAT (13H S(NUCLEONS)=,F10.3,3X,11H S(DELTAS)=,F10.3/)          P-N01540
  148 FORMAT (1H1,25H # OF TRANSPARENT EVENTS=,I4//52H THE AVERAGES ARE P-N01550
     -CALCULATED ON THE REMAINING EVENTS/)                              P-N01560
  149 FORMAT (1H , 'COMPOSITION OF THE REMNANT',/                       P-N01570
     -'    N      P     D-    D0    D+    D++'/6F7.2/)                  P-N01580
  150 FORMAT (1H ,    'PROB. PARTICLE #1 KEEPS ITS CHARGE=',F6.2)       P-N01590
  151 FORMAT (1H1,16H   (3) REMNANTS ///)                               P-N01600
  152 FORMAT (1H ,27H  <E*> (AFTER DELTA DECAY)=,F8.2,4H MEV/)          P-N01610
  153 FORMAT (1H ,' EJECTILE KINETIC ENERGY SPECTRUM, STEP= ',F8.2,     P-N01620
     -  4H MEV/)                                                        P-N01630
  154 FORMAT (16H   NEUTRONS     ///44H DISTRIBUTION OF THE # OF EJECTILP-N01640
     -ES(0,1,2..)//)                                                    P-N01650
  155 FORMAT (16H   PROTONS      ///44H DISTRIBUTION OF THE # OF EJECTILP-N01660
     -ES(0,1,2..)//)                                                    P-N01670
  156 FORMAT (1H ,'DISTRIBUTION OF THE REMNANTS IN THE CHARGE LOSS(H)-  P-N01680
     - NEUTRON LOSS(V) PLANE (0,1,2,...)'/'ATTENTION : ZERO OF THE SCALEHE-01950
     -S CORRESPONDS TO THE NINTH ENTRY'//)
C                                                                       P-N01700
CCC   READING OF THE DATA                                               P-N01710
C                                                                       P-N01720
C-------EXPLANATION OF SOME PHYSICAL QUANTITIES------------------------ P-N01730
C           (BASICALLY INPUT DATA, WHEN RUN AS A PROGRAM)
C                                                                       P-N01740
C     IA1=MASS NUMBER OF THE INCIDENT ION                               P-N01750
C                                                                       P-N01760
C     IA2=MASS NUMBER OF THE TARGET                                     P-N01770
C     IZ1=ATOMIC NUMBER OF THE PROJECTILE                               P-N01780
C                                                                       P-N01790
C     IZ2=ATOMIC NUMBER OF THE TARGET                                   P-N01800
C     R01=RADIUS PARAMETER OF THE PROJECTILE                            P-N01810
C           R01 SHOULD BE PUT TO 1.000                                  P-N01820
C     R02=RADIUS PARAMETER OF THE TARGET                                P-N01830
C     ADIF1=DIFFUSENESS OF THE PROJECTILE                               P-N01840
C           ADIF1 SHOULD BE PUT TO 1.0000                               P-N01850
C     ADIF2=DIFFUSENESS OF THE TARGET                                   P-N01860
C                                                                       P-N01870
C     TLAB = INCIDENT ENERGY  (IN MEV)                                  P-N01880
C                                                                       P-N01890
C     K1=1 REFERENCE FRAME=C.M. OF THE COVERING MATTER                  P-N01900
C     K1=2 REFERENCE FRAME=C.M. FRAME OF THE INCIDENT ION AND ITS       P-N01910
C          INTERCEPT                                                    P-N01920
C     K1=3 REFERENCE FRAME=C.M. FRAME OF A N-N SYSTEM WITH THE SAME     P-N01930
C          KINEMATICS AS THE TWO IONS                                   P-N01940
C     K1=4 REFERENCE FRAME=C.M. FRAME FOR THE TOTAL  SYSTEM             P-N01950
C     K1=5 REFERENCE FRAME=LAB SYSTEM                                   P-N01960
C              K1 SHOULD BE PUT TO 5 IN THIS VERSION                    P-N01970
C     K2=0 RELATIVISTIC KINEMATICS                                      P-N01980
C     K2=1 NON-RELATIVISTIC KINEMATICS(ABANDONNED IN THIS VERSION)      P-N01990
C     K3=0  DELTAS ARE PRODUCED                                         P-N02000
C     K3=1  NO DELTA PRODUCTION                                         P-N02010
C     K4=0 THE DELTA IS GIVEN A VANISHING LIFETIME                      P-N02020
C     K4=1 THE DELTA HAS A VERY LARGE LIFETIME                          P-N02030
C     K4=2 THE DELTA HAS A EXPONENTIALLY RANDOM LIFETIME                P-N02040
C     K5=0 NO DELTA-NUCLEON,DELTA-DELTA INTERACTIONS                    P-N02050
C     K5=1 DELTA-NUCLEON=DELTA-DELTA=NUCLEON-NUCLEON ELASTIC X-SECTION  P-N02060
C     K6=0 NO ANGULAR MOMENTUM CONSERVATION                             P-N02070
C     K6=1 ANGULAR MOMENTUM CONSERVATION                                P-N02080
C                                                                       P-N02090
C     B=IMPACT PARAMETER                                                P-N02100
C                                                                       P-N02110
C     RBL(PBL) IS THE RADIUS IN REAL(MOMENTUM) SPACE OF THE VOLUME      P-N02120
C         ON WHICH THE NUCLEONS ARE COUNTED TO EVALUATE THE PAULI       P-N02130
C         BLOCKING FACTORS                                              P-N02140
C     RECOMMENDED VALUES ARE 2 FM AND 200 MEV/C RESPECTIVELY            P-N02150
C                                                                       P-N02160
C     NRUN= NUMBER OF  RUNS (IRRELEVANT HERE)                           P-N02170
C                                                                       P-N02180
C     NTM=NUMBER OF INTERMEDIATE TIMES AT WHICH THE SPATIAL AND MOMENTUMP-N02190
C         DISTRIBUTIONS ARE STORED (NOT RELEVANT HERE)                  P-N02200
C                                                                       P-N02210
C     TEM(IT)=VALUES OF THE INTERMEDIATE TIMES                          P-N02220
C                                                                       P-N02230
C     V0=DEPTH OF THE NUCLEON POTENTIAL                                 P-N02240
C     V1=DEPTH OF THE DELTA POTENTIAL                                   P-N02250
C         IN THIS VERSION V1 SHOULD BE EQUAL TO V0                      P-N02260
C         (ATTENTION! A  V0 > 0 CORRESPONDS TO AN ATTRACTIVE POTENTIAL) P-N02270
C                                                                       P-N02280
C     NCASE=NUMBER OF CELLS NECESSARY TO DEPICT THE SPATIAL DISTRIBUTIONP-N02290
C         IN THE X AND Z DIRECTIONS                                     P-N02300
C                                                                       P-N02310
C     XP,XG,ZP,ZG=LIMITS OF THE BOX IN THE X AND Z DIRECTIONS           P-N02320
C     IN ORDER TO HAVE THE SAME SCALE IN THE TWO DIRECTIONS ,PUT XG-XP= P-N02330
C         0.9*(ZG-ZP)                                                   P-N02340
C     DY=DIMENSION OF THE BOX IN THE Y DIRECTION                        P-N02350
C                                                                       P-N02360
C     RAP1,RAP2,PP1,PP2=LIMITS OF THE BOX IN THE RAPIDITY-P(PERPENDICU- P-N02370
C          -LAR) SPACE                                                  P-N02380
C     THE BOX IS DIVIDED INTO 15 CELLS ALONG P(PERP) AND 40 CELLS ALONG P-N02390
C         THE RAPIDITY                                                  P-N02400
C                                                                       P-N02410
C     EN,SN=AVERAGE ENERGY CARRIED AWAY BY A NUCLEON,SEPARATION ENERGY  P-N02420
C                                                                       P-N02430
C       (NTM,NCASE,XP,...,SN are not used here)
C-----------------------------------------------------------------------P-N02440
C 	if(ibert.eq.0)then
c      write(6,*)'*************** VERSION INCL 4.1************************'
c      write(6,*)'* stopping time and potential can be changed     *******'
c      write(6,*)'***********      input of first random numbers   *******'
c      write(6,*)'* bimpact is output instead of sepa as in INCL 3.0 *****'
c      write(6,*)'* implementation of surface W.S.          4/2000  ******'
c      write(6,*)'* interaction only with "participants"    4/2000  ******'
c      write(6,*)'* CDPP:Coherent dynamical Pauli Principle 5/2001  ******'
c      write(6,*)'* Paris momentum density for the deuteron 4/2001  ******'
c      write(6,*)'* ND-NN (*3) and phase space Delta width  4/2001  ******'
c      write(6,*)'************* INCL 4.0 -> INCL4.1 ********************'
c      write(6,*)'* No lower cut on pi-N interaction	   2/2002  ****'       
c      write(6,*)'* Init of first avatars for participants  2/2002  ****'     
c      write(6,*)'* Output of excit energy for absorption   2/2002  ****'
c      write(6,*)'********************************************************'
c      write(6,*)' '
c	endif


   91 continue
c deutons
      beproj=0.
      ia2=f(1)
      iz2=f(2)
      r02=1.12
      kindf7=int(f(7)+0.1)

      if(kindf7.lt.0) kindf7=kindf7-1

      if(kindf7.eq.-12) then
         ia_be = 12
         iz_be = 6
C         call mglms(ia_be, iz_be, 3, bind_be)
         bind_be = 92.16
         rms_be = 2.36
         pms_be = 100.0
      endif

      if(kindf7.gt.0) then
         ia1=ia1t(kindf7)
         iz1=iz1t(kindf7)
         fmpinc=fmpinct(kindf7)
         rms1=rms1t(kindf7)
         pf1=pf1t(kindf7)
      else
C         ia1=abs(kindf7)
C         ia_be=ia1
         ia1 = ia_be
         iz1=iz_be
C         write(6,*)'ia1 = ',ia1
C         write(6,*)'iz1 = ',iz1
         fmpinc=ia1*fmp - bind_be
         rms1=rms_be
         pf1=pms_be
      endif
      if(idebug.eq.1) then
          write(6,*) 'rms1 = ',rms1
          write(6,*) 'pf1 = ',pf1
      endif
c      if (kindf7.le.2) then
c        ia1=1
c        iz1=2-kindf7
c        fmpinc=fmp
c       else
c        ia1=0
c        iz1=4-kindf7
c        fmpinc=fmpi
c      endif
c deutons
      TLAB=F(3)
c     READ 13,K1,K2,K3,K4,K5,K6                                         P-N02520
      K1=5                                                              P-N02530
      K2=0                                                              P-N02540
      k3=0
      k4=2
      k5=1
      k6=0
c material number:      
      IMAT=F(9)+0.5
c espace de phases test (R et P) pour PAULI: 
C Valeur recommandee par J.C. V-test=0.589 h**3:
      rbl=2.
      pbl=200.
           
C Valeur pour avoir V-test=h**3 avec 200 MeV/c:
C      rbl=2.386
C      pbl=200.
      
C Valeur pour avoir V-test=h**3 avec 2fm:
C      rbl=2.0
C      pbl=238.6

C Valeur pour avoir V-test=2.38 h**3 (avec pbl=200)
      rbl=3.1848

C Cette valeur donne vol =0.5 dans pauli: f=0 ou 1, pas de seconde chance.
C      rbl = 1.8935
C      pbl=200.
   
      XRR=RBL
      XRR2=XRR*XRR
      PR=PBL
      PR2=PR*PR
      
c     READ 12,RBL,PBL                                                   P-N02560
c surface 
C      ADIF=0.455		
c      adif=0.0
c fin surface 
c     READ 13,NTM                                                       P-N02590
c     READ 12,(TEM(IT),IT=1,NTM)                                        P-N02600
      ntm=1
      tem(1)=100000.
c   temfin (time at which the INC is stopped), tmax5 defined after chosing b
 
c   on conserve pour memoire
      nrun=1
      ncase=25
      xp=10.
      xg=-10.
      zp=-12.5
      zg=12.5
      dy=2.
      rap1=-3.
      rap2=3.
      pp1=0.
      pp2=1200.
      n2=ncase*ncase
      en=20.
      sn=8.
c     READ 12,V0,V1                                                     P-N02610
c a way to change v0=f(5) f(5) not used here
      V0=f(5)
      V1=V0                                                             P-N02620
      RATHR=0.
      RAMASS=0.

C                                                                       P-N02780
CCC   CONSTANTS AND DERIVED DATA                                        P-N02790
C                                                                       P-N02800

      PF=1.37*HC                                                        P-N02820
      TF=SQRT(PF*PF+FMP*FMP)-FMP
      G0=115.                                                           P-N02830
      SQ0=4000.                                                         P-N02840
      TH=0.                                                             P-N02850
      PM2=FMP*FMP                                                       P-N02860
      IA=IA1+IA2                                                        P-N02870
      ITM=NTM+1                                                         P-N02880
      A1=IA1                                                            P-N02890
      A2=IA2                                                            P-N02900
      R2=R02*A2**0.33333333                                             P-N02910
C parametres moyens de densite de la cible (fermi 2 parametres)
C      R0 = (2.745E-4*A2+1.063)*A2**0.33333333
C		IF (NOSURF.LE.0) THEN	!Flag sans surface (=1)
C      ADIF = 1.63E-4*A2+0.510
C      		ELSE
C      ADIF = 0.
C      		END IF
C      	IF(f(7).EQ.2.) PBEAM=SQRT(f(3)*(f(3)+1879.13126))   !neutron
C	IF(f(7).EQ.6.) PBEAM=SQRT(f(3)*(f(3)+2.*1874.35))   !deuton
C      RMAXWS = R0 + XFOISA*ADIF
C parametres moyens de densite de la cible (fermi 2 parametres)
      IF (IA2.GE.28) THEN
      	R0 = (2.745E-4*IA2+1.063)*IA2**0.33333333
      	ADIF = 1.63E-4*IA2+0.510
        RMAXWS = R0 + XFOISA*ADIF
      ELSE IF(IA2.GE.19) THEN
      	R0 = R_light_nuc(IA2)
      	ADIF = a_light_nuc(IA2)
        RMAXWS = R0 + XFOISA*ADIF
      ELSE IF(IA2.GE.6) THEN
      	R0 = 1.581*a_light_nuc(IA2)
     s            *(2.+5.*R_light_nuc(IA2))/(2.+3.*R_light_nuc(IA2))
      	ADIF = a_light_nuc(IA2)
	RMAXWS = 5.5 + 0.3*(IA2-6.)/12.
      ELSE IF(IA2.GE.2) THEN
      	IF(IA2.EQ.2) THEN
		R0=rms1t(6)
		PF=pfln(6)
	        TF=tfln(6)
		V0=vnuc(6)
	ENDIF
      	IF(IA2.EQ.3.AND.IZ2.EQ.1) THEN
		R0=rms1t(7)
		PF=pfln(7)
	        TF=tfln(7)
		V0=vnuc(7)
	ENDIF
      	IF(IA2.EQ.3.AND.IZ2.EQ.2) THEN
		R0=rms1t(8)
		PF=pf1t(8)
	        TF=tfln(8)
		V0=vnuc(8)
	ENDIF
      	IF(IA2.EQ.4) THEN
		R0=rms1t(9)
		PF=pf1t(9)
	        TF=tfln(9)
		V0=vnuc(9)
	ENDIF
	V1=V0
	ADIF=0.57735*R0
	RMAXWS=R0+2.5
      END IF
      IF(NOSURF.GT.0) THEN
      	ADIF=0.
	RMAXWS=R0
      END IF
      DRWS = RMAXWS/29.
C on impose le rayon du noyau:********************
C   ....voir coherence FUNCTION WSAX(R) et DERIVWSAX(R)
 
		R2=R0
                
                
C*************************************************
      
      BEV=IZ2*1.44/(1.2*R2)                                             P-N02920
      TNR=TLAB                                                          P-N02930
      BINC=SQRT(TNR*TNR+2.*TLAB*FMPINC)/(TNR+FMPINC)                    P-N02940
      GINC=1./SQRT(1.-BINC*BINC)                                        P-N02950
      PINC=FMPINC*BINC*GINC                                             P-N02970
      DPP=PF/6.                                                         P-N02980
      DHE=25.                                                           P-N02990
      DH1=50.                                                           P-N03000
      DH2=50.                                                           P-N03010
      DH3=PINC/50.                                                      P-N03020
      DH4=(FMP+TLAB)/50.                                                P-N03030
      DH5=TLAB/50.                                                      P-N03040
      
      DO I=1,IA
            JPARTICIP(I)=0
      END DO 
c
c     skip initialisations
c
      NEX=0                                                             P-N03700
      PLOSX=0.                                                          P-N03710
      PLOSY=0.                                                          P-N03720
      PLOSZ=0.                                                          P-N03730
      PLOST=0.                                                          P-N03740
      NPRIM=0                                                           P-N03750
      EFIN1=0.                                                          P-N03760
      EFIN=0.                                                           P-N03770
      EPIONS=0.                                                         P-N03780
      EFRUN=0.                                                          P-N03790
      EREMS=0.                                                          P-N03800
      YPREM1=0.                                                         HE-04200
      YPREM2=0.                                                         HE-04210
      YPREM3=0.                                                         HE-04220
      NTRAN=0                                                           P-N03810
      AVERL=0.0                                                         P-N03820
      NAVL1=0                                                           P-N03830
      NAVL2=0                                                           P-N03840
      NAVL3=0                                                           P-N03850
          Iavat=0
C                                                                       P-N03860
CCC   GENERATION OF THE INITIAL DISTRIBUTION IN THE REST FRAME OF THE IOP-N03870
C                                                                       P-N03880
c surface
C      bmax=R2+2.2*adif
C      r2i=(r2-2.2*adif)
C      r2s=(r2+2.2*adif)
C **********
      
      bmax=rmaxws		! maximum extension of the nucleus ( W.S.)

C*********************************************************************           
c fin surface (que faire aux pions ?)
C      if (kindf7.gt.2) bmax=r2+2.2
C      if (kindf7.gt.2) bmax=bmax	! A.B. (avec W.S., idem les nucleons)
c deutons CV 22/01/2001
      if (kindf7.lt.6.and.kindf7.gt.0)  then
       bmax=bmax     ! comme Alain
      else
         beproj=fmpinc-ia1*fmp
         bmax=rmaxws+rms1
      endif   
c deutons     
       IRUN=0                                                            P-N03890
      call ribm(al,ial)
      B=sqrt(al)*bmax
      bred=b/r2
      BIMPACT=B

      IF(NOSURF.NE.-2) THEN    ! La suite, c'est la version temps AVANT 2001

           IF(NOSURF.LE.0) THEN
C**********************************      
      tnor=70./30.6349         
C********************************** 
           ELSE
      tnor=1.
           ENDIF               
      IF(NOSURF.EQ.0) bred = 0.
                                            
c       write (6,*) bred
       if (kindf7.le.2.and.kindf7.gt.0) then
        if (tlab.lt.400.) then
          cb0=6.86-0.0035*tlab
          eb0=0.32-0.00005*tlab
        else
          if (tlab.lt.1000.) then
            cb0=5.23+0.000575*tlab
            eb0=0.32-0.00005*tlab
          else
            cb0=5.73+0.00007*tlab
            eb0=0.283-0.000013*tlab
          endif
        endif
        temfin=1.25*cb0/amax1(1.,0.854+.438*bred)
     -     *ia2**(eb0/amax1(1.,.941+.177*bred))
        temfin = temfin*tnor
C**********************************            
      else
       if (kindf7.lt.6.and.kindf7.gt.0)  then      
c here for pions:
       temfin=30.*(float(IA2)/208.)**0.25*(1.-0.2*bred)*(1.-tlab/1250.)
C correction for pions in the case NOSURF=0 or -1 (A.B., C.V. 2/2002)
       temfin=temfin*tnor
       else
c deutons
         tlabu=tlab/ia1
         if (tlabu.le.400) then    
          coeffb0=-0.0035*tlabu+6.86
          expob0=-0.00005*tlabu+0.32
         else
          if (tlabu.le.1000) then    
             coeffb0=0.000575*tlabu+5.23
             expob0=-0.00005*tlabu+0.32
          else
c             if (tlab.le.2000) then    
              coeffb0=0.00007*tlabu+5.73
              expob0=-0.000013*tlabu+0.283
          endif
         endif
         if (bred.le.0.33333) then
          xc=1.
          xe=1.
         else
          xc=0.438*bred+0.854
          xe=0.177*bred+0.941
         endif
         temfin=1.25*(coeffb0/xc)*ia2**(expob0/xe)        
C**********************************
C Same renormalisation of time for p,n and composit particles.
         temfin = temfin*tnor
       endif                 
      endif
      
      ELSE        ! Ici,NOSURF=-2 c'est la fonction temps 2001 (avec surface).
	   IF(kindf7.GE.3.AND.kindf7.LE.5) THEN
c here for pions (arbitrary multiplied by 2 for surface) (A.B., C.V.) 2/2002:
C       temfin=60.*(float(IA2)/208.)**0.25*(1.-0.2*bred)*(1.-tlab/1250.)
C modified in april 2003 (more reasonable but not yet checked!)
        temfin=25.5*ia2**0.16  ! Pb208->60fm/c
       	   ELSE
c here for other hadrons
         temfin = 29.8*ia2**0.16  ! Pb208->70fm/c
	   ENDIF
      ENDIF
c   
c       write (6,*) 'temfin,bred,b',temfin,bred,b,tnor
c deutons
C**********************************            
c a way to change stopping time f(6) not used here
	factemp=f(6)
C ATTENTION !!! 30/04/2001 scaling time is now a multiplication factor
     	temfin=temfin*factemp
	
      exi=0.
      nbquit=0
      iqe=0
      idecf=0
      TMAX5=temfin+0.1
      NPION=0
      EFER=0.                                                           P-N03920
c deutons
c      IF (IA1.GT.1) GO TO 77                                            P-N03930
      IF (IA1.GT.1) GO TO 7    
c deutons
      IF (IA1.EQ.1) THEN
        NESC(1)=0                                                       P-N03940
        IND2(1)=2*IZ1-1                                                 P-N03950
        IND1(1)=0                                                       P-N03960
        X1(1)=0.                                                        P-N03970
        X2(1)=0.                                                        P-N03980
        X3(1)=0.                                                        P-N03990
        P1(1)=0.                                                        P-N04000
        P2(1)=0.                                                        P-N04010
        P3(1)=0.                                                        P-N04020
        hel(1)=0.
        JPARTICIP(1)=1
       ELSE
        NPION=1
        IPI(1)=8-2*KINDF7
        Y1(1)=0.
        Y2(1)=0.
        Y3(1)=0.
        Q1(1)=0.
        Q2(1)=0.
        Q3(1)=0.
        Q4(1)=FMPI
      ENDIF
c deutons
      go to 9

 7    continue
      if(kindf7.eq.6) then
C     Deuteron in:
         call rgauss(xga)
         x1(1) = xga*rms1*0.57735
         call rgauss(xga)
         x2(1) = xga*rms1*0.57735
         call rgauss(xga)
         x3(1) = xga*rms1*0.57735
         
         x1(2) = -x1(1)
         x2(2) = -x2(1)
         x3(2) = -x3(1)
C     Deuteron density from Paris potential in q space:
         call ribm(xq,iy10)
         qdeut=splineab(xq)*197.3289
         call ribm(u,iy11)
         cstet=u*2-1
         sitet=sqrt(1.0-cstet**2)
         call ribm(v,iy12)
         phi=2.0*3.141592654*v

         p1(1)=qdeut*sitet*cos(phi)
         p2(1)=qdeut*sitet*sin(phi)
         p3(1)=qdeut*cstet
         eps(1)=w(p1(1),p2(1),p3(1),fmp)

         p1(2)=-p1(1)
         p2(2)=-p2(1)
         p3(2)=-p3(1)
         eps(2)=eps(1)

         jparticip(1) = 1
         hel(1) = 0
         nesc(1) = 0
         ind2(1) = 1
         ind1(1) = 0
         jparticip(2) = 1
         hel(2) = 0
         nesc(2) = 0
         ind2(2) = -1
         ind1(2) = 0
      else
C     Composite heavier than deuteron:
         s1t1=0.0
         s2t1=0.0
         s3t1=0.0
         sp1t1=0.0
         sp2t1=0.0
         sp3t1=0.0
         do 8 i=1,ia1
            hel(i) = 0
            nesc(i) = 0
            ind2(i) = 1
            if(i.gt.iz1) ind2(i) = -1
            ind1(i) = 0
            jparticip(i) = 1
            call rgauss(xga)
            x1(i) = xga*rms1*0.57735
            s1t1 = s1t1 + x1(i)
            call rgauss(xga)
            x2(i) = xga*rms1*0.57735
            s2t1 = s2t1 + x2(i)
            call rgauss(xga)
            x3(i) = xga*rms1*0.57735
            s3t1 = s3t1 + x3(i)

C     Density of composite as a gaussian in q space:
            call rgauss(xga)
            p1(i) = xga*pf1*0.57735
            call rgauss(xga)
            p2(i) = xga*pf1*0.57735
            call rgauss(xga)
            p3(i) = xga*pf1*0.57735
            eps(i) = w(p1(i),p2(i),p3(i),fmp)
            if(idebug.eq.1) then
               write(6,*) 'nucleon = ', i
               write(6,*) 'x1 = ', x1(i)
               write(6,*) 'x2 = ', x2(i)
               write(6,*) 'x3 = ', x3(i)
               write(6,*) 'p1 = ', p1(i)
               write(6,*) 'p2 = ', p2(i)
               write(6,*) 'p3 = ', p3(i)
               write(6,*) 'eps = ', eps(i)
            endif
            sp1t1 = sp1t1 + p1(i)
            sp2t1 = sp2t1 + p2(i)
            sp3t1 = sp3t1 + p3(i)
 8       continue

         s1t1 = s1t1/ia1
         s2t1 = s2t1/ia1
         s3t1 = s3t1/ia1
         sp1t1 = sp1t1/ia1
         sp2t1 = sp2t1/ia1
         sp3t1 = sp3t1/ia1
         do i=1,ia1
            x1(i) = x1(i) - s1t1
            x2(i) = x2(i) - s2t1
            x3(i) = x3(i) - s3t1
            p1(i) = p1(i) - sp1t1
            p2(i) = p2(i) - sp2t1
            p3(i) = p3(i) - sp3t1
            eps(i) = w(p1(i),p2(i),p3(i),fmp)
         enddo
      endif
    9 continue
c deutons
C**********************************************************
C Target preparation for 1<A<5 (with sum of momentum =0)
      IF(IA2.GE.2.AND.IA2.LE.4) THEN
1633  s1t1=0.
      s2t1=0.
      s3t1=0.
      sp1t1=0.
      sp2t1=0.
      sp3t1=0.
      EFER=0.
      DO 1330 I=IA1+1,IA-1
      IND2(I)=1
      NESC(I)=0
      IF (I.GT.IZ2+IA1) IND2(I)=-1
      DO 633 J=1,7
        call ribm(t(j),ial)
  633 CONTINUE
      T(2)=-1.+2.*T(2)                                              
      T(3)=6.283185*T(3)                                        
      T(5)=-1.+2.*T(5)                                       
      T(6)=6.283185*T(6)                                      
      T1=T(2)                                            
      T2=SQRT(1.-T1*T1)                                       
      T3=COS(T(3))                                             
      T4=SIN(T(3))                                              
      T5=T(5)                                                  
      T6=SQRT(1.-T5*T5)                                         
      T7=COS(T(6))                                             
      T8=SIN(T(6))                                             
	     IF (NOSURF.EQ.1) THEN
      X=R2*T(1)**0.33333333                                      
      Y=PF*T(4)**0.33333333                                    
	     ELSE
c surface..W.S.: impulsion (sphere dure), puis R(q)
       t33=t(7)**0.33333333
C      r44=(ri4+t(7)*(rs4-ri4))**0.25
C      x33=t(4)**0.33333333
C      x=x33*r44
       Y=PF*t33
       RR=FLIN(t33)
       X=RR*T(4)**0.33333333       
c fin surface on a redefini x et y
	     END IF
      X1(I)=X*T2*T3                                                 
      X2(I)=X*T2*T4                                             
      X3(I)=X*T1                                                    
      s1t1=s1t1+x1(i)
      s2t1=s2t1+x2(i)
      s3t1=s3t1+x3(i)
      P1(I)=Y*T6*T7                                              
      P2(I)=Y*T6*T8                                                  
      P3(I)=Y*T5                                                   
      sp1t1=sp1t1+p1(i)
      sp2t1=sp2t1+p2(i)
      sp3t1=sp3t1+p3(i)
      IND1(I)=0                                                       
      EPS(I)=W(P1(I),P2(I),P3(I),FMP)                             
      hel(i)=0.
      EFER=EFER+EPS(I)-FMP                                        
      
1330  CONTINUE 
      	hel(ia)=0 
      	nesc(ia)=0
      	ind2(ia)=-1
      	ind1(ia)=0
      	x1(ia)=-s1t1
      	x2(ia)=-s2t1
      	x3(ia)=-s3t1
      	p1(ia)=-sp1t1
      	p2(ia)=-sp2t1
      	p3(ia)=-sp3t1
      p_mod = SQRT(P1(IA)**2+P2(IA)**2+P3(IA)**2)
      IF(p_mod.GT.(PF+0.05)) GO TO 1633
        eps(ia)=w(p1(ia),p2(ia),p3(ia),fmp)
        EFER=EFER+EPS(ia)-FMP                                        
C for rho(r),rho(q) checking
C      DO I=IA1+1,IA
C      r_dist = SQRT(X1(I)*X1(I)+X2(I)*X2(I)+X3(I)*X3(I))
C      CALL HFILL(3,r_dist,0.,1.)
C      p_mod = SQRT(P1(I)**2+P2(I)**2+P3(I)**2)
C      CALL HFILL(4,p_mod,0.,1.)
C      ENDDO
      
      
      END IF       !(IA2.GE.2.AND.IA2.LE.4)
C******************************************************
C Target preparation for A>4
      IF(IA2.GT.4) THEN
        X1_target=0.
        X2_target=0.
        X3_target=0.
      DO 1 I=IA1+1,IA                                                   P-N04030
      NESC(I)=0                                                         P-N04040
      IND2(I)=1                                                         P-N04050
      IF (I.GT.IZ2+IA1) IND2(I)=-1                                      P-N04060
c surface ajout de t(7) surface.f avait do 6 j=2,7 ici 1,6 ?
    5 DO 6 J=1,7                                                        P-N04070
      call ribm(t(j),ial)                                               P-N04110
    6 continue                                                          P-N04120
      T(2)=-1.+2.*T(2)                                                  P-N04130
      T(3)=6.283185*T(3)                                                P-N04140
      T(5)=-1.+2.*T(5)                                                  P-N04150
      T(6)=6.283185*T(6)                                                P-N04160
      T1=T(2)                                                           P-N04170
      T2=SQRT(1.-T1*T1)                                                 P-N04180
      T3=COS(T(3))                                                      P-N04190
      T4=SIN(T(3))                                                      P-N04200
      T5=T(5)                                                           P-N04210
      T6=SQRT(1.-T5*T5)                                                 P-N04220
      T7=COS(T(6))                                                      P-N04230
      T8=SIN(T(6))                                                      P-N04240
	     IF (NOSURF.EQ.1) THEN
      X=R2*T(1)**0.33333333                                             P-N04250
      Y=PF*T(4)**0.33333333                                             P-N04260
	     ELSE
c surface..W.S.: impulsion (sphere dure), puis R(q)
       t33=t(7)**0.33333333
C      r44=(ri4+t(7)*(rs4-ri4))**0.25
C      x33=t(4)**0.33333333
C      x=x33*r44
       Y=PF*t33
       RR=FLIN(t33)
       X=RR*T(4)**0.33333333       
c fin surface on a redefini x et y
	     END IF
      X1(I)=X*T2*T3                                                     P-N04270
      X2(I)=X*T2*T4                                                     P-N04280
      X3(I)=X*T1                                                        P-N04290
      P1(I)=Y*T6*T7                                                     P-N04300
      P2(I)=Y*T6*T8                                                     P-N04310
      P3(I)=Y*T5                                                        P-N04320
	X1_target=X1_target + X1(I)
	X2_target=X2_target + X2(I)
	X3_target=X3_target + X3(I)
    4 IND1(I)=0                                                         P-N04330
      EPS(I)=W(P1(I),P2(I),P3(I),FMP)                                   P-N04340
      hel(i)=0.
      EFER=EFER+EPS(I)-FMP                                              P-N04350
C for rho(r) checking
C      r_dist = SQRT(X1(I)*X1(I)+X2(I)*X2(I)+X3(I)*X3(I))
C      CALL HFILL(3,r_dist,0.,1.)
    1 CONTINUE                                                          P-N04360
	X1_target=X1_target/IA2
	X2_target=X2_target/IA2
	X3_target=X3_target/IA2
	
      END IF                !(IA2.GT.4)

      EFRUN=EFRUN+EFER                                                  P-N04370
C                                                                       P-N04380
CCC   LOCATION OF incident particle AT POINT (B,Z)                      P-N04390
C                                                                       P-N04400
      R22=R2*R2                                                         P-N04410
C correction 24 janv 99 A.B.
C      Z=R22-B*B                                                        P-N04420
C      Z=r2s*r2s - B*B
       Z=bmax*bmax - B*B       ! for the Wood-Saxon density...
       
      IF (Z.LT.0.) Z=0.
   20 Z=SQRT(Z)                                                         P-N04460
C Random azimuthal direction of the impact parameter (sept 99)

      if (kindf7.le.2.and.kindf7.gt.0) then
      CALL RIBM(TBID,IY14)
      TBID=TBID*6.283185
        X1(1)=X1(1)+B*COS(TBID)                                        
        X2(1)=X2(1)+B*SIN(TBID)
        X3(1)=X3(1)-Z
C pour le ntuple des avatars:
                IF(Kveux.EQ.1) THEN
		R1_in(1)=X1(1)
		R1_in(2)=X2(1)
		R1_in(3)=X3(1)
		ENDIF
       else
        if (kindf7.lt.6.and.kindf7.gt.0) then ! pour les pions on laisse
      CALL RIBM(TBID,IY14)
      TBID=TBID*6.283185
        y1(1)=y1(1)+b*COS(TBID)
        y2(1)=y2(1)+b*SIN(TBID)
        y3(1)=y3(1)-z
	else
c deutons
         nmiss=0.
         xlengm=1000.
         do 21 i=1,ia1
         x3(i)=x3(i)/ginc
         zai2=RMAXWS*RMAXWS-(b+x1(i))**2-x2(i)**2
         if (zai2.lt.0.) go to 22
         ztu=-sqrt(zai2)
c r22 remplace par RMAXWS*RMAXWS et r2 par RMAXWS CV correct ?
         za_i=2.*RMAXWS+ztu
         xleng=za_i-x3(i)
         if (xleng.gt.xlengm) go to 21
         ilm=i
         xlengm=xleng
         ztouch=ztu
         go to 21
   22    nmiss=nmiss+1      
   21    continue
         if (nmiss.eq.ia1) then
             nopart=-1
           return
         else
          zshif=x3(ilm)-ztouch
      CALL RIBM(TBID,IY14)
      TBID=TBID*6.283185
          do i=1,ia1
          xxx=x1(i)+b
          x1(i)=xxx*cos(TBID)-x2(i)*sin(TBID)
          x2(i)=xxx*sin(TBID)+x2(i)*cos(TBID)
          x3(i)=x3(i)-zshif
c           x1(i)=x1(i)+b
c           x3(i)=x3(i)-zshif
          enddo
        if (abs(x1(ilm)**2+x2(ilm)**2+x3(ilm)**2-RMAXWS*RMAXWS).gt.0.01)
     -      write (6,*) 'wrong position'
         endif
        endif         
      endif  
      
C     for rho(r),rho(q) checking
      DO I=IA1+1,IA
         r_dist = SQRT(X1(I)*X1(I)+X2(I)*X2(I)+X3(I)*X3(I))
         CALL HFILL(3,r_dist,0.,1.)
         p_mod = SQRT(P1(I)**2+P2(I)**2+P3(I)**2)
         CALL HFILL(4,p_mod,0.,1.)
         CALL HFILL(21,r_dist,p_mod,1.)
      ENDDO

c deutons
c      endif

C Initial momentum for all type of incident particles:
      XL1=B*PINC*SIN(TBID)                                  
      XL2=-B*PINC*COS(TBID)                                          
      XL3=0.                                           

C                                                                       P-N04490
CCC   TRANSCRIPTION IN THE GENERAL FRAME OF REFERENCE                   P-N04500
CCC        (HERE,=LAB FRAME)                                            P-N04510
C                                                                       P-N04520
   32 BE=0.                                                             P-N04530
      GE=1.                                                             P-N04540
      B1=(BINC-BE)/(1.-BE*BINC)                                         P-N04550
      B2=-BE                                                            P-N04560
      G1=1./SQRT(1.-B1*B1)                                              P-N04570
      G2=1.                                                             P-N04580
c deutons
C Here for nucleons
      if (kindf7.le.2.and.kindf7.gt.0) then
        EPS(1)=G1*FMP+V0                                                P-N04590
        P3(1)=SQRT(EPS(1)**2-FMP**2)
      else
C Here for pions
       if (kindf7.lt.6.and.kindf7.gt.0) then
        q4(1)=g1*fmpi
        q3(1)=b1*q4(1)
       else
C Here for composit projectiles:
C    The kinetic energy is below the threshold. Put all
C    Fermi momentum to 0... projectile nucleons not on shell!
	energie_in=tlab+fmpinc
      IF((energie_in).LE.(ia1*FMP)) THEN
      	DO i=1,ia1
		eps(i)=energie_in/ia1
		P1(i)=0.
		P2(i)=0.
		P3(i)=PINC/ia1
	ENDDO
	GO TO 1871
      ENDIF
C    Here the composit is above threshold
      DO i=1,ia1	!save E,P in the composit rest frame
      	eps_c(i)=eps(i)
	p3_c(i)=p3(i)
      ENDDO
	
       nbtest=ia1-1
       if(kindf7.EQ.6) nbtest=2
       IFLAG=0
1870        sueps=0.

       IFLAG=IFLAG+1
       
        do i=1,ia1
        tte=eps_c(i)
        eps(i)=g1*(eps_c(i)+b1*p3_c(i))
        p3(i)=g1*(b1*tte+p3_c(i))
        sueps=sueps+eps(i)
        enddo
	
        cobe=(tlab+fmpinc)/sueps
	
C OFF shell problem for incident clusters (A.B. 2/2002)
C	IFLAGRE=0
			IF(IFLAG.EQ.nbtest) THEN !Too much..all momentum to 0
				DO klm=1,ia1
				eps_c(klm)=FMP
				P1(klm)=0.
				P2(klm)=0.
				P3_c(klm)=0
				ENDDO
			GO TO 1870
			ENDIF
        do i=1,ia1
	ARG = (cobe*eps(i))**2-pm2
		IF (ARG.LE.0.) THEN	! put maximum momentum to 0. 
			i_emax=1	!find maximum
			ener_max=eps(1)
			DO klm=2,ia1	
				IF(eps(klm).GT.ener_max) THEN
				ener_max=eps(klm)
				i_emax=klm
				ENDIF
			ENDDO
			eps_c(i_emax)=FMP
			P1(i_emax)=0.
			P2(i_emax)=0.
			P3_c(i_emax)=0.
			
			IF(i_emax.EQ.ia1) THEN	!circular permut if the last one
			epsv=eps_c(ia1)		!permutation circulaire
			p1v=p1(ia1)
			p2v=p2(ia1)
			p3v=p3_c(ia1)
			  DO k=ia1-1,1,-1
			    eps_c(k+1)=eps_c(k)
			    p1(k+1)=p1(k)
			    p2(k+1)=p2(k)
			    p3_c(k+1)=p3_c(k)
			  ENDDO
			  eps_c(1)=epsv
			  p1(1)=p1v
			  p2(1)=p2v
			  p3_c(1)=p3v		!fin permut.
			  
			ENDIF
			sp1t1=0.		!re-compute the last one 
      			sp2t1=0.
      			sp3t1=0.
			DO j=1,ia1-1
      			  sp1t1=sp1t1+p1(j)
      			  sp2t1=sp2t1+p2(j)
      			  sp3t1=sp3t1+p3_c(j)
			ENDDO
      			p1(ia1)=-sp1t1
      			p2(ia1)=-sp2t1
      			p3_c(ia1)=-sp3t1
      			eps_c(ia1)=w(p1(ia1),p2(ia1),p3_c(ia1),fmp)	
		
			GO TO 1870      !..and boost all of them.
		ENDIF	
        enddo			
			
        do i=1,ia1
	ARG = (cobe*eps(i))**2-pm2
        comom=sqrt(ARG/(eps(i)**2-pm2))
        p1(i)=comom*p1(i)
        p2(i)=comom*p2(i)
        p3(i)=comom*p3(i)
        eps(i)=eps(i)*cobe
        if (abs(am(p1(i),p2(i),p3(i),eps(i))-fmp).gt.0.01)
     -             write (6,*) 'wrong correction',i                   
        enddo	
        eps(ilm)=eps(ilm)+v0  
	            
       endif
      endif
      
1871  CONTINUE
C for rho(r),rho(q) checking
C      DO I=1,IA1
C      r_dist = SQRT(X1(I)*X1(I)+X2(I)*X2(I)+X3(I)*X3(I))
C      CALL HFILL(3,r_dist,0.,1.)
C      p_mod = SQRT(P1(I)**2+P2(I)**2+P3_c(I)**2)
C      CALL HFILL(4,p_mod,0.,1.)
C      ENDDO
c deutons
C Write informations at the first call to INCL
C      if(ibert.ne.0) GO TO 6000                                         P-N04610
C      PRINT 102,K1,K2,K3,K4,K5,K6
C      PRINT 107,LINE
C      WRITE(6,*) 'R*P cell for pauli stat: ',RBL,PBL
C      WRITE(6,*) 'Fermi momentum: ',PF
C      PRINT 108                                                         P-N04630
C      PRINT 106,B1,B2,G1,G2                                             P-N04640
C      PRINT 107,LINE                                                    P-N04650
C      PRINT 107,LINE                                                    P-N04660
C 6000 CONTINUE                                                          P-N04670
C      	ibert=1
	
	
C                                                                       P-N04680
CCC   EVALUATION OF THE TIMES T(A,B)                                    P-N04690
C                                                                       P-N04700
      K=0                                                               P-N04710
      KCOL=0                                                            P-N04720
      if (kindf7.le.2.and.kindf7.gt.0) then
c     DO 40 I=2,IA                                                      P-N04730
c modif S.Vuillier tient compte propagation projectile,1e collision
c imposee pour lui (c'est une maniere de faire!!)
      DO 40 I=1,IA
      TREF=REF(X1(I),X2(I),X3(I),P1(I),P2(I),P3(I),EPS(I),R22)          P-N04740
      IF (TREF.GT.TMAX5) GO TO 45                                       P-N04750
      K=K+1                                                             P-N04760
      CROIS(K)=TREF                                                     P-N04770
      IND(K)=I                                                          P-N04780
      JND(K)=-1                                                         P-N04790
   45 I1=I-1                                                            P-N04800
      if (i.eq.1) go to 40
C ICI on ne calcule que les interactions NN impliquant le projectile !!! (2/02)
                          IF (I1.GT.IA1) I1=IA1
C**********************************************************************               
      DO 41 J=1,I1                                                      P-N04810
c    no collisions before the first collision of the incident particle
c      do 41 j=1,1
      CALL TIME (I,J)                                                   P-N04820
      IF (TA.LT.0.) GO TO 41                                            P-N04830
      IF(TA.GT.TMAX5) GO TO 41                                          P-N04840
      EIJ=AM(P1(I)+P1(J),P2(I)+P2(J),P3(I)+P3(J),EPS(I)+EPS(J))         P-N04850
      IF (EIJ.LT.1925.) GO TO 41                                        P-N04860
      ISOS=IND2(I)+IND2(J)                                              P-N04870
      IF (31.*RAB2.GT.STOT(EIJ,0,ISOS)) GO TO 41                        P-N04880
      K=K+1                                                             P-N04890
      IF (J.EQ.1) KCOL=KCOL+1                                           P-N04900
      CROIS(K)=TA                                                       P-N04910
      IND(K)=I                                                          P-N04920
      JND(K)=J                                                          P-N04930
   41 CONTINUE                                                          P-N04940
   40 CONTINUE                                                          P-N04950
      else
c deutons
      if (kindf7.lt.6.and.kindf7.gt.0) then
C Here for incoming pions:
      do i=ia1+1,ia
      TREF=REF(X1(I),X2(I),X3(I),P1(I),P2(I),P3(I),EPS(I),R22)
       IF (TREF.LT.TMAX5) then
         K=K+1
         CROIS(K)=TREF
         IND(K)=I
         JND(K)=-1
       ENDIF
      enddo
      call new2(y1(1),y2(1),y3(1),q1(1),q2(1),q3(1),q4(1),1,0)
C Modif A.B. 21/06/2002: Should check at least one valid collision
C    with incoming pion.
C      kcol=1
       IF(K.NE.0) KCOL=1
c      endif
c deutons
      else
      
      do 38 i=1,ia1
      nesc(i)=1
      if (i.ne.ilm) go to 36
      TREF=REF(X1(I),X2(I),X3(I),P1(I),P2(I),P3(I),EPS(I),R22)
      nesc(i)=0
      npproj(i)=0
      go to 37
   36 T1=X1(i)*P1(i)+X2(i)*P2(i)+X3(i)*P3(i)                      
      T2=P1(i)*P1(i)+P2(i)*P2(i)+P3(i)*P3(i)               
      T3=T1/T2                                                          
      T4=X1(i)*X1(i)+X2(i)*X2(i)+X3(i)*X3(i)                          
c incoming nucleons enter potential at maximum radius (modif. 13/06/01)
      T5=T3*T3+(RMAXWS*RMAXWS-T4)/T2                                                                                             
      IF (T5.LT.0.) GO TO 38                   
      TREF=(-T3-SQRT(T5))*EPS(I)  
      IF (TREF.GT.TMAX5) GO TO 38
      npproj(i)=1
   37 K=K+1                                                            
      CROIS(K)=TREF                                                     
      IND(K)=I                                                         
      JND(K)=-1
   38 continue
      kcol=1
      
      DO  39 I=IA1+1,IA
      npproj(i)=0
      TREF=REF(X1(I),X2(I),X3(I),P1(I),P2(I),P3(I),EPS(I),R22)
      IF (TREF.LT.TMAX5) then
       K=K+1
       CROIS(K)=TREF
       IND(K)=I
       JND(K)=-1
       endif
      CALL TIME (I,ILM)                                               
      IF (TA.LT.0.) GO TO 39                                            
      IF (TA.GT.TMAX5) GO TO 39                                          
      EIJ=AM(P1(I)+P1(ILM),P2(I)+P2(ILM),P3(I)+P3(ILM),EPS(I)+EPS(ILM))
      IF (EIJ.LT.1925.) GO TO 39                                        
      ISOS=IND2(I)+IND2(ILM)                                             
      IF (31.*RAB2.GT.STOT(EIJ,0,ISOS)) GO TO 39                        
      K=K+1                                                             
      KCOL=KCOL+1                                           
      CROIS(K)=TA                                                       
      IND(K)=I                                                          
      JND(K)=ILM
   39 CONTINUE 
      endif
      endif                 
c deutons
      IF (KCOL.NE.0) GO TO 48                                           P-N04960
      nopart=-1                                                         P-N04970
c pour eviter renvoi des resultats du run precedent CV 7/7/98
       iarem=ia2
       izrem=iz2
       esrem=0.
       erecrem=0.
c fin ajout CV
      return
C                                                                       P-N05000
CCC   INITIALIZATION  AT THE BEGINNING OF THE RUN                       P-N05010
C                                                                       P-N05020
   48 TIMI=0.                                                           P-N05030
      TIM=0.                                                            P-N05040
      NCOL=0                                                            P-N05050
C compteur des collisions a deux corps (call collis acceptes par pauli)
      NCOL_2C=0
      
C---  NPION=0                                                           P-N05060
      MRNN=0                                                            P-N05070
      MRND=0                                                            P-N05080
      MRDD=0                                                            P-N05090
      MRDN=0                                                            P-N05100
      MRDP=0                                                            P-N05110
      MRPD=0                                                            P-N05120
      MCDD=0                                                            P-N05130
      MPAUL2=0                                                          P-N05140
      MPAUL1=0                                                          P-N05150
C      ITCH=IZ2                                                          P-N05160
C Approx. (not considering escaping protons of incident clusters) 11/03 A.B.
      ITCH = IZ1 + IZ2 -1
      DO 47 I=1,IA                                                      P-N05170
      TLG(I)=0.                                                         P-N05180
   47 NC(I)=0                                                           P-N05190
      ITT=1                                                             P-N05200
C Tableau des energies a l'initialisation (avatar.hbk)
        IF(Kveux.EQ.1) THEN
	DO i=ia1+1,IA
		EPSd(i)=EPS(i)
	ENDDO
	IFLAG20=0
	IFLAG40=0
	IFLAG60=0
	ENDIF

C                                                                       P-N05210
CCC   SEARCH FOR THE SMALLEST POSITIVE T(A,B)                           P-N05220
C                                                                       P-N05230
C Pour tests, =0 interdit les  reflexions avant un avatar du projectile,
C =1 comme avant (reflexions autorisees). (A.B. 3/2002)
C                irst_avatar=0
                irst_avatar=1
   
  449 NEXT=1                                                            P-N05240
      INDIC(NEXT)=1                                                     P-N05250
   44 IF(NEXT.EQ.0) GO TO 449                                           P-N05260
      IDEP=INDIC(NEXT)+1                                                P-N05270
      TAU=CROIS(IDEP-1)                                                 P-N05280
      IF(IDEP.GT.K) GO TO 448                                           P-N05290
      DO 42 I=IDEP,K                                                    P-N05300
      IF (CROIS(I).GT.TAU) GO TO 42                                     P-N05310
      TAU=CROIS(I)                                                      P-N05320
      NEXT=NEXT+1                                                       P-N05330
      INDIC(NEXT)=I                                                     P-N05340
   42 CONTINUE                                                          P-N05350
  448 IMIN=INDIC(NEXT)                                                  P-N05360
      L1=IND(IMIN)                                                      P-N05370
      L2=JND(IMIN)                                                      P-N05380
C TEST le 20/3/2003: tue sinon le dernier avatar?
C      K=K-1                                                             P-N05390
      IF (K.EQ.0) GO TO 230                                             P-N05400
      K=K-1                                                             P-N05390
      NEXT=NEXT-1                                                       P-N05410
c     IF (IMIN.GE.K) GO TO 46                                           P-N05420
c correction S.Vuillier 25/1/96 decalage temps correct
      IF (IMIN.GT.K) GO TO 46
      DO 43 I=IMIN,K                                                    P-N05430
      CROIS(I)=CROIS(I+1)                                               P-N05440
      IND(I)=IND(I+1)                                                   P-N05450
   43 JND(I)=JND(I+1)                                                   P-N05460
   46 TIM=TIMI+TAU                                                      P-N05470
C Tableau des energies a t=20,40,60 fm/c (avatar.hbk)
        IF(Kveux.EQ.1) THEN
	IF(IFLAG20.EQ.0.AND.TIM.GE.20.) THEN
		IFLAG20=1
		DO i=1,IA
			IF(NESC(i).EQ.0) then
			IF(JPARTICIP(i).EQ.1) then
				EPS2(i)=EPS(i)
			ELSE
				EPS2(i)=0.
			ENDIF
			ELSE
			EPS2(i)=0.
			ENDIF
		ENDDO
	ENDIF
	IF(IFLAG40.EQ.0.AND.TIM.GE.40.) THEN
		IFLAG40=1
		DO i=1,IA
			IF(NESC(i).EQ.0) then
			IF(JPARTICIP(i).EQ.1) then
				EPS4(i)=EPS(i)
			ELSE
				EPS4(i)=0.
			ENDIF
			ELSE
			EPS4(i)=0.
			ENDIF
		ENDDO
	ENDIF
	IF(IFLAG60.EQ.0.AND.TIM.GE.60.) THEN
		IFLAG60=1
		DO i=1,IA
			IF(NESC(i).EQ.0) then
			IF(JPARTICIP(i).EQ.1) then
				EPS6(i)=EPS(i)
			ELSE
				EPS6(i)=0.
			ENDIF
			ELSE
			EPS6(i)=0.
			ENDIF
		ENDDO
	ENDIF

        ENDIF
C                                                                       P-N05480
CCC   COUNTING THE POSITIONS AND THE MOMENTA AT INTERMEDIATE TIMES      P-N05490
C                                                                       P-N05500
c  645 IF (TIM.LT.TEM(ITT)) GO TO 49                                    P-N05510
  645 CONTINUE
C Modif: Pas de reflexions avant au moins un avatar du (des) nucleon incident
C   Celui-ci ne peut etre qu'une collision NN (ou piN)
C           	WRITE(6,*) 'L1,L2',L1,L2,jparticip(L1),irst_avatar
					
               IF((irst_avatar.EQ.0).AND.(L2.EQ.-1)) GO TO 44
	       
		irst_avatar=irst_avatar+1
		
C		WRITE(6,*) 'avatar',irst_avatar				  
C               WRITE(6,*) 'L1,L2=',L1,L2
  
      if (tim.lt.temfin) go to 49

       
      go to 255
   49 CONTINUE                                                          P-N06330
C      IF (K.EQ.0) GO TO 312                                            P-N06340
       if (k.eq.0) go to 255
C L1 va a la surface du noyau:
      IF (L2.EQ.-1) GO TO 220                                           P-N06350
      IF(K4-1) 803,803,830                                              P-N06360
C L1 est un delta: 
  830 IF(L2.EQ.0) GO TO 220                                             P-N06370
C Interaction pi(L1-IA)-nucleon(L2)
      IF(L1.GT.IA) GO TO 801                                            P-N06380
  803 CONTINUE                                                          P-N06390

C Pas de collision entre 2 non participants:
      IF(JPARTICIP(L1).EQ.0.AND.JPARTICIP(L2).EQ.0) GO TO 44

C                                                                       P-N06400
CCC   PARAMETERS FOR THE NEXT COLLIDING PAIR                            P-N06410
C                                                                       P-N06420
      T(10)=EPS(L1)+EPS(L2)                                             P-N06430
      T0=1./T(10)                                                       P-N06440
      B1=(P1(L1)+P1(L2))*T0                                             P-N06450
      B2=(P2(L1)+P2(L2))*T0                                             P-N06460
      B3=(P3(L1)+P3(L2))*T0                                             P-N06470
      S=(1.-B1*B1-B2*B2-B3*B3)*T(10)*T(10)                              P-N06480
      SQ=SQRT(S)                                                        P-N06490
c      if (io.eq.342) write (6,*) 'sq=',sq
      IF(SQ.LT.1925.5) GO TO 44                                         P-N06500
      TA=TAU/EPS(L1)                                                    P-N06510
      X1L1=X1(L1)+P1(L1)*TA                                             P-N06520
      X2L1=X2(L1)+P2(L1)*TA                                             P-N06530
      X3L1=X3(L1)+P3(L1)*TA                                             P-N06540
      TA=TAU/EPS(L2)                                                    P-N06550
      X1L2=X1(L2)+P1(L2)*TA                                             P-N06560
      X2L2=X2(L2)+P2(L2)*TA                                             P-N06570
      X3L2=X3(L2)+P3(L2)*TA                                             P-N06580
C                                                                       P-N06590
CCC   TEST ON THE MINIMUM DISTANCE OF APPROACH                          P-N06600
C                                                                       P-N06610
      T(11)=X1L1-X1L2                                                   P-N06620
      T(12)=X2L1-X2L2                                                   P-N06630
      T(13)=X3L1-X3L2                                                   P-N06640
      T(14)=T(11)*T(11)+T(12)*T(12)+T(13)*T(13)                         P-N06650
      T(15)=B1*T(11)+B2*T(12)+B3*T(13)                                  P-N06660
      T(16)=B1*B1+B2*B2+B3*B3                                           P-N06670
      BB2=T(14)+T(15)*T(15)/(1.-T(16))                                  P-N06680
      IF (K3.EQ.1) GO TO 260                                            P-N06690
      IF (K4.EQ.0) GO TO 260                                            P-N06700
      MG=IND1(L1)+IND1(L2)                                              P-N06710
      ISOS=IND2(L1)+IND2(L2)                                            P-N06720
      IF (MG.NE.1) GO TO 260                                            P-N06730
      LDEL=L2                                                           P-N06740
      IF(MG-IND1(L1).EQ.0) LDEL=L1                                      P-N06750
      XX10=SQRT(EPS(LDEL)**2-P1(LDEL)**2-P2(LDEL)**2-P3(LDEL)**2)       P-N06760
      ISA=IND2(LDEL)                                                    P-N06770
      BMAX2=STOT(SQ,MG,ISOS)/31.415926                                  P-N06780
      IF (K5.EQ.0.AND.MG.NE.0) BMAX2=BMAX2-SEL(SQ,MG,ISOS)/31.415926    P-N06790
      GO TO 261                                                         P-N06800
  260 BMAX2=STOT(SQ,MG,ISOS)/31.41592                                   P-N06810
  261 IF (BB2.LT.BMAX2) GO TO 220                                       P-N06820
      IF (K.EQ.0) GO TO 230                                             P-N06830
      GO TO 44                                                          P-N06840
C                                                                       P-N06850
CCC   EVALUATION OF THE POSITIONS AT TIME=TIM                           P-N06860
C                                                                       P-N06870
  220 CONTINUE                                                          P-N06880
c      if (io.eq.342.or.io.eq.701) write (6,*) l1,l2, ind1(l1),
c     - ind1(l2),ind2(l1),ind2(l2), eps(l1)+eps(l2)
      TIMI=TIM                                                          P-N06890
        IF(Kveux.EQ.1) THEN
         Iavat=Iavat+1
	 TIMEavat(Iavat)=TIM
	 L1avat(Iavat)=L1
	 L2avat(Iavat)=L2
	 ENERGYavat(Iavat)=SQ
	 IF(L1.LE.ia) THEN 
	 	JPARTL1(Iavat)=JPARTICIP(L1)
	 ELSE
	 	JPARTL1(Iavat)=0
	 ENDIF
	 IF(L2.GT.0) THEN
	 	JPARTL2(Iavat)=JPARTICIP(L2)
	 ELSE
	 	JPARTL2(Iavat)=0
	 ENDIF
	ENDIF
C Gel des nucleons non participants sur le premier avatar (NN)=(L1,1)      
      IF (irst_avatar.EQ.1) THEN
      	   DO i=1,L1,L1-1
C	        WRITE(6,*) 'L1,L2,i=',L1,L2,i
      		TA=TAU/EPS(I)                                                
      		X1(I)=X1(I)+P1(I)*TA                                      
      		X2(I)=X2(I)+P2(I)*TA                                    
      		X3(I)=X3(I)+P3(I)*TA
           ENDDO
	   DO i=1,K
	   	CROIS(i)=CROIS(i)+TAU
	   ENDDO                                   
      ELSE
         
      DO 51 I=1,IA                                                      P-N06900
      TA=TAU/EPS(I)                                                     P-N06910
      X1(I)=X1(I)+P1(I)*TA                                              P-N06920
      X2(I)=X2(I)+P2(I)*TA                                              P-N06930
   51 X3(I)=X3(I)+P3(I)*TA                                              P-N06940
      ENDIF

      IF(NPION.EQ.0) GO TO 840                                          P-N06950
      DO 804 I=1,NPION                                                  P-N06960
      TA=TAU/Q4(I)                                                      P-N06970
      Y1(I)=Y1(I)+Q1(I)*TA                                              P-N06980
      Y2(I)=Y2(I)+Q2(I)*TA                                              P-N06990
  804 Y3(I)=Y3(I)+Q3(I)*TA                                              P-N07000
  840 IF(L2.EQ.0) GO TO 805                                             P-N07010
  
C coupure sur l'energie cinetique maximale **********************************
C      EPSMAX=0.
C      DO klm=1,ia
C      	  IF(eps(klm).GT.EPSMAX.AND.nesc(klm).EQ.0) EPSMAX=EPS(klm)
C      ENDDO
C            
C      IF (EPSMAX.LT.(FMP+f(5)+10.)) GO TO 255
C****************************************************************************

C Reflexions sur le potentiel, sortie eventuelle de la particule:
      IF (L2.EQ.-1) GO TO 600                                           P-N07020
      
      IF(L1.GT.IA) GO TO 831                                            P-N07030
C                                                                       P-N07040
CCC   COLLISION OF PARTICLES L1 AND L2                                  P-N07050
C                                                                       P-N07060
      ICH1=IND1(L1)                                                     P-N07070
      ICH2=IND1(L2)                                                     P-N07080
      ICH3=IND2(L1)                                                     P-N07090
      ICH4=IND2(L2)                                                     P-N07100
      AML1=SQRT(EPS(L1)**2-P1(L1)**2-P2(L1)**2-P3(L1)**2)               P-N07120
      AML2=SQRT(EPS(L2)**2-P1(L2)**2-P2(L2)**2-P3(L2)**2)               P-N07130
      GL1=EPS(L1)/AML1                                                  P-N07140
      GL2=EPS(L2)/AML2                                                  P-N07150
C++++++++++++++ L-CONSERVATION +++++++++++++++++++++++++++++++++++++++
      IF (K6.EQ.1) THEN                                                 P-N07160
      T(31)=(AML1*X1(L1)+AML2*X1(L2))/(AML1+AML2)                       P-N07170
      T(32)=(AML1*X2(L1)+AML2*X2(L2))/(AML1+AML2)                       P-N07180
      T(33)=(AML1*X3(L1)+AML2*X3(L2))/(AML1+AML2)                       P-N07190
      TT31 =X1(L1)-X1(L2)                                               P-N07200
      TT32 =X2(L1)-X2(L2)                                               P-N07210
      TT33 =X3(L1)-X3(L2)                                               P-N07220
      T(34)=(AML2*P1(L1)-AML1*P1(L2))/(AML1+AML2)                       P-N07230
      T(35)=(AML2*P2(L1)-AML1*P2(L2))/(AML1+AML2)                       P-N07240
      T(36)=(AML2*P3(L1)-AML1*P3(L2))/(AML1+AML2)                       P-N07250
      TT34=P1(L1)+P1(L2)                                                P-N07260
      TT35=P2(L1)+P2(L2)                                                P-N07270
      TT36=P3(L1)+P3(L2)                                                P-N07280
      ENDIF                                                             P-N07290
C++++++++++++++++ L-CONSERVATION ++++++++++++++                         P-N07300
      T(21)=P1(L1)                                                      P-N07310
      T(22)=P2(L1)                                                      P-N07320
      T(23)=P3(L1)                                                      P-N07330
      T(24)=EPS(L1)                                                     P-N07340
      T(25)=P1(L2)                                                      P-N07350
      T(26)=P2(L2)                                                      P-N07360
      T(27)=P3(L2)                                                      P-N07370
      T(28)=EPS(L2)                                                     P-N07380
C Info delta ou nucleon:
	      IF(Kveux.EQ.1) THEN
		DEL1avat(Iavat)=IND1(L1)
                DEL2avat(Iavat)=IND1(L2)                                               
	      ENDIF                                               
      CALL LOREN(P1(L1),P2(L1),P3(L1),-B1,-B2,-B3,EPS(L1))              P-N07390
      CALL LOREN(P1(L2),P2(L2),P3(L2),-B1,-B2,-B3,EPS(L2))              P-N07400
c      if (io.eq.342.or.io.eq.701) write (6,*) l1,l2, ind1(l1),
c     - ind1(l2),ind2(l1),ind2(l2), eps(l1),eps(l2)
      CALL COLLIS(P1(L1),P2(L1),P3(L1),EPS(L1),P1(L2),P2(L2),P3(L2),EPS(P-N07410
     -L2),T(12),T(13),T(14),T(15),NP,IP,K2,K3,K4,K5,IND1(L1),IND1(L2),INP-N07420
     -D2(L1),IND2(L2))                                                  P-N07430
c      if (io.eq.342.or.io.eq.701) write (6,*) l1,l2, ind1(l1),
c     - ind1(l2),ind2(l1),ind2(l2), eps(l1),eps(l2)
      CALL LOREN(P1(L1),P2(L1),P3(L1),B1,B2,B3,EPS(L1))                 P-N07440
      CALL LOREN(P1(L2),P2(L2),P3(L2),B1,B2,B3,EPS(L2))                 P-N07450

      IF (IND1(L1).EQ.1) GO TO 243                                      P-N07460
      CALL PAULI(L1,RBL,PBL,XBL1)                                       P-N07470
      call ribm(rndm,iy11) 
      IF (RNDM.GT.1.-XBL1) GO TO 248                                    P-N07530
  243 IF (IND1(L2).EQ.1) GO TO 241                                      P-N07540
      CALL PAULI(L2,RBL,PBL,XBL2)                                       P-N07550
      call ribm(rndm,iy11) 
      IF (RNDM.GT.1.-XBL2) GO TO 248                                    P-N07610
      GO TO 241                                                         P-N07620
  248 MPAUL1=MPAUL1+1                                                   P-N07630
        IF(Kveux.EQ.1) Bloc_Paul(Iavat)=1
C Restitution de L1 et L2 si rejet de la col. par Pauli:
      P1(L1)=T(21)                                                      P-N07640
      P2(L1)=T(22)                                                      P-N07650
      P3(L1)=T(23)                                                      P-N07660
      EPS(L1)=T(24)                                                     P-N07670
      P1(L2)=T(25)                                                      P-N07680
      P2(L2)=T(26)                                                      P-N07690
      P3(L2)=T(27)                                                      P-N07700
      EPS(L2)=T(28)                                                     P-N07710
      IND1(L1)=ICH1                                                     P-N07720
      IND1(L2)=ICH2                                                     P-N07730
      IND2(L1)=ICH3                                                     P-N07740
      IND2(L2)=ICH4                                                     P-N07750
      IF (K.EQ.0) GO TO 230                                             P-N07760
      DO 284 I=1,K                                                      P-N07770
  284 CROIS(I)=CROIS(I)-TAU                                             P-N07780
C      GO TO 449                                                         P-N07790
C Pour le temps de calcul (A.B. 02/2002)
				GO TO 44

241   CONTINUE
      
C La premiere collision a deux corps ne peut pas baisser l'energie
C    du nucleon de recul (bloque par pauli dans un noyau cible froid).
C    (Ici, toujours L2 < L1)
       NCOL_2C=NCOL_2C+1
       IF(NCOL_2C.EQ.1) THEN
       do icomp=1,ia1
C Test on the first collision modified 4/07/2001 for direct AND exchange.
         IF(icomp.EQ.L1.OR.icomp.EQ.L2) THEN
           Xavant=MIN(T(24),T(28))
           Xapres=MIN(EPS(L1),EPS(L2))             
           IF(Xapres.LE.Xavant) THEN
        	IF(Kveux.EQ.1) Bloc_CDPP(Iavat)=1
	   GO TO 248
	   ENDIF
         ENDIF
       enddo
C pour le ntuple des avatars:
        	IF(Kveux.EQ.1) THEN
		R1_first_avat(1)=X1(1)
		R1_first_avat(2)=X2(1)
		R1_first_avat(3)=X3(1)
		ENDIF
       ELSE
C Les collisions suivantes ne penvent conduire a un noyau de A nucleons
C  sous l'energie de Fermi et dans une config. d'energie inferieure a
C  EFER-(IA2-NBALTTF)*TF).
       EGS=0.
       NBALTTF=0
       DO I=1,IA
         IF(NESC(I).EQ.0) THEN
	    IF(SQRT(P1(I)**2+P2(I)**2+P3(I)**2).LT.PF) THEN
	       NBALTTF=NBALTTF+1
	       EGS=EGS+EPS(I)-FMP
	    ENDIF
	 ENDIF
       END DO
       IF(EGS.LT.(EFER-(IA2-NBALTTF)*TF)) THEN
        IF(Kveux.EQ.1) Bloc_CDPP(Iavat)=1 
       GO TO 248       
       ENDIF
       ENDIF
       
        IF(Kveux.EQ.1) THEN
         Bloc_CDPP(Iavat)=0
         Bloc_Paul(Iavat)=0
        ENDIF
      
      JPARTICIP(L1)=1
      JPARTICIP(L2)=1
      
		IF (NOSURF.LE.0) THEN
c surface
      pppp=sqrt(p1(l1)**2+p2(l1)**2+p3(l1)**2)
      rrrr=sqrt(x1(l1)**2+x2(l1)**2+x3(l1)**2)
      if (pppp.le.pf) then
	XV=pppp/pf
	rcorr=FLIN(XV)
      if (rrrr.gt.rcorr) then
      x1(l1)=x1(l1)*rcorr/rrrr
      x2(l1)=x2(l1)*rcorr/rrrr
      x3(l1)=x3(l1)*rcorr/rrrr
      endif
      endif
      pppp=sqrt(p1(l2)**2+p2(l2)**2+p3(l2)**2)
      rrrr=sqrt(x1(l2)**2+x2(l2)**2+x3(l2)**2)
      if (pppp.le.pf) then
	XV=pppp/pf
	rcorr=FLIN(XV)
      if (rrrr.gt.rcorr) then
      x1(l2)=x1(l2)*rcorr/rrrr
      x2(l2)=x2(l2)*rcorr/rrrr
      x3(l2)=x3(l2)*rcorr/rrrr
      endif
      endif
		END IF
      IF (NP.EQ.0) GO TO 240                                            P-N07800
      NPION=NPION+1                                                     P-N07810
      CALL LOREN(T(12),T(13),T(14),B1,B2,B3,T(15))                      P-N07820
      Q1(NPION)=T(12)                                                   P-N07830
      Q2(NPION)=T(13)                                                   P-N07840
      Q3(NPION)=T(14)                                                   P-N07850
      Q4(NPION)=T(15)                                                   P-N07860
  240 NCOL=NCOL+1                                                       P-N07870
      IF (L2.NE.1) GO TO 870                                            P-N07880
      
C critere pour la leading particle: AVANT impulsion longitudinale max L=1
C Change en fevrier 2002: leading part. = energie totale max (L=1)
      IF (P3(L2).GT.P3(L1)) GO TO 870                                   P-N07890
C***************************************************************
C ATTENTION, Il faut mieux faire et selectionner la plus grande energie
C   des particules participantes (JPARTICIP()=1) et dans le noyau (NESC()=0)!
C      IF (EPS(L2).GT.EPS(L1)) GO TO 870
C***************************************************************      
      XR1=P1(L1)                                                        P-N07900
      XR2=P2(L1)                                                        P-N07910
      XR3=P3(L1)                                                        P-N07920
      XR4=EPS(L1)                                                       P-N07930
      XR5=X1(L1)                                                        P-N07940
      XR6=X2(L1)                                                        P-N07950
      XR7=X3(L1)                                                        P-N07960
      XR8=GL1                                                           P-N07970
      IXR1=IND1(L1)                                                     P-N07980
      IXR2=IND2(L1)                                                     P-N07990
      IXR3=ICH1                                                         P-N08000
      P1(L1)=P1(L2)                                                     P-N08010
      P2(L1)=P2(L2)                                                     P-N08020
      P3(L1)=P3(L2)                                                     P-N08030
      EPS(L1)=EPS(L2)                                                   P-N08040
      X1(L1)=X1(L2)                                                     P-N08050
      X2(L1)=X2(L2)                                                     P-N08060
      X3(L1)=X3(L2)                                                     P-N08070
      GL1=GL2                                                           P-N08080
      IND1(L1)=IND1(L2)                                                 P-N08090
      IND2(L1)=IND2(L2)                                                 P-N08100
      ICH1=ICH2                                                         P-N08110
      P1(L2)=XR1                                                        P-N08120
      P2(L2)=XR2                                                        P-N08130
      P3(L2)=XR3                                                        P-N08140
      EPS(L2)=XR4                                                       P-N08150
      X1(L2)=XR5                                                        P-N08160
      X2(L2)=XR6                                                        P-N08170
      X3(L2)=XR7                                                        P-N08180
      GL2=XR8                                                           P-N08190
      IND1(L2)=IXR1                                                     P-N08200
      IND2(L2)=IXR2                                                     P-N08210
      ICH2=IXR3                                                         P-N08220
      IF(ICH1+ICH2-IND1(L1)-IND1(L2).NE.0.OR.ICH1+ICH2.NE.1) GO TO 870  P-N08230
      IF (K.EQ.0) GO TO 870                                             P-N08240
      DO 871 I=1,K                                                      P-N08250
      IF (IND(I).NE.1) GO TO 872                                        P-N08260
      IF (JND(I).NE.0) GO TO 872                                        P-N08270
      IND(I)=L1                                                         P-N08280
      GO TO 870                                                         P-N08290
  872 IF (IND(I).NE.L1) GO TO 871                                       P-N08300
      IF (JND(I).NE.0) GO TO 871                                        P-N08310
      IND(I)=1                                                          P-N08320
      GO TO 870                                                         P-N08330
  871 CONTINUE                                                          P-N08340
  870 TLG(L1)=TH*EPS(L1)/SQRT(EPS(L1)**2-P1(L1)**2-P2(L1)**2-P3(L1)**2) P-N08350
      TLG(L2)=TH*EPS(L2)/SQRT(EPS(L2)**2-P1(L2)**2-P2(L2)**2-P3(L2)**2) P-N08360
      NC(L1)=NC(L1)+1                                                   P-N08370
      NC(L2)=NC(L2)+1                                                   P-N08380
      LED=0                                                             P-N08390
      IF(ICH1+ICH2-IND1(L1)-IND1(L2)) 500,501,502                       P-N08400
  500 MRND=MRND+1                                                       P-N08410
      GO TO 503                                                         P-N08420
  501 IF(ICH1+ICH2-1) 504,505,506                                       P-N08430
  504 MRNN=MRNN+1                                                       P-N08440
      GO TO 503                                                         P-N08450
  505 MRDD=MRDD+1                                                       P-N08460
      LED=1                                                             P-N08470
      GO TO 503                                                         P-N08480
  506 MCDD=MCDD+1                                                       P-N08490
      LED=1                                                             P-N08500
      GO TO 503                                                         P-N08510
  502 MRDN=MRDN+1                                                       P-N08520
  503 CONTINUE                                                          P-N08530
  242 CONTINUE                                                          P-N08540
C*******************************************************************
C Test d'arret sur max des energies des particules dans le noyau
C +10Mev=parametre d'energie au dessus du potentiel f(5)
C      IF (eps(1).LT.(FMP+f(5)+10.)) GO TO 255
C*******************************************************************

C                                                                       P-N08550
CCC   REEVALUATION OF THE TIMES T(A,B) FOR (A OR B)=(L1 OR L2)          P-N08560
C                                                                       P-N08570
C+++++++++++ L-CONSERVATION ++++++++++++++++++++++++++++++++++          P-N08580
      IF (K6.EQ.1) THEN                                                 P-N08590
      AML1=AM(P1(L1),P2(L1),P3(L1),EPS(L1))                             P-N08600
      AML2=AM(P1(L2),P2(L2),P3(L2),EPS(L2))                             P-N08610
      T(37)=(AML2*P1(L1)-AML1*P1(L2))/(AML1+AML2)                       P-N08620
      T(38)=(AML2*P2(L1)-AML1*P2(L2))/(AML1+AML2)                       P-N08630
      T(39)=(AML2*P3(L1)-AML1*P3(L2))/(AML1+AML2)                       P-N08640
      T(40)=SQRT(T(34)*T(34)+T(35)*T(35)+T(36)*T(36))                   P-N08650
      T(41)=SQRT(T(37)*T(37)+T(38)*T(38)+T(39)*T(39))                   P-N08660
      RHOPI=TT31*T(34)+TT32*T(35)+TT33*T(36)                            P-N08670
      T(43)=TT31-RHOPI*T(34)/T(40)**2                                   P-N08680
      T(44)=TT32-RHOPI*T(35)/T(40)**2                                   P-N08690
      T(45)=TT33-RHOPI*T(36)/T(40)**2                                   P-N08700
      T(46)=SQRT(T(43)*T(43)+T(44)*T(44)+T(45)*T(45))                   P-N08710
      T(43)=T(43)/T(46)                                                 P-N08720
      T(44)=T(44)/T(46)                                                 P-N08730
      T(45)=T(45)/T(46)                                                 P-N08740
      CIF=(T(34)*T(37)+T(35)*T(38)+T(36)*T(39))/T(40)/T(41)             P-N08750
c trouble with forward scattering 22/3/95
      if(abs(cif).gt.1.)cif=sign(1.,cif)
      SIF=SQRT(1.-CIF*CIF)                                              P-N08760
      T(37)=(T(34)*CIF/T(40)+T(43)*SIF)*T(41)                           P-N08770
      T(38)=(T(35)*CIF/T(40)+T(44)*SIF)*T(41)                           P-N08780
      T(39)=(T(36)*CIF/T(40)+T(45)*SIF)*T(41)                           P-N08790
      TRI=SQRT(TT31*TT31+TT32*TT32+TT33*TT33)                           P-N08800
      CCHI=RHOPI/TRI/T(40)                                              P-N08810
      SCHI=SQRT(1.-CCHI*CCHI)                                           P-N08820
      C1=CIF*CCHI-SIF*SCHI                                              P-N08830
      C2=SIF*CCHI+CIF*SCHI                                              P-N08840
      TT31=(C1*T(34)/T(40)+C2*T(43))*TRI*T(40)/T(41)                    P-N08850
      TT32=(C1*T(35)/T(40)+C2*T(44))*TRI*T(40)/T(41)                    P-N08860
      TT33=(C1*T(36)/T(40)+C2*T(45))*TRI*T(40)/T(41)                    P-N08870
      X1(L1)=T(31)+AML2*TT31/(AML1+AML2)                                P-N08880
      X2(L1)=T(32)+AML2*TT32/(AML1+AML2)                                P-N08890
      X3(L1)=T(33)+AML2*TT33/(AML1+AML2)                                P-N08900
      X1(L2)=T(31)-AML1*TT31/(AML1+AML2)                                P-N08910
      X2(L2)=T(32)-AML1*TT32/(AML1+AML2)                                P-N08920
      X3(L2)=T(33)-AML1*TT33/(AML1+AML2)                                P-N08930
      P1(L1)=AML1*TT34/(AML1+AML2)+T(37)                                P-N08940
      P2(L1)=AML1*TT35/(AML1+AML2)+T(38)                                P-N08950
      P3(L1)=AML1*TT36/(AML1+AML2)+T(39)                                P-N08960
      EPS(L1)=W(P1(L1),P2(L1),P3(L1),AML1)                              P-N08970
      P1(L2)=AML2*TT34/(AML1+AML2)-T(37)                                P-N08980
      P2(L2)=AML2*TT35/(AML1+AML2)-T(38)                                P-N08990
      P3(L2)=AML2*TT36/(AML1+AML2)-T(39)                                P-N09000
      EPS(L2)=W(P1(L2),P2(L2),P3(L2),AML2)                              P-N09010
      ENDIF                                                             P-N09020
C+++++++++++++++ L-CONSERVATION ++++++++++++++++++++++++++++++++++++    P-N09030
c      if (io.eq.342.or.io.eq.701) write (6,*) l1,l2, ind1(l1),
c     -ind1(l2),ind2(l1),ind2(l2), eps(l1)+eps(l2)
      IF (K.EQ.0) GO TO 514                                             P-N09040
      KD=0                                                              P-N09050
      CCR=TAU                                                           P-N09060
      DO 50 I=1,K                                                       P-N09070
      I20=I-KD                                                          P-N09080
      IF (K4.NE.2.OR.LED.NE.1) GO TO 512                                P-N09090
      IF (JND(I).EQ.0) GO TO 511                                        P-N09100
  512 CONTINUE                                                          P-N09110
      IF (IND(I).EQ.L1) GO TO 52                                        P-N09120
      IF (IND(I).EQ.L2) GO TO 52                                        P-N09130
      IF (JND(I).EQ.L2) GO TO 52                                        P-N09140
      IF (JND(I).EQ.L1) GO TO 52                                        P-N09150
      GO TO 513                                                         P-N09160
  511 IF (IND(I).EQ.L1.AND.IND1(L1).EQ.1) CROIS(I)=(CROIS(I)-CCR)*EPS(L1P-N09170
     -)/SQRT(EPS(L1)**2-P1(L1)**2-P2(L1)**2-P3(L1)**2)/GL1+CCR          P-N09180
      IF (IND(I).EQ.L2.AND.IND1(L2).EQ.1) CROIS(I)=(CROIS(I)-CCR)*EPS(L2P-N09190
     -)/SQRT(EPS(L2)**2-P1(L2)**2-P2(L2)**2-P3(L2)**2)/GL2+CCR          P-N09200
  513 CONTINUE                                                          P-N09210
      CROIS(I20)=CROIS(I)-CCR                                           P-N09220
      IND(I20)=IND(I)                                                   P-N09230
      JND(I20)=JND(I)                                                   P-N09240
      GO TO 50                                                          P-N09250
   52 KD=KD+1                                                           P-N09260
   50 CONTINUE                                                          P-N09270
      K=K-KD                                                            P-N09280
  514 CALL NEWT(L1,L2)                                                  P-N09290
      TREF=REF(X1(L1),X2(L1),X3(L1),P1(L1),P2(L1),P3(L1),EPS(L1),R22)   P-N09300
      IF (TREF.GT.TMAX5) GO TO 515                                      P-N09310
      K=K+1                                                             P-N09320
      CROIS(K)=TREF                                                     P-N09330
      IND(K)=L1                                                         P-N09340
      JND(K)=-1                                                         P-N09350
  515 TREF=REF(X1(L2),X2(L2),X3(L2),P1(L2),P2(L2),P3(L2),EPS(L2),R22)   P-N09360
      IF (TREF.GT.TMAX5) GO TO 516                                      P-N09370
      K=K+1                                                             P-N09380
      CROIS(K)=TREF                                                     P-N09390
      IND(K)=L2                                                         P-N09400
      JND(K)=-1                                                         P-N09410
  516 IF (K4.EQ.2) GO TO 848                                            P-N09420
      IF (K) 230,230,449                                                P-N09430
C                                                                       P-N09440
  848 IF (NPION.EQ.0) GO TO 844                                         P-N09450
      IF (IND1(L1).EQ.1) GO TO 843                                      P-N09460
      DO 842 K20=1,NPION                                                P-N09470
  842 CALL NEW3(Y1(K20),Y2(K20),Y3(K20),Q1(K20),Q2(K20),Q3(K20),Q4(K20),P-N09480
     -K20,L1)                                                           P-N09490
  843 IF (IND1(L2).EQ.1) GO TO 844                                      P-N09500
      DO 847 K20=1,NPION                                                P-N09510
  847 CALL NEW3(Y1(K20),Y2(K20),Y3(K20),Q1(K20),Q2(K20),Q3(K20),Q4(K20),P-N09520
     -K20,L2)                                                           P-N09530
  844 CONTINUE                                                          P-N09540
C                                                                       P-N09550
      IF(IND1(L1)+IND1(L2).LE.ICH1+ICH2) GO TO 849                      P-N09560
      IF(IND1(L1)-ICH1.NE.1) GO TO 820                                  P-N09570
      LNEW=L1                                                           P-N09580
      GO TO 821                                                         P-N09590
  820 IF(IND1(L2)-ICH2.NE.1) GO TO 849                                  P-N09600
      LNEW=L2                                                           P-N09610
  821 call ribm(rndm,iy17)
C largeur variable du delta (Phase Space Factor introduced 4/2001)
C (180.**3 = 5832000.)
      AMLNEW=SQRT(EPS(LNEW)**2-P1(LNEW)**2-P2(LNEW)**2-P3(LNEW)**2) 
      GEFF=EPS(LNEW)/AMLNEW
      QQQ=SQRT((AMLNEW**2-(FMP+FMPI)**2)*(AMLNEW**2-(FMP-FMPI)**2))
     s     /(2.*AMLNEW)
      PSF=QQQ**3/(QQQ**3+5832000.)
      TDEL=-HC/(G0*PSF)*ALOG(RNDM)*GEFF                                  
      
      IF (TDEL.GT.TMAX5) GO TO 849                                      P-N09700
      K=K+1                                                             P-N09710
      CROIS(K)=TDEL                                                     P-N09720
      IND(K)=LNEW                                                       P-N09730
      JND(K)=0                                                          P-N09740
  849 IF (K.EQ.0) GO TO 230                                             P-N09750
      GO TO 449                                                         P-N09760
C                                                                       P-N09770
CCC   DECAY OF THE DELTA PARTICLE                                       P-N09780
C                                                                       P-N09790
  805 NPION=NPION+1                                                     P-N09800
      ICHD=IND2(L1)                                                     P-N09810
      T(31)=P1(L1)                                                      P-N09820
      T(32)=P2(L1)                                                      P-N09830
      T(33)=P3(L1)                                                      P-N09840
      T(34)=EPS(L1)                                                     P-N09850
      var_ab=EPS(L1)**2-P1(L1)**2-P2(L1)**2-P3(L1)**2
      YM(NPION)=0.
      IF(var_ab.GT.0.)
     s                 YM(NPION)=SQRT(var_ab)                           P-N09860
        IF(Kveux.EQ.1) THEN
            DEL1avat(Iavat)=IND1(L1)
	    ENERGYavat(Iavat)=YM(NPION)
	ENDIF
  800 CALL DECAY2(P1(L1),P2(L1),P3(L1),EPS(L1),Q1(NPION),Q2(NPION),Q3(NPP-N09870
     -ION),Q4(NPION),YM(NPION),FMP,FMPI,hel(l1))                        P-N09880

C DECAY
      IF (IND2(L1)*IND2(L1).EQ.9) GO TO 806                             P-N10210
      call ribm(rndm,ial)
 
      IF (RNDM.LT.0.333333333) GO TO 837                                P-N10260
      IPI(NPION)=0                                                      P-N10270
      GO TO 809                                                         P-N10280
  837 IPI(NPION)=IND2(L1)*2                                             P-N10290
      IND2(L1)=-IND2(L1)                                                P-N10300
      GO TO 809                                                         P-N10310
  806 IND2(L1)=IND2(L1)/3                                               P-N10320
      IPI(NPION)=2*IND2(L1)                                             P-N10330
  809 CONTINUE                                                          P-N10340
      IND1(L1)=0                                                        P-N10350
      TLG(L1)=0.                                                        P-N10360
C escape ?
      IF (NESC(L1).GT.0) GO TO 850                                      P-N09890
      ITESTE=0
      CALL PAULI(L1,RBL,PBL,XPB)                                        P-N09900
      call ribm(rndm,iy11) 
C Pauli blocking?
      IF (RNDM.LE.XPB) GO TO 1848                                       P-N09960
C Le decay ne peut conduire a un noyau de A nucleons
C  sous l'energie de Fermi et dans une config. d'energie inferieure a
C  EFER-(IA2-NBALTTF)*TF).
       EGS=0.
       NBALTTF=0
       ITESTE=1
       DO I=1,IA
         IF(NESC(I).EQ.0) THEN
	    IF(SQRT(P1(I)**2+P2(I)**2+P3(I)**2).LT.PF) THEN
	       NBALTTF=NBALTTF+1
	       EGS=EGS+EPS(I)-FMP
	    ENDIF
	 ENDIF
       END DO
       IF(EGS.GE.(EFER-(IA2-NBALTTF)*TF)) GO TO 850
C ATTENTION, logique negative!!! Liberer le goto si on veut supprimer la
C   sequence precedente (CDPP sur Delta-> pi N)
C        GO TO 850
	  IF(Kveux.EQ.1) Bloc_CDPP(Iavat)=1
C It is blocked!      
1848      MPAUL2=MPAUL2+1                                               P-N09970
	  IF(Kveux.EQ.1) Bloc_Paul(Iavat)=1
C largeur variable du delta (Phase Space Factor introduced 4/2001)
C (180.**3 = 5832000.)
      QQQ=SQRT((YM(NPION)**2-(FMP+FMPI)**2)
     s   *(YM(NPION)**2-(FMP-FMPI)**2))/(2.*YM(NPION))
      PSF=QQQ**3/(QQQ**3+5832000.)
      TDEL=HC*T(34)/(G0*PSF*YM(NPION))                                 
      
      IF (ITESTE.EQ.0) TDEL=TDEL*XPB/(1.000001-XPB)                     P-N09990
      IF (TDEL.GT.TMAX5) GO TO 853                                      P-N10000
      K=K+1                                                             P-N10010
      CROIS(K)=TDEL                                                     P-N10020
      IND(K)=L1                                                         P-N10030
      JND(K)=0                                                          P-N10040
  853 P1(L1)=T(31)                                                      P-N10050
      P2(L1)=T(32)                                                      P-N10060
      P3(L1)=T(33)                                                      P-N10070
      EPS(L1)=T(34)                                                     P-N10080
      IND1(L1)=1                                                        P-N10090
      IND2(L1)=ICHD                                                     P-N10100
      NPION=NPION-1                                                     P-N10110
      IF (K.EQ.0) GO TO 230                                             P-N10120
      DO 854 I=1,K                                                      P-N10130
  854 CROIS(I)=CROIS(I)-TAU                                             P-N10140
      GO TO 449                                                         P-N10150

C Valid decay of the Delta
850   CONTINUE
        IF(Kveux.EQ.1) THEN
	 Bloc_Paul(Iavat)=0
	 Bloc_CDPP(Iavat)=0
	ENDIF
		IF (NOSURF.LE.0) THEN
c surface
      pppp=sqrt(p1(l1)**2+p2(l1)**2+p3(l1)**2)
      rrrr=sqrt(x1(l1)**2+x2(l1)**2+x3(l1)**2)
      if (pppp.le.pf) then
C      rcorr=(ri4+pppp*(rs4-ri4)/pf)**0.25
	XV=pppp/pf
	rcorr=FLIN(XV)
      if (rrrr.gt.rcorr) then
c      write(6,*)'ligne 1392', io,rrrr,rcorr
      x1(l1)=x1(l1)*rcorr/rrrr
      x2(l1)=x2(l1)*rcorr/rrrr
      x3(l1)=x3(l1)*rcorr/rrrr
      endif
      endif

		END IF
      NCOL=NCOL+1                                                       P-N10160
      MRDP=MRDP+1                                                       P-N10170
      Y1(NPION)=X1(L1)                                                  P-N10180
      Y2(NPION)=X2(L1)                                                  P-N10190
      Y3(NPION)=X3(L1)                                                  P-N10200

      IF (K.EQ.0) GO TO 4047                                            P-N10370
      KD=0                                                              P-N10380
      CCR=TAU                                                           P-N10390
      DO 700 I=1,K                                                      P-N10400
      I20=I-KD                                                          P-N10410
      IF(IND(I).EQ.L1) GO TO 812                                        P-N10420
      IF(JND(I).EQ.L1) GO TO 812                                        P-N10430
      CROIS(I20)=CROIS(I)-CCR                                           P-N10440
      IND(I20)=IND(I)                                                   P-N10450
      JND(I20)=JND(I)                                                   P-N10460
      GO TO 700                                                         P-N10470
  812 KD=KD+1                                                           P-N10480
  700 CONTINUE                                                          P-N10490
      K=K-KD                                                            P-N10500
 4047 IF (NESC(L1).NE.0) GO TO 845                                      P-N10510
      CALL NEW1(L1)                                                     P-N10520
      K=K+1                                                             P-N10530
      CROIS(K)=REF(X1(L1),X2(L1),X3(L1),P1(L1),P2(L1),P3(L1),EPS(L1),R22P-N10540
     -)                                                                 P-N10550
      IND(K)=L1                                                         P-N10560
      JND(K)=-1                                                         P-N10570
C                                                                       P-N10580
  814 IF (NPION.LE.1) GO TO 845                                         P-N10590
      N20=NPION-1                                                       P-N10600
      DO 846 K20=1,N20                                                  P-N10610
  846 CALL NEW3(Y1(K20),Y2(K20),Y3(K20),Q1(K20),Q2(K20),Q3(K20),Q4(K20),P-N10620
     -K20,L1)                                                           P-N10630
  845 CALL NEW2(Y1(NPION),Y2(NPION),Y3(NPION),Q1(NPION),Q2(NPION),Q3(NPIP-N10640
     -ON),Q4(NPION),NPION,L1)                                           P-N10650
      IF(K.EQ.0) GO TO 230                                              P-N10660
      GO TO 449                                                         P-N10670
C                                                                       P-N10680
CCC   PION-NUCLEON COLLISION                                            P-N10690
C                                                                       P-N10700
  801 LP=L1-IA                                                          P-N10710
      DIS1=X1(L2)-Y1(LP)+(P1(L2)/EPS(L2)-Q1(LP)/Q4(LP))*TAU             P-N10720
      DIS2=X2(L2)-Y2(LP)+(P2(L2)/EPS(L2)-Q2(LP)/Q4(LP))*TAU             P-N10730
      DIS3=X3(L2)-Y3(LP)+(P3(L2)/EPS(L2)-Q3(LP)/Q4(LP))*TAU             P-N10740
      DIST=DIS1*DIS1+DIS2*DIS2+DIS3*DIS3                                P-N10750
      T(10)=EPS(L2)+Q4(LP)                                              P-N10760
      T0=1./T(10)                                                       P-N10770
      B1=(P1(L2)+Q1(LP))*T0                                             P-N10780
      B2=(P2(L2)+Q2(LP))*T0                                             P-N10790
      B3=(P3(L2)+Q3(LP))*T0                                             P-N10800
      S=(1.-B1*B1-B2*B2-B3*B3)*T(10)*T(10)                              P-N10810
      SQ=SQRT(S)                                                        P-N10820
      CG=4+IND2(L2)*IPI(LP)                                             P-N10830
C      IF(SQ.LT.1117..OR.SQ.GT.3000.) GO TO 832                          P-N10840
      IF(SQ.GT.3000.) GO TO 832                
      IF(31.41592*DIST.GT.SPN(SQ)*CG/6.) GO TO 832                      P-N10850
      GO TO 220                                                         P-N10860
  832 IF (K.EQ.0) GO TO 230                                             P-N10870
      GO TO 44                                                          P-N10880
  831 call ribm(rndm,iy19)
      GEFF=T(10)/SQ                                                     P-N10940
      GG=G0                                                             P-N10950
      IF (SQ.GT.1500.) GG=200.                                          P-N10960
C largeur variable du delta (Phase Space Factor introduced 4/2001)
C (180.**3 = 5832000.)
      QQQ=SQRT((SQ**2-(FMP+FMPI)**2)*(SQ**2-(FMP-FMPI)**2))
     s     /(2.*SQ)
      PSF=QQQ**3/(QQQ**3+5832000.)
      TDEL=-HC/(GG*PSF)*ALOG(RNDM)*GEFF         
      
      IND1(L2)=1                                                        P-N10980
      IND2(L2)=IND2(L2)+IPI(LP)                                         P-N10990
      NC(L2)=NC(L2)+1                                                   P-N11000
      EPS(L2)=T(10)                                                     P-N11010
      P1(L2)=P1(L2)+Q1(LP)                                              P-N11020
      P2(L2)=P2(L2)+Q2(LP)                                              P-N11030
      P3(L2)=P3(L2)+Q3(LP)                                              P-N11040
C Ce nucleon (ici delta) devient un participant:
      JPARTICIP(L2)=1

		IF (NOSURF.LE.0) THEN
c surface
      pppp=sqrt(p1(l2)**2+p2(l2)**2+p3(l2)**2)
      rrrr=sqrt(x1(l2)**2+x2(l2)**2+x3(l2)**2)
      if (pppp.le.pf) then
C      rcorr=(ri4+pppp*(rs4-ri4)/pf)**0.25
	XV=pppp/pf
	rcorr=FLIN(XV)
      if (rrrr.gt.rcorr) then
c      write(6,*) 'ligne 1489',io,rrrr,rcorr
      x1(l2)=x1(l2)*rcorr/rrrr
      x2(l2)=x2(l2)*rcorr/rrrr
      x3(l2)=x3(l2)*rcorr/rrrr
      endif
      endif
c fin surface
		END IF
C                                                                       P-N11050
COMM  COMMENT :DIFFERENCE WITH STANDARD CASCADE :                       P-N11060
COMM     THE DELTA IS LOCATED AT THE NUCLEON SITE TO AVOID PROBLEMS     P-N11070
COMM     WITH THE REFLEXION ON THE WALL                                 P-N11080
C                                                                       P-N11090
      IF(LP.EQ.NPION) GO TO 701                                         P-N11100
      LP1=LP+1                                                          P-N11110
      DO 702 I10=LP1,NPION                                              P-N11120
      IPI(I10-1)=IPI(I10)                                               P-N11130
      YM(I10-1)=YM(I10)                                                 P-N11140
      Q1(I10-1)=Q1(I10)                                                 P-N11150
      Q2(I10-1)=Q2(I10)                                                 P-N11160
      Q3(I10-1)=Q3(I10)                                                 P-N11170
      Q4(I10-1)=Q4(I10)                                                 P-N11180
      Y1(I10-1)=Y1(I10)                                                 P-N11190
      Y2(I10-1)=Y2(I10)                                                 P-N11200
  702 Y3(I10-1)=Y3(I10)                                                 P-N11210
  701 NPION=NPION-1                                                     P-N11220
      NCOL=NCOL+1                                                       P-N11230
      MRPD=MRPD+1                                                       P-N11240
      IF (K.EQ.0) GO TO 704                                             P-N11250
      KD=0                                                              P-N11260
      CCR=TAU                                                           P-N11270
      DO 810 I=1,K                                                      P-N11280
      I20=I-KD                                                          P-N11290
      IF (IND(I).EQ.L1) GO TO 59                                        P-N11300
      IF (IND(I).EQ.L2) GO TO 59                                        P-N11310
      IF (JND(I).EQ.L2) GO TO 59                                        P-N11320
      IF (JND(I).EQ.L1) GO TO 59                                        P-N11330
      CROIS(I20)=CROIS(I)-CCR                                           P-N11340
      IND(I20)=IND(I)                                                   P-N11350
      JND(I20)=JND(I)                                                   P-N11360
      GO TO 810                                                         P-N11370
   59 KD=KD+1                                                           P-N11380
  810 CONTINUE                                                          P-N11390
      K=K-KD                                                            P-N11400
      DO 703 I10=1,K                                                    P-N11410
      IF(IND(I10).LE.L1) GO TO 703                                      P-N11420
      IND(I10)=IND(I10)-1                                               P-N11430
  703 CONTINUE                                                          P-N11440
  704 CALL NEW1(L2)                                                     P-N11450
      IF(TDEL.GT.TMAX5) GO TO 447                                       P-N11460
      K=K+1                                                             P-N11470
      CROIS(K)=TDEL                                                     P-N11480
      IND(K)=L2                                                         P-N11490
      JND(K)=0                                                          P-N11500
  447 K=K+1                                                             P-N11510
      CROIS(K)=REF(X1(L2),X2(L2),X3(L2),P1(L2),P2(L2),P3(L2),EPS(L2),R22P-N11520
     -)                                                                 P-N11530
      IND(K)=L2                                                         P-N11540
      JND(K)=-1                                                         P-N11550
      GO TO 449                                                         P-N11560
C                                                                       P-N11570
CCC   REFLECTION ON OR TRANSMISSION THROUGH THE POTENTIAL WALL          P-N11580
C                                                                       P-N11590

  600 CONTINUE
c deutons pas bien compris ici CV ?
      if (npproj(l1).eq.0) go to 608
      if (ind1(l1).ne.0) write (6,*) 'wrong reentering particle'
      if (x1(l1)*p1(l1)+x2(l1)*p2(l1)+x3(l1)*p3(l1).gt.0.)
     -write (6,*) 'wrong reentering particle'
      var_ab=p1(l1)**2+p2(l1)**2+p3(l1)**2
            GPSG=0.
      IF (var_ab.GT.0.)
     s      GPSG=SQRT(((EPS(L1)+v0)**2-PM2)/var_ab)    
      P1(L1)=GPSG*P1(L1)                                                
      P2(L1)=GPSG*P2(L1)                                               
      P3(L1)=GPSG*P3(L1)                                               
      EPS(L1)=EPS(L1)+v0
      npproj(l1)=0
      nesc(l1)=0
c      write (6,*) 'rentree'
c      go to 607    
c reevaluation of the times tab after entrance of 2nd,..nucleon 
c of the projectile (goto 602 instead of 607 modif. 13/06/01)
      go to 602    
c deutons
C Pour un non participant la transmission est impossible:
608   CONTINUE   
        IF(Kveux.EQ.1) THEN
		DEL1avat(Iavat)=IND1(L1)
		ENERGYavat(Iavat)=EPS(L1)-FMP
	ENDIF
      IF(JPARTICIP(L1).EQ.0) GO TO 601
  	IF(Kveux.EQ.1) GO_OUT(Iavat)=1
      IF (IND1(L1).EQ.0) GO TO 605                                      P-N11600
      FM=AM(P1(L1),P2(L1),P3(L1),EPS(L1))                               P-N11610
      POT=V1                                                            P-N11620
      GO TO 606                                                         P-N11630
  605 FM=FMP                                                            P-N11640
      POT=V0                                                            P-N11650
  606 TP=BARR(EPS(L1)-FM,IND2(L1),ITCH,R2,V0)                           P-N11660
	IF(Kveux.EQ.1) ENERGYavat(Iavat)=EPS(L1)-FM
      call ribm(rndm,iy11)                                              P-N11710


C      enerbid=EPS(L1)-FM
C      IF(enerbid.GE.45.) THEN
C      IF(enerbid.LE.65..AND.IND2(L1).GT.0) THEN
C      WRITE(6,*) 'ener,BC,R2,V0,POT',enerbid,TP,R2,V0,POT
C      ENDIF
C      ENDIF
      
      IF (RNDM.GT.TP) GO TO 601                                         P-N11720
	
C-------------------------------------------------------------------
C Ici la particule L1 s'echappe du noyau:
      NESC(L1)=1                                                        P-N11730
      nbquit=nbquit+1
      ITCH=ITCH-(1+IND2(L1))/2                                          P-N11740
      var_ab=p1(l1)**2+p2(l1)**2+p3(l1)**2
             GPSG=0.
      IF(var_ab.GT.0.)
     s       GPSG=SQRT(((EPS(L1)-POT)**2-FM*FM)/var_ab)                 P-N11750
      P1(L1)=GPSG*P1(L1)                                                P-N11760
      P2(L1)=GPSG*P2(L1)                                                P-N11770
      P3(L1)=GPSG*P3(L1)                                                P-N11780
      EPS(L1)=EPS(L1)-POT                                               P-N11790

C Comptage des particules hors du noyau (7/6/2002):
C (remnant minimum=1 nucleon)
	IF(nbquit.GE.(IA-1)) GO TO 255
      GO TO 602                                                         P-N11800
C Here no transmission possible
  601 PSPR=X1(L1)*P1(L1)+X2(L1)*P2(L1)+X3(L1)*P3(L1)                    P-N11810
        IF(Kveux.EQ.1) GO_OUT(Iavat)=0
C Surface: Modif A.B. pour tenir compte du rayon variable du noyau.
C        (X2COUR remplace R22 le rayon**2 fixe du noyau)
      X2COUR=X1(L1)**2+X2(L1)**2+X3(L1)**2
      P1(L1)=P1(L1)-2.*X1(L1)*PSPR/X2COUR        
      P2(L1)=P2(L1)-2.*X2(L1)*PSPR/X2COUR             
      P3(L1)=P3(L1)-2.*X3(L1)*PSPR/X2COUR
C Fin modif surface A.B.             
  602 IF (K.EQ.0) GO TO 607                                             P-N11850
      KD=0                                                              P-N11860
      CCR=TAU                                                           P-N11870
      DO 603 I=1,K                                                      P-N11880
      I20=I-KD                                                          P-N11890
      IF(JND(I).EQ.L1) GO TO 604                                        P-N11900
      IF(IND(I).EQ.L1.AND.JND(I).NE.0) GO TO 604                        P-N11910
      CROIS(I20)=CROIS(I)-CCR                                           P-N11920
      IND(I20)=IND(I)                                                   P-N11930
      JND(I20)=JND(I)                                                   P-N11940
      GO TO 603                                                         P-N11950
  604 KD=KD+1                                                           P-N11960
  603 CONTINUE                                                          P-N11970
      K=K-KD                                                            P-N11980
  607 IF (NESC(L1).EQ.1) GO TO 613                                      P-N11990
      CALL NEW1(L1)                                                     P-N12000
      TREF=REF(X1(L1),X2(L1),X3(L1),P1(L1),P2(L1),P3(L1),EPS(L1),R22)   P-N12010
      IF (TREF.GT.TMAX5) GO TO 615                                      P-N12020
      K=K+1                                                             P-N12030
      CROIS(K)=TREF                                                     P-N12040
      IND(K)=L1                                                         P-N12050
      JND(K)=-1                                                         P-N12060
  615 IF (NPION.EQ.0) GO TO 613                                         P-N12070
      IF (IND1(L1).EQ.1) GO TO 613                                      P-N12080
      DO 614 K20=1,NPION                                                P-N12090
  614 CALL NEW3(Y1(K20),Y2(K20),Y3(K20),Q1(K20),Q2(K20),Q3(K20),Q4(K20),P-N12100
     -K20,L1)                                                           P-N12110
  613 IF (K.EQ.0) GO TO 230                                             P-N12120
      GO TO 449                                                         P-N12130
C                                                                       P-N12140
CCC   DECAY OF THE SURVIVING DELTAS                                     P-N12150
C                                                                       P-N12160
  230 CONTINUE                                                          P-N12170
  255 CONTINUE                                                          P-N12200
      IF (K3.EQ.1) GO TO 256                                            P-N12210
      IF (K4.EQ.0) GO TO 256                                            P-N12220    
      NPIDIR=NPION                                                      P-N12230
      DO 257 I=1,IA                                                     P-N12240
      IF (IND1(I).EQ.0) GO TO 257                                       P-N12250
      NPION=NPION+1                                                     P-N12260
      var_ab=EPS(I)**2-P1(I)**2-P2(I)**2-P3(I)**2
            YM(NPION)=0.
      IF(var_ab.GT.0.)
     s      YM(NPION)=SQRT(var_ab)                                      P-N12270
      xy1=p1(i)
      xy2=p2(i)
      xy3=p3(i)
      xye=eps(i)
        IF(Kveux.EQ.1) THEN
                Iavat=Iavat+1
		TIMEavat(Iavat)=TIM
		L1avat(Iavat)=I
		L2avat(Iavat)=-2
		ENERGYavat(Iavat)=YM(NPION)
		Bloc_Paul(Iavat)=0
		Bloc_CDPP(Iavat)=0
		DEL1avat(Iavat)=IND1(L1)
	        JPARTL1(Iavat)=1
	        JPARTL2(Iavat)=0
	ENDIF
      CALL DECAY2(P1(I),P2(I),P3(I),EPS(I),Q1(NPION),Q2(NPION),Q3(NPION)P-N12280
     -,Q4(NPION),YM(NPION),FMP,FMPI,hel(i))                             P-N12290
      if(nesc(i).eq.0) idecf=1
		IF (NOSURF.LE.0) THEN
c surface
      if (nesc(i).eq.0.) then
      pppp=sqrt(p1(i)**2+p2(i)**2+p3(i)**2)
      rrrr=sqrt(x1(i)**2+x2(i)**2+x3(i)**2)
      if (pppp.le.pf) then
	XV=pppp/pf
	rcorr=FLIN(XV)
      if (rrrr.gt.rcorr) then
      x1(i)=x1(i)*rcorr/rrrr
      x2(i)=x2(i)*rcorr/rrrr
      x3(i)=x3(i)*rcorr/rrrr
      endif
      endif
      endif
c fin surface
		END IF
      IF (IND2( I)*IND2( I).EQ.9) GO TO 280                             P-N12300      
      call ribm(rndm,ial)
 
      IF (RNDM*3..LT.1.) GO TO 283                                      P-N12350
      IPI(NPION)=0                                                      P-N12360
      GO TO 285                                                         P-N12370
  283 IPI(NPION)=IND2( I)*2                                             P-N12380
      IND2( I)=-IND2(I)                                                 P-N12390
      GO TO 285                                                         P-N12400
  280 IND2( I)=IND2( I)/3                                               P-N12410
      IPI(NPION)=2*IND2( I)                                             P-N12420
  285 Y1(NPION)=X1(I)                                                   P-N12430
      Y2(NPION)=X2(I)                                                   P-N12440
      Y3(NPION)=X3(I)                                                   P-N12450
  257 CONTINUE                                                          P-N12460
        IF(Kveux.EQ.1) THEN
		Bavat=B
		NOPARTavat=NOPART
		NCOLavat=NCOL
		NB_AVAT=Iavat
	ENDIF
C                                                                       P-N12470
CCC   FINAL PROPERTIES OF THE INCOMING NUCLEON AND OF THE REMNANT       P-N12480
C     BEFORE EVAPORATION                                                P-N12490
C                                                                       P-N12500
C Tableau des energies a la fin (avatar.hbk)
        IF(Kveux.EQ.1) THEN
	DO i=1,IA
		IF(NESC(i).EQ.0) then
		IF(JPARTICIP(i).EQ.1) then
			EPSf(i)=EPS(i)
		ELSE
			EPSf(i)=0.
		ENDIF
		ELSE
		EPSf(i)=0.
		ENDIF
	ENDDO
	ENDIF

  256 ELEAD=0.                                                          P-N12510
      LEAD=0                                                            P-N12520
      NPX=0                                                             P-N12530
      EREM=0.                                                           P-N12540
      IZREM=0                                                           P-N12550
      INREM=0                                                           P-N12560
      iarem=0
      RCM1=0.                                                           P-N12570
      RCM2=0.                                                           P-N12580
      RCM3=0.                                                           P-N12590
      PREM1=0.                                                          P-N12600
      PREM2=0.                                                          P-N12610
      PREM3=0.                                                          P-N12620
      POUT1=0.                                                          HE-13640
      POUT2=0.                                                          HE-13650
      POUT3=0.                                                          HE-13660
      EOUT =0.                                                          HE-13670
      cmultn=0.
      if (kindf7.le.2.and.kindf7.gt.0) then
          if (ncol.eq.0.or.nc(1).eq.0) then
           go to 9100
          endif
      else
        if (kindf7.le.5.and.kindf7.gt.0) then
           if (ncol.eq.0) then 
            go to 9100
           endif
        else
C Ici faisceau composite: modif A.B. 2/2002 pour tous les composites:
	nsum_col=0
	DO i=1,ia1
		nsum_col = nsum_col + NC(i)
	ENDDO        
	    if (ncol.eq.0.or.nsum_col.eq.0) then
C            if (ncol.eq.0.or.(nc(1).eq.0.and.nc(2).eq.0)) then

c       if (ncol.eq.0) then
c        write (6,*)'ncol=', ncol,nc(1),timi,b
            go to 9100
            endif
        endif
      endif
      go to 9101
c pour eviter renvoi des resultats du run precedent CV 20/11/98
 9100  iarem=ia2
       izrem=iz2
       esrem=0.
       erecrem=0.
       nopart=-1
c fin ajout CV
        return

 9101 nopart=0
      ekout=0.

      DO 258 I=1,IA                                                     P-N12660
c      IF (EPS(I).LT.ELEAD) GO TO 259                                   P-N12670
c      LEAD=I                                                           P-N12680
c      ELEAD=EPS(I)                                                     P-N12690
  259 IF (NESC(I).EQ.0) GO TO 254                                       P-N12700
      XL1=XL1-X2(I)*P3(I)+X3(I)*P2(I)                                   P-N12710
      XL2=XL2-X3(I)*P1(I)+X1(I)*P3(I)                                   P-N12720
      XL3=XL3-X1(I)*P2(I)+X2(I)*P1(I)                                   P-N12730
c   ici ajout de pout CV le 5/7/95
        POUT1=POUT1+P1(I)                                               HE-13740
        POUT2=POUT2+P2(I)                                               HE-13750
        POUT3=POUT3+P3(I)                                               HE-13760
c      write(6,*)'eout nucleon',eout
       EOUT =EOUT+EPS(I)-fmp                                            HE-13770
c      write(6,*)'eout nucleon 2eme',eout
c
c      KKE=INT((EPS(I)-FMP)/DHE)+1
      IC33=(IND2(I)+3)/2
 
 
      nopart=nopart+1
      kind(nopart)=3-ic33
C     Spectators of composite projectiles (7/2006, AB)
C      if(npproj(i).eq.1) kind(nopart) = -kind(nopart)

      ep(nopart)=eps(i)-fmp
      bmass(nopart)=fmp
      ekout=ekout+ep(nopart)
      ptotl=sqrt(p1(i)**2+p2(i)**2+p3(i)**2)
      alpha(nopart)=p1(i)/ptotl
      beta(nopart)=p2(i)/ptotl
      gam(nopart)=p3(i)/ptotl
 
      go to 258
  254 T(4)=X1(I)*X1(I)+X2(I)*X2(I)+X3(I)*X3(I)                          P-N12780
      EREM=EREM+EPS(I)-FMP                                              P-N12800
      RCM1=RCM1+X1(I)                                                   P-N12840
      RCM2=RCM2+X2(I)                                                   P-N12850
      RCM3=RCM3+X3(I)                                                   P-N12860
      PREM1=PREM1+P1(I)                                                 P-N12870
      PREM2=PREM2+P2(I)                                                 P-N12880
      PREM3=PREM3+P3(I)                                                 P-N12890
c     IF (IND2(I).EQ.-1) INREM=INREM+1                                  P-N12900
c     IF (IND2(I).EQ. 1) IZREM=IZREM+1                                  P-N12910
	izrem=izrem+(1+ind2(i))/2
	iarem=iarem+1
c      IRR=INT(5.*T(4)/R22)+1                                           P-N12920
c      PABS=SQRT(P1(I)*P1(I)+P2(I)*P2(I)+P3(I)*P3(I))                   P-N12930
c      CTH=P3(I)/PABS                                                   P-N12940
c      IPZ=INT((CTH+1.)*3.)+1
c      IPTR=INT(PABS/DPP)+1                                             P-N12960
c      IF (IRR.GT.10.OR.IPZ.GT.8.OR.IPTR.GT.20) GO TO 258               P-N12970
c      NENTR(IRR,IPZ,IPTR,10)=NENTR(IRR,IPZ,IPTR,10)+1                  P-N12980
  258 CONTINUE                                                          P-N12990
c    correction pions 21/3/95 JC
      ichpion=0
      if(npion.ne.0)then
      do ipion=1,npion
      pout1=pout1+q1(ipion)
      pout2=pout2+q2(ipion)
      pout3=pout3+q3(ipion)
      eout=eout+q4(ipion)
      xl1=xl1-y2(ipion)*q3(ipion)+y3(ipion)*q2(ipion)
      xl2=xl2-y3(ipion)*q1(ipion)+y1(ipion)*q3(ipion)
      xl3=xl3-y1(ipion)*q2(ipion)+y2(ipion)*q1(ipion)
      ichpion=ichpion+ipi(ipion)/2
      nopart=nopart+1
      kind(nopart)=4-ipi(ipion)/2
      ptotl=sqrt(q1(ipion)**2+q2(ipion)**2+q3(ipion)**2)
      ep(nopart)=q4(ipion)-fmpi
      bmass(nopart)=fmpi
      ekout=ekout+ep(nopart)
      alpha(nopart)=q1(ipion)/ptotl
      beta(nopart)=q2(ipion)/ptotl
      gam(nopart)=q3(ipion)/ptotl
      enddo
      endif
c  fin correction pions sur impulsion et moment angulaire et charge
c ici ajout de PFREM CV le 5/7/95
      PFREM1=-POUT1                                                     HE-14070
      PFREM2=-POUT2                                                     HE-14080
      PFREM3=PINC-POUT3                                                 HE-14090
c
	inrem=iarem-izrem
      IEJP=IZ2-IZREM                                                    P-N13000
      IEJN=IA2-INREM-IZ2                                                P-N13010
      IREM=INREM+IZREM                                                  P-N13020

C Intrinsic momentum of the remnant (A.B. 05/2001): 
C    momentum of projectile minus momentum of all outgoing particles
C                           minus angular momentum of the remnant computed
C                                  from the sum of all inside nucleons.
C      XL1=XL1-RCM2/IREM*PREM3+RCM3/IREM*PREM2
C      XL2=XL2-RCM3/IREM*PREM1+RCM1/IREM*PREM3
C      XL3=XL3-RCM1/IREM*PREM2+RCM2/IREM*PREM1
C                           Here the remnant momentum is Pin - SIG(pout),
C   and the distance with respect to the barycenter of the actual target
      XL1=XL1-(RCM2/IREM - X2_target)*PFREM3
     s       +(RCM3/IREM - X3_target)*PFREM2
      XL2=XL2-(RCM3/IREM - X3_target)*PFREM1
     s       +(RCM1/IREM - X1_target)*PFREM3
      XL3=XL3-(RCM1/IREM - X1_target)*PFREM2
     s       +(RCM2/IREM - X2_target)*PFREM1
      L=SQRT(XL1*XL1+XL2*XL2+XL3*XL3)/HC + 0.5
      
            
c      AVERL=AVERL+L
c      NAVL1=NAVL1+INT(XL1/HC)
c      NAVL2=NAVL2+INT(XL2/HC)
c      NAVL3=NAVL3+INT(XL3/HC)
      IEJ=IA2-IREM                                                      P-N13110
c      EFIN=EFIN+EPS(1)                                                 P-N13120
c      EFIN1=EFIN1+ELEAD
c      IF (LEAD.NE.1) NEX=NEX+1
c      NCLE=NC(LEAD)+1
c      IHC(NCLE)=IHC(NCLE)+1
c      PLOSX=PLOSX+P1(LEAD)
c      PLOSY=PLOSY+P2(LEAD)
c      PLOSZ=PLOSZ+P3(LEAD)
c      PLOST=PLOST+P1(LEAD)**2+P2(LEAD)**2
c      NPRIM=NPRIM+NC(LEAD)
c      MH1=INT(P1(LEAD)/DH1+25.)
c      IF (MH1.LT.1.OR.MH1.GT.50) GO TO 620
c      IHF1(MH1)=IHF1(MH1)+1
c  620 MH2=INT(P2(LEAD)/DH2+25.)
c      IF (MH2.LT.1.OR.MH2.GT.50) GO TO 621
c      IHF2(MH2)=IHF2(MH2)+1
c  621 MH3=INT(P3(LEAD)/DH3)+1
c      IF(MH3.LT.1.OR.MH3.GT.50) GO TO 622
c      IHF3(MH3)=IHF3(MH3)+1
c  622 MH4=INT(EPS(LEAD)/DH4)+1
c      IF(MH4.LT.1.OR.MH4.GT.50) GO TO 623
c      IHF4(MH4)=IHF4(MH4)+1
  623 CONTINUE
c       IF (IEJP.GE.0.AND.IEJN.GE.0) THEN
c      IHP(2,IEJP+1)=IHP(2,IEJP+1)+1
c      IHP(1,IEJN+1)=IHP(1,IEJN+1)+1
c       ENDIF
c      if(iejn.gt.-9.and.iejp.gt.-9)then
c       if(iejnn.le.91.and.iejpp.le.91)
c     -IHREM(IEJN+9,IEJP+9)=IHREM(IEJN+9,IEJP+9)+1
c      endif
      EH5=EREM-(FLOAT(IREM)/A2)**1.666667*EFER                          P-N13370
      sepa=(ia2-irem)*(v0-tf)
      eh6=eh5
c deutons ajout beproj ?????? On retire beproj (18/06/2002 AB CV)
C      eh5=erem-efer-beproj+(ia2-irem)*tf
      eh5=erem-efer+(ia2-irem)*tf
      IF (EH5.LT.0.) EH5=0.00000001                                     P-N13380

       
      
c      eh5=erem-efer+(ia2-irem)*tf
      xlab=tlab-eout-eh5-sepa
c      if (abs (xlab).gt.1.) write(6,*)'xlab',xlab,eout,eh5,sepa, io
c      sepa=(ia2-irem)*v0-(1.-(FLOAT(IREM)/A2)**1.666667)*EFER
C      IF (EH5.LT.0.) EH5=0.00000001                                     P-N13380
       ecoreh5=0.
       if (iqe.eq.1) then
        eh5=0.00000001
        else
        if (eh5.lt.0.) then
         if (npion.eq.0) then
           nopart=-1
           return
         else
           ecoreh5=-eh5
           eh5=0.000000001
         endif
        endif
       endif
      if (idecf.ne.0.and.eh5.lt.0.) then
       ecoreh5=-eh5
       eh5=0.000001
      endif 
c      if (iqe.eq.1) write(6,*) 'iqe=1',' io ',io
c      if (iqe.eq.1) eh5=0.00000001
c      EREMS=EREMS+EH5
c       r111=sqrt(x1(1)**2+x2(1)**2+x3(1)**2)
c      if (nopart.eq.0) print *,'no part',nc(1),eps(1),r111,r2
      iarem=irem
      pfreml2=pfrem1**2+pfrem2**2+pfrem3**2
      if (pfreml2.gt.1.0e-12) then
      pfreml=sqrt(pfreml2)
      alrem=pfrem1/pfreml
      berem=pfrem2/pfreml
      garem=pfrem3/pfreml
      else
      alrem=0.
      berem=0.
      garem=1.0
      endif
      erecrem=pfreml2/(SQRT(pfreml2+(fmp*iarem)**2) + fmp*iarem)
C      erecrem=pfreml2/(2.*fmp*iarem)
      IF(iarem.EQ.1) erecrem = erecrem + eh5 
c correction recul
      erecg=erecrem+ecoreh5
C Correction energie d'excitation pour une absorption (A.B., C.V. 2/2002)
      ESREM=eh5
      
      if (ekout.lt.0.001) return
c       esrem=esrem-erecrem
c on ote l'instruction precedente car esrem toujours nulle 14/9/99

       if (erecg.gt.0.25) then
           fffc=(ekout-erecg)/ekout
           if (fffc.lt.0.0) fffc=0.
c           nopart=-100
c           return
c           else
           do ipart=1,nopart
           ep(ipart)=ep(ipart)*fffc
           enddo
c           endif
c           eout=eout-erecrem
       endif

c Modif BOUDARD juillet 99 (il faut tenir compte de la renormalisation
C   des energies pour les impulsions.)
           pfrem1=0.
           pfrem2=0.
           pfrem3=pinc
         do ipart=1,nopart
           xmodp=sqrt(ep(ipart)*(2.*bmass(ipart)+ep(ipart)))
           pfrem1=pfrem1-alpha(ipart)*xmodp
           pfrem2=pfrem2-beta(ipart)*xmodp
           pfrem3=pfrem3-gam(ipart)*xmodp
         enddo
c Fin modif A.B.
      pfreml2=pfrem1**2+pfrem2**2+pfrem3**2
      erecrem=pfreml2/(SQRT(pfreml2+(fmp*iarem)**2) + fmp*iarem)
C      erecrem=pfreml2/(2.*fmp*iarem)

      if (pfreml2.gt.1.0e-12) then
      pfreml=sqrt(pfreml2)
      alrem=pfrem1/pfreml
      berem=pfrem2/pfreml
      garem=pfrem3/pfreml
      else
      alrem=0.
      berem=0.
      garem=1.0
      endif
C Fin  modif A.B. pour incl3 
      esrem=eh5     
C If the remnant is a nucleon, NO excitation energy
      IF(IAREM.EQ.1) esrem=0.
	  return
 
   77 continue                                                          P-N17750
        print *,' Ia1 > 1 ! '
 	  return
      END                                                               P-N17760

c
C COLLIS                                                                P-N17770
      SUBROUTINE COLLIS(P1,P2,P3,E1,POUT11,POUT12,POUT13,EOUT1,Q1,Q2,Q3,P-N17780
     -Q4,NP,IP,K2,K3,K4,K5,M1,M2,IS1,IS2)                               P-N17790
      COMMON/BL6/XX10,ISA                                               P-N17800
      COMMON/BL8/RATHR,RAMASS                                           P-N17810
      common/hazard/ial,IY1,IY2,IY3,IY4,IY5,IY6,IY7,IY8,IY9,IY10,
     s               IY11,IY12,IY13,IY14,IY15,IY16,IY17,IY18,IY19
      common/bl9/hel(300),l1,l2
      DIMENSION EX(3),EY(3),EZ(3),QQ(3)                                 P-N17820
      DATA XM,XM2,XMDEL,EZERO,XPI/938.2796,8.8037E5,1232.,1876.6,138./  P-N17830
C      DATA IY1,IY2,IY3,IY4,IY5,IY6,IY7,IY8,IY10,IY11,IY12,IY13/         P-N17840
C     1 12345,22345,32345,42345,52345,62345,72345,82345,34567,47059,21033P-N17850
C     1 12345,22345,32345,42345,52345,62345,72345,82341,34567,47059,21033P-N17850
C     2,32835/                                                           P-N17860
C      data iy9/15637/
      PCM(E,A,C)=0.5*SQRT((E**2-(A+C)**2)*(E**2-(A-C)**2))/E            P-N17870
      ISO=IS1+IS2                                                       P-N17880
      NP=0                                                              P-N17890
      PSQ=P1*P1+P2*P2+P3*P3                                             P-N17900
      PNORM=SQRT(PSQ)                                                   P-N17910
      ECM=E1+EOUT1                                                      P-N17920
      IF(ECM.LT.1925.) GO TO 160                                        P-N17930
      IF (K3.EQ.1) GO TO 17                                             P-N17940
      IF(ECM.LT.2065.+RATHR) GO TO 17                                   P-N17950
      IF(K4-1) 18,10,20                                                 P-N17960
   10 IF (M1+M2-1) 19,14,13                                             P-N17970
   19 IF (ECM-2170.4-RAMASS) 17,17,18                                   P-N17980
   20 IF(M1+M2-1) 18,14,13                                              P-N17990
C                                                                       P-N18000
CCC   TEST ON THE RECOMBINATION POSSIBILITY FOR THE N-DELTA SYSTEM      P-N18010
C                                                                       P-N18020
   14 IF (K5.EQ.0) GO TO 170                                            P-N18030
      call ribm(rndm,iy11)
      S1=SEL(ECM,1,ISO)                                                 P-N18090
  206 IF (M1.EQ.0) GO TO 200                                            P-N18100
      XX10=SQRT(E1*E1-PSQ)                                              P-N18110
      ISA=IS1                                                           P-N18120
      GO TO 201                                                         P-N18130
  200 XX10=SQRT(EOUT1**2-PSQ)                                           P-N18140
      ISA=IS2                                                           P-N18150
  201 CONTINUE                                                          P-N18160
      S=S1+SREC(ECM,XX10,ISO,ISA)                                       P-N18170
      A=(S-S1)/S                                                        P-N18180
      IF (RNDM-A) 170,170,17                                            P-N18190
C                                                                       P-N18200
CCC   TEST FOR THE BEHAVIOUR OF THE DELTA-DELTA SYSTEM                  P-N18210
C                                                                       P-N18220
   13 IF (K5.EQ.0) GO TO 160                                            P-N18230
      GO TO 17                                                          P-N18240
C                                                                       P-N18250
CCC   TEST ON THE INELASTICITY                                          P-N18260
C                                                                       P-N18270
   18 call ribm(rndm,iy1)
 
      S=SEL(ECM,0,ISO)                                                  P-N18330
      A=SPRO(ECM,ISO)                                                   P-N18340
      A=S/(S+A)                                                         P-N18350
      IF(RNDM.GT.A) GO TO 100                                           P-N18360
C                                                                       P-N18370
CCC   ELASTIC SCATTERING                                                P-N18380
CCC      FIT OF THE B PARAMETER IN THE DIFFERENTIAL X-SECTION:          P-N18390
CCC       TAKEN FROM :J.C.,D.L'HOTE, J.VANDERMEULEN,NIMB111(1996)215    P-N18400
CCC       FOR PN :IMPROVEMENT OF THE BACKWARD SCATTERING ACCORDING
CCC          J.C ET AL, PRC56(1997)2431
C                                                                       P-N18430
   17 CONTINUE
      PL=0.5*ECM*SQRT(ECM**2-4.*XM2)/XM                                 P-N18440
      X=0.001*PL                                                        P-N18450
      IF (ISO.EQ.0) GO TO 80                                            P-N18460
   83 IF (PL.GT.2000.) GO TO 81                                         P-N18470
      X=X*X                                                             P-N18480
      X=X**4                                                            P-N18490
      B=5.5E-6*X/(7.7+X)                                                P-N18500
      GO TO 82                                                          P-N18510
   81 B=(5.34+0.67*(X-2.))*1.E-6                                        P-N18520
      GO TO 82                                                          P-N18530
C   80 IF (PL.GT.1600.) GO TO 83                                         P-N18540
   80 IF (PL.LT.800.) THEN   
       B=(7.16-1.63*X)*1.E-6                                            P-N18550
       B=B/(1.+EXP(-(X-0.45)/0.05))                                      P-N18560
      ELSE
       IF (PL.LT.1100.) THEN
         B=(9.87-4.88*X)*1.E-6
       ELSE
         B=(3.68+0.76*X)*1.E-6
       ENDIF
      ENDIF
   82 BTMAX=4.*PSQ*B                                                    P-N18570
      Z=EXP(-BTMAX)                                                     P-N18580
      call ribm(rndm,iy2)
      ranres=rndm 
      Y=1.-RNDM*(1.-Z)                                                  P-N18640
      T=ALOG(Y)/B                                                       P-N18650
      IEXPI=0

      IF (M1+M2.EQ.0.AND.ISO.EQ.0) THEN                                        

      APT=1.                                                       
      IF (PL.GT.800.)THEN
         APT=(800./PL)**2                                          
         CPT=AMAX1(6.23*EXP(-1.79*X),0.3)
         ALPHAC=100.*1.E-6
         AAA=(1+APT)*(1-EXP(-BTMAX))/B
         ARGU=PSQ*ALPHAC
            IF (ARGU.GE.8) THEN
               ARGU=0.
            ELSE ARGU=EXP(-4.*ARGU)
            ENDIF
         AAC=CPT*(1.-ARGU)/ALPHAC
c         write (6,*) aaa,aac
c          AAC=0.
         FRACPN=AAA/(AAC+AAA)
c          fracpn=0.
         call ribm(rndm,iy8)

         IF (RNDM.GT.FRACPN) THEN 
           Z=EXP(-4.*PSQ*ALPHAC)
           IEXPI=1
c             write (6,*) 'iexpi=', iexpi
           Y=1.-RANRES*(1.-Z)                                                 
           T=ALOG(Y)/ALPHAC             
         ENDIF        
      ENDIF
      ENDIF

      CTET=1.+0.5*T/PSQ                                                 P-N18660
      if(abs(ctet).gt.1.)ctet=sign(1.,ctet)      
      STET=SQRT(1.-CTET**2)                                             P-N18670
      call ribm(rndm,iy3)
 
      FI=6.2832*RNDM                                                    P-N18730
      CFI=COS(FI)                                                       P-N18740
      SFI=SIN(FI)                                                       P-N18750
      XX=P1*P1+P2*P2                                                    P-N18760
      ZZ=P3*P3                                                          P-N18770
      IF (XX.LT.ZZ*1.E-8) GO TO 30                                      P-N18780
      YN=SQRT(XX)                                                       P-N18790
      ZN=YN*PNORM                                                       P-N18800
      EZ(1)=P1/PNORM                                                    P-N18810
      EZ(2)=P2/PNORM                                                    P-N18820
      EZ(3)=P3/PNORM                                                    P-N18830
      EX(1)=P2/YN                                                       P-N18840
      EX(2)=-P1/YN                                                      P-N18850
      EX(3)=0.                                                          P-N18860
      EY(1)=P1*P3/ZN                                                    P-N18870
      EY(2)=P2*P3/ZN                                                    P-N18880
      EY(3)=-XX/ZN                                                      P-N18890
      P1=(EX(1)*CFI*STET+EY(1)*SFI*STET+EZ(1)*CTET)*PNORM               P-N18900
      P2=(EX(2)*CFI*STET+EY(2)*SFI*STET+EZ(2)*CTET)*PNORM               P-N18910
      P3=(EX(3)*CFI*STET+EY(3)*SFI*STET+EZ(3)*CTET)*PNORM               P-N18920
      GO TO 31                                                          P-N18930
   30 P1=P3*CFI*STET                                                    P-N18940
      P2=P3*SFI*STET                                                    P-N18950
      P3=P3*CTET                                                        P-N18960
   31 POUT11=-P1                                                        P-N18970
      POUT12=-P2                                                        P-N18980
      POUT13=-P3                                                        P-N18990
C
C%%%  BACKWARD SCATTERING ACCORDING THE PARAMETRIZATION OF REF
C%%%     PRC56(1997)1
C
      IF (M1+M2.EQ.1) GO TO 133                                         P-N19000
      IF (ISO.NE.0) GO TO 133
      call ribm(rndm,iy8)
 
      APT=1.                                                            P-N19050
      IF (PL.GT.800.)THEN
         APT=(800./PL)**2                                               P-N19060
      ENDIF
      IF (IEXPI.EQ.1.OR.RNDM.GT.1./(1.+APT)) THEN
      II=IS1                                                            P-N19080
      IS1=IS2                                                           P-N19090
      IS2=II                                                            P-N19100
      ENDIF
  133 CONTINUE                                                          P-N19110
      RETURN                                                            P-N19120
C                                                                       P-N19130
CCC   DELTA PRODUCTION                                                  P-N19140
C                                                                       P-N19150
C        THE PRODUCTION IS NOT ISOTROPIC IN THIS VERSION                P-N19160
C        IT HAS THE SAME EXP(B*T) STRUCTURE AS THE NN ELASTIC SCATTERINGP-N19170
C          (FORMULA 2.3 OF J.CUGNON ET AL, NUCL PHYS A352(1981)505)     P-N19180
C        PARAMETRIZATION OF B TAKEN FROM REF. PRC56(1997)2431
C                                                                       P-N19190
  100 IF (K4.NE.1) GO TO 101                                            P-N19200
      XMDEL=1232.+RAMASS                                                P-N19210
      GO TO 103                                                         P-N19220
  101 call ribm(rndm,iy10)
 
      Y=TAN(3.1415926*(RNDM-0.5))                                       P-N19280
      X=1232.+0.5*130.*Y+RAMASS                                         P-N19290
      IF (X.LT.XM+XPI+2.) GO TO 101                                     P-N19300
      IF (ECM.LT.X+XM+1.) GO TO 101                                     P-N19310
C
C%%%   GENERATION OF THE DELTA MASS WITH THE PENETRATION FACTOR
C          (SEE PRC56(1997)2431)
      Y=ECM**2
      Q2=(Y-1076.**2)*(Y-800.**2)/Y/4.                                 
      Q3=(SQRT(Q2))**3                                                  
      F3MAX=Q3/(Q3+180.0**3)                                       
      Y=X**2
      Q2=(Y-1076.**2)*(Y-800.**2)/Y/4.                                 
      Q3=(SQRT(Q2))**3                                                  
      F3=Q3/(Q3+180.0**3)
      call ribm(rndm,iy11) 
      if (rndm.gt.f3/f3max) go to 101                          
      XMDEL=X                                                           P-N19320
  103 PIN=PNORM                                                         P-N19330
      PNORM=PCM(ECM,XM,XMDEL)                                           P-N19340
      if (pnorm.le.0) pnorm=0.000001
      INDEX=0                                                           P-N19350
      INDEX2=0                                                          P-N19360
      call ribm(rndm,iy4)
 
      IF (RNDM.LT.0.5) INDEX=1                                          P-N19420
      IF (ISO.EQ.0) THEN
       call ribm(rndm,iy5)
 
       IF (RNDM.LT.0.5) INDEX2=1 
      ENDIF
      call ribm(rndm,iy6)

c      X=(ECM-EZERO)*0.001                                               P-N19480
c      X=(3.65*X)**6                                                     P-N19490
c      B=6.*X/(1.+X)*1.E-6                                               P-N19500
      x=0.001*0.5*ECM*SQRT(ECM**2-4.*XM2)/XM
      if(x.lt.1.4)then
	  b=(5.287/(1.+exp((1.3-x)/0.05)))*1.e-6
	else
	  b=(4.65+0.706*(x-1.4))*1.e-6
      endif
      XKH=2.*B*PIN*PNORM                                                P-N19510
      CTET=1.+ALOG(1.-RNDM*(1.-EXP(-2.*XKH)))/XKH                       P-N19520
      if(abs(ctet).gt.1.)ctet=sign(1.,ctet)
      STET=SQRT(1.-CTET**2)                                             P-N19530
      call ribm(rndm,iy7)
 
      FI=6.2832*RNDM                                                    P-N19590
      CFI=COS(FI)                                                       P-N19600
      SFI=SIN(FI)                                                       P-N19610
c   delta production: correction of the angular distribution 02/09/02

      XX=P1*P1+P2*P2                                                    P-N18760
      ZZ=P3*P3                                                          P-N18770
      IF (XX.LT.ZZ*1.E-8) GO TO 230                                     P-N18780
      YN=SQRT(XX)                                                       P-N18790
      ZN=YN*PIN                                                         P-N18800
      EZ(1)=P1/PIN                                                      P-N18810
      EZ(2)=P2/PIN                                                      P-N18820
      EZ(3)=P3/PIN                                                      P-N18830
      EX(1)=P2/YN                                                       P-N18840
      EX(2)=-P1/YN                                                      P-N18850
      EX(3)=0.                                                          P-N18860
      EY(1)=P1*P3/ZN                                                    P-N18870
      EY(2)=P2*P3/ZN                                                    P-N18880
      EY(3)=-XX/ZN                                                      P-N18890
      XP1=(EX(1)*CFI*STET+EY(1)*SFI*STET+EZ(1)*CTET)*PNORM              P-N18900
      XP2=(EX(2)*CFI*STET+EY(2)*SFI*STET+EZ(2)*CTET)*PNORM              P-N18910
      XP3=(EX(3)*CFI*STET+EY(3)*SFI*STET+EZ(3)*CTET)*PNORM              P-N18920
      GO TO 231
  230 CONTINUE            
      XP1=PNORM*STET*CFI                                                P-N19620
      XP2=PNORM*STET*SFI                                                P-N19630
      XP3=PNORM*CTET                                                    P-N19640
  231 CONTINUE
c     end of correction angular distribution of delta production
      E3=SQRT(XP1*XP1+XP2*XP2+XP3*XP3+XM*XM)                            P-N19650
      IF(K4.NE.0) GO TO 161                                             P-N19660
C                                                                       P-N19670
CCC   DECAY OF THE DELTA PARTICLE (K4=0)                                P-N19680
C                                                                       P-N19690
      NP=1                                                              P-N19700
      IP=0                                                              P-N19710
      QQ(1)=XP1                                                         P-N19720
      QQ(2)=XP2                                                         P-N19730
      QQ(3)=XP3                                                         P-N19740
      QQ4=SQRT(XP1*XP1+XP2*XP2+XP3*XP3+XMDEL*XMDEL)                     P-N19750
      heli=ctet**2
      CALL DECAY2(QQ(1),QQ(2),QQ(3),QQ4,Q1,Q2,Q3,Q4,XMDEL,XM,XPI,heli)  P-N19760
      IF (INDEX.EQ.0) GO TO 149                                         P-N19770
      P1=QQ(1)                                                          P-N19780
      P2=QQ(2)                                                          P-N19790
      P3=QQ(3)                                                          P-N19800
      POUT11=-XP1                                                       P-N19810
      POUT12=-XP2                                                       P-N19820
      POUT13=-XP3                                                       P-N19830
      EOUT1=E3                                                          P-N19840
      GO TO 153                                                         P-N19850
  149 POUT11=QQ(1)                                                      P-N19860
      POUT12=QQ(2)                                                      P-N19870
      POUT13=QQ(3)                                                      P-N19880
      EOUT1=E1                                                          P-N19890
      P1=-XP1                                                           P-N19900
      P2=-XP2                                                           P-N19910
      P3=-XP3                                                           P-N19920
      E1=E3                                                             P-N19930
  153 IF (ISO.EQ.0) GO TO 150                                           P-N19940
      IF (RNDM.GT.0.333333) GO TO 151                                   P-N19950
      IS1=-IS1                                                          P-N19960
      IP=-2*IS1                                                         P-N19970
  151 RETURN                                                            P-N19980
  150 IF (INDEX.EQ.1) RETURN                                            P-N19990
      IF (RNDM.LT.0.5) GO TO 152                                        P-N20000
      IS1=1                                                             P-N20010
      IS2=1                                                             P-N20020
      IP=-2                                                             P-N20030
      RETURN                                                            P-N20040
  152 IS1=-1                                                            P-N20050
      IS2=-1                                                            P-N20060
      IP=2                                                              P-N20070
      RETURN                                                            P-N20080
  160 POUT11=-P1                                                        P-N20090
      POUT12=-P2                                                        P-N20100
      POUT13=-P3                                                        P-N20110
      RETURN                                                            P-N20120

C                                                                       P-N20130
CCC   LONG-LIVED DELTA                                                  P-N20140
C                                                                       P-N20150
  161 IF (INDEX.EQ.1) GO TO 162                                         P-N20160
      P1=XP1                                                            P-N20170
      P2=XP2                                                            P-N20180
      P3=XP3                                                            P-N20190
c      E1=E3                                                             P-N20200
c      EOUT1=ECM-E1                                                      P-N20210
       eout1=e3
       e1=ecm-eout1
c      M2=1                                                              P-N20220
       m1=1

      GO TO 163                                                         P-N20230
  162 P1=-XP1                                                           P-N20240
      P2=-XP2                                                           P-N20250
      P3=-XP3                                                           P-N20260
      EOUT1=E3                                                          P-N20270
      E1=ECM-EOUT1                                                      P-N20280
      M1=1                                                              P-N20290
C  163 IF (ISO.EQ.0) GO TO 160                                           P-N20300
C%%%  SYMMETRIZATION OF CHARGES IN pn -> N DELTA
C%%%     THE TEST ON "INDEX" ABOVE SYMETRIZES THE EXCITATION OF ONE 
C%%%     OF THE NUCLEONS WITH RESPECT TO THE DELTA EXCITATION
C%%%     (SEE NOTE 16/10/97)
C     
  163 IF (ISO.EQ.0) THEN
        IF (INDEX2.EQ.1) THEN
          ISI=IS1
          IS1=IS2
          IS2=ISI
        ENDIF
      GO TO 160
      ENDIF
      hel(l1)=ctet**2
      call ribm(rndm,iy8)
c      call ribm(rndm,iy9)
 
      IF (RNDM.LT.0.25) GO TO 160                                       P-N20350
      IS1=3*IS1*M1-(1-M1)*IS1                                           P-N20360
      IS2=3*IS2*M2-(1-M2)*IS2                                           P-N20370
      GO TO 160                                                         P-N20380
C                                                                       P-N20390
CCC  RECOMBINATION PROCESS                                              P-N20400
C                                                                       P-N20410
  170 PNORM=PCM(ECM,XM,XM)                                              P-N20420
      call ribm(rndm,iy12)
 
      CTET=-1.+2.*RNDM                                                  P-N20480
      if(abs(ctet).gt.1.)ctet=sign(1.,ctet)
      STET=SQRT(1.-CTET*CTET)                                           P-N20490
      call ribm(rndm,iy13)
 
      FI=6.2832*RNDM                                                    P-N20550
      CFI=COS(FI)                                                       P-N20560
      SFI=SIN(FI)                                                       P-N20570
      P1=PNORM*STET*CFI                                                 P-N20580
      P2=PNORM*STET*SFI                                                 P-N20590
      P3=PNORM*CTET                                                     P-N20600
      M1=0                                                              P-N20610
      M2=0                                                              P-N20620
      E1=SQRT(P1*P1+P2*P2+P3*P3+XM*XM)                                  P-N20630
      EOUT1=ECM-E1                                                      P-N20640
      IF (ISO.EQ.0) GO TO 160                                           P-N20650
      IS1=ISO/2                                                         P-N20660
      IS2=ISO/2                                                         P-N20670
      GO TO 160                                                         P-N20680
      END                                                               P-N20690 

C DECAY                                                                 P-N20700
      SUBROUTINE DECAY1(P1,P2,P3,WP,Q1,Q2,Q3,WQ,XI,X1,X2)               P-N20710
C                                                                       P-N20720
CCC   THIS ROUTINE DESCRIBES THE ISOTROPIC DECAY OF A PARTICLE OF MASS  P-N20730
CCC     XI IN 2 PARTICLES OF MASSES X1,X2                               P-N20740
CCC   IN THE INPUT, P1,P2,P3 IS THE MOMENTUM OF PARTICLE XI             P-N20750
CCC   IN THE OUTPUT, P1,P2,P3 IS THE MOMENTUM OF PARTICLE X1 , WHILE    P-N20760
CCC     Q1,Q2,Q3 IS THE MOMENTUM OF PARTICLE X2                         P-N20770
C                                                                       P-N20780
       common/hazard/ial,IY1,IY2,IY3,IY4,IY5,IY6,IY7,IY8,IY9,IY10,
     s               IY11,IY12,IY13,IY14,IY15,IY16,IY17,IY18,IY19

C      DATA IY8,IY9/82345,92345/                                         P-N20790
      PCM(E,A,C)=0.5*SQRT((E**2-(A+C)**2)*(E**2-(A-C)**2))/E            P-N20800
      XE=WP                                                             P-N20810
      B1=P1/XE                                                          P-N20820
      B2=P2/XE                                                          P-N20830
      B3=P3/XE                                                          P-N20840
      XQ=PCM(XI,X1,X2)                                                  P-N20850
      call ribm(rndm,iy8)

      CTET=-1.+2.*RNDM                                                  P-N20910
      if(abs(ctet).gt.1.)ctet=sign(1.,ctet)
      STET=SQRT(1.-CTET**2)                                             P-N20920
      call ribm(rndm,iy9)

      FI=6.2832*RNDM                                                    P-N20980
      CFI=COS(FI)                                                       P-N20990
      SFI=SIN(FI)                                                       P-N21000
      Q1=XQ*STET*CFI                                                    P-N21010
      Q2=XQ*STET*SFI                                                    P-N21020
      Q3=XQ*CTET                                                        P-N21030
      W1=Q1*Q1+Q2*Q2+Q3*Q3                                              P-N21040
      WQ=SQRT(W1+X2*X2)                                                 P-N21050
      P1=-Q1                                                            P-N21060
      P2=-Q2                                                            P-N21070
      P3=-Q3                                                            P-N21080
      WP=SQRT(W1+X1*X1)                                                 P-N21090
      CALL LOREN(Q1,Q2,Q3,B1,B2,B3,WQ)                                  P-N21100
      CALL LOREN(P1,P2,P3,B1,B2,B3,WP)                                  P-N21110
      RETURN                                                            P-N21120
      END                                                               P-N21130
    
C DECAY2                                                               
      SUBROUTINE DECAY2(P1,P2,P3,WP,Q1,Q2,Q3,WQ,XI,X1,X2,hel)            
C                                                                       
CCC   THIS ROUTINE DESCRIBES THE ANISOTROPIC DECAY OF A PARTICLE OF MASS
CCC     XI INTO 2 PARTICLES OF MASSES X1,X2                             
CCC   THE ANISOTROPY IS SUPPOSED TO FOLLOW A 1+3*hel*(cos(theta))**2
CCC   LAW WITH RESPECT TO THE DIRECTION OF THE INCOMING PARTICLE
CCC   IN THE INPUT, P1,P2,P3 IS THE MOMENTUM OF PARTICLE XI             P-N20750
CCC   IN THE OUTPUT, P1,P2,P3 IS THE MOMENTUM OF PARTICLE X1 , WHILE    P-N20760
CCC     Q1,Q2,Q3 IS THE MOMENTUM OF PARTICLE X2                         P-N20770
C                                                                       P-N20780
       common/hazard/ial,IY1,IY2,IY3,IY4,IY5,IY6,IY7,IY8,IY9,IY10,
     s               IY11,IY12,IY13,IY14,IY15,IY16,IY17,IY18,IY19

C      DATA IY8,IY9,IY10/82345,92345,45681/                           
      PCM(E,A,C)=0.5*SQRT((E**2-(A+C)**2)*(E**2-(A-C)**2))/E            P-N20800
      XE=WP                                                             P-N20810
      B1=P1/XE                                                          P-N20820
      B2=P2/XE                                                          P-N20830
      B3=P3/XE                                                          
      XQ=PCM(XI,X1,X2)                                                  
  100 call ribm(rndm,iy8)
      CTET=-1.+2.*RNDM                                                  P-N20910
      if(abs(ctet).gt.1.)ctet=sign(1.,ctet)
      STET=SQRT(1.-CTET**2)                                             P-N20920
      call ribm(rndm,iy10)
      if (rndm.gt.((1.+3.*hel*ctet**2)/(1.+3.*hel))) go to 100
      call ribm(rndm,iy9)      
      FI=6.2832*RNDM                                                    P-N20980
      CFI=COS(FI)                                                       P-N20990
      SFI=SIN(FI)
      beta=sqrt(b1*b1+b2*b2+b3*b3)
      if (beta.lt.1.0e-10) go to 101
      sal=sqrt(b1**2+b2**2)/beta
      cal=b3/beta
      if (sal.lt.1.0e-6) go to 101
      t1=ctet+cal*stet*sfi/sal
      t2=stet/sal                                                       
      q1=xq*(b1*t1+b2*t2*cfi)/beta
      q2=xq*(b2*t1-b1*t2*cfi)/beta
      q3=xq*(b3*t1/beta-t2*sfi)
      go to 102
  101 Q1=XQ*STET*CFI                                                    
      Q2=XQ*STET*SFI                                                    
      Q3=XQ*CTET
  102 hel=0.                                                       
      W1=Q1*Q1+Q2*Q2+Q3*Q3                                              P-N21040
      WQ=SQRT(W1+X2*X2)                                                 P-N21050
      P1=-Q1                                                            P-N21060
      P2=-Q2                                                            P-N21070
      P3=-Q3                                                            P-N21080
      WP=SQRT(W1+X1*X1)                                                 P-N21090
      CALL LOREN(Q1,Q2,Q3,B1,B2,B3,WQ)                                  P-N21100
      CALL LOREN(P1,P2,P3,B1,B2,B3,WP)                                  P-N21110
      RETURN                                                            P-N21120
      END                                                               P-N21130
C SEL                                                                   P-N21140
      FUNCTION SEL(E,M,I)                                               P-N21150
C         FIT BY J.VANDERMEULEN                                         P-N21160
C         LOW ENRGY FIT OF J.C., D. L'HOTE, J. VDM, NIM B111(1996)215
C         I=2,0,-2  FOR PP,PN,NN                                        P-N21170
C         M=0,1,2 FOR NUCLEON-NUCLEON,NUCLEON-DELTA,DELTA,DELTA         P-N21180
C                                                                       P-N21190
      SCALE=1.                                                          P-N21200
      PLAB=E*SQRT(E*E-3.52E6)/1876.6                                    P-N21210
      P1=0.001*PLAB                                                     P-N21220
      IF(PLAB.GT.2000.) GO TO 13                                        P-N21230
      IF (M-1) 1,3,3                                                    P-N21240
    1 IF (I.EQ.0) GO TO 2                                               P-N21250
    3 IF (PLAB.LT.800.) GO TO 4                                         P-N21260
      IF (PLAB.GT.2000.) GO TO 13                                       P-N21270
      SEL=(1250./(50.+P1)-4.*(P1-1.3)**2)*SCALE                         P-N21280
      RETURN
    4 IF (PLAB.LT.440.) THEN                                            P-N21290
       SEL=34.*(P1/0.4)**(-2.104)
      ELSE
       SEL=(23.5+1000.*(P1-0.7)**4)*SCALE                               P-N21300
      ENDIF
      RETURN                                                            P-N21310
   13 SEL=77./(P1+1.5)*SCALE                                            P-N21320
      RETURN                                                            P-N21330
    2 IF (PLAB.LT.800.) GO TO 11                                        P-N21340
      IF (PLAB.GT.2000.) GO TO 13                                       P-N21350
      SEL=31./SQRT(P1)*SCALE                                            P-N21360
      RETURN                                                            P-N21370
   11 IF (PLAB.LT.450.) THEN
       ALP=ALOG(P1)
       SEL=6.3555*EXP(-3.2481*ALP-0.377*ALP*ALP)
      ELSE
       SEL=(33.+196.*SQRT(ABS(P1-0.95)**5))*SCALE                       P-N21380
      ENDIF
      RETURN                                                            P-N21390
      END                                                               P-N21400
C STOT                                                                  P-N21410
      FUNCTION STOT(E,M,I)                                              P-N21420
C        TOTAL CROSS-SECTIONS                                           P-N21430
C         I=2,0,-2  FOR PP,PN,NN                                        P-N21440
C         M=0,1,2 FOR NUCLEON-NUCLEON,NUCLEON-DELTA,DELTA,DELTA         P-N21450
C                                                                       P-N21460
      COMMON/BL6/XX10,ISA                                               P-N21470
      IF (M-1) 1,2,3                                                    P-N21480
    1 SINE=SPRO(E,I)                                                    P-N21490
      GO TO 4                                                           P-N21500
    2 SINE=SREC(E,XX10,I,ISA)                                           P-N21510
      GO TO 4                                                           P-N21520
    3 SINE=0.                                                           P-N21530
    4 STOT=SINE+SEL(E,M,I)                                              P-N21540
      call hfill(301, STOT, 0.0, 1.0)
      if (m.eq.0) then
         if (i.eq.2) then
            call hfill(305, E, STOT, 1.0)
         endif
         if (i.eq.0) then
            call hfill(306, E, STOT, 1.0)
         endif
         if (i.eq.-2) then
            call hfill(307, E, STOT, 1.0)
         endif
      endif
      if (m.eq.1) then
         call hfill(303, E, STOT, 1.0)
      endif
      if (m.eq.2) then
         call hfill(304, E, STOT, 1.0)
      endif
      RETURN                                                            P-N21550
      END                                                               P-N21560
C SREC                                                                  P-N21570
      FUNCTION SREC(E,D,I,ISA)                                          P-N21580
      common/kind/kindf7
      ein=e
      IF (I*I.EQ.16) GO TO 2                                            P-N21590
      IF(E.GT.938.3+D) GO TO 1                                          P-N21600
    2 SREC=0.                                                           P-N21610
      RETURN                                                            P-N21620
    1 IF(E.LT.938.3+D+2.) E=938.3+D+2.                                  P-N21630
      S=E*E                                                             P-N21640
      X=(S-3.523E6)/(S-(938.3+D)**2)                                    P-N21650
      Y=S/(S-(D-938.3)**2)                                              P-N21660
      SREC=0.5*X*Y*SPRO(E,I)                                            P-N21670
      SREC=SREC*(32.+I*I*(ISA*ISA-5))/64.                               P-N21680
      SREC=SREC/(1.+0.25*I*I)                                           P-N21690
c++  modification for pion-induced cascade (see JC and MC LEMAIRE,NPA489(88)781
C      if (kindf7.gt.2) srec=3.*srec
      srec=3.*srec  !pi absorption increased also for internal pions (7/3/01)
      e=ein
      RETURN                                                            P-N21700
      END                                                               P-N21710
C SPRO                                                                  P-N21720
      FUNCTION SPRO(E,I)                                                P-N21730
C         DELTA PRODUCTION CROSS-SECTIONS                               P-N21740
C         FIT BY J.VANDERMEULEN                                         P-N21750
C         I=2,0,-2  FOR PP,PN,NN                                        P-N21760
C                                                                       P-N21770
      COMMON/BL8/RATHR,RAMASS  
      SCALI=1.                                                          P-N21780
      EE=E-RATHR                                                        P-N21790
      IF (EE*EE-3.53E6.LT.0.) GO TO 22                                  P-N21800
      PLAB=EE*SQRT(EE*EE-3.52E6)/1876.6                                 P-N21810
      P1=0.001*PLAB                                                     P-N21820
      IF (PLAB.GT.800.) GO TO 1                                         P-N21830
   22 SPRO=0.                                                           P-N21840
      RETURN                                                            P-N21850
    1 IF (I*I.EQ.4) GO TO 10                                            P-N21860
      IF (PLAB.LT.2000.) GO TO 2                                        P-N21870
      SPRO=(42.-77./(P1+1.5))*SCALI                                     P-N21880
      RETURN                                                            P-N21890
    2 IF (PLAB.LT.1000.) GO TO 3                                        P-N21900
      SPRO=(24.2+8.9*P1-31.1/SQRT(P1))*SCALI                            P-N21910
      RETURN                                                            P-N21920
    3 SPRO=(33.+196.*SQRT(ABS(P1-0.95)**5)-31.1/SQRT(P1))*SCALI         P-N21930
      RETURN                                                            P-N21940
   10 IF (PLAB.LT.2000.) GO TO 11                                       P-N21950
      SPRO=(41.+(60.*P1-54.)*EXP(-1.2*P1)-77./(P1+1.5))*SCALI           P-N21960
      RETURN                                                            P-N21970
   11 IF (PLAB.LT.1500.) GO TO 12                                       P-N21980
      SPRO=41.+60.*(P1-0.9)*EXP(-1.2*P1)-1250./(P1+50.)+4.*(P1-1.3)**2  P-N21990
      SPRO=SPRO*SCALI                                                   P-N22000
      RETURN                                                            P-N22010
   12 SPRO=23.5+24.6/(1.+EXP(-10.*P1+12.))-1250./(P1+50.)+4.*(P1-1.3)**2P-N22020
      SPRO=SPRO*SCALI                                                   P-N22030
      RETURN                                                            P-N22040
      END                                                               P-N22050
C SPN                                                                   P-N22060
      FUNCTION SPN(X)                                                   P-N22070
C                                                                       P-N22080
C      SIGMA(PI+ + P) IN THE (3,3) REGION                               P-N22090
C       NEW FIT BY J.VANDERMEULEN + CONSTANT VALUE ABOVE THE (3,3)      P-N22100
C                                                RESONANCE              P-N22110
C                                                                       P-N22120
      COMMON/BL8/RATHR,RAMASS                                           REL21800
      Y=X*X                                                             P-N22130
      Q2=(Y-1076.**2)*(Y-800.**2)/Y/4.                                  P-N22140
      IF (Q2.GT.0.) GO TO 1                                             P-N22150
      SPN=0.                                                            P-N22160
      RETURN                                                            P-N22170
    1 Q3=(SQRT(Q2))**3                                                  P-N22180
      F3=Q3/(Q3+180.0**3)                                               P-N22190
      SPN=326.5/(((X-1215.-RAMASS)*2./110.)**2+1.)                      P-N22200
      SPN=SPN*(1.-5.*RAMASS/1215.)
      SPN=SPN*F3
C      SPN=SPN*3     ! multiplication by 3 suppressed  june/2001
      RETURN                                                            P-N22210
      END                                                               P-N22220
C TIME 
      SUBROUTINE TIME (I,J)                                             P-N22240
      DIMENSION T(10)                                                   P-N22250
      COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TAP-N22260
      COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N22270
      T(1)=P1(I)/EPS(I)-P1(J)/EPS(J)                                    P-N22280
      T(2)=P2(I)/EPS(I)-P2(J)/EPS(J)                                    P-N22290
      T(3)=P3(I)/EPS(I)-P3(J)/EPS(J)                                    P-N22300
      T(4)=X1(I)-X1(J)                                                  P-N22310
      T(5)=X2(I)-X2(J)                                                  P-N22320
      T(6)=X3(I)-X3(J)                                                  P-N22330
      T(7)=T(1)*T(4)+T(2)*T(5)+T(3)*T(6)                                P-N22340
      T(10)=T(1)*T(1)+T(2)*T(2)+T(3)*T(3)                               P-N22350
        IF(T(10).LE.1.e-10) THEN
		TA=100000
        ELSE	
                TA=-T(7)/T(10)                                          P-N22360
        ENDIF
      RAB2=T(4)*T(4)+T(5)*T(5)+T(6)*T(6)+TA*T(7)                        P-N22370
      RETURN                                                            P-N22380
      END                                                               P-N22390
C NEWT                                                                  P-N22400
      SUBROUTINE NEWT(L1,L2)                                            P-N22410
      DIMENSION IND(20000),JND(20000)                                   P-N22420
      COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TAP-N22430
      COMMON/BL2/CROIS(19900),K,IND,JND                                 P-N22440
      COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N22450
      COMMON/BL4/TMAX5                                                  P-N22460
      COMMON/BL5/TLG(300),NESC(300)                                     P-N22470
      COMMON/BL6/XX10,ISA                                               P-N22480

      EXTERNAL TIME

      AM(A,B,C,D)=SQRT(D*D-A*A-B*B-C*C)                                 P-N22490
      IA=IA1+IA2                                                        P-N22500
      DO 52 I=1,IA                                                      P-N22510
      IF (NESC(I).NE.0) GO TO 52                                        P-N22520
      IF (I-L2) 53,52,54                                                P-N22530
   53 IG=L2                                                             P-N22540
      ID=I                                                              P-N22550
      KG=L1                                                             P-N22560
      KD=I                                                              P-N22570
      GO TO 55                                                          P-N22580
   54 IF (I-L1) 56,52,57                                                P-N22590
   56 KG=L1                                                             P-N22600
      KD=I                                                              P-N22610
      GO TO 58                                                          P-N22620
   57 KG=I                                                              P-N22630
      KD=L1                                                             P-N22640
   58 IG=I                                                              P-N22650
      ID=L2                                                             P-N22660
   55 CALL TIME(IG,ID)                                                  P-N22670
      IF (TA.LT.0.) GO TO 50                                            P-N22680
      IF(TA.GT.TMAX5) GO TO 50                                          P-N22690
      IF (TA.LT.TLG(L2)) GO TO 50                                       P-N22700
      IF (IND1(IG)+IND1(ID).GT.0) GO TO 60                              P-N22710
      E=AM(P1(IG)+P1(ID),P2(IG)+P2(ID),P3(IG)+P3(ID),EPS(IG)+EPS(ID))   P-N22720
      IF (E.LT.1925.) GO TO 50                                          P-N22730
      IY=IND1(IG)+IND1(ID)                                              P-N22740
      IF (IY.NE.1) GO TO 61                                             P-N22750
      IX=IG*IND1(IG)+ID*IND1(ID)                                        P-N22760
      XX10=AM(P1(IX),P2(IX),P3(IX),EPS(IX))                             P-N22770
      ISA=IND2(IX)                                                      P-N22780
   61 IF (31.*RAB2.GT.STOT(E,IY,IND2(IG)+IND2(ID))) GO TO 50            P-N22790
   60 CONTINUE                                                          P-N22800
      K=K+1                                                             P-N22810
      CROIS(K)=TA                                                       P-N22820
      IND(K)=IG                                                         P-N22830
      JND(K)=ID                                                         P-N22840
   50 CALL TIME (KG,KD)                                                 P-N22850
      IF (TA.LT.0.) GO TO 52                                            P-N22860
      IF(TA.GT.TMAX5) GO TO 52                                          P-N22870
      IF (TA.LT.TLG(L1)) GO TO 52                                       P-N22880
      IF (IND1(KG)+IND1(KD).GT.0) GO TO 62                              P-N22890
      E=AM(P1(KG)+P1(KD),P2(KG)+P2(KD),P3(KG)+P3(KD),EPS(KG)+EPS(KD))   P-N22900
      IF (E.LT.1925.) GO TO 52                                          P-N22910
      IY=IND1(KG)+IND1(KD)                                              P-N22920
      IF (IY.NE.1) GO TO 63                                             P-N22930
      IX=KG*IND1(KG)+KD*IND1(KD)                                        P-N22940
      XX10=AM(P1(IX),P2(IX),P3(IX),EPS(IX))                             P-N22950
      ISA=IND2(IX)                                                      P-N22960
   63 IF (31.*RAB2.GT.STOT(E,IY,IND2(KG)+IND2(KD))) GO TO 52            P-N22970
   62 CONTINUE                                                          P-N22980
      K=K+1                                                             P-N22990
      CROIS(K)=TA                                                       P-N23000
      IND(K)=KG                                                         P-N23010
      JND(K)=KD                                                         P-N23020
   52 CONTINUE                                                          P-N23030
      RETURN                                                            P-N23040
      END                                                               P-N23050
C   NEW1                                                                P-N23060
      SUBROUTINE NEW1(L1)                                               P-N23070
      DIMENSION IND(20000),JND(20000)                                   P-N23080
      COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TAP-N23090
      COMMON/BL2/CROIS(19900),K,IND,JND                                 P-N23100
      COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N23110
      COMMON/BL4/TMAX5                                                  P-N23120
      COMMON/BL5/TLG(300),NESC(300)                                     P-N23130
      COMMON/BL6/XX10,ISA                                               P-N23140

      EXTERNAL TIME

      AM(A,B,C,D)=SQRT(D*D-A*A-B*B-C*C)                                 P-N23150
      IA=IA1+IA2                                                        P-N23160
      DO 52 I=1,IA                                                      P-N23170
      IF (NESC(I).NE.0) GO TO 52                                        P-N23180
      IF(I-L1) 53,52,54                                                 P-N23190
   53 CALL TIME(I,L1)                                                   P-N23200
      IF(TA.LT.0.) GO TO 52                                             P-N23210
      IF(TA.GT.TMAX5) GO TO 52                                          P-N23220
      IF (IND1(I)+IND1(L1).GT.0) GO TO 60                               P-N23230
      E=AM(P1(I)+P1(L1),P2(I)+P2(L1),P3(I)+P3(L1),EPS(I)+EPS(L1))       P-N23240
      IF (E.LT.1925.) GO TO 52                                          P-N23250
      IY=IND1(I)+IND1(L1)                                               P-N23260
      IF (IY.NE.1) GO TO 61                                             P-N23270
      IX=I*IND1(I)+L1*IND1(L1)                                          P-N23280
      XX10=AM(P1(IX),P2(IX),P3(IX),EPS(IX))                             P-N23290
      ISA=IND2(IX)                                                      P-N23300
   61 IF (31.*RAB2.GT.STOT(E,IY,IND2(I)+IND2(L1))) GO TO 52             P-N23310
   60 CONTINUE                                                          P-N23320
      K=K+1                                                             P-N23330
      CROIS(K)=TA                                                       P-N23340
      IND(K)=L1                                                         P-N23350
      JND(K)=I                                                          P-N23360
      GO TO 52                                                          P-N23370
   54 CALL TIME(I,L1)                                                   P-N23380
      IF(TA.LT.0.) GO TO 52                                             P-N23390
      IF(TA.GT.TMAX5) GO TO 52                                          P-N23400
      IF (IND1(I)+IND1(L1).GT.0) GO TO 70                               P-N23410
      E=AM(P1(I)+P1(L1),P2(I)+P2(L1),P3(I)+P3(L1),EPS(I)+EPS(L1))       P-N23420
      IF (E.LT.1925.) GO TO 52                                          P-N23430
      IY=IND1(I)+IND1(L1)                                               P-N23440
      IF (IY.NE.1) GO TO 71                                             P-N23450
      IX=I*IND1(I)+L1*IND1(L1)                                          P-N23460
      XX10=AM(P1(IX),P2(IX),P3(IX),EPS(IX))                             P-N23470
      ISA=IND2(IX)                                                      P-N23480
   71 IF (31.*RAB2.GT.STOT(E,IY,IND2(I)+IND2(L1))) GO TO 52             P-N23490
   70 CONTINUE                                                          P-N23500
      K=K+1                                                             P-N23510
      CROIS(K)=TA                                                       P-N23520
      IND(K)=I                                                          P-N23530
      JND(K)=L1                                                         P-N23540
   52 CONTINUE                                                          P-N23550
      RETURN                                                            P-N23560
      END                                                               P-N23570
C   NEW2                                                                P-N23580
      SUBROUTINE NEW2(Y1,Y2,Y3,Q1,Q2,Q3,Q4,NPION,L1)                    P-N23590
      DIMENSION IND(20000),JND(20000)                                   P-N23600
      COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TAP-N23610
      COMMON/BL2/CROIS(19900),K,IND,JND                                 P-N23620
      COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N23630
      COMMON/BL4/TMAX5                                                  P-N23640
      COMMON/BL5/TLG(300),NESC(300)                                     P-N23650
      DIMENSION T(10)                                                   P-N23660
      IA=IA1+IA2                                                        P-N23670
      DO 52 I=1,IA                                                      P-N23680
      IF (NESC(I).NE.0) GO TO 52                                        P-N23690
      IF(I.EQ.L1) GO TO 52                                              P-N23700
      IF(IND1(I).EQ.1) GO TO 52                                         P-N23710
      T(1)=P1(I)/EPS(I)-Q1/Q4                                           P-N23720
      T(2)=P2(I)/EPS(I)-Q2/Q4                                           P-N23730
      T(3)=P3(I)/EPS(I)-Q3/Q4                                           P-N23740
      T(4)=X1(I)-Y1                                                     P-N23750
      T(5)=X2(I)-Y2                                                     P-N23760
      T(6)=X3(I)-Y3                                                     P-N23770
      T(7)=T(1)*T(4)+T(2)*T(5)+T(3)*T(6)                                P-N23780
      IF(T(7).GT.0.) GO TO 52                                           P-N23790
      T(10)=T(1)*T(1)+T(2)*T(2)+T(3)*T(3)                               P-N23800
      TA=-T(7)/T(10)                                                    P-N23810
      IF(TA.GT.TMAX5) GO TO 52                                          P-N23820
      XX2=T(4)*T(4)+T(5)*T(5)+T(6)*T(6)+TA*T(7)                         P-N23830
      E=SQRT((EPS(I)+Q4)**2-(P1(I)+Q1)**2-(P2(I)+Q2)**2-(P3(I)+Q3)**2)  P-N23840
      IF (31.*XX2.GT.SPN(E)) GO TO 52                                   P-N23850
      K=K+1                                                             P-N23860
      CROIS(K)=TA                                                       P-N23870
      IND(K)=IA+NPION                                                   P-N23880
      JND(K)=I                                                          P-N23890
   52 CONTINUE                                                          P-N23900
      RETURN                                                            P-N23910
      END                                                               P-N23920
C LOREN                                                                 P-N23930
      SUBROUTINE LOREN(Q1,Q2,Q3,B1,B2,B3,E)                             P-N23940
C                                                                       P-N23950
CCC   TRANSFORMS MOMENTUM Q AND ENERGY E FROM A FRAME MOVING WITH       P-N23960
CCC   VELOCITY BETA                                                     P-N23970
C                                                                       P-N23980
      BB2=B1*B1+B2*B2+B3*B3                                             P-N23990
      BQ=B1*Q1+B2*Q2+B3*Q3                                              P-N24000
      GAM2=1./(1.-BB2)                                                  P-N24010
      GAM=SQRT(GAM2)                                                    P-N24020
      C=GAM2/(GAM+1.)                                                   P-N24030
      G=C*BQ+GAM*E                                                      P-N24040
      E=GAM*(E+BQ)                                                      P-N24050
      Q1=Q1+B1*G                                                        P-N24060
      Q2=Q2+B2*G                                                        P-N24070
      Q3=Q3+B3*G                                                        P-N24080
      RETURN                                                            P-N24090
      END                                                               P-N24100
C NEW3                                                                  P-N24110
      SUBROUTINE NEW3(Y1,Y2,Y3,Q1,Q2,Q3,Q4,NPION,L1)                    P-N24120
      DIMENSION IND(20000),JND(20000)                                   P-N24130
      COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TAP-N24140
      COMMON/BL2/CROIS(19900),K,IND,JND                                 P-N24150
      COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N24160
      COMMON/BL4/TMAX5                                                  P-N24170
      COMMON/BL5/TLG(300),NESC(300)                                     P-N24180
      DIMENSION T(10)                                                   P-N24190
c  skips collisions with escaped particle 22/3/95
      if(nesc(l1).gt.0)go to 52
      T(1)=P1(L1)/EPS(L1)-Q1/Q4                                         P-N24200
      T(2)=P2(L1)/EPS(L1)-Q2/Q4                                         P-N24210
      T(3)=P3(L1)/EPS(L1)-Q3/Q4                                         P-N24220
      T(4)=X1(L1)-Y1                                                    P-N24230
      T(5)=X2(L1)-Y2                                                    P-N24240
      T(6)=X3(L1)-Y3                                                    P-N24250
      T(7)=T(1)*T(4)+T(2)*T(5)+T(3)*T(6)                                P-N24260
      IF(T(7).GT.0.) GO TO 52                                           P-N24270
      T(10)=T(1)*T(1)+T(2)*T(2)+T(3)*T(3)                               P-N24280
      TA=-T(7)/T(10)                                                    P-N24290
      IF(TA.GT.TMAX5) GO TO 52                                          P-N24300
      IF (TA.LT.TLG(L1)) GO TO 52                                       P-N24310
      XX2=T(4)*T(4)+T(5)*T(5)+T(6)*T(6)+TA*T(7)                         P-N24320
      E=SQRT((EPS(L1)+Q4)**2-(P1(L1)+Q1)**2-(P2(L1)+Q2)**2-(P3(L1)+Q3)  P-N24330
     -**2)                                                              P-N24340
      IF (31.*XX2.GT.SPN(E)) GO TO 52                                   P-N24350
      K=K+1                                                             P-N24360
      CROIS(K)=TA                                                       P-N24370
      IA=IA1+IA2                                                        P-N24380
      IND(K)=IA+NPION                                                   P-N24390
      JND(K)=L1                                                         P-N24400
   52 RETURN                                                            P-N24410
      END                                                               P-N24420
C BARR                                                                  P-N24430
      FUNCTION BARR(E,IZ,IZN,R,V0)                                      P-N24440
C                                                                       P-N24450
CCC   BARR=TRANSMISSION PROBABILITY FOR A NUCLEON OF KINETIC ENERGY     P-N24460
CCC         E ON THE EDGE OF THE WELL OF DEPTH V0 (NR APPROXIMATION)    P-N24470
CCC   IZ IS THE ISOSPIN OF THE NUCLEON,IZN THE INSTANTENEOUS CHARGE     P-N24480
CCC    OF THE NUCLEUS AND R IS THE TARGET RADIUS                        P-N24490
C                                                                       P-N24500
      IF (E.GT.V0) GO TO 1                                              P-N24510
      BARR=0.                                                           P-N24520
      RETURN                                                            P-N24530
    1 X=SQRT(E*(E-V0))                                                  P-N24540
      BARR=4.*X/(E+E-V0+X+X)                                            P-N24550
      IF (IZ.GT.0) GO TO 2                                              P-N24560
      RETURN                                                            P-N24570
    2 B=IZN*1.44/R                                                      P-N24580
      PX=SQRT((E-V0)/B)                                                 P-N24590
      IF (PX.LT.1.) GO TO 3                                             P-N24600
      RETURN                                                            P-N24610
    3 G=IZN/137.03*SQRT(2.*938.3/(E-V0))*(ACOS(PX)-PX*SQRT(1.-PX*PX))   P-N24620
      IF (G.GT.35.)THEN                                                 P-N24630
        BARR=0.                                                         P-N24640
       ELSE                                                             P-N24650
        BARR=BARR*EXP(-2.*G)                                            P-N24660
      ENDIF                                                             P-N24670
      RETURN                                                            P-N24680
      END                                                               P-N24690
C REF                                                                   P-N24700
      FUNCTION REF(X1,X2,X3,P1,P2,P3,E,R2)                              P-N24710
      common/bl10/ri4,rs4,r2i,r2s,PDUMMY
      COMMON/SAXW/ XX(30,500),YY(30,500),SS(30,500),NBPINTER,IMAT
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
      data pf,pf2,pf3 /270.33936,73083.4,19756627./
C surface : modif de REF                                                                       P-N24720
CCC   REF=TIME NECESSARY FOR A NUCLEON TO REACH THE SURFACE             P-N24730
C                                                                       P-N24740
      T2=P1*P1+P2*P2+P3*P3                                              P-N24760
      p=sqrt(t2)
      r=R2
		IF (NOSURF.LE.0) THEN
C Modif pour W.S.:
	XV=p/pf
        if(t2.le.pf2) then
           r=FLIN(XV)
           call hfill(202, XV, 0.0, 1.0);
           call hfill(203, r, 0.0, 1.0);
           if(xv.gt.1.0) then
              call hfill(204, r, 0.0, 1.0);
           endif
        else
           r=rmaxws
        endif
	r=r*r
		END IF		
   21 t4=x1*x1+x2*x2+x3*x3
      if (t4.gt.r) go to 2
      T1=X1*P1+X2*P2+X3*P3                                              P-N24750
      T3=T1/T2                                                          P-N24770
    4 T5=T3*T3+(R-T4)/T2   
      if (t5.gt.0) go to 1
      ref=10000.
      return
    1 REF=(-T3+SQRT(T5))*E                                              P-N24830
      return
    2 s=sqrt(r*0.99/t4)
      x1=x1*s
      x2=x2*s
      x3=x3*s
      go to 21
      END   
C                                                                       
      SUBROUTINE PAULI(L,XR,PR,F)                                       P-N24870
C                                                                       P-N24880
C     THIS SUBROUTINE CALCULATES THE OCCUPATION IN PHASE SPACE AROUND   P-N24890
C         NUCLEON  L  , BY COUNTING THE PARTICLES IN A VOLUME AROUND L  P-N24900
C     THE VOLUME IS THE PRODUCT OF A SPHERE OF RADIUS XR IN R-SPACE BY  P-N24910
C         A SPHERE OF RADIUS PR IN MOMENTUM SPACE                       P-N24920
C     AVERAGE IS TAKEN ON THE SPIN ONLY                                 P-N24930
C                                                                       P-N24940
      COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TAP-N24950
      COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N24960
      COMMON/BL5/TLG(300),NESC(300)                                     P-N24970
      COMMON/SAXW/ XX(30,500),YY(30,500),SS(30,500),NBPINTER,IMAT
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
      IF (NPAULSTR.EQ.2) then
      f=0.
      return 	
      endif
      IF (NPAULSTR.EQ.1) THEN

c Pauli strict
       f=0.
       pmod=sqrt(p1(l)**2+p2(l)**2+p3(l)**2)
       if (pmod.lt.270.) f=1.
       return

      ELSE

      XR2=XR*XR                                                         P-N24980
      PR2=PR*PR                                                         P-N24990
      VOL=(4.*3.1415926/3.)**2*(XR*PR/(2.*3.1415926*197.13))**3         P-N25000
      RS=SQRT(X1(L)*X1(L)+X2(L)*X2(L)+X3(L)*X3(L))                      P-N25010
		IF (NOSURF.LE.0) THEN
C modifs A.B.: R2 -> RMAXWS pour la densite en W.S.
	rdeq=RMAXWS
		ELSE
	rdeq=R0
		END IF

C      IF (RS-XR.LE.R2) GO TO 3                                         P-N25020
      IF (RS-XR.LE.rdeq) GO TO 3
      F=0.                                                              P-N25030
      RETURN                                                            P-N25040
C    3 IF (RS+XR.GT.R2) VOL=VOL*0.5*(R2-RS+XR)/XR                       P-N25050
    3 IF (RS+XR.GT.rdeq) VOL=VOL*0.5*(rdeq-RS+XR)/XR
    2 IA=IA1+IA2                                                        P-N25060
      NL=0                                                              P-N25070
      DO 1 I=1,IA                                                       P-N25080
      IF (NESC(I).GT.0) GO TO 1                                         P-N25090
      IF (IND1(I).GT.0) GO TO 1                                         P-N25100
      IF (IND2(I).NE.IND2(L)) GO TO 1                                   P-N25110
      DX2=(X1(L)-X1(I))**2+(X2(L)-X2(I))**2+(X3(L)-X3(I))**2            P-N25120
      IF (DX2.GT.XR2) GO TO 1                                           P-N25130
      DP2=(P1(L)-P1(I))**2+(P2(L)-P2(I))**2+(P3(L)-P3(I))**2            P-N25140
      IF (DP2.GT.PR2) GO TO 1                                           P-N25150
      NL=NL+1                                                           P-N25160
    1 CONTINUE                                                          P-N25170
      F=(NL-1)/VOL/2.                                                   P-N25180
      IF (F.GT.1.) F=1.                                                 P-N25190

C Coupure imposee sur l'impulsion:
C       pmod=sqrt(p1(l)**2+p2(l)**2+p3(l)**2)
C       if (pmod.lt.150.) F=1.

      RETURN                                                            P-N25200

      ENDIF

      END                                                               P-N25210

C
CCC   IBM STANDARD RANDOM NUMBER
C
      SUBROUTINE RIBM(RNDM,IAL)
      common/debug/idebug
      if(idebug.eq.1) then
         rndm = ranecu(ial)
      else
      ial=ial*65539
      if(ial)1187,1188,1188
1187  ial=ial+2147483647+1
1188  al=ial
      rndm=al*0.4656613E-9
      endif
      return
      end
c-------------------------------------------------------------------------------
C
CCC   GAUSSIAN RANDOM NUMBER GENERATOR
C      
      subroutine rgauss(xg)
       common/hazard/ial,IY1,IY2,IY3,IY4,IY5,IY6,IY7,IY8,IY9,IY10,
     s               IY11,IY12,IY13,IY14,IY15,IY16,IY17,IY18,IY19
C      DATA IAL/25643/,IBL/45287/
C      ibl remplace par IY11
    1 xg=0.
      do 6 j=1,12
      call ribm(xxg,ial)
    6 xg=xg+xxg
      xg=xg-6.
      if (xg*xg.gt.9) go to 1
C      call ribm(xshuf,ibl)
      call ribm(xshuf,IY11)
      if (xshuf.gt.0.5) call ribm(xxg,ial)
    2 return
      end
c-------------------------------------------------------------------------------
      SUBROUTINE FORCE_ABS(IPROJO,AT,ZT,EP,BMAX,PT,PROBA)
      REAL*8 AP,ZP,A,Z,E 
      A=AT
      Z=ZT
      E=EP
      SIG_G=31.41592654*BMAX*BMAX
      	IF(IPROJO.EQ.1) THEN
      		AP=1.d0
      		ZP=1.d0
      	ELSE
      		AP=1.d0
      		ZP=0.d0
      	END IF
	
	CALL xabs2(ZP,AP,Z,A,E,sig_exp)	
	CALL sig_reac(IPROJO,E,A,sig_incl)
      
C      WRITE(6,*) 'SIG_G, SIG_EXP, SIG_INCL',SIG_G, SIG_EXP, SIG_INCL
      PROBA=(sig_exp-PT*sig_incl)/(PT*(sig_g-sig_incl))
      
      IF(PROBA.LE.0.) PROBA=0.
      IF(PROBA.GT.1.) PROBA=1.
      	
      RETURN
      END

      SUBROUTINE xabs2(zp,ap,zt,at,ep,sig)	                                   
c  NASA subroutine                                                      
      implicit double precision (a-h,o-z)
      REAL*4 sig                               
      save                                                              
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,  
     2 dp2th=dp2/dp3)                                                   
c absoprption xsec revised version rkt-97/5                             
c neutron data from barashenkov                                         
c this gives absorption xsec for given zp,ap,zt,at,e (mev/nucleon)      
c arguement changed to mev; then e=ep/ap mev/nucleon                    
c can be used for neutrons also.                                        
c     this has coulomb as ours                                          
      e=ep/ap                                                           
      xq1=zp                                                            
      xq2=ap	                                                           
      xq3=zt                                                            
      xq4=at	                                                           
c                                                                       
c nucleon-nucleon inelastc xsec not included here                       
c                                                                       
      if (nint(ap*at).eq.1.or.nint(zp+zt).eq.1) then                    
c     if (nint(zp*zt).eq.1.or.nint(zp+zt).eq.1) then                    
        sig=dp0                                                        
        return                                                          
      endif                                                             
      rp=radi_us(ap)                                                     
      rt=radi_us(at)                                                     
      vp=(dp1+dpth)*dppi*rp**3                                          
      vt=(dp1+dpth)*dppi*rt**3                                          
      density=dph*((ap/vp)+(at/vt))                                     
      const=1.75d0*density/8.824728d-02                                 
c     if (zt.lt.zp) then                                                
      if (zt.lt.zp.or.(zt.eq.zp.and.at.lt.ap)) then                     
        xzt=zp                                                          
        xat=ap                                                          
        zp=zt                                                           
        ap=at                                                           
        zt=xzt                                                          
        at=xat                                                          
      endif                                                             
      if (nint(ap).eq.1) const=2.05d0                                   
      if (nint(zp).eq.2.and.nint(ap).eq.4) then                         
        const1=2.77-at*8.0d-03+(at*at)*1.8d-05                          
      endif                                                             
      if (nint(zp).eq.3) const=const/3.                                 
      t1=40.d0                                                          
      if (nint(zp).eq.0) then                                           
        if (nint(at).ge.11.and.nint(at).lt.40) t1=30.d0                 
        if (nint(zt).eq.14) t1=35.                                      
        if (nint(zt).eq.26) t1=30.                                      
      endif                                                             
      gcm=(ap*(dp1+e/938.d0)+at)/(ap**2+at**2+dp2*ap*(e+938.d0)*at/938  
     1 .d0)**dph                                                        
      bcm=sqrt(dp1-dp1/gcm**2)                                          
      plab=ap*sqrt(dp2*938.d0*e+e*e)                                    
      ecmp=gcm*(e+938.d0)*ap-bcm*gcm*plab-ap*938.d0                     
      ecmt=gcm*938.d0*at-at*938.d0                                      
      rela=ecmp+ecmt                                                    
      ecm=rela                                                          
      if (ecm.lt.0.1d0*rela) ecm=0.1d0*rela                             
      rm=(197.32/137.01)*zp*zt/ecm                                      
      bigr=rp+rt+1.2*(ap**dpth+at**dpth)/(ecm**dpth)                    
      bigb=1.44*zp*zt/bigr                                              
      if (nint(zp).eq.1.and.nint(at).gt.56) bigb=0.90*bigb              
      if (nint(ap).gt.56.and.nint(zt).eq.1) bigb=0.90*bigb              
      if (nint(ap).eq.1.and.nint(at).eq.12) bigb=3.5*bigb               
      if (nint(ap).eq.1) then                                           
        if (nint(at).le.16.and.nint(at).ge.13) then                     
          bigb=(at/7.)*bigb                                             
        endif                                                           
        if (nint(zt).eq.12) bigb=1.8*bigb                               
        if (nint(zt).eq.14) bigb=1.4*bigb                               
        if (nint(zt).eq.20) bigb=1.3*bigb                               
      endif                                                             
      if (nint(ap).eq.1.and.nint(at).lt.4) bigb=21.0*bigb               
      if (nint(ap).lt.4.and.nint(at).eq.1) bigb=21.0*bigb               
      if (nint(ap).eq.1.and.nint(at).eq.4) bigb=27.0*bigb               
      if (nint(ap).eq.4.and.nint(at).eq.1) bigb=27.0*bigb               
      if (nint(zp).eq.0.or.nint(zt).eq.0) bigb=dp0	                     
      xsec=dp10*dppi*bigr*bigr*(dp1-bigb/ecm)			                        
      xm=1.                                                             
      if (nint(zp).eq.0) then                                           
        if (nint(at).lt.200) then                                       
          x1=2.83-3.1d-02*at+1.7d-04*at*at                              
          if (x1.le.1) x1=1.	                                           
          sl=dp1                                                        
          if (nint(at).eq.12) sl=1.6                                    
          if (nint(at).lt.12) sl=0.6                                    
          xm=(1-x1*exp(-e/(sl*x1)))                                     
          if (e.lt.20) print *, 'e,xm=',e,xm                            
        else                                                            
          xm=(1-0.3*exp(-(e-1)/15))*(1-exp(-(e-0.9))) !7/11/96          
        endif                                                           
      endif                                                             
      if (nint(zp).eq.2.and.nint(ap).eq.4) then                         
        const=const1-0.8/(1+exp((250.d0-e)/75.d0))                      
      endif                                                             
      if (nint(zp).eq.1.and.nint(ap).eq.1) then                         
        if (nint(at).gt.45) then                                        
          t1=40.d0+at/dp3                                               
        endif                                                           
        if (nint(at).lt.4) then                                         
          t1=55                                                         
        endif                                                           
        const=2.05d0-0.05/(1+exp((250.d0-e)/75.d0))                     
        if (nint(at).lt.4) const=1.7                                    
        if (nint(zt).eq.12) then                                        
          t1=40.d0                                                      
          const=2.05d0-dp3/(1.+texp((e-20.d0)/dp10))                    
        endif                                                           
        if (nint(zt).eq.14) then                                        
          t1=40.d0                                                      
          const=2.05d0-1.75d0/(1.+texp((e-20.d0)/dp10))                 
        endif                                                           
        if (nint(zt).eq.18) then                                        
          t1=40.d0                                                      
          const=2.05d0-dp2/(1.+texp((e-20.d0)/dp10))                    
        endif                                                           
        if (nint(zt).eq.20) then                                        
          t1=40.d0                                                      
          const=2.05d0-dp1/(1.+texp((e-40.d0)/dp10))                    
          const=const-0.25d0/(1+exp((250.d0-e)/75.d0))                  
        endif                                                           
        if (nint(zt).ge.35) then                                        
          phst=(nint(zt)-35.d0)/260.d0	                                 
          const=const-phst/(1+exp((250.d0-e)/75.d0))                    
        endif                                                           
      endif                                                             
      if (nint(zp).eq.0.and.nint(ap).eq.1) then                         
        const=2*(0.134457/density)                                      
        if (nint(at).gt.140.and.nint(at).lt.200) const=const-1.5d0*(at- 
     1 dp2*zt)/at                                                       
        if (nint(at).lt.60) const=const-1.5d0*(at-dp2*zt)/at            
        if (nint(at).le.40) const=const+0.25d0/(dp1+texp(-(170.d0-e)/   
     1 100.d0))                                                         
        if (nint(zt).gt.82) const=const-zt/(at-zt)                      
        if (nint(zt).ge.82) const=const-dp2/(1.+texp((e-20.d0)/20.d0))  
        if (nint(zt).le.20.and.nint(zt).ge.10) const=const-dp1/(dp1+texp
     1 ((e-20.d0)/dp10))                                                
      endif                                                             
      ce=const*(1.-exp(-e/t1))-0.292*exp(-e/792)*cos(0.229*e**0.453)    
      term1=(at*ap)**dpth/(at**dpth+ap**dpth)                           
      delta=1.615*term1-0.873d0*ce                                      
      delta=delta+0.140*term1/ecm**dpth                                 
      delta=delta+0.794d0*(at-dp2*zt)*zp/(at*ap)                        
      delta=-delta                                                      
      beta=1.                                                           
      twxsec=dp10*dppi*1.26d0*1.26d0*beta*(0.873d0*ap**dpth+0.873d0*at* 
     1 *dpth-delta)**2                                                  
      sig=twxsec*(dp1-bigb/ecm)*xm                                     
      if (sig.lt.dp0) sig=dp0
c      	write(6,*)ep,sig,ecm,t1,const,xm                                         
      zp=xq1                                                            
      ap=xq2                                                            
      zt=xq3                                                            
      at=xq4                                                            
      return                                                            
      end                                                               
      function radi_us (a)                                               
c  NASA subroutine                                                      
      implicit double precision (a-h,o-z)                               
      save                                                              
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,  
     2 dp2th=dp2/dp3)                                                   
      dimension na(23), rms(23)                                         
      data na /1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,2
     1 5,26/                                                            
      data rms /0.85d0,2.095d0,1.976d0,1.671d0,2.57d0,2.41d0,2.519d0,2.4
     1 5d0,2.42d0,2.471d0,2.440d0,2.58d0,2.611d0,2.730d0,2.662d0,2.727d0
     2 ,2.900d0,3.040d0,2.969d0,2.94d0,3.075d0,3.11d0,3.06/             
      fact=sqrt(dp5/dp3)                                                
      ia=a+0.4d0                                                        
      radi_us=fact*(0.84d0*a**dpth+0.55d0)                               
      do 20 i=1,23                                                      
      if (ia.eq.na(i)) go to 10                                         
      go to 20                                                          
   10 radi_us=fact*rms(i)                                                
   20 continue                                                          
      return                                                            
      end                                                               

      function texp (x)                                                 
c  NASA subroutine                                                      
      implicit double precision (a-h,o-z)                               
      save                                                              
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,  
     2 dp2th=dp2/dp3)                                                   
c                                                                       
c   purpose                                                             
c   to eliminate over/under flow of cpu if exp is used                  
c                                                                       
c   result = texp(x)                                                    
      if (x.lt.-80.d0) x=-80.d0                                         
      if (x.gt.80.d0) x=80.d0                                           
      texp=exp(x)                                                       
      end                                                               

      SUBROUTINE sig_reac(iprojo,e,a,sig)
C parametrisation de la section efficace de reaction calculee par INCL4.1
c iprojo=1 proton incident, iprojo=2, neutron incident).
C Entre Al et U, entre 10 et 100 MeV protons, 20 et 100 MeV neutrons.
C Bon ordre de grandeur pour les noyaux legers (C, O ...), tres faux
C a energie sup a 100 MeV.
      implicit double precision (a-h,o-z)
      REAL*4 sig                               
      save 
                                                                   
      DIMENSION coefp(3,4),coefn(3,4),Apow(3),Epow(5),
     s          coef2p(3,5),coef2n(3,5)
      DATA coefp/-5.926d-9,6.8945d-6,-6.098d-6,
     s           2.1544d-6,-1.848d-3,-5.982d-4,
     s           -2.59d-4,0.17595d+0,1.1741d+0,
     s           1.1504d-3,2.8281d+0,-28.730d+0/ 
      DATA coefn/1.6105d-9,3.3985d-7,1.4678d-5,
     s           -5.35d-7,-3.465d-4,-0.01633d+0,
     s           -4.755d-6,0.07608d+0,2.8135d+0,
     s           -3.622d-3,3.5924d+0,-38.294d+0/ 
      DATA coef2p/6.8108d-9,-2.163d-7,2.1898d-6,
     s           -2.187d-6,7.8331d-5,-7.164d-4,
     s           2.3651d-4,-9.690d-3,0.076424d+0,
     s           -9.195d-3,0.5030d+0,-2.4979d+0,
     s           -0.01087d+0,2.6494d+0,-2.5173d+0/ 
         
      IF(a.GE.27.d0) THEN      
      	ii=3        
      	jj=4
      ELSE
        ii=3
	jj=5
      ENDIF
      
      DO i=1,ii
      	Apow(i)=a**(ii-i)
      END DO
      
      DO j=1,jj
        Epow(j)=e**(jj-j)
      END DO
	      
      res=0.d0
      IF(a.GE.27.d0) THEN 
           
      IF(IPROJO.EQ.1) THEN
	      DO i=1,ii
	      	DO j=1,jj
	      res = res + coefp(i,j)*Apow(i)*Epow(j)
C      WRITE(6,*) 'i,j:',i,j,coefp(i,j)
	      	END DO
	      END DO
      ELSE
	      DO i=1,ii
	      	DO j=1,jj
	      res = res + coefn(i,j)*Apow(i)*Epow(j)
C      WRITE(6,*) 'i,j:',i,j,coefn(i,j)
	      	END DO
	      END DO
      ENDIF
      
      ELSE
      	
C      IF(IPROJO.EQ.1) THEN
	      DO i=1,ii
	      	DO j=1,jj
	      res = res + coef2p(i,j)*Apow(i)*Epow(j)
C      WRITE(6,*) 'i,j:',i,j,coefp(i,j)
	      	END DO
	      END DO
      
      ENDIF
      
      sig = res
      
      RETURN
      END

      SUBROUTINE COULOMB_TRANSM(E,FM1,Z1,FM2,Z2,PROBA)      
C calcul du coulombien dans LAHET (proba de transmission ou
C d'entree dans le potentiel nucleaire).
      REAL*8 ETA,RHO,clmb1
      C2 = .00516
      C3 = .007165
      uma=938.0
           
      ECM=E*FM2/(FM1+FM2)
      fm= FM1*FM2*uma/(FM1+FM2)
      R1=0.
      IF(FM1.GE.2.) R1=1.2*FM1**0.33333333
      R2=1.2*FM2**0.33333333
      ETA = C2*Z1*Z2*SQRT(fm/ECM)
      RHO = C3*(R1+R2)*SQRT(fm*ECM)
      proba = clmb1(RHO,ETA,ml)
      
      RETURN
      END

      function clmb1 (ro,eta,ml)                                        
      implicit double precision (a-h,o-z)                               
      save                                                              
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,  
     2 dp2th=dp2/dp3)                                                   
c                                                                       
c     ETA = C2*Z1*Z2*SQRT(M/E)                                          
c     RHO = C3*(R1+R2)*SQRT(ME)                                         
c     M = reduced mass in MeV                                           
c     E = kinetic energy in C of M                                      
c     C2 = .00516 = fine structure constant / sqrt(2)                   
c     C3 = .007165 = sqrt(2)/(hbar-c)                                   
c                                                                       
      parameter (prm1=69.06d0)                                          
      parameter (ln0=81, lt0=21)                                        
      parameter (ln1=61, lt1=61)                                        
      dimension psi0(ln0), trans0(lt0,ln0), x0(lt0), f0(ln0), psi1(ln1),
     1 trans1(lt1,ln1), x1(lt1), f1(ln1)                                
      data pi /3.14159d0/                                               
      data c0 /.11225d0/, c1 /dph/, gamma /0.5772157d0/, s3 /0.2020569d0
     1 /, s4 /0.08232323d0/                                             
      y=dp2*eta                                                         
      psi=ro*y                                                          
      if (ro.gt.y) then                                                 
        if (psi.gt.dp4.and.psi.lt.50.d0) then                           
          prob=clmb2(ro,eta,dumm)                                       
        else                                                            
          x=exp(log(eta)/6.d0)                                          
          x3=x*x*x                                                      
          temp=c0+c1*x3                                                 
          temp=temp+ro*x                                                
          arg=dp1-y*x/temp                                              
          prob=sqrt(arg)                                                
        endif                                                           
        ml=0                                                            
      else                                                              
        x=ro/y                                                          
        if (psi.le.psi0(1)) then                                        
          t=min(pi*y,prm1)                                              
          cx=t/(exp(t)-dp1)                                             
          t1=cos(ro)*(dp1-.75d0*psi**2+dp5*x*psi**2)-dph*psi*ro*sin(ro) 
          t2=dp1+dph*psi*(dp1-x/6.d0)                                   
          if (eta.gt.dp1) then                                          
            t3=log(psi)+dp2*gamma-dp1+dp1/(12.d0*eta**2)+dp1/(12.d1*eta*
     1      *4)                                                         
          else                                                          
            t3=log(dp2*ro)+gamma-dp1/(dp1+eta**2)+s3*eta**2+s4*eta**4   
          endif                                                         
          g=t1+psi*t2*t3                                                
          f=cx*ro*t2                                                    
          prob=cx/(g**2+f**2)                                           
          ml=3                                                          
        elseif (psi.le.psi0(ln0)) then                                  
          if (x.le.x0(1)) then                                          
            temp=log(psi/psi0(1))                                       
            j0=1+int(temp/delp0)                                        
            j0=min(max(j0,1),ln0-1)                                     
            temp=temp-(j0-1)*delp0                                      
            t=f0(j0)+(f0(j0+1)-f0(j0))*temp/delp0                       
            xk=x*sqrt(psi)                                              
            prob=(dp1+3.33d-1*x+3.d-1*xk+1.d-1*xk**2)*exp(t)            
            t=min(pi*y,prm1)                                            
            cx=t/(exp(t)-dp1)                                           
            prob=cx/prob**2                                             
            ml=1                                                        
          else                                                          
            temp1=log(x/x0(1))                                          
            i0=1+int(temp1/delx0)                                       
            i0=min(max(i0,1),lt0-1)                                     
            temp1=temp1-(i0-1)*delx0                                    
            temp2=log(psi/psi0(1))                                      
            j0=1+int(temp2/delp0)                                       
            j0=min(max(j0,1),ln0-1)                                     
            temp2=temp2-(j0-1)*delp0                                    
            t1=trans0(i0,j0)+(trans0(i0+1,j0)-trans0(i0,j0))*temp1/delx0
            t2=trans0(i0,j0+1)+(trans0(i0+1,j0+1)-trans0(i0,j0+1))*temp1
     1      /delx0                                                      
            t=t1+(t2-t1)*temp2/delp0                                    
            ml=2                                                        
            prob=exp(t)                                                 
          endif                                                         
        elseif (psi.le.psi1(ln1)) then                                  
          if (x.le.x1(1)) then                                          
            temp=log(psi/psi1(1))                                       
            j0=1+int(temp/delp1)                                        
            j0=min(max(j0,1),ln1-1)                                     
            temp=temp-(j0-1)*delp1                                      
            t=f1(j0)+(f1(j0+1)-f1(j0))*temp/delp1                       
            xk=x*sqrt(psi)                                              
            prob=(dp1+3.33d-1*x+3.d-1*xk+1.d-1*xk**2)*exp(t)            
            t=min(pi*y,prm1)                                            
            cx=t/(exp(t)-dp1)                                           
            prob=cx/prob**2                                             
            ml=1                                                        
          else                                                          
            temp1=log(x/x1(1))                                          
            i0=1+int(temp1/delx1)                                       
            i0=min(max(i0,1),lt1-1)                                     
            temp1=temp1-(i0-1)*delx1                                    
            temp2=log(psi/psi1(1))                                      
            j0=1+int(temp2/delp1)                                       
            j0=min(max(j0,1),ln1-1)                                     
            temp2=temp2-(j0-1)*delp1                                    
            t1=trans1(i0,j0)+(trans1(i0+1,j0)-trans1(i0,j0))*temp1/delx1
            t2=trans1(i0,j0+1)+(trans1(i0+1,j0+1)-trans1(i0,j0+1))*temp1
     1      /delx1                                                      
            t=t1+(t2-t1)*temp2/delp1                                    
            ml=2                                                        
            prob=exp(t)                                                 
          endif                                                         
        else                                                            
          prob=clmb2(ro,eta,dumm)                                       
          ml=4                                                          
        endif                                                           
      endif                                                             
      clmb1=prob                                                        
      return                                                            
      entry rdtab(lun)                                                  
      read (lun) delp0,delx0                                            
      read (lun) psi0                                                   
      read (lun) x0                                                     
      read (lun) f0                                                     
      read (lun) trans0                                                 
      read (lun) delp1,delx1                                            
      read (lun) psi1                                                   
      read (lun) x1                                                     
      read (lun) f1                                                     
      read (lun) trans1                                                 
      return                                                            
      end                                                               
      function clmb2 (ro,eta,t1)                                        
      implicit double precision (a-h,o-z)                               
      save                                                              
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,  
     2 dp2th=dp2/dp3)                                                   
      parameter (ln=101)                                                
      dimension t0(ln)                                                  
      data t0 /dp0,.1083d0,.1369d0,.1572d0,.1736d0,.1876d0,.2d0,.2113d0,
     1 .2216d0,.2312d0,.2403d0,.2489d0,.2571d0,.265d0,.2725d0,.2798d0,.2
     2 869d0,.2938d0,.3006d0,.3071d0,.3136d0,.3199d0,.3261d0,.3322d0,.33
     3 82d0,.3442d0,.3499d0,.3557d0,.3615d0,.3672d0,.3729d0,.3785d0,.384
     4 1d0,.3897d0,.3952d0,.4008d0,.4063d0,.4118d0,.4173d0,.4228d0,.4283
     5 d0,.4338d0,.4393d0,.4448d0,.4504d0,.4559d0,.4615d0,.4671d0,.4728d
     6 0,.4784d0,.4841d0,.4899d0,.4957d0,.5015d0,.5074d0,.5133d0,.5193d0
     7 ,.5253d0,.5315d0,.5376d0,.5439d0,.5503d0,.5567d0,.5632d0,.5698d0,
     8 .5765d0,.5833d0,.5903d0,.5973d0,.6045d0,.6118d0,.6193d0,.6269d0,.
     9 6346d0,.6426d0,.6507d0,.659d0,.6675d0,.6763d0,.6853d0,.6945d0,.70
     $ 4d0,.7139d0,.724d0,.7345d0,.7453d0,.7566d0,.7683d0,.7805d0,.7932d
     $ 0,.8065d0,.8205d0,.8352d0,.8508d0,.8673d0,.8849d0,.9038d0,.9243d0
     $ ,.9467d0,.9715d0,dp1/                                            
      data x1 /1.d-2/, xi /1.d2/                                        
      x=dp1/(dp1+sqrt(dph*ro*eta))                                      
      if (x.lt.x1) then                                                 
        temp=t0(2)*(x/x1)**(dpth)                                       
      else                                                              
        i=xi*x                                                          
        i=i+1                                       
        i=max(min(i,ln-1),2)                                            
        temp=t0(i)+(t0(i+1)-t0(i))*(x-i*x1)/x1                          
      endif                                                             
      t1=dp1-temp                                                       
      prob=dp1-dp2*t1*eta/ro                                            
      clmb2=max(prob,dp0)                                               
      return                                                            
      end                                                               
C------------------------------------------------------------------------------
      SUBROUTINE FORCE_ABSOR(nopart,F,IAREM,IZREM,ESREM,ERECREM,
     s  ALREM,BEREM,GAREM,JREM) 
C Absorption forcee pour p (10-100 MeV) et n (20-100MeV)

      DIMENSION F(15)
      REAL*4 ia1,iz1
      
      COMMON/hazard/ial,IY(19)
C Dialogue with INCL for nucleus density and parameters.
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
C RMS espace R, espace P, Fermi momentum and energy for light gauss nuc.      
      COMMON/light_gaus_nuc/rms1t(9),pf1t(9),pfln(9),tfln(9),vnuc(9)
      
     
     
      IF(nopart.NE.-1) RETURN
      
      FMP=938.2796
      ia2=f(1)
      SEP=6.8309
      IF(ia2.GT.4) GOTO 10
      	   IF(ia2.EQ.2) itg=6
	   IF(ia2.EQ.3.AND.f(2).EQ.1) itg=7
	   IF(ia2.EQ.3.AND.f(2).EQ.2) itg=8
	   IF(ia2.EQ.4) itg=9
           SEP=VNUC(itg)-tfln(itg)
10    CONTINUE 
     
      IF(f(3).GE.10..AND.f(3).LE.100.) THEN
      	    IF(f(7).EQ.1.)THEN
            
            ia1=1.
            iz1=1.
            fmpinc=938.2796
            PBEAM2=F(3)*(F(3)+2.*fmpinc)
	    BMAXT=BMAX
	    CALL COULOMB_TRANSM(f(3),ia1,iz1,f(1),f(2),PROBA_TRANS)
	    CALL FORCE_ABS(1,F(1),F(2),F(3),BMAXT,PROBA_TRANS,PROBA)
	    CALL RIBM(ALEA,IY(5))
	    IF(ALEA.GT.PROBA) RETURN
		IAREM=f(1)+ia1
		IZREM=f(2)+iz1
C	    	ESREM=f(3)+SEP

		del=SQRT(((f(1)+1.)*fmpinc+f(3))**2-PBEAM2)
		ERECREM=PBEAM2/((f(1)+1.)*fmpinc+f(3) + del)
                
                ESREM=f(3)+SEP -ERECREM
		
		ALREM=0.00001
		BEREM=0.0
		GAREM=0.99999
		BIMPACT=0.
		JREM=0.
		nopart = 0
	    	RETURN
	    ELSE IF(f(7).EQ.2.AND.f(3).GE.20.) THEN
      
      	    ia1=1.
      	    iz1=0.
            fmpinc=938.2796
            PBEAM2=F(3)*(F(3)+2.*fmpinc)
	    BMAXT=BMAX
	    CALL COULOMB_TRANSM(f(3),ia1,iz1,f(1),f(2),PROBA_TRANS)
	    CALL FORCE_ABS(2,F(1),F(2),F(3),BMAXT,PROBA_TRANS,PROBA)
	    CALL RIBM(ALEA,IY(5))
	    IF(ALEA.GT.PROBA) RETURN
		IAREM=f(1)+ia1
		IZREM=f(2)
C	    	ESREM=f(3)+SEP

		del=SQRT(((f(1)+1.)*fmpinc+f(3))**2-PBEAM2)
		ERECREM=PBEAM2/((f(1)+1.)*fmpinc+f(3) + del)
		          
                ESREM=f(3)+SEP -ERECREM

		ALREM=0.00001
		BEREM=0.0
		GAREM=0.99999
		BIMPACT=0.
		JREM=0.
	    	nopart = 0
	    	RETURN
	    END IF
      END IF
      
      RETURN
      END      	

