C********************************************************************** 
C                                                                       
C      PROGRAM Thin_Target_Monte-Carlo_Spallation_Code (TT_MSC)
C
C				INCL4.2
C
C**********************************************************************
C	      The aim of this code is to provide a self consistent MAIN
C	to be used for a spallation calculation in the case of a thin 
C	target (one interaction, no transport) starting from the intra
C	nuclear cascade of Liège (INCL), and coupled with various
C	evaporation and fission codes.
C	      The idea is to provide a tool for physicists to
C 	check various models and options compared with observables
C	obtained in the one interaction limit, and/or to make a good
C	choice before using the proposed physics models in complex 
C	geometries treated by transport codes (LAHET, MCNPX, HERMES,
C	GEANT).
C	TT_MSC 	is developped at Saclay by A. BOUDARD and C. VOLANT.
C	INCL	is developped at Liège by J. CUGNON
C	KHS	is developped at GSI by K.H. SCHMIDT,
C			and at Santiago de Compostella by J. BENLLIURE.
C	GEM     is developped in Tokyo by S. FURIHATA 
C		(For each code the list of contributors is much larger,
C		but the names should be understood also as possible
C		contacts)
C	The Physics case is defined by:
C		a projectile (proton, neutron, deuteron, triton, 3He, 4He,
C				pi+, pi- or pi0) and its kinetic energy.
C		a target nucleus at rest. All inputs and outputs are
C			expressed in the "target at rest" system.
C		a number of shots (number of projectiles sent in a rather
C			large "geometrical" cross section).
C	Outputs are given in a NTUPLE to be treated by PAW. 
C		For each projectile-target interaction, all produced particle
C		is recorded (nature, energy, momentum, theta, phi direction)
C		as well as more global variables as impact parameter 
C		or intermediate quantities as
C		excitation, intrinsic spin, momentum of the remnant
C		nucleus at the end of the intra nuclear cascade 
C		or nature of the eventual fissioning nucleus.
C		
C		Global numbers are also given for normalisation to an 
C		absolute cross section for any choice of particles
C		or correlation between them.  
C**********************************************************************     
C                                                                       
C     INPUT:  (in a file of name cu42.in)
C		path and name of a file (.out) to receive comments and results
C		Primary seed for random gene., Choice of evapo:
C						1=GSI KHS model
C						2=GEM model
C		Number of projectiles, print for event #, seeds given 0/1=Y/N
C			(It is possible to print the status of the 20 seeds
C			on a FORTRAN test, and then to give them here to go
C			directly to a bugged event for Monte-Carlo debugging)
C		type of projectile, total kinetic energy (MeV)
C		           1.00 FOR PROTON
C 		           2.00 FOR NEUTRON
C 		           3.00 FOR PI+
C 		           4.00 FOR PI0
C 		           5.00 FOR PI-
C		           6.00 FOR DEUTERON
C		           7.00 FOR TRITON
C		           8.00 FOR HE3
C		           9.00 FOR HE4
C		Target mass number and charge (A and Z) 
C		Nuclear potential, scaling factor of the standard cascade
C			 stoping time, NOSURF, XFOISA, NPAULSTR
C			Standard values: 45. 1. -2  8. 0			                        
C                                 NOSURF=1 Sharp density (hard sphere), 
C                                 NOSURF=0 Wood-Saxon density, stopping time "70" 
C                                                without B (impact) dependence.
C                                 NOSURF=-1 Wood-Saxon density, stopping time "70" 
C                                                with B (impact) dependence
C                                 NOSURF=-2 Wood-Saxon density, official stopping time
C						of the INCL4 version: 70*(At/208)**0.16 
C                                 
C                                XFOISA    Rmaxws = R0 + XFOISA*Diffuseness
C					   (Bmax = Rmaxws for pions and nucleons
C					   Bmax = Rmaxws + rms1t (data) for composits)
C         			 NPAULSTR  0=Pauli statistic
C					   1=Pauli strict 
C					   2=No Pauli blocking
C		Path to found table files (for KHS evapo-fission)
C		Path and name of a file to receive the ntuple of results
C
C The evaporation-fission step can be choosen with choice_evap between
C		KHS_V3		(1) and default
C		GEM		(2)		(08/2002)
C               Input parameters for KHS are explicit in:
C				abla_v3p.f, subroutine init_evapora
C		Input parameters for GEM are in the local file "in":
C				(init_AB.f, subroutine init_ab read(in)) 
C*************************************************************************************
C
C     OUTPUT:     in the .OUT file, characteristics of the run,
C			number of fissions, transparents (no nuclear interaction),
C			nuclear absorption,
C			Geometrical cross section used and total reaction cross
C			section found,
C			number of neutrons produced by nuclear interaction in
C				various energy bins,
C			execution time.
C		  in the NTUPLE:             
******************************************************************
*      1   * I*4  *             * VAR_NTP  * MASSINI	     A of the remnant
*      2   * I*4  *             * VAR_NTP  * MZINI	     Z    "        "
*      3   * R*4  *             * VAR_NTP  * EXINI	     Excit energy " "
*      4   * I*4  *             * VAR_NTP  * MULNCASC	     Cascade n multip.
*      5   * I*4  *             * VAR_NTP  * MULNEVAP	     Evapo   "      "
*      6   * I*4  *             * VAR_NTP  * MULNTOT	     Total   "      "
*      7   * R*4  *             * VAR_NTP  * BIMPACT	     Impact parameter
*      8   * I*4  *             * VAR_NTP  * JREMN	     Remnant Intrinsic Spin
*      9   * I*4  *             * VAR_NTP  * KFIS	     Fission 1/0=Y/N
*     10   * R*4  *             * VAR_NTP  * ESTFIS		Excit energy at fis
*     11   * I*4  *             * VAR_NTP  * IZFIS		Z of fiss nucleus
*     12   * I*4  *             * VAR_NTP  * IAFIS		A of "          "
*     13   * I*4  *[0,250]      * VAR_NTP  * NTRACK		Number of particles
*     14   * I*4  *             * VAR_NTP  * ITYP(NTRACK)	emitted in cascade (0)
*                                                                       or evapo   (1)
*     15   * I*4  *             * VAR_NTP  * AVV(NTRACK)		A (-1 for pions)
*     16   * I*4  *             * VAR_NTP  * ZVV(NTRACK)		Z
*     17   * R*4  *             * VAR_NTP  * ENERJ(NTRACK)		kinetic energy
*     18   * R*4  *             * VAR_NTP  * PLAB(NTRACK)		momentum
*     19   * R*4  *             * VAR_NTP  * TETLAB(NTRACK)		Theta (deg)
*     20   * R*4  *             * VAR_NTP  * PHILAB(NTRACK)		Phi   (deg)
******************************************************************
C   Additional informations:
C                                                                       
c  This version has a unique random generator. One event can be reproduced
C  if ialview and IYV(1...19) are printed during the test and then given
C  as input for initialisation of the common hazard.
C
C With INCL4.2: This version takes into account the coulomb entrance
C               barrier (from LAHET) as an additional source of transparent.
C               Below 100 MeV, absorption is increased so that the experiment
C               reaction cross section (from LAHET-TRIPATHY) is obtained.
C The evaporation-fission step can be choosen with choice_evap between
C		KHS_V3		(1) and default
C		GEM		(2)		(08/2002)
C
C 9/2002, the nuclear mass is taken from KHS model. Solves some incoherences
C    in the kinematics computed in the MAIN. Still problems with GEM.
C    MGLMS(A,Z,0,EL) est la fonction de masse cohérente avec KHS evapo-fis.
C    Attention aux parametres, ici 0=OPTSHP, NO microscopic corrections!
C
C 11/2003 Corrections proposées par John Hendricks (MCNPX): précisions
C	numériques pour le calcul des angles (translab)
C *********************************************************************** 
C *********************************************************************** 
C
	PROGRAM  TT_MSC
	
C *****************   Definition des COMMON *****************************
cjcd
      common /inout/ in, io, itty, iscrt                                
cjcd
  
C Dialogue with INCL for nucleus density and parameters.
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
      
C Input for INCL, to be known at the CALL PNUC (SUPPRIMER icoup!!)
      COMMON/CALINCL/FINPUT(15),icoup
      
C input: should contain Z and A for NBMAT nucleus considered in this problem 
      INTEGER ZMAT,AMAT     
      COMMON/MAT/ZMAT(500),AMAT(500),BMAX_GEO(6,500),NBMAT
      	   
C ial generateur pour le cascade (et les IY pour eviter les correlations)
      DIMENSION IY(19),IYV(19)
      COMMON/hazard/ial,IY

      REAL*8 AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,R_0,R_P,R_T,PI,   
     s       BFPRO,SNPRO,SPPRO,SHELL                                    
      INTEGER*4 IMAX,IRNDM
      COMMON /ABLAMAIN/ AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,R_0,   
     s                  R_P,R_T,IMAX,IRNDM,PI,BFPRO,SNPRO,SPPRO,SHELL   
      
          
      INTEGER*4 IFIS,OPTSHP,OPTXFIS,OPTLES,OPTCOL               
      REAL*8 AKAP,BET,HOMEGA,KOEFF
      COMMON /FISS/ AKAP,BET,HOMEGA,KOEFF,IFIS,OPTSHP,                   
     s       OPTXFIS,OPTLES,OPTCOL
                                           
      REAL*8 AV,AS,AK
      INTEGER*4 OPTAFAN                                                   
      COMMON /ALD/AV,AS,AK,OPTAFAN   
      
      REAL*8 EFA                                                                                                 
      COMMON /FB/ EFA(0:100,0:160)
                                            
      INTEGER I
      COMMON/NEVENT/I
      
      REAL*8 ME(100),RR(100)
      INTEGER AA(100),ZZ(100),NPMAX
      COMMON /MER/AA,ZZ,ME,RR,NPMAX                                              
      
	parameter (max=250)                                                                       
	REAL*4 EXINI,ENERJ,BIMPACT,PLAB,TETLAB,PHILAB,ESTFIS      
	INTEGER AVV,ZVV,JREMN,KFIS,IZFIS,IAFIS      
        COMMON/VAR_NTP/MASSINI,MZINI,EXINI,MULNCASC,MULNEVAP,
     +MULNTOT,BIMPACT,JREMN,KFIS,ESTFIS,IZFIS,IAFIS,NTRACK,
     +ITYPCASC(max),AVV(max),ZVV(max),ENERJ(max),PLAB(max),
     +TETLAB(max),PHILAB(max)
     
      REAL*4 Bavat,TIME,ENERGY,EPSd,EPS2,EPS4,EPS6,EPSf
      INTEGER Bloc_Paul,Bloc_CDPP,GO_OUT,avm,DEL1,DEL2
      PARAMETER (avm=1000)
      COMMON/VAR_AVAT/Kveux,Bavat,NOPART,NCOL,
     s R1_in(3),R1_first_avat(3),
     s EPSd(250),EPS2(250),EPS4(250),EPS6(250),EPSf(250),      
     s NB_AVAT,
     s TIME(avm),L1(avm),L2(avm),JPARTL1(avm),JPARTL2(avm),
     s DEL1(avm),DEL2(avm),ENERGY(avm),Bloc_Paul(avm),
     s Bloc_CDPP(avm),GO_OUT(avm)

      REAL*4 acv,zpcv,pcv,xcv,ycv,zcv
      COMMON/volant/acv(200),zpcv(200),pcv(200),xcv(200),
     s              ycv(200),zcv(200),iv 

C ****************** Fin de definition des COMMON ***********************

C ****************** Definition Variables locales ***********************	
      DIMENSION f(15),kind(300),ep(300),alpha(300),beta2(300),gam(300)
      
      REAL*8 MTOTA,MALPHA1,MALPHA2
      REAL*8 PLEVA,PTEVA,PXEVA,PYEVA,PEVA,FEE
      REAL*8 FFPLEVA1,FFPXEVA1,FFPYEVA1,FFPLEVA2,FFPXEVA2,FFPYEVA2 
      INTEGER*4 AFP,ZFP,INUM                     
       
      REAL*4 rndm                                                  
      REAL*8 ZF,AF,ZPRF,APRF,EE,EL,FMP,FMN         
      REAL*8 JPRF                       
      REAL*8 ZF1,AF1,ZF2,AF2
      REAL*8 TREM,REMMASS,GAMREM,ETREM,CSREM(3)
      REAL*8 MASSEF,MASSE1,MASSE2,EF,B,T1,P1,T2,P2,UMA,MELEC
      REAL*8 CTET1,PHI1,R(3,3),PLAB1,GAM1,ETA1,CSDIR1(3)
      REAL*8 CTET2,PHI2,PLAB2,GAM2,ETA2,CSDIR2(3)
      REAL*8 MCOREM,GAMFIS,ETFIS,MPROJO
      REAL*8 S1x,S1y,S1z,S2x,S2y,S2z,PXMOY,PYMOY,PZMOY
      REAL*8 fmcv,e_evapo,m_frag,m_cour
                                                                               
      INTEGER*4 J,FF,INTTYPE,FF1,FTYPE1,FF2,FTYPE2
      
      INTEGER choice_evap
      REAL*4 BIMPAC
      
      REAL muln_t,muln_2,muln_20,muln_max
 
      REAL*4 MASSDEF
      CHARACTER*80 stringout,RACINE
      
      DIMENSION PFIS_REM(3),PF1_REM(3),PF2_REM(3),PFIS_TRAV(3)
      REAL*8 MASSE,EPF1_IN,EPF2_IN,EPF1_OUT,EPF2_OUT
      
C Pour GEM:      
      REAL*8 ERECI,BET0I(3),WT,EKINF(300),BET1F(300,3)
      DIMENSION IZF_AB(300),IAF_AB(300)
            
C ***************************  DATA  ************************************                                                
      DATA UMA,MELEC/931.4942,0.511/
      DATA FMP,FMN/938.27231,939.56563/
C generateurs secondaires:
      DATA (IY(ihaz),ihaz=1,19)/
     s 21033,17563,27563,33657,43657,56375,66375,77365,87365,
     s 12345,22345,32345,42345,52345,62345,72345,82345,34567,47059/


C FILE des donnees ecrites en dur pour WORSHOP-DEBUG
	OPEN(5,file='cu42.in',status='old')
      READ(5,*) stringout	! output file
      	OPEN(6,file=stringout)
	
C ****** Time start *********************

          CALL getime(0)

	
C sortie (WRITE) pour INIT_INCL
	io=6	

C ****** Petit commentaire ***************
      WRITE(6,*) ' '
      WRITE(6,*) 'Thin_Target_Monte-Carlo_Spallation_Code (TT_MSC)'
      WRITE(6,*) ' '
	
	
C ****** Lecture de donnees specifiques au test du programme *****
C Premiere graine, choix d'evapo (1=KHS, 2=GEM)
      READ(5,*)ial,choice_evap
C icoup=number of tried shots	
C ievtest= number of the event for which we want a print of the common hazard 
C inithaz=1 pour initialisation du common hazard.
      READ(5,*)icoup,ievtest,inithaz
      		INIT_GRAINE = 1		!init of seeds in INIT_INCL
      IF (inithaz.EQ.1) THEN
      	INIT_GRAINE = 0      
        READ(5,*) ial,(IY(ihaz),ihaz=1,19)
      ENDIF
C	write(6,*)' on est parti pour ',icoup,' runs',inithaz,ievtest
C	READ(5,*) f(7),f(3)
C	write(6,*)' type projectile, energie projectile ',f(7),f(3)
c      READ(5,*)icoup
      write(6,*)' on est parti pour ',icoup,' runs'
      READ(5,*) finput(7),finput(3)
      write(6,*)' type projectile, energie projectile ',
     s                                         finput(7),finput(3)
      READ(5,*) finput(1),finput(2)
      write(6,*)'Masse, Charge cible ',finput(1),finput(2)
      ZMAT(1)=finput(2)+0.5
      AMAT(1)=finput(1)+0.5
      NBMAT=1      
      
C ******  INIT of INCL ********************
C --------------   Lecture des donnees File 5 ----------------------------
      
      READ(5,*) finput(5),finput(6),NOSURF,XFOISA,NPAULSTR
C racine du chemin pour trouver les tables (.tab,.dat)
      READ(5,*) RACINE
      
          CALL INIT_INCL(INIT_GRAINE)
	  
	  DO i=1,15
	     F(i)=FINPUT(i)
	  ENDDO
	  F(9)=1.

C ******  INIT of EVAPORATION *************
      IF (choice_evap.EQ.2) THEN
C 		******  INIT of EVAPORATION GEM *********
	  WRITE(6,*) '*********************************************'
          WRITE(6,*) ' PARAMETERS of GEM (fission-evapo)'
          WRITE(6,*) ' '
C ::::Commented out because init_ab not found!
C          		CALL init_ab(RACINE) 
      ELSE	  
C 		******  INIT of EVAPORATION KHS *********	  
	  WRITE(6,*) '*********************************************'
          WRITE(6,*) ' PARAMETERS of KHSv3p (fission-evapo)'
          WRITE(6,*) ' '
	  		CALL INIT_EVAPORA(RACINE) 
      ENDIF                                     
C ******  INIT of mass tables *************

          CALL INIPACE(RACINE)
          
C ******  INIT of PAWC for NTUPLE *********
	  WRITE(6,*) '*********************************************'
          WRITE(6,*) ' Output file for NTUPLE:'
         
          CALL HBINIT(Je_veux)
          
C ************************************************************************
C --------------  INIT of counters and variables. -----------------------
	ntrans=0
	nabs=0
	nretir=0
        nfis = 0
        muln_t=0.
        muln_2=0.
	ener_2=0.
        muln_20=0.
	ener_20=0.
        muln_max=0.
	ener_max=0.
	nbevhbk=0
	
	impulse1=0    ! conservation d'impulsion fausse sans fission
	impulse2=0
	impulse3=0
	imp_f1=0      ! conservation d'impulsion fausse avec fission
	imp_f2=0
        imp_f3=0
	
        S1x=0.
        S2x=0.
        S1y=0.
        S2y=0.
        S1z=0.
        S2z=0.
 
	
	IF(f(7).EQ.3.) THEN
	MPROJO=139.56995
	PBEAM=SQRT(f(3)*(f(3)+2.*MPROJO))  !pi+
	AP=0.D0
	ZP=1.D0
	ENDIF	  
	IF(f(7).EQ.4.) THEN
	MPROJO=134.9764
	PBEAM=SQRT(f(3)*(f(3)+2.*MPROJO))  !pi0
	AP=0.D0
	ZP=0.D0
	ENDIF
	IF(f(7).EQ.5.) THEN
	MPROJO=139.56995
	PBEAM=SQRT(f(3)*(f(3)+2.*MPROJO))  !pi-
	AP=0.D0
	ZP=-1.D0
	ENDIF
	
C Coulomb en entree seulement pour les particules ci-dessous
	IF(f(7).EQ.1.) THEN
	MPROJO=938.27231
	PBEAM=SQRT(f(3)*(f(3)+2.*MPROJO))   !proton
	AP=1.D0
	ZP=1.D0
	ENDIF
	IF(f(7).EQ.2.) THEN
	MPROJO=939.56563
	PBEAM=SQRT(f(3)*(f(3)+2.*MPROJO))   !neutron
	AP=1.D0
	ZP=0.D0
	ENDIF
	IF(f(7).EQ.6.) THEN
	MPROJO=1875.61276
	PBEAM=SQRT(f(3)*(f(3)+2.*MPROJO))  !deuton
	AP=2.D0
	ZP=1.D0
	ENDIF
	IF(f(7).EQ.7.) THEN
	MPROJO=2808.95
	PBEAM=SQRT(f(3)*(f(3)+2.*MPROJO))  !triton
	AP=3.D0
	ZP=1.D0
	ENDIF
	IF(f(7).EQ.8.) THEN
	MPROJO=2808.42
	PBEAM=SQRT(f(3)*(f(3)+2.*MPROJO))  !3He
	AP=3.D0
	ZP=2.D0
	ENDIF
	IF(f(7).EQ.9.) THEN
	MPROJO=3727.42
	PBEAM=SQRT(f(3)*(f(3)+2.*MPROJO))  !4He
	AP=4.D0
	ZP=2.D0
	ENDIF
        IF(f(7).EQ.-12.) THEN
           MPROJO=12.0*939.56563
           PBEAM=SQRT(f(3)*(f(3)+2.*MPROJO)) !12C
           AP=12.D0
           ZP=6.D0
        ENDIF
	       
	AT=f(1)
	ZT=f(2)
	EAP=f(3)
	
	f(4)=0.		!Seuil sortie proton
	f(8)=0.		!Seuil sortie neutron
	
C calculation of the transmission probability through the coulomb barrier
		APRO=AP
		ZPRO=ZP
		PROBA_TRANS=1.
	IF((f(7).EQ.1.).OR.(f(7).GE.6.)) THEN
	  CALL COULOMB_TRANSM(f(3),APRO,ZPRO,f(1),f(2),PROBA_TRANS)
	ENDIF
	  	ntrans_coul=icoup*(1.-PROBA_TRANS)
		icoup=icoup-ntrans_coul
          WRITE(6,*) 'icoup, proba_trans',icoup,proba_trans
C ************************************************************************
C                       Main loop on the number of shot
C ************************************************************************
      DO 100 I = 1,icoup
        INUM=I		!nouvel argument de EVAPO (V3)

        IF(i.EQ.ievtest) THEN
	    WRITE(6,*) 'Le tir selectionne commence: i=',i
	ENDIF
        
990	CONTINUE	!Retirage si NOPART= -100 (void event)
            
C stockage des generateurs pour reproduire les bugs...
	   ialview = ial
              	DO ihaz=1,19
              	   IYV(ihaz)=IY(ihaz)
              	ENDDO
              	
C          WRITE(6,*)'Event',i

      MTOTA = 0.D0		!Counters of A masses emitted as alpha                             
      MALPHA1 = 0.D0                            
      MALPHA2 = 0.D0
                                  
      AP1 = IDNINT(AP)                                                  
      ZP1 = IDNINT(ZP)                                                  


C*********************** Appel de la CASCADE ****************************            
	IBERT=1
	IF(I.EQ.1)IBERT=0
	IBERT=I
			
      CALL PNU(IBERT,F,NOPART,KIND,EP,ALPHA,BETA2,GAM,IZREM,
     #         IAREM,ESREM,ERECREM,ALREM,BEREM,GAREM,BIMPAC,JREM)

C Pour evapo directe (sans calcul d'INCL a basse energie)
c	NOPART=0
c	IAREM=f(1)+1
c	IZREM=f(2)+1
c	ESREM=f(3)+7.
c	ERECREM=PBEAM**2/(2.*938.2796*IAREM)
c	ALREM=0.00001
c	BEREM=0.0
c	GAREM=0.99999
c	BIMPACT=0.
c	JREM=0.





c***ouput***
C  NOPART=-1 PSEUDO REACTION (VOID EVENT)
C          0 ABSORPTION
C         >0 TRUE EVENT, = NUMBER OF PARTICLES EMITTED (EXCLUDING THE REMNANT)
C  FOR N=1,NOPART:
C  KIND(N)= TYPE OF PARTICLES (SAME CONVENTION AS FOR F(7), BUT IN INTEGERS)
C  EP(N)=  KINETIC ENERGY
C  ALPHA(N),BETA(N),GAM(N)= DIRECTION COSINES
C  IZREM= Z (REMNANT)
C  IAREM= A (REMNANT)
C  ESREM= EXCITATION ENERGY OF THE REMNANT
C  ERECREM= RECOIL ENERGY OF THE REMNANT
C  ALREM,BEREM,GAREM=DIRECTION COSINES OF THE REMNANT
c************************************************************************
C cascade imposee pour test p+Pb->p(5deg)+Pb(E*=200MeV)
C************************************************************************
C		NOPART=1
C		KIND(1)=1
C		EP(1)=799.835
C		ALPHA(1)=0.08716
C		BETA2(1)=0.
C		GAM(1)=0.99619
C		IZREM=82
C		IAREM=208
C		ESREM=200.
C		ERECREM=0.18870
C		ALREM=-0.47101
C		BEREM=0.
C		GAREM=0.88213
C		BIMPACT=2.

C	WRITE(6,*) 'cascade, NOPART, Z,A,E*,T remnant:', 
C     s              NOPART,IZREM,IAREM,ESREM,ERECREM
C      DO ipr=1,NOPART
C      WRITE(6,*) ipr,KIND(ipr),EP(ipr),ALPHA(ipr),BETA2(ipr),GAM(ipr)
C      ENDDO

C Absorption forcée pour p (10-100 MeV) et n (20-100MeV)
      CALL FORCE_ABSOR(nopart,F,IAREM,IZREM,ESREM,ERECREM,
     s  ALREM,BEREM,GAREM,JREM)

	IF(nopart.eq.-1)THEN
		ntrans=ntrans+1
		go to 99	! tirage non exploite, au suivant!
	ENDIF
	
	IF(nopart.eq.0)THEN
		nabs=nabs+1


C		go to 99	! tirage non exploite, au suivant!
C Et si! Il faut evaporer une absorption... (Important a basse energie 
C                                           en dessous de 200 MeV)
  	ENDIF
	
	IF(nopart.eq.-100)THEN
		nretir=nretir+1 
		go to 990	! Retirage
	ENDIF



C ----------------- OK, valuable cascade, we can continue --------------

      ZPRF = DBLE(IZREM)	!NUCLEAR CHARGE OF THE PREFRAGMENT
      APRF = DBLE(IAREM)	!MASS NUMBER OF THE PREFRAGMENT
      EE = DBLE(ESREM)		!EXCITATION ENERGY OF THE PREFRAGMENT
      JPRF=0.			!ANGULAR MOMENTUM OF THE PREFRAGMENT
C                                                                       
C MEAN ANGULAR MOMENTUM OF PREFRAGMENT                                  
C                                                                       
       JPRF = 0.165D0 * AT**(2.D0/3.D0) *                               
     &       APRF*(AT - APRF)/(AT - 1.D0)                               
       IF (JPRF.LT.0)JPRF = 0.D0
C    check M.de Jong, Ignatyuk, Schmidt Nuc.Phys A 613, pg442, 7th line
       JPRF = DSQRT(2*JPRF)

C      CALL HFILL(1,REAL(JREM),0.,1.)
C      CALL HFILL(2,SNGL(JPRF),0.,1.)
      
      JPRF=JREM
      JREMN=JREM		!JREM copie dans le NTUPLE


      NUMPI=0			!Compteurs de pions, neutrons protons
      MULTN=0
      MULTP=0

C  Ecriture dans le NTUPLE des particules de cascade (sauf REMNANT)      
      NTRACK=nopart		!Nombre de particules pour ce tir
      MASSINI=iarem
      MZINI=izrem
      EXINI=esrem
      BIMPACT=BIMPAC
C Three ways to compute the mass of the remnant: 
C		-From the output of the cascade and the canonic mass
C		-From energy balance (input - all emitted energies)
C		-Following the approximations of the Cugnon code (ESREM...)
      MCOREM = MPROJO +f(3) 
     s   +pace2(DBLE(f(1)),DBLE(f(2))) +f(1)*UMA -f(2)*MELEC

      PXBIL=0.
      PYBIL=0.
      PZBIL=0.         
      DO 900,j=1,nopart
	ITYPCASC(j)=1
c kind(): 1=proton, 2=neutron, 3=pi+, 4=pi0, 5=pi -	 
	if(kind(j).eq.1) then
		AVV(j)=1
		ZVV(j)=1
C		PLAB(j)=SQRT(EP(j)*(EP(j)+1876.54462))     !OK
		    PLAB(j)=SQRT(EP(j)*(EP(j)+1876.5592))  !Cugnon
		MULTP=MULTP+1
		MCOREM=MCOREM -EP(j) -938.27231
	endif
	if(kind(j).eq.2) then
		AVV(j)=1
		ZVV(j)=0
C		PLAB(j)=SQRT(EP(j)*(EP(j)+1879.13126))	   !OK	
		    PLAB(j)=SQRT(EP(j)*(EP(j)+1876.5592))  !Cugnon
		MULTN=MULTN+1
		MCOREM=MCOREM -EP(j) -939.56563
	endif	
	if(kind(j).eq.3) then
		AVV(j)=-1
		ZVV(j)=1
C		PLAB(j)=SQRT(EP(j)*(EP(j)+279.1399))       !OK
	  	    PLAB(j)=SQRT(EP(j)*(EP(j)+276.0))      !Cugnon
		NUMPI=NUMPI+1
		MCOREM=MCOREM -EP(j) -139.56995
	endif
	if(kind(j).eq.4) then
		AVV(j)=-1
		ZVV(j)=0
C		PLAB(j)=SQRT(EP(j)*(EP(j)+269.9528))       !OK
		    PLAB(j)=SQRT(EP(j)*(EP(j)+276.0))      !Cugnon
		NUMPI=NUMPI+1
		MCOREM=MCOREM -EP(j) -134.9764
	endif
	if(kind(j).eq.5) then
		AVV(j)=-1
		ZVV(j)=-1
C		PLAB(j)=SQRT(EP(j)*(EP(j)+279.1399))       !OK
		    PLAB(j)=SQRT(EP(j)*(EP(j)+276.0))      !Cugnon
		NUMPI=NUMPI+1
		MCOREM=MCOREM -EP(j) -139.56995
	endif
	if(kind(j).eq.6) then
		AVV(j)=2
		ZVV(j)=1
C		PLAB(j)=SQRT(EP(j)*(EP(j)+279.1399))       !OK
		    PLAB(j)=SQRT(EP(j)*(EP(j)+2.*1874.34))      !Cugnon
		NUMPI=NUMPI+1
		MCOREM=MCOREM -EP(j) -1874.34
	endif
	if(kind(j).eq.7) then
		AVV(j)=3
		ZVV(j)=1
C		PLAB(j)=SQRT(EP(j)*(EP(j)+279.1399))       !OK
		    PLAB(j)=SQRT(EP(j)*(EP(j)+2.*2806.359))      !Cugnon
		NUMPI=NUMPI+1
		MCOREM=MCOREM -EP(j) -2806.359
	endif
	if(kind(j).eq.8) then
		AVV(j)=3
		ZVV(j)=2
C		PLAB(j)=SQRT(EP(j)*(EP(j)+279.1399))       !OK
		    PLAB(j)=SQRT(EP(j)*(EP(j)+2.*2807.119))      !Cugnon
		NUMPI=NUMPI+1
		MCOREM=MCOREM -EP(j) -2807.119
	endif
	if(kind(j).eq.9) then
		AVV(j)=4
		ZVV(j)=2
C		PLAB(j)=SQRT(EP(j)*(EP(j)+279.1399))       !OK
		    PLAB(j)=SQRT(EP(j)*(EP(j)+2.*3724.818))      !Cugnon
		NUMPI=NUMPI+1
		MCOREM=MCOREM -EP(j) -3724.818
	endif
	ENERJ(j)=EP(j)
	TETLAB(j)=180.*ACOS(GAM(j))/3.141592654
	PHILAB(j)=180.*ATAN2(BETA2(j),ALPHA(j))/3.141592654
	
	PXBIL = PXBIL+PLAB(j)*ALPHA(j)
	PYBIL = PYBIL+PLAB(j)*BETA2(j)
	PZBIL = PZBIL+PLAB(j)*GAM(j)
	
C	WRITE(6,*) 'cascade, ntpl:',AVV(j),ZVV(j)
C     s             ,ENERJ(j),PLAB(j),TETLAB(j),PHILAB(j)

900   CONTINUE

C calcul de la masse (impulsion) du remnant coherente avec la conservation d'energie:
C        MCOREM = MCOREM -ERECREM -ESREM			!OK
C        PCOREM=SQRT(ERECREM*(ERECREM+2.*(MCOREM+ESREM)))	!OK
        	PCOREM=SQRT(ERECREM*(ERECREM +2.*938.2796*IAREM))  !Cugnon
C        	PCOREM=SQRT(ERECREM*2.*938.2796*IAREM)  !Cugnon
        	MCOREM=938.2796*IAREM 				!Cugnon
		
C Note: Il faut negliger l'energie d'excitation (ESREM) pour que le bilan 
C d'impulsion soit correct a la sortie de la cascade.....et prendre la
C masse MCOREM comme ci-dessus (fausse de ~1GeV par rapport aux tables...)        
        PXREM=PCOREM*ALREM
        PYREM=PCOREM*BEREM
        PZREM=PCOREM*GAREM
        
        PXBIL=PXBIL+PXREM
        PYBIL=PYBIL+PYREM
        PZBIL=PZBIL+PZREM
	
        
C        WRITE(6,*)'pxyzbil,pxyzcorem',PXBIL,PYBIL,PZBIL,
C     s                                PXREM,PYREM,PZREM,MCOREM,ESREM

C pour test, somme des impulsions
C        S1x=S1x+PXBIL
C        S1y=S1y+PYBIL
C        S1z=S1z+PZBIL
C        S2x=S2x+PXBIL**2
C        S2y=S2y+PYBIL**2
C        S2z=S2z+PZBIL**2
        
        
C Checkind longitudinal momentum conservation at this step.
C... and the transvers 9/2002
        IF((ABS(PZBIL-PBEAM).GT.5.)
     s                    .OR.(SQRT(PXBIL**2+PYBIL**2).GE.3.)) THEN
           WRITE(6,*) 'Bad momentum conservation after INCL:'
           WRITE(6,*) 'I,ERECREM,ESREM,IAREM,ALREM,BEREM,GAREM',
     s                 I,ERECREM,ESREM,IAREM,ALREM,BEREM,GAREM,PZBIL
           WRITE(6,*) 'kind,ep,a,b,g',
     s     (KIND(k),EP(k),ALPHA(k),BETA2(k),GAM(k),k=1,NOPART)    

        ENDIF
       
C --------------------------------------------------

	iv=0		!Init du compteur des part evaporees
	FF=-1		!Indice de fission (1) ou d'evaporation (0)
	                !initialise a -1 pour les cas sans fission ni evapo
	KFIS=0		!Drapeau de fission copie dans le NTUPLE
	ESTFIS=0.
	IZFIS=0
	IAFIS=0
	
C      IF ((ZPRF.GT.0).AND.(APRF.GT.0).AND.(ESREM.GT.0.01)) THEN

C **************************************************************************
C                            EVAPORATION step                                       
C **************************************************************************                            
      IF (choice_evap.EQ.2) THEN
C **************************  EVAPORATION GEM *********
          EE=EE/APRF 
      ERECI=DBLE(ERECREM) 
      BET0I(1)=DBLE(ALREM)
      BET0I(2)=DBLE(BEREM)
      BET0I(3)=DBLE(GAREM)
      WT=1.
      NOCASI=INUM
      
      
      	                                 
C      CALL GEM(IAREM,IZREM,EE,ERECI,BET0I,WT,NOCASI,
C     s  IZF_AB,IAF_AB,EKINF,BET1F,NBPART)
      
C 
C GEM donne NBPART particules (<300) dans les tableaux:
C    IZF_AB    Z
C    IAF_AB    A
C    EKINF     Kinetic energy
C    BET1F(300,3)  Cosinus directeurs                                                                      
C 

C Copie des resultats dans le NTUPLE
      NTRACK=NBPART+NOPART
      DO i33=1,NBPART
         AVV(NOPART+i33)=IAF_AB(i33)
	 ZVV(NOPART+i33)=IZF_AB(i33)
	 ITYPCASC(NOPART+i33)=0
	 ENERJ(NOPART+i33)=EKINF(i33)
	 TETLAB(NOPART+i33)=180.*ACOS(BET1F(i33,3))/3.141592654                                                                      
	 PHILAB(NOPART+i33)=180.*ATAN2(BET1F(i33,2),BET1F(i33,1))
     s	 /3.141592654
	 ALOC=IAF_AB(i33)
	 ZLOC=IZF_AB(i33)
         AVVMASS=pace2(DBLE(ALOC),DBLE(ZLOC)) 
     s                             + ALOC*UMA - ZLOC*MELEC
         PLAB(NOPART+i33)=SQRT(EKINF(i33)*(EKINF(i33)+2.*AVVMASS))	 
      ENDDO
      
      iv=NBPART
C *************************  FIN EVPORATION GEM ***************      
      ELSE	  
C *************************  EVAPORATION KHS *********	  

C$$$      IF(ESREM.GE.1.E-3) THEN	       
                                                
C$$$      CALL EVAPORA(ZPRF,APRF,EE,JPRF,ZF,AF,MTOTA,PLEVA,PXEVA,           
C$$$     &             PYEVA,FF,INTTYPE,INUM)
       
C$$$      ELSE
C$$$         FF=0
C$$$	 ZF=ZPRF
C$$$	 AF=APRF
C$$$	 PXEVA=PXREM
C$$$	 PYEVA=PYREM
C$$$	 PLEVA=PZREM
C$$$      ENDIF

     
C$$$C                                                                       
C$$$C AFP,ZFP is the final fragment if no fission occurs (FF=0)                   
C$$$C In case of fission (FF=1) it is the nucleus that undergoes fission.          
C$$$C                                                                       
C$$$           ZFP = IDNINT(ZF)                                                  
C$$$           AFP = IDNINT(AF)
                                              

C$$$      IF (FF.EQ.1) THEN   
C$$$C ---------------------  Here, a FISSION occures --------------------------
C$$$C                                                                       
C$$$C FEE: (EE) energy of fissioning nucleus ABOVE the fission barrier.          
C$$$C                                                                       
C$$$             nfis = nfis +1
C$$$             FEE = EE                                                                                                    
C$$$	KFIS=1		!Drapeau de fission copie dans le NTUPLE
C$$$C
C$$$C  calcul des impulsions des particules evaporees (avant fission) 
C$$$C                dans le systeme labo:
C$$$c
C$$$      TREM = DBLE(ERECREM)
C$$$C      REMMASS = pace2(APRF,ZPRF) + APRF*UMA - ZPRF*MELEC	!Canonic
C$$$C      REMMASS = MCOREM  + DBLE(ESREM)				!OK
C$$$C      REMMASS = MCOREM						!Cugnon
C$$$C      GAMREM = (REMMASS + TREM)/REMMASS
C$$$C      ETREM = DSQRT(TREM*(TREM + 2.*REMMASS))/REMMASS

C$$$C This is not treated as accurately as for the non fission case for which
C$$$C the remnant mass is computed to satisfy the energy conservation 
C$$$C of evaporated particles. But it is not bad and more canonical!      
C$$$      REMMASS = pace2(APRF,ZPRF) + APRF*UMA - ZPRF*MELEC+DBLE(ESREM) !canonic
C$$$C Essais avec la masse de KHS (9/2002):
C$$$	   CALL MGLMS(APRF,ZPRF,0,EL)
C$$$      REMMASS = ZPRF*FMP + (APRF-ZPRF)*FMN + EL + DBLE(ESREM)
       

C$$$      GAMREM=DSQRT(PCOREM**2+REMMASS**2)/REMMASS
C$$$      ETREM=PCOREM/REMMASS
      
C$$$      CSREM(1)=ALREM
C$$$      CSREM(2)=BEREM
C$$$      CSREM(3)=GAREM
      
C$$$C Pour Vérif Remnant = evapo(Pre fission) + Noyau_fissionant (système  Remnant)
C$$$	Bil_E=0.
C$$$	Bil_Px=0.
C$$$	Bil_Py=0.
C$$$	Bil_Pz=0.
C$$$      DO iloc=1,iv
C$$$        CALL MGLMS(DBLE(acv(iloc)),DBLE(zpcv(iloc)),0,EL)
C$$$	MASSE = zpcv(iloc)*FMP + (acv(iloc)-zpcv(iloc))*FMN + EL
C$$$	Bil_E = Bil_E + DSQRT(pcv(iloc)**2 + MASSE**2)
C$$$	Bil_Px = Bil_Px	+ pcv(iloc)*xcv(iloc)
C$$$	Bil_Py = Bil_Py	+ pcv(iloc)*ycv(iloc)
C$$$	Bil_Pz = Bil_Pz	+ pcv(iloc)*zcv(iloc)
C$$$      ENDDO
C$$$C Ce bilan (impulsion nulle) est parfait. (Bil_Px=Bil_Px+PXEVA....)

C$$$        NDEC = 1
C$$$      IF(iv.NE.0) THEN
C$$$        CALL TRANSLAB(GAMREM,ETREM,CSREM,NOPART,NDEC)
C$$$      ENDIF          
C$$$	NBPEVAP = iv	!Nombre de particules d'evaporation traitees
C$$$C                                                                       
C$$$C Now calculation of the fission fragment distribution including                  
C$$$C evaporation from the fragments.                                           
C$$$C                                   

C$$$C Distribution of the fission fragments:
                                                                       
C$$$         CALL FISSION_DISTRI(SNGL(AF),SNGL(ZF),SNGL(EE),AFF1,           
C$$$     &        ZFF1,EFF1,AFF2,ZFF2,EFF2)
C$$$C verif des A et Z decimaux:
C$$$      NA_F=AF+0.5
C$$$      NZ_F=ZF+0.5
C$$$      	IZFIS=NZ_F	!Copie dans le NTUPLE
C$$$      	IAFIS=NA_F
C$$$      NA_PF1=AFF1+0.5
C$$$      NZ_PF1=ZFF1+0.5     
C$$$      NA_PF2=AFF2+0.5
C$$$      NZ_PF2=ZFF2+0.5
C$$$      IF((NA_F.NE.(NA_PF1+NA_PF2)).OR.(NZ_F.NE.(NZ_PF1+NZ_PF2)))
C$$$     s    THEN
C$$$	  WRITE(6,*) 'Problemes arrondis dans la fission'
C$$$	  WRITE(6,*) 'AF,ZF,AFF1,ZFF1,AFF2,ZFF2',
C$$$     s                AF,ZF,AFF1,ZFF1,AFF2,ZFF2    
C$$$	  WRITE(6,*) 'A,Z,A1,Z1,A2,Z2 integer',
C$$$     s                NA_F,NZ_F,NA_PF1,NZ_PF1,NA_PF2,NZ_PF2 
C$$$      ENDIF   
C$$$C Calcul de l'impulsion des PF dans le syteme noyau de fission:
C$$$           Kboud = IDNINT(ZF)                                                  
C$$$           Jboud = IDNINT(AF-ZF)                                             
C$$$           EF = EFA(Kboud,Jboud)	!barriere de fission
C$$$           ESTFIS=EE+EF   		!Copie dans le NTUPLE   
           
     
C$$$C           MASSEF = pace2(AF,ZF)
C$$$C      	   MASSEF = MASSEF + AF*UMA - ZF*MELEC + EE + EF
C$$$C           MASSE1 = pace2(DBLE(AFF1),DBLE(ZFF1))
C$$$C      	   MASSE1 = MASSE1 + AFF1*UMA - ZFF1*MELEC + EFF1
C$$$C           MASSE2 = pace2(DBLE(AFF2),DBLE(ZFF2))
C$$$C      	   MASSE2 = MASSE2 + AFF2*UMA - ZFF2*MELEC + EFF2
C$$$C        WRITE(6,*) 'MASSEF,MASSE1,MASSE2',MASSEF,MASSE1,MASSE2
C$$$C MGLMS est la fonction de masse cohérente avec KHS evapo-fis.
C$$$C   Attention aux parametres, ici 0=OPTSHP, NO microscopic correct. 
C$$$	   CALL MGLMS(AF,ZF,0,EL)
C$$$	   MASSEF = ZF*FMP + (AF-ZF)*FMN + EL + EE + EF
C$$$	   CALL MGLMS(DBLE(AFF1),DBLE(ZFF1),0,EL)
C$$$      	   MASSE1 = ZFF1*FMP + (AFF1-ZFF1)*FMN + EL + EFF1
C$$$	   CALL MGLMS(DBLE(AFF2),DBLE(ZFF2),0,EL)
C$$$      	   MASSE2 = ZFF2*FMP + (AFF2-ZFF2)*FMN + EL + EFF2
C$$$C        WRITE(6,*) 'MASSEF,MASSE1,MASSE2',MASSEF,MASSE1,MASSE2	   
C$$$           B = MASSEF-MASSE1-MASSE2
C$$$	   IF(B.LT.0.) THEN
C$$$	   	B=0.
C$$$		WRITE(6,*) 'anomalie dans la fission:', 
C$$$     s     inum,AF,ZF,massef,AFF1,ZFF1,masse1,AFF2,ZFF2,masse2
C$$$	   ENDIF
C$$$           T1=B*(B+2.*MASSE2)/(2.*MASSEF)
C$$$           P1 = DSQRT(T1*(T1+2.*MASSE1))
           
C$$$           CALL RIBM(rndm,IY(14))
C$$$           CTET1 = 2.*rndm-1.
C$$$           CALL RIBM(rndm,IY(10))
C$$$           PHI1 = rndm*2.*3.141592654
           
C$$$C ----Coefs de la transformation de Lorentz (noyau de fission -> Remnant) 
C$$$      PEVA = PXEVA**2+PYEVA**2+PLEVA**2
C$$$      GAMFIS = DSQRT(MASSEF**2 + PEVA)/MASSEF
C$$$      PEVA = DSQRT(PEVA)
C$$$      ETFIS = PEVA/MASSEF
      
C$$$C ----Matrice de rotation (noyau de fission -> Remnant)
C$$$      SITET = 0.
C$$$      IF(PEVA.GE.1.E-4)SITET = SQRT(PXEVA**2+PYEVA**2)/PEVA
C$$$      IF(SITET.GT.1.E-5)THEN
C$$$        CSTET = PLEVA/PEVA
C$$$        SIPHI = PYEVA/(SITET*PEVA)
C$$$        CSPHI = PXEVA/(SITET*PEVA)
	
C$$$	R(1,1) = CSTET*CSPHI
C$$$	R(1,2) = -SIPHI
C$$$	R(1,3) = SITET*CSPHI
C$$$	R(2,1) = CSTET*SIPHI
C$$$	R(2,2) = CSPHI
C$$$	R(2,3) = SITET*SIPHI
C$$$	R(3,1) = -SITET
C$$$	R(3,2) = 0.
C$$$	R(3,3) = CSTET
C$$$      ELSE
C$$$734	R(1,1) = 1.
C$$$	R(1,2) = 0.
C$$$	R(1,3) = 0.
C$$$	R(2,1) = 0.
C$$$	R(2,2) = 1.
C$$$	R(2,3) = 0.
C$$$	R(3,1) = 0.
C$$$	R(3,2) = 0.
C$$$	R(3,3) = 1.
C$$$      ENDIF
	           
C$$$c test de verif:                                      
C$$$         IF( (ZFF1.LE.0.D0).OR.(AFF1.LE.0.D0).OR.(AFF1.LT.ZFF1)) THEN   
C$$$                        WRITE(6,*) ZF,AF,EE,ZFF1,AFF1                                
C$$$         ELSE
                                                                    
C$$$C ---------------------- PF1 will evaporate 
C$$$         EPF1_IN=DBLE(EFF1)
C$$$	 EPF1_OUT=EPF1_IN
C$$$         CALL EVAPORA(DBLE(ZFF1),DBLE(AFF1),EPF1_OUT,0.D0,            
C$$$     &   ZF1,AF1,MALPHA1,FFPLEVA1,FFPXEVA1,FFPYEVA1,FF1,FTYPE1,INUM)
     
C$$$C On ajoute le fragment:
C$$$         iv = iv +1
C$$$         acv(iv) = AF1
C$$$         zpcv(iv) = ZF1        
C$$$         PEVA = DSQRT(FFPXEVA1**2+FFPYEVA1**2+FFPLEVA1**2)
C$$$         pcv(iv) = PEVA
C$$$	IF(PEVA.GT.0.001) THEN
C$$$         xcv(iv) = FFPXEVA1/PEVA
C$$$         ycv(iv) = FFPYEVA1/PEVA
C$$$         zcv(iv) = FFPLEVA1/PEVA 
C$$$        ELSE
C$$$         xcv(iv)=1.
C$$$         ycv(iv)=0.
C$$$         zcv(iv)=0.
C$$$        END IF
	        
C$$$C Pour Vérif evapo de PF1 dans le systeme du Noyau Fissionant
C$$$	Bil1_E=0.
C$$$	Bil1_Px=0.
C$$$	Bil1_Py=0.
C$$$	Bil1_Pz=0.
C$$$      DO iloc=NBPEVAP+1,iv
C$$$	CALL MGLMS(DBLE(acv(iloc)),DBLE(zpcv(iloc)),0,EL)
C$$$      	MASSE = zpcv(iloc)*FMP + (acv(iloc)-zpcv(iloc))*FMN + EL 
C$$$	Bil1_E = Bil1_E + DSQRT(pcv(iloc)**2 + MASSE**2)
C$$$	Bil1_Px = Bil1_Px	+ pcv(iloc)*xcv(iloc)
C$$$	Bil1_Py = Bil1_Py	+ pcv(iloc)*ycv(iloc)
C$$$	Bil1_Pz = Bil1_Pz	+ pcv(iloc)*zcv(iloc)
C$$$      ENDDO
       
C$$$C ----Calcul des cosinus directeurs de PF1 dans le Remnant et calcul
C$$$c des coefs pour la transformation de Lorentz Systeme PF --> Systeme Remnant
	
C$$$        CALL TRANSLABPF(MASSE1,T1,P1,CTET1,PHI1,GAMFIS,ETFIS,R,
C$$$     s   PLAB1,GAM1,ETA1,CSDIR1)
     
C$$$C
C$$$C  calcul des impulsions des particules evaporees dans le systeme Remnant:
C$$$c
C$$$         CALL TRANSLAB(GAM1,ETA1,CSDIR1,NOPART,NBPEVAP+1)
C$$$         MEMIV = NBPEVAP+1		!Memoires pour la future transformation
C$$$         MEMPAW = NOPART		!remnant->labo pour PF1 ET PF2.
C$$$         LMI_PF1 = NOPART + NBPEVAP+1	!indices min et max dans /VAR_NTP/
C$$$	 LMA_PF1 = NOPART + iv		! des particules issues de PF1
C$$$	 NBPEVAP = iv	!Nombre de particules d'evaporation traitees

C$$$         END IF
C$$$C --------------------- End of PF1 calculation

C$$$c test de verif:                                                                                                         
C$$$         IF( (ZFF2.LE.0.D0).OR.(AFF2.LE.0.D0).OR.(AFF2.LT.ZFF2)) THEN   
C$$$           		WRITE(6,*) ZF,AF,EE,ZFF2,AFF2                                
C$$$         ELSE                                                           
                                                                    
C$$$C ---------------------- PF2 will evaporate 
C$$$         EPF2_IN=DBLE(EFF2)
C$$$	 EPF2_OUT=EPF2_IN
C$$$         CALL EVAPORA(DBLE(ZFF2),DBLE(AFF2),EPF2_OUT,0.D0,            
C$$$     &   ZF2,AF2,MALPHA2,FFPLEVA2,FFPXEVA2,FFPYEVA2,FF2,FTYPE2,INUM)        
     
C$$$C On ajoute le fragment:
C$$$         iv = iv +1
C$$$         acv(iv) = AF2
C$$$         zpcv(iv) = ZF2        
C$$$         PEVA = DSQRT(FFPXEVA2**2+FFPYEVA2**2+FFPLEVA2**2)
C$$$         pcv(iv) = PEVA
C$$$	IF(PEVA.GT.0.001) THEN
C$$$         xcv(iv) = FFPXEVA2/PEVA
C$$$         ycv(iv) = FFPYEVA2/PEVA
C$$$         zcv(iv) = FFPLEVA2/PEVA 
C$$$        ELSE
C$$$         xcv(iv)=1.
C$$$         ycv(iv)=0.
C$$$         zcv(iv)=0.
C$$$        END IF        
C$$$C Pour Vérif evapo de PF1 dans le systeme du Noyau Fissionant
C$$$	Bil2_E=0.
C$$$	Bil2_Px=0.
C$$$	Bil2_Py=0.
C$$$	Bil2_Pz=0.
C$$$      DO iloc=NBPEVAP+1,iv
C$$$	CALL MGLMS(DBLE(acv(iloc)),DBLE(zpcv(iloc)),0,EL)
C$$$      	MASSE = zpcv(iloc)*FMP + (acv(iloc)-zpcv(iloc))*FMN + EL 
C$$$	Bil2_E = Bil2_E + DSQRT(pcv(iloc)**2 + MASSE**2)
C$$$	Bil2_Px = Bil2_Px	+ pcv(iloc)*xcv(iloc)
C$$$	Bil2_Py = Bil2_Py	+ pcv(iloc)*ycv(iloc)
C$$$	Bil2_Pz = Bil2_Pz	+ pcv(iloc)*zcv(iloc)
C$$$      ENDDO

C$$$C ----Calcul des cosinus directeurs de PF2 dans le Remnant et calcul
C$$$c des coefs pour la transformation de Lorentz Systeme PF --> Systeme Remnant
C$$$     	T2 = B - T1
C$$$     	CTET2 = -CTET1
C$$$     	PHI2 = DMOD(PHI1+3.141592654,6.283185308D+00)
C$$$     	P2 = DSQRT(T2*(T2+2.*MASSE2))
     
C$$$        CALL TRANSLABPF(MASSE2,T2,P2,CTET2,PHI2,GAMFIS,ETFIS,R,
C$$$     s   PLAB2,GAM2,ETA2,CSDIR2)

C$$$C
C$$$C  calcul des impulsions des particules evaporees dans le systeme Remnant:
C$$$c
C$$$        CALL TRANSLAB(GAM2,ETA2,CSDIR2,NOPART,NBPEVAP+1 )
C$$$         LMI_PF2 = NOPART + NBPEVAP+1	!indices min et max dans /VAR_NTP/
C$$$	 LMA_PF2 = NOPART + iv		! des particules issues de PF2

C$$$        END IF
C$$$C --------------------- End of PF2 calculation

C$$$C Pour vérifications: calculs du noyau fissionant et des PF dans 
C$$$C    le systeme du remnant.
C$$$      DO iloc=1,3
C$$$        PFIS_REM(iloc)=0.
C$$$      ENDDO 
C$$$      CALL LOR_AB(GAMFIS,ETFIS,MASSEF,PFIS_REM,EFIS_REM,PFIS_TRAV)
C$$$      CALL ROT_AB(R,PFIS_TRAV,PFIS_REM)
C$$$C      WRITE(6,*) (PFIS_REM(iloc),iloc=1,3)
C$$$C      WRITE(6,*) CTET1,CTET2,PHI1,PHI2,P1,P2
      
C$$$      STET1=SQRT(1.-CTET1**2)
C$$$      PF1_REM(1)=P1*STET1*COS(PHI1)
C$$$      PF1_REM(2)=P1*STET1*SIN(PHI1)
C$$$      PF1_REM(3)=P1*CTET1
C$$$      CALL LOR_AB(GAMFIS,ETFIS,MASSE1+T1,PF1_REM,E1_REM,PFIS_TRAV)
C$$$      CALL ROT_AB(R,PFIS_TRAV,PF1_REM)  

C$$$      STET2=SQRT(1.-CTET2**2)
C$$$      PF2_REM(1)=P2*STET2*COS(PHI2)
C$$$      PF2_REM(2)=P2*STET2*SIN(PHI2)
C$$$      PF2_REM(3)=P2*CTET2
C$$$      CALL LOR_AB(GAMFIS,ETFIS,MASSE2+T2,PF2_REM,E2_REM,PFIS_TRAV)
C$$$      CALL ROT_AB(R,PFIS_TRAV,PF2_REM)
C$$$C Verif 0: Remnant = evapo_pre_fission + Noyau Fissionant
C$$$	Bil_E = REMMASS - EFIS_REM - Bil_E
C$$$	Bil_Px = Bil_Px + PFIS_REM(1)  
C$$$	Bil_Py = Bil_Py + PFIS_REM(2)  
C$$$	Bil_Pz = Bil_Pz + PFIS_REM(3)  
C$$$C Verif 1: noyau fissionant = PF1 + PF2 dans le systeme remnant
C$$$	Bilan_E = EFIS_REM - E1_REM - E2_REM
C$$$	Bilan_PX = PFIS_REM(1) - PF1_REM(1) - PF2_REM(1)
C$$$	Bilan_PY = PFIS_REM(2) - PF1_REM(2) - PF2_REM(2)
C$$$	Bilan_PZ = PFIS_REM(3) - PF1_REM(3) - PF2_REM(3)
C$$$C Verif 2: PF1 et PF2 egaux a toutes leurs particules evaporees
C$$$C   (Systeme remnant)
C$$$      IF((LMA_PF1-LMI_PF1).NE.0) THEN
C$$$	Bil_E_PF1 = E1_REM - EPF1_OUT
C$$$	Bil_PX_PF1 = PF1_REM(1) 
C$$$	Bil_PY_PF1 = PF1_REM(2) 
C$$$	Bil_PZ_PF1 = PF1_REM(3)
C$$$        DO ipf1=LMI_PF1,LMA_PF1
C$$$		 Bil_E_PF1 = Bil_E_PF1 
C$$$     s		 - (PLAB(ipf1)**2 + ENERJ(ipf1)**2)/(2.*ENERJ(ipf1))
C$$$		 CST = COS(TETLAB(ipf1)/57.2957795)
C$$$		 SST = SIN(TETLAB(ipf1)/57.2957795)
C$$$		 CSF = COS(PHILAB(ipf1)/57.2957795)
C$$$		 SSF = SIN(PHILAB(ipf1)/57.2957795)
C$$$		 Bil_PX_PF1 = Bil_PX_PF1 - PLAB(ipf1)*SST*CSF
C$$$		 Bil_PY_PF1 = Bil_PY_PF1 - PLAB(ipf1)*SST*SSF
C$$$		 Bil_PZ_PF1 = Bil_PZ_PF1 - PLAB(ipf1)*CST		 
C$$$        ENDDO
C$$$	ENDIF
	 
C$$$      IF((LMA_PF2-LMI_PF2).NE.0) THEN
C$$$	Bil_E_PF2 =  E2_REM - EPF2_OUT
C$$$	Bil_PX_PF2 = PF2_REM(1) 
C$$$	Bil_PY_PF2 = PF2_REM(2) 
C$$$	Bil_PZ_PF2 = PF2_REM(3)
C$$$        DO ipf2=LMI_PF2,LMA_PF2
C$$$		 Bil_E_PF2 = Bil_E_PF2 
C$$$     s		 - (PLAB(ipf2)**2 + ENERJ(ipf2)**2)/(2.*ENERJ(ipf2))
C$$$		 CST = COS(TETLAB(ipf2)/57.2957795)
C$$$		 SST = SIN(TETLAB(ipf2)/57.2957795)
C$$$		 CSF = COS(PHILAB(ipf2)/57.2957795)
C$$$		 SSF = SIN(PHILAB(ipf2)/57.2957795)
C$$$		 Bil_PX_PF2 = Bil_PX_PF2 - PLAB(ipf2)*SST*CSF
C$$$		 Bil_PY_PF2 = Bil_PY_PF2 - PLAB(ipf2)*SST*SSF
C$$$		 Bil_PZ_PF2 = Bil_PZ_PF2 - PLAB(ipf2)*CST		 
C$$$        ENDDO
C$$$	ENDIF 
C$$$C
C$$$C ---- Transformation systeme Remnant -> systeme labo. (evapo de PF1 ET PF2)
C$$$C
C$$$	CALL TRANSLAB(GAMREM,ETREM,CSREM,MEMPAW,MEMIV)
	
C$$$C *******************  END of fission calculations ************************

C$$$      ELSE
       
C$$$C ************************ Evapo sans fission *****************************
C$$$C Here, FF=0, --> Evapo sans fission, on ajoute le fragment:
C$$$C *************************************************************************
C$$$         iv = iv +1
C$$$         acv(iv) = AF
C$$$         zpcv(iv) = ZF
C$$$         PEVA = DSQRT(PXEVA**2+PYEVA**2+PLEVA**2)
C$$$         pcv(iv) = PEVA
C$$$	IF(PEVA.GT.0.001) THEN
C$$$         xcv(iv) = PXEVA/PEVA
C$$$         ycv(iv) = PYEVA/PEVA
C$$$         zcv(iv) = PLEVA/PEVA        
C$$$        ELSE
C$$$         xcv(iv)=1.
C$$$         ycv(iv)=0.
C$$$         zcv(iv)=0.
C$$$        END IF        
	
C$$$C
C$$$C  calcul des impulsions des particules evaporees dans le systeme labo:
C$$$c
C$$$      TREM = DBLE(ERECREM)
C$$$C      REMMASS = pace2(APRF,ZPRF) + APRF*UMA - ZPRF*MELEC	!Canonic
C$$$C      REMMASS = MCOREM  + DBLE(ESREM)				!OK
C$$$      REMMASS = MCOREM						!Cugnon
C$$$C      GAMREM = (REMMASS + TREM)/REMMASS			!OK
C$$$C      ETREM = DSQRT(TREM*(TREM + 2.*REMMASS))/REMMASS		!OK
C$$$      CSREM(1)=ALREM
C$$$      CSREM(2)=BEREM
C$$$      CSREM(3)=GAREM
      

C$$$      e_evapo=0.
C$$$      DO j=1,iv
C$$$        CALL MGLMS(DBLE(acv(j)),DBLE(zpcv(j)),0,EL)
C$$$	fmcv = zpcv(j)*FMP + (acv(j)-zpcv(j))*FMN + EL
C$$$	e_evapo = e_evapo + DSQRT(pcv(j)**2 + fmcv**2)
C$$$      ENDDO
      
C$$$C Redefinition pour conservation d'impulsion!!!
C$$$C   this mass obtained by energy balance is very close to the
C$$$C   mass of the remnant computed by pace2 + excitation energy (EE). (OK)      
C$$$      REMMASS = e_evapo
      
C$$$      GAMREM=DSQRT(PCOREM**2+REMMASS**2)/REMMASS
C$$$      ETREM=PCOREM/REMMASS
      
C$$$        CALL TRANSLAB(GAMREM,ETREM,CSREM,NOPART,1)
                  
C$$$C End of the (FISSION - NO FISSION) condition (FF=1 or 0)                                          
C$$$      END IF 
C *********************** FIN de l'EVAPO KHS ******************** 
                                                          
      ENDIF 	!choix de l'evaporation (KHS-GEM)  
                                        
C                                                                       
C *************************************************************************                                                                      
C                         FIN DE L'EVAPORATION 
C *************************************************************************

C      ELSE	! Evaporation impossible
      	NTRACK=NTRACK+1		! on recopie le remnant dans le ntuple
      	ITYPCASC(NTRACK)=1
      	AVV(NTRACK)=IAREM
      	ZVV(NTRACK)=IZREM
      	PLAB(NTRACK)=PCOREM
      	ENERJ(NTRACK)=SQRT(PCOREM**2+MCOREM**2)-MCOREM
	TETLAB(NTRACK)=180.*ACOS(GAREM)/3.141592654
	PHILAB(NTRACK)=180.*ATAN2(BEREM,ALREM)/3.141592654

C      END IF	! Fin du test evapo possible
      
      	  
C count of n and p during the evaporation step, number of particles ...          
	MULTEN=0
	MULTEP=0
	
        IF(iv.gt.0) THEN
	
      IF (choice_evap.EQ.2) THEN
        FF=0				!Pour les bilans en impuls
        DO kcv=NOPART+1,NOPART+iv  	!GEM   
      if((AVV(kcv).eq.1).and.(ZVV(kcv).eq.0))MULTEN=MULTEN+1
      if((AVV(kcv).eq.1).and.(ZVV(kcv).eq.1))MULTEP=MULTEP+1
 	ENDDO
      ELSE        
        DO kcv=1,iv			!KHS
 	   NTRACK=NTRACK+1    !crucial...add the evaporated particles  	   
      if((AVV(NTRACK).eq.1).and.(ZVV(NTRACK).eq.0))MULTEN=MULTEN+1
      if((AVV(NTRACK).eq.1).and.(ZVV(NTRACK).eq.1))MULTEP=MULTEP+1
 	ENDDO
      ENDIF
       	
 	ENDIF
	  
          IZTOT = 0
          IATOT = 0
          

      multeen=MULTN+MULTEN
      MULNCASC=MULTN
      MULNEVAP=MULTEN
      MULNTOT=multeen     
      multeep=MULTP+MULTEP      

C ---- Various checks on conservation numbers -------------------
      
      Iabil=0
      Izbil=0
      PXBIL=0.
      PYBIL=0.
      PZBIL=0.
      ebilan=0.
      DO ib=1,NTRACK
        Izbil=Izbil+ZVV(ib)
        if(AVV(ib).gt.0) Iabil=Iabil+AVV(ib)
        ebilan=ebilan+ENERJ(ib)
        if((AVV(ib).lt.0).AND.(ZVV(ib).EQ.0)) ebilan=ebilan+134.9764
        if((AVV(ib).lt.0).AND.(ZVV(ib).NE.0)) ebilan=ebilan+139.56995
        TET=3.141592654*TETLAB(ib)/180.
        PHI=3.141592654*PHILAB(ib)/180.
        CSTET=COS(TET)
        SITET=SIN(TET)
        CSPHI=COS(PHI)
        SIPHI=SIN(PHI)
	PXBIL = PXBIL + PLAB(ib)*SITET*CSPHI
	PYBIL = PYBIL + PLAB(ib)*SITET*SIPHI
	PZBIL = PZBIL + PLAB(ib)*CSTET
      ENDDO
	
C conservation of momentum (po=pz, px=py=0.)      
C        WRITE(6,*) 'i,p0,pz,px,py',P0,PZBIL,PXBIL,PYBIL

	P0 = PBEAM
	IF(FF.EQ.0) THEN
C Tests sans fission
		IF ((ABS(P0-PZBIL).GE.0.01)
     s               .OR.(ABS(PXBIL).GE.0.01)
     s               .OR.(ABS(PYBIL).GE.0.01)) THEN
                                impulse1=impulse1+1
	 	ENDIF
		IF((ABS(P0-PZBIL).GE.2.)
     s               .OR.(ABS(PXBIL).GE.1.)
     s               .OR.(ABS(PYBIL).GE.1.)) THEN
                                impulse2=impulse2+1 
		ENDIF		
        ENDIF
	IF(FF.EQ.1) THEN
C Tests avec fission
		IF ((ABS(P0-PZBIL).GE.0.01)
     s               .OR.(ABS(PXBIL).GE.0.01)
     s               .OR.(ABS(PYBIL).GE.0.01)) THEN
                                imp_f1=imp_f1+1
	 	ENDIF
		IF((ABS(P0-PZBIL).GE.2.)
     s               .OR.(ABS(PXBIL).GE.1.)
     s               .OR.(ABS(PYBIL).GE.1.)) THEN
                                imp_f2=imp_f2+1 
		ENDIF		
	ENDIF
		IF ((ABS(P0-PZBIL).GE.5.)
     s               .OR.(ABS(PXBIL).GE.3.0)
     s               .OR.(ABS(PYBIL).GE.3.0)) THEN
                                IF(FF.EQ.0) 
     s			            impulse3=impulse3+1
     				IF(FF.EQ.1)
     s				    imp_f3=imp_f3+1
        IF(choice_evap.EQ.2) GO TO 496  !Pas d'impression pour GEM
C Pour KHS, impression si trop mauvaise conservation de P (PL>5, PT>3):
		WRITE(6,*) '****  Bad momentum conservation  i: ',I	
		WRITE(6,*) 'p0,pz,px,py,FF',P0,PZBIL,PXBIL,PYBIL,FF
	        WRITE(6,*) 'ial,IYV:',ialview,(IYV(ihaz),ihaz=1,19)
		WRITE(6,*) 'ZFF1,AFF1,EFF1,ZF1,AF1,MALPHA1,FFPLEVA1'
     s                     ,',FFPXEVA1,FFPYEVA1,FF1,FTYPE1'
		WRITE(6,*) ZFF1,AFF1,EFF1,ZF1,AF1,MALPHA1,FFPLEVA1
     s                    ,FFPXEVA1,FFPYEVA1,FF1,FTYPE1
		WRITE(6,*) ZFF2,AFF2,EFF2,ZF2,AF2,MALPHA2,FFPLEVA2
     s                    ,FFPXEVA2,FFPYEVA2,FF2,FTYPE2
                WRITE(6,*)'MASSE1,T1,P1,CTET1,PHI1,GAMFIS,ETFIS,R,'
     s   ,'PLAB1,GAM1,ETA1,CSDIR1'
		WRITE(6,*) MASSE1,T1,P1,CTET1,PHI1,GAMFIS,ETFIS,R,
     s   PLAB1,GAM1,ETA1,CSDIR1
		WRITE(6,*) MASSE2,T2,P2,CTET2,PHI2,GAMFIS,ETFIS,R,
     s   PLAB2,GAM2,ETA2,CSDIR2
        IF (FF.EQ.1) THEN
      WRITE(6,*) 'Bilans de fission (système remnant) i:',I
      WRITE(6,*) '	 Remnant,      Bilan (Evapo pré fis)'
      WRITE(6,*) 'E  ',REMMASS,Bil_E 
      WRITE(6,*) 'Px ','          0',Bil_Px 
      WRITE(6,*) 'Py ','          0',Bil_Py 
      WRITE(6,*) 'Pz ','          0',Bil_Pz 
      WRITE(6,*) '       N Fis,    PF1,    PF2,    Bilan'
      WRITE(6,*) 'E  ',EFIS_REM, E1_REM, E2_REM, Bilan_E 
      WRITE(6,*) 'Px ',PFIS_REM(1),PF1_REM(1),PF2_REM(1),Bilan_PX
      WRITE(6,*) 'Py ',PFIS_REM(2),PF1_REM(2),PF2_REM(2),Bilan_PY
      WRITE(6,*) 'Pz ',PFIS_REM(3),PF1_REM(3),PF2_REM(3),Bilan_PZ
      WRITE(6,*) 'NOPART,...',NOPART,LMI_PF1,LMA_PF1,LMI_PF2,LMA_PF2
      WRITE(6,*) '       PF1,    Evapo PF1      (Dans son systeme)'
      WRITE(6,*) 'E  ',MASSE1,Bil1_E 
      WRITE(6,*) 'Px ','          0',Bil1_Px 
      WRITE(6,*) 'Py ','          0',Bil1_Py 
      WRITE(6,*) 'Pz ','          0',Bil1_Pz 
      WRITE(6,*) '       PF1,    PF1-Evapo PF1'
      WRITE(6,*) 'E  ',E1_REM,Bil_E_PF1
      WRITE(6,*) 'Px ',PF1_REM(1),Bil_PX_PF1
      WRITE(6,*) 'Py ',PF1_REM(2),Bil_PY_PF1
      WRITE(6,*) 'Pz ',PF1_REM(3),Bil_PZ_PF1
      WRITE(6,*) '       PF2,    Evapo PF2      (Dans son systeme)'
      WRITE(6,*) 'E  ',MASSE2,Bil2_E 
      WRITE(6,*) 'Px ','          0',Bil2_Px 
      WRITE(6,*) 'Py ','          0',Bil2_Py 
      WRITE(6,*) 'Pz ','          0',Bil2_Pz 
      WRITE(6,*) '       PF2,    PF2-Evapo PF2'
      WRITE(6,*) 'E  ',E2_REM,Bil_E_PF2
      WRITE(6,*) 'Px ',PF2_REM(1),Bil_PX_PF2
      WRITE(6,*) 'Py ',PF2_REM(2),Bil_PY_PF2
      WRITE(6,*) 'Pz ',PF2_REM(3),Bil_PZ_PF2
      WRITE(6,*) 'Energie des Résidus de fis. ',
     s  '(apres évapo) :',EPF1_OUT,EPF2_out
        ENDIF

		
		CALL PRINTCOM
C		STOP
		ENDIF


496     CONTINUE        
C Testing permanently A and Z conservation for the Ntuple:
	NATOT = NINT(f(1)+AP)
	NZTOT = NINT(f(2)+ZP)	
	IF((NATOT.NE.Iabil).OR.(NZTOT.NE.Izbil)) THEN
	WRITE(6,*) '**** Bad mass or charge conservation  i: ',I
     	write(6,*)'BILAN:A,Z,E,ERECREM',Iabil,Izbil,ebilan,erecrem
        WRITE(6,*) 'Run,nopart=',i,nopart
	WRITE(6,*) 'ial,IYV:',ialview,(IYV(ihaz),ihaz=1,19)
		WRITE(6,*) 'i,p0,pz,px,py',I,P0,PZBIL,PXBIL,PYBIL
		WRITE(6,*) 'ZFF1,AFF1,EFF1,ZF1,AF1,MALPHA1,FFPLEVA1
     s                    ,FFPXEVA1,FFPYEVA1,FF1,FTYPE1'
		WRITE(6,*) ZFF1,AFF1,EFF1,ZF1,AF1,MALPHA1,FFPLEVA1
     s                    ,FFPXEVA1,FFPYEVA1,FF1,FTYPE1
		WRITE(6,*) ZFF2,AFF2,EFF2,ZF2,AF2,MALPHA2,FFPLEVA2
     s                    ,FFPXEVA2,FFPYEVA2,FF2,FTYPE2
                WRITE(6,*)'MASSE1,T1,P1,CTET1,PHI1,GAMFIS,ETFIS,R,'
     s   ,'PLAB1,GAM1,ETA1,CSDIR1'
		WRITE(6,*) MASSE1,T1,P1,CTET1,PHI1,GAMFIS,ETFIS,R,
     s   PLAB1,GAM1,ETA1,CSDIR1
		WRITE(6,*) MASSE2,T2,P2,CTET2,PHI2,GAMFIS,ETFIS,R,
     s   PLAB2,GAM2,ETA2,CSDIR2
		
		CALL PRINTCOM
C        STOP
        ENDIF
        
        IF(NTRACK.GT.0) THEN
		nbevhbk=nbevhbk+1
C		IF((inum.GE.92800).AND.(inum.LE.113300))
C     s  		WRITE(6,*) 'iboucle, nbevhbk:',inum,nbevhbk

C Evenement entre dans le NTUPLE:     
		CALL HFNTB(101,'VAR_NTP')
		IF(Je_veux.EQ.1) CALL HFNTB(102,'VAR_AVAT')
	ENDIF

C Verif de reproduction d'un evenement par reinitialisation du common hazard
c    (impression de l'evenement numero ievtest).
	IF(I.EQ.ievtest) THEN
c	if(inum.ge.1)then
c	write(6,*)'Numero',inum
	     WRITE(6,*) 'ial,IYV:',ialview,(IYV(ihaz),ihaz=1,19)
	     WRITE(6,*) MASSINI,MZINI,EXINI,MULNCASC,MULNEVAP
     	     WRITE(6,*) MULNTOT,BIMPACT,NTRACK
     		DO itest=1,NTRACK
                  WRITE(6,*) itest,ITYPCASC(itest),AVV(itest),
     s             ZVV(itest),ENERJ(itest),PLAB(itest),
     s             TETLAB(itest),PHILAB(itest)
                ENDDO
	ENDIF

C Imressions pour vérif sur un test a choisir:
C        IF (EXINI.GE.58.) THEN
C	WRITE(6,*) 'inum,IZREM,IAREM,MCOREM',inum,IZREM,IAREM,MCOREM
C        WRITE(6,*) 'E*,T_REM,A,B,G',ESREM,ERECREM,ALREM,BEREM,GAREM
C		CALL PRINTCOM
C	ENDIF
	
C calcul des multiplicites de neutrons:
	DO ib=1,NTRACK
	  IF((AVV(ib).EQ.1).AND.(ZVV(ib).EQ.0)) THEN
	    IF(ENERJ(ib).GT.20.) THEN
	    	muln_max=muln_max+1.
		ener_max=ener_max+ENERJ(ib)
	    ELSEIF(ENERJ(ib).GT.2) THEN
	    	muln_20=muln_20+1.
		ener_20=ener_20+ENERJ(ib)
	    ELSEIF(ENERJ(ib).GE.0.) THEN
	    	muln_2=muln_2+1.
		ener_2=ener_2+ENERJ(ib)
	    ELSE	    
		    WRITE(6,*) 'IL Y A UN DEFAUT!!!'
		    WRITE(6,*) 'tir numero:',INUM
	            WRITE(6,*) 'ial,IYV:',ialview,(IYV(ihaz),ihaz=1,19)
		    WRITE(6,*) MASSINI,MZINI,EXINI,MULNCASC,MULNEVAP
		    WRITE(6,*) MULNTOT,BIMPACT,NTRACK
	      DO ikl=1,NTRACK
	      WRITE(6,*)ikl,ITYPCASC(ikl),AVV(ikl),ZVV(ikl),ENERJ(ikl)
	      ENDDO
	    ENDIF
	  ENDIF
	END DO
	
C ************************************************************************
C                     Fin de la boucle sur les tirages
C ************************************************************************
99	CONTINUE
100     CONTINUE

C Print de verif
C	nbons=icoup-ntrans-nabs-nretir
C	PXMOY=S1x/nbons
C	PYMOY=S1y/nbons
C	PZMOY=S1z/nbons
C	SIGX=DSQRT(S2x/nbons-PXMOY**2)
C	SIGY=DSQRT(S2y/nbons-PYMOY**2)
C	SIGZ=DSQRT(DABS(S2z/nbons-PZMOY**2))
C	WRITE(6,*) 'pz,sigpz,px,sigpx,py,sigpy',
C    s  PZMOY,SIGZ,PXMOY,SIGX,PYMOY,SIGY

        f_cross_sect=(icoup-ntrans)
	f_cross_sect=f_cross_sect/(icoup+ntrans_coul)
	

        WRITE(6,*) ' '
	WRITE(6,*) 'END OF SHOTS'
        WRITE(6,*) '        ',nfis,' fissions'
	WRITE(6,*) '        ',ntrans,' transp nuclear   ',
     s                        ntrans_coul,' transp coulomb'
	WRITE(6,*) '        ',nabs,' absorptions'
	WRITE(6,*) '        ',nretir,' retirages cascade (nopart=-100)'
	WRITE(6,*) 'Geometrical cross section (mb) ',31.4159*BMAX**2
	WRITE(6,*) 'React. cross section (mb): ',
     s                      f_cross_sect*31.4159*BMAX**2
        WRITE(6,*) ' '
	WRITE(6,*) ' Bad momentum conservation without fission:',
     s             impulse1,impulse2,impulse3
	WRITE(6,*) ' Bad momentum conservation with fission:',
     s             imp_f1,imp_f2,imp_f3
        WRITE(6,*)' (number of evts above respect. 0.01, 2 and 5 MeV/c)'
        WRITE(6,*) ' '
	
        WRITE(6,*)' last random',ial
        fmuln=muln_2/(icoup-ntrans)
	ener_2=ener_2/(icoup-ntrans)
        err=SQRT(muln_2)/(icoup-ntrans)
        WRITE(6,*) ' '
        WRITE(6,1516) fmuln,err,ener_2
1516   FORMAT( F8.3,'+/-',F8.3,' neutrons (0-2 MeV) per interaction',
     s  /,10x,F8.3,' MeV = Mean energy carried by ALL of them')
        fmuln=muln_20/(icoup-ntrans)
	ener_20=ener_20/(icoup-ntrans)
        err=SQRT(muln_20)/(icoup-ntrans)
        WRITE(6,1517) fmuln,err,ener_20
1517   FORMAT( F8.3,'+/-',F8.3,' neutrons (2-20 MeV) per interaction',
     s  /,10x,F8.3,' MeV = Mean energy carried by ALL of them')
        fmuln=muln_max/(icoup-ntrans)
	ener_max=ener_max/(icoup-ntrans)
        err=SQRT(muln_max)/(icoup-ntrans)
        WRITE(6,1518) fmuln,err,ener_max
1518   FORMAT( F8.3,'+/-',F8.3,' neutrons (above 20 MeV)',
     s  ' per interaction',
     s  /,10x,F8.3,' MeV = Mean energy carried by ALL of them')
        WRITE(6,*) ' '
      
        call hbclose(Je_veux)

C ******* Time stop ****************

      CALL getime(1)
        
      STOP
      END
C*****************************************************************************          
      subroutine hbinit(Je_veux)                                                         
      implicit real*8 (a-h,o-z)
                                                                                     
      character*80 string
      parameter(lrecl=4096)
      parameter (max=250)
      COMMON/PAWC/W(6000000)
	real*4 EXINI,ENERJ,BIMPACT,PLAB,TETLAB,PHILAB,ESTFIS
	integer AVV,ZVV,JREMN,KFIS,IZFIS,IAFIS      
      COMMON/VAR_NTP/MASSINI,MZINI,EXINI,MULNCASC,MULNEVAP,
     +MULNTOT,BIMPACT,JREMN,KFIS,ESTFIS,IZFIS,IAFIS,NTRACK,
     +ITYPCASC(max),AVV(max),ZVV(max),ENERJ(max),PLAB(max),
     +TETLAB(max),PHILAB(max)
      
      REAL*4 Bavat,TIME,ENERGY,EPSd,EPS2,EPS4,EPS6,EPSf
      REAL*4 R1_in,R1_first_avat
      INTEGER Bloc_Paul,Bloc_CDPP,GO_OUT,avm,DEL1,DEL2
      PARAMETER (avm=1000)
      COMMON/VAR_AVAT/Kveux,Bavat,NOPART,NCOL,
     s R1_in(3),R1_first_avat(3),
     s EPSd(250),EPS2(250),EPS4(250),EPS6(250),EPSf(250),      
     s NB_AVAT,
     s TIME(avm),L1(avm),L2(avm),JPARTL1(avm),JPARTL2(avm),
     s DEL1(avm),DEL2(avm),ENERGY(avm),Bloc_Paul(avm),
     s Bloc_CDPP(avm),GO_OUT(avm)
       
      call hlimit(6000000)
      READ(5,*)string
c      string='/home/crash2/il2_crash/volant/kheinz/r2bis45.hbook'
      write(6,*)string
      id=11
      call hropen(id,'NTP',string,'N',lrecl,ierr)
      call hbnt(101,'Bidon test',' ')

      call hbname(101,'VAR_NTP',
     +MASSINI,'MASSINI,MZINI,EXINI,MULNCASC,MULNEVAP,
     +MULNTOT,BIMPACT,JREMN,KFIS,ESTFIS,IZFIS,IAFIS,NTRACK[0,250],
     +'//'ITYP(NTRACK),AVV(NTRACK):I*4,ZVV(NTRACK):I*4,
     +ENERJ(NTRACK),PLAB(NTRACK),TETLAB(NTRACK),PHILAB(NTRACK)')
     
C      CALL HBOOK1(1,'J_rem_INCL',100,0.,100.,0.)
C      CALL HBOOK1(2,'J_rem_KHS',100,0.,100.,0.)

C Pour verif de la densite du noyau cible:
C	CALL HBOOK1(3,'r_distrib',900,0.,15.,0.)
C	CALL HBOOK1(4,'p_distrib',900,0.,300.,0.)

C Second NTUPLE optionnel pour etude des avatars en fonction du temps 
C de la cascade.

      READ(5,*) Je_veux,string
           Kveux=Je_veux	!Pour transmission a INCL
      IF (Je_veux.NE.1) RETURN
      WRITE(6,*) ' Ntuple for avatars study:' 
      write(6,*) string
      
      id=12
      call hropen(id,'AVATAR',string,'N',lrecl,ierr)
      call hbnt(102,'Avatars test',' ')

c$$$      call hbname(102,'VAR_AVAT',
c$$$     sBavat,'BIMPACT,NOPART,NCOL,R1IN(3),R1FIRSTAVAT(3),
c$$$     sEPSD(250),EPS2(250),EPS4(250),EPS6(250),EPSF(250),
c$$$     sNB_AVAT[0,1000],
c$$$     sTIME(NB_AVAT),L1(NB_AVAT),L2(NB_AVAT),
c$$$     sJPARL1(NB_AVAT),JPARTL2(NB_AVAT),
c$$$     sDEL1(NB_AVAT):I*4,DEL2(NB_AVAT):I*4,
c$$$     sENERGY(NB_AVAT),Bloc_Paul(NB_AVAT):I*4,Bloc_CDPP(NB_AVAT):I*4,
c$$$     sGO_OUT(NB_AVAT):I*4')
      

      return                                                                    
      end                                                                       
C***************************************************************************
      subroutine hbclose(Je_veux)                                              
      implicit real*8 (a-h,o-z)                                                 
                                                          
      call hrout(0,icycle,' ')
      call hrend('NTP')
      IF(Je_veux.EQ.1) CALL HREND('AVATAR')                                               
      return                                                                    
      end                                
C ------------------------------------------------------------------
      SUBROUTINE PRINTCOM
c Impression des common
	parameter (max=250)                                                                       
	real*4 EXINI,ENERJ,BIMPACT,PLAB,TETLAB,PHILAB,ESTFIS
	integer AVV,ZVV,JREMN,KFIS,IZFIS,IAFIS
        common/VAR_NTP/MASSINI,MZINI,EXINI,MULNCASC,MULNEVAP,
     +MULNTOT,BIMPACT,JREMN,KFIS,ESTFIS,IZFIS,IAFIS,NTRACK,
     +ITYPCASC(max),AVV(max),ZVV(max),ENERJ(max),PLAB(max),
     +TETLAB(max),PHILAB(max)
	
	WRITE(6,*) 'MASSINI,MZINI,EXINI,MULNCASC,MULNEVAP,
     +MULNTOT,BIMPACT,NTRACK',MASSINI,MZINI,EXINI,MULNCASC,MULNEVAP,
     +MULNTOT,BIMPACT,NTRACK
        WRITE(6,*) 'ITYPCASC,AVV,ZVV,ENERJ,PLAB,TETLAB,PHILAB'
     	DO i=1,NTRACK
     	WRITE(6,*) ITYPCASC(i),AVV(i),ZVV(i),ENERJ(i),PLAB(i),
     +TETLAB(i),PHILAB(i)
        ENDDO
      RETURN
      END
C ------------------------------------------------------------------
      SUBROUTINE TRANSLAB(GAMREM,ETREM,CSREM,NOPART,NDEC)
c Ce subroutine transforme dans un repere 1 les impulsions pcv des 
c particules acv, zcv et de cosinus directeurs xcv, ycv, zcv calculees 
c dans un repere 2.    
c La transformation de lorentz est definie par GAMREM (gamma) et
c ETREM (eta). La direction  du repere 2 dans 1 est donnees par les 
c cosinus directeurs ALREM,BEREM,GAREM (axe oz du repere 2).
c L'axe oy(2) est fixe par le produit vectoriel oz(1)*oz(2).
c Le calcul est fait pour les particules de NDEC a iv du common volant.
C Resultats dans le NTUPLE (common VAR_NTP) decale de NOPART (cascade).
    
      REAL*8  GAMREM,ETREM,ER,PLABI(3),PLABF(3),R(3,3)
      real*8  MASSE,PTRAV2,CSREM(3),UMA,MELEC,EL
      real*4 acv,zpcv,pcv,xcv,ycv,zcv
      common/volant/acv(200),zpcv(200),pcv(200),xcv(200),
     s              ycv(200),zcv(200),iv
      
	parameter (max=250)                                                                       
	real*4 EXINI,ENERJ,BIMPACT,PLAB,TETLAB,PHILAB,ESTFIS
	integer AVV,ZVV,JREMN,KFIS,IZFIS,IAFIS
        common/VAR_NTP/MASSINI,MZINI,EXINI,MULNCASC,MULNEVAP,
     +MULNTOT,BIMPACT,JREMN,KFIS,ESTFIS,IZFIS,IAFIS,NTRACK,
     +ITYPCASC(max),AVV(max),ZVV(max),ENERJ(max),PLAB(max),
     +TETLAB(max),PHILAB(max)
      
      DATA UMA,MELEC/931.4942,0.511/

C Matrice de rotation dans le labo:
        SITET = SQRT(CSREM(1)**2+CSREM(2)**2)
        IF(SITET.GT.1.E-6)THEN
        CSTET = CSREM(3)
        SIPHI = CSREM(2)/SITET
        CSPHI = CSREM(1)/SITET	

	R(1,1) = CSTET*CSPHI
	R(1,2) = -SIPHI
	R(1,3) = SITET*CSPHI
	R(2,1) = CSTET*SIPHI
	R(2,2) = CSPHI
	R(2,3) = SITET*SIPHI
	R(3,1) = -SITET
	R(3,2) = 0.
	R(3,3) = CSTET
	ELSE

	R(1,1) = 1.
	R(1,2) = 0.
	R(1,3) = 0.
	R(2,1) = 0.
	R(2,2) = 1.
	R(2,3) = 0.
	R(3,1) = 0.
	R(3,2) = 0.
	R(3,3) = 1.
	ENDIF
	
     
      DO i=NDEC,iv
        intp = i + NOPART
        AVV(intp) = NINT(acv(i))
        ZVV(intp) = NINT(zpcv(i))
        ITYPCASC(intp) = 0    
C Transformation de Lorentz Remnan --> Labo:
	IF (AVV(intp).EQ.-1) THEN
		MASSE=138.00	!Cugnon
C		IF (AVV(intp).EQ.1)  MASSE=938.2796	!Cugnon
C		IF (AVV(intp).EQ.4)  MASSE=3727.42	!OK
	ELSE
           CALL MGLMS(DBLE(acv(i)),DBLE(zpcv(i)),0,EL)
	   MASSE = zpcv(i)*938.27 + (acv(i)- zpcv(i))*939.56 + EL
	END IF
	
	ER = DSQRT(pcv(i)**2 + MASSE**2)
	
	PLABI(1) = pcv(i)*xcv(i)
	PLABI(2) = pcv(i)*ycv(i)
	PLABI(3) = ER*ETREM + GAMREM*pcv(i)*zcv(i)
	 
	
	PTRAV2 = PLABI(1)**2 +PLABI(2)**2 +PLABI(3)**2
	PLAB(intp) = DSQRT(PTRAV2) 
	ENERJ(intp) = DSQRT(PTRAV2 + MASSE**2) - MASSE
	
C Rotation dans le labo:
	DO j=1,3
	    PLABF(j) = 0.
		DO k=1,3
		    PLABF(j) = PLABF(j) + R(j,k)*PLABI(k)
		END DO
	END DO
C impulsions dans le nouveau systeme copiees dans /volant/
	pcv(i) = PLAB(intp)
        ptrav2=sqrt(plabf(1)**2+plabf(2)**2+plabf(3)**2)
      	IF(ptrav2.GE.1.e-6) THEN
	xcv(i) = PLABF(1)/ptrav2
	ycv(i) = PLABF(2)/ptrav2
	zcv(i) = PLABF(3)/ptrav2
	ELSE
	xcv(i)=1.
	ycv(i)=0.
	zcv(i)=0.
	ENDIF
c impulsions dans le nouveau systeme copiees dans /VAR_NTP/	
	IF(PLAB(intp).GE.1.e-6) THEN

	bidon=PLABF(3)/PLAB(intp)
	if(bidon.gt.1)bidon=1.
	if(bidon.lt.-1.)bidon=-1.
	TETLAB(intp) = ACOS(bidon)
          SITET = SIN(TETLAB(intp))
          PHILAB(intp) = ATAN2(PLABF(2),PLABF(1))        
          TETLAB(intp) = TETLAB(intp)*57.2957795
          PHILAB(intp) = PHILAB(intp)*57.2957795
	ELSE
	  TETLAB(intp) = 90.
	  PHILAB(intp) = 0.
	ENDIF
      END DO
      
      RETURN
      END
C-------------------------------------------------------------------------
      SUBROUTINE TRANSLABPF(MASSE1,T1,P1,CTET1,PHI1,GAMREM,ETREM,R,
     s   PLAB1,GAM1,ETA1,CSDIR)
C Calcul de l'impulsion du PF (PLAB1, cos directeurs CSDIR(3)) dans le
C systeme remnant et des coefs de Lorentz GAM1,ETA1 de passage  
c du systeme PF --> systeme remnant.
c 
C Input: MASSE1, T1 (energie cinetique), CTET1,PHI1 (cosTHETA et PHI)
C                    (le PF dans le systeme du Noyau de Fission (NF)).
C	 GAMREM,ETREM les coefs de Lorentz systeme NF --> syst remnant, 
C        R(3,3) la matrice de rotation systeme NF--> systeme remnant.
C
C      
     	REAL*8 MASSE1,T1,P1,CTET1,PHI1,GAMREM,ETREM,R(3,3),
     s   PLAB1,GAM1,ETA1,CSDIR(3),ER,SITET,PLABI(3),PLABF(3)
     
	ER = T1 + MASSE1
	
	SITET = DSQRT(1.-CTET1**2)
	
C ----Transformation de Lorentz Noyau fissionnant --> Remnant:	
	PLABI(1) = P1*SITET*COS(PHI1)
	PLABI(2) = P1*SITET*SIN(PHI1)
	PLABI(3) = ER*ETREM + GAMREM*P1*CTET1
	
C ----Rotation du syst Noyaut Fissionant vers syst remnant:
	DO j=1,3
	    PLABF(j) = 0.
		DO k=1,3
		    PLABF(j) = PLABF(j) + R(j,k)*PLABI(k)
		END DO
	END DO
C ----Cosinus directeurs et coefs de la transf de Lorentz dans le
c     nouveau systeme:	
        PLAB1 = PLABF(1)**2+PLABF(2)**2+PLABF(3)**2
        GAM1 = DSQRT(MASSE1**2 + PLAB1)/MASSE1
        PLAB1 = DSQRT(PLAB1)
        ETA1 = PLAB1/MASSE1
        
	IF(PLAB1.LE.1.E-6) THEN
	   CSDIR(1)=0.
	   CSDIR(2)=0.
	   CSDIR(3)=1.
	   ELSE   
        DO i=1,3
           CSDIR(i) = PLABF(i)/PLAB1
        END DO
	ENDIF
        
        RETURN
        END
C-------------------------------------------------------------------------
      subroutine getime(ic)
      dimension jd(3),jt(3),ti(2)                                       
      common/cput/cpa

      if(ic.eq.0) then
       cpa=etime(ti)
       call idate(jd)                                                    
       call itime(jt)                                                    
       write(*,10)jd(2),jd(1),jd(3),jt
   10 format(20hCalculation starts :,2x,i2.2,1h/,i2.2,1h/,i4.4,1x,i2.2
     &      ,1h:,i2.2,1h:,i2.2)       
      else
       cpe=etime(ti)
       call idate(jd)                                                    
       call itime(jt)                                                    
       write(*,20)jd(2),jd(1),jd(3),jt
   20 format(20hCalculation ends   :,2x,i2.2,1h/,i2.2,1h/,i4.4,1x,i2.2
     &      ,1h:,i2.2,1h:,i2.2)       
       write(*,30)cpe-cpa
   30 format(20hExecution time is  :,f10.2,4h sec)       
      endif
      return
      end
C------------------------------------------------------------------------
      SUBROUTINE LOR_AB(GAM,ETA,Ein,Pin,Eout,Pout)
C  Transformation de lorentz brute pour vérifs.
C	P(3) = P_longitudinal (transformé)
C	P(1) et P(2) = P_transvers (non transformés)
      DIMENSION Pin(3),Pout(3)
      REAL*8 GAM,ETA,Ein

      Pout(1) = Pin(1)
      Pout(2) = Pin(2) 
      Eout = GAM*Ein + ETA*Pin(3)
      Pout(3) = + ETA*Ein + GAM*Pin(3)
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE ROT_AB(R,Pin,Pout)
C  Rotation d'un vecteur
      DIMENSION Pin(3),Pout(3)
      REAL*8 R(3,3)
      
      DO i=1,3
      Pout(i) = 0.
      	DO j=1,3
		Pout(i) = Pout(i) + R(i,j)*Pin(j)
	ENDDO
      ENDDO
      
      RETURN
      END         
 

