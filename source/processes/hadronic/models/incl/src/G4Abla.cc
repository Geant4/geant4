//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4Abla.cc,v 1.1 2007-05-23 10:25:37 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include "G4Abla.hh"
#include <time.h>

#include "G4InclAblaHardcodedData.hh"

G4Abla::G4Abla()
{

}

G4Abla::G4Abla(G4Hazard *hazard, G4Volant *volant)
{
  volant = volant;
  hazard = hazard;
  
  pace = (G4Pace*) malloc(sizeof(G4Pace));
  ald = (G4Ald*) malloc(sizeof(G4Ald));
  ablamain = (G4Ablamain*) malloc(sizeof(G4Ablamain));
  emdpar = (G4Emdpar*) malloc(sizeof(G4Emdpar));
  eenuc = (G4Eenuc*) malloc(sizeof(G4Eenuc));
  ec2sub = (G4Ec2sub*) malloc(sizeof(G4Ec2sub));
  ecld = (G4Ecld*) malloc(sizeof(G4Ecld)); 
  fb = (G4Fb*) malloc(sizeof(G4Fb));
  fiss = (G4Fiss*) malloc(sizeof(G4Fiss));
  opt = (G4Opt*) malloc(sizeof(G4Opt));

}

G4Abla::~G4Abla()
{

}

// Evaporation code
void G4Abla::initEvapora()
{
//      6	C *******************************************************************                                                                      
//      7	C                                                                       
//      8	C      SUBROUTINE ABLAINIT(STATUS,TSTAT,NAME,FPATH)
//      9	C********************************************************************                      
//     10	
//     11	      SUBROUTINE INIT_EVAPORA(RACINE)
//     12	                      
//     13	C********************************************************************                                                                       
//     14	C     ON INPUT:  INPUT PARAMETERS FROM FILE                             
//     15	C---------------------------------------------------------------------  
//     16	C     ON OUTPUT:                                                        
//     17	C     STATUS - FLAG FOR END OF INPUT FILE                               
//     18	C     TSTAT  - FLAG FOR NTUPLE-OUTPUT                                   
//     19	C     NAME   - NAME FOR ISOTOPIC PRODUCTION CROSS SECTION FILES         
//     20	C     FPATH  - PATH FOR  "          "        "        "    "            
//     21	C---------------------------------------------------------------------
//     22	C
//     23	C     Modification 5-january-2000 by KHS and BJ
//     24	C
//     25	C     New treatment of dissipation. 
//     26	C     See report of Beatriz Jurado, Jan. 2000
//     27	C
//     28	C---------------------------------------------------------------------
//     29	C   
//     30	C     MODIFICATION 6-aug-1999 by JB and MVR
//     31	C
//     32	C Some problems arised from an uncorrect evaluation of the fission barrier
//     33	C { 1) shell correction ( ECGNZ(J,K) ) was not subctracted (20-jul-99)
//     34	C   2) fiss. barrier EF was calc. before ang. mom. correct. (6-aug-99) }
//     35	C
//     36	C---------------------------------------------------------------------  
//     37	C     PROJECTILE AND TARGET PARAMETERS + CROSS SECTIONS                 
//     38	C     COMMON /ABLAMAIN/ AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,       
//     39	C                       R_0,R_P,R_T, IMAX,IRNDM,PI,                     
//     40	C                       BFPRO,SNPRO,SPPRO,SHELL                         
//     41	C                                                                       
//     42	C     AP,ZP,AT,ZT   - PROJECTILE AND TARGET MASSES                      
//     43	C     EAP,BETA      - BEAM ENERGY PER NUCLEON, V/C                      
//     44	C     BMAXNUC       - MAX. IMPACT PARAMETER FOR NUCL. REAC.             
//     45	C     CRTOT,CRNUC   - TOTAL AND NUCLEAR REACTION CROSS SECTION          
//     46	C     R_0,R_P,R_T,  - RADIUS PARAMETER, PROJECTILE+ TARGET RADII        
//     47	C     IMAX,IRNDM,PI - MAXIMUM NUMBER OF EVENTS, DUMMY, 3.141...         
//     48	C     BFPRO         - FISSION BARRIER OF THE PROJECTILE                 
//     49	C     SNPRO         - NEUTRON SEPARATION ENERGY OF THE PROJECTILE       
//     50	C     SPPRO         - PROTON    "           "   "    "   "              
//     51	C     SHELL         - GROUND STATE SHELL CORRECTION                     
//     52	C---------------------------------------------------------------------  
//     53	C                                                                       
//     54	C     ENERGIES WIDTHS AND CROSS SECTIONS FOR EM EXCITATION              
//     55	C     COMMON /EMDPAR/ EGDR,EGQR,FWHMGDR,FWHMGQR,CREMDE1,CREMDE2,        
//     56	C                     AE1,BE1,CE1,AE2,BE2,CE2,SR1,SR2,XR                
//     57	C                                                                       
//     58	C     EGDR,EGQR       - MEAN ENERGY OF GDR AND GQR                      
//     59	C     FWHMGDR,FWHMGQR - FWHM OF GDR, GQR                                
//     60	C     CREMDE1,CREMDE2 - EM CROSS SECTION FOR E1 AND E2                  
//     61	C     AE1,BE1,CE1     - ARRAYS TO CALCULATE                             
//     62	C     AE2,BE2,CE2     - THE EXCITATION ENERGY AFTER E.M. EXC.           
//     63	C     SR1,SR2,XR      - WITH MONTE CARLO                                
//     64	C---------------------------------------------------------------------  
//     65	C                                                                       
//     66	C     DEFORMATIONS AND G.S. SHELL EFFECTS                               
//     67	C     COMMON /ECLD/   ECGNZ,ECFNZ,VGSLD,ALPHA                           
//     68	C                                                                       
//     69	C     ECGNZ - GROUND STATE SHELL CORR. FRLDM FOR A SPHERICAL G.S.       
//     70	C     ECFNZ - SHELL CORRECTION FOR THE SADDLE POINT (NOW: == 0)         
//     71	C     VGSLD - DIFFERENCE BETWEEN DEFORMED G.S. AND LDM VALUE            
//     72	C     ALPHA - ALPHA GROUND STATE DEFORMATION (THIS IS NOT BETA2!)       
//     73	C             BETA2 = SQRT(5/(4PI)) * ALPHA                             
//     74	C---------------------------------------------------------------------  
//     75	C                                                                       
//     76	C     ARRAYS FOR EXCITATION ENERGY BY STATISTICAL HOLE ENERY MODEL      
//     77	C     COMMON /EENUC/  SHE, XHE                                          
//     78	C                                                                       
//     79	C     SHE, XHE - ARRAYS TO CALCULATE THE EXC. ENERGY AFTER              
//     80	C                ABRASION BY THE STATISTICAL HOLE ENERGY MODEL          
//     81	C---------------------------------------------------------------------  
//     82	C                                                                       
//     83	C     G.S. SHELL EFFECT                                                 
//     84	C     COMMON /EC2SUB/ ECNZ                                              
//     85	C                                                                       
//     86	C     ECNZ G.S. SHELL EFFECT FOR THE MASSES (IDENTICAL TO ECGNZ)        
//     87	C---------------------------------------------------------------------  
//     88	C                                                                       
//     89	C     OPTIONS AND PARAMETERS FOR FISSION CHANNEL                        
//     90	C     COMMON /FISS/    AKAP,BET,HOMEGA,KOEFF,IFIS,                       
//     91	C                            OPTSHP,OPTXFIS,OPTLES,OPTCOL               
//     92	C                                                                       
//     93	C     AKAP   - HBAR**2/(2* MN * R_0**2) = 10 MEV                        
//     94	C     BET    - REDUCED NUCLEAR FRICTION COEFFICIENT IN (10**21 S**-1)   
//     95	C     HOMEGA - CURVATURE OF THE FISSION BARRIER = 1 MEV                 
//     96	C     KOEFF  - COEFFICIENT FOR THE LD FISSION BARRIER == 1.0            
//     97	C     IFIS   - 0/1 FISSION CHANNEL OFF/ON                               
//     98	C     OPTSHP - INTEGER SWITCH FOR SHELL CORRECTION IN MASSES/ENERGY     
//     99	C            = 0 NO MICROSCOPIC CORRECTIONS IN MASSES AND ENERGY        
//    100	C            = 1 SHELL ,  NO PAIRING                                    
//    101	C            = 2 PAIRING, NO SHELL                                      
//    102	C            = 3 SHELL AND PAIRING                                      
//    103	C     OPTCOL - 0/1 COLLECTIVE ENHANCEMENT SWITCHED ON/OFF               
//    104	C     OPTXFIS- 0,1,2 FOR MYERS & SWIATECKI, DAHLINGER, ANDREYEV         
//    105	C              FISSILITY PARAMETER.                                     
//    106	C     OPTLES - CONSTANT TEMPERATURE LEVEL DENSITY FOR A,Z > TH-224      
//    107	C     OPTCOL - 0/1 COLLECTIVE ENHANCEMENT OFF/ON                        
//    108	C---------------------------------------------------------------------  
//    109	C                                                                       
//    110	C     OPTIONS                                                           
//    111	C     COMMON /OPT/    OPTEMD,OPTCHA,EEFAC                               
//    112	C                                                                       
//    113	C     OPTEMD - 0/1  NO EMD / INCL. EMD                                  
//    114	C     OPTCHA - 0/1  0 GDR / 1 HYPERGEOMETRICAL PREFRAGMENT-CHARGE-DIST. 
//    115	C              ***  RECOMMENDED IS OPTCHA = 1 ***                       
//    116	C     EEFAC  - EXCITATION ENERGY FACTOR, 2.0 RECOMMENDED                
//    117	C---------------------------------------------------------------------  
//    118	C                                                                       
//    119	C     FISSION BARRIERS                                                  
//    120	C     COMMON /FB/     EFA                                               
//    121	C     EFA    - ARRAY OF FISSION BARRIERS                                
//    122	C---------------------------------------------------------------------  
//    123	C                                                                       
//    124	C     LEVEL DENSITY PARAMETERS                                          
//    125	C     COMMON /ALD/    AV,AS,AK,OPTAFAN                                  
//    126	C     AV,AS,AK - VOLUME,SURFACE,CURVATURE DEPENDENCE OF THE             
//    127	C                LEVEL DENSITY PARAMETER                                
//    128	C     OPTAFAN - 0/1  AF/AN >=1 OR AF/AN ==1                             
//    129	C               RECOMMENDED IS OPTAFAN = 0                              
//    130	C---------------------------------------------------------------------  
//    131	C   ____________________________________________________________________
//    132	C  /                                                                    
//    133	C  /  INITIALIZES PARAMETERS IN COMMON /ABRAMAIN/, /EMDPAR/, /ECLD/ ... 
//    134	C  /  PROJECTILE PARAMETERS, EMD PARAMETERS, SHELL CORRECTION TABLES.   
//    135	C  /  CALCULATES MAXIMUM IMPACT PARAMETER FOR NUCLEAR COLLISIONS AND    
//    136	C  /  TOTAL GEOMETRICAL CROSS SECTION + EMD CROSS SECTIONS              
//    137	C   ____________________________________________________________________
//    138	C                                                                       
//    139	C                                                                       
//    201	C                                                                       
//    202	C---------- SET INPUT VALUES                                            
//    203	C                                                                       
//    204	C *** INPUT FROM UNIT 10 IN THE FOLLOWING SEQUENCE !                    
//    205	C     AP1 =    INTEGER  !                                               
//    206	C     ZP1 =    INTEGER  !                                               
//    207	C     AT1 =    INTEGER  !                                               
//    208	C     ZT1 =    INTEGER  !                                               
//    209	C     EAP =    REAL     !                                               
//    210	C     IMAX =   INTEGER  !                                               
//    211	C     IFIS =   INTEGER SWITCH FOR FISSION                               
//    212	C     OPTSHP = INTEGER SWITCH FOR SHELL CORRECTION IN MASSES/ENERGY     
//    213	C            =0 NO MICROSCOPIC CORRECTIONS IN MASSES AND ENERGY         
//    214	C            =1 SHELL , NO PAIRING CORRECTION                           
//    215	C            =2 PAIRING, NO SHELL CORRECTION                            
//    216	C            =3 SHELL AND PAIRING CORRECTION IN MASSES AND ENERGY       
//    217	C     OPTEMD =0,1  0 NO EMD, 1 INCL. EMD                                
//    218	C               ELECTROMAGNETIC DISSOZIATION IS CALCULATED AS WELL.     
//    219	C     OPTCHA =0,1  0 GDR- , 1 HYPERGEOMETRICAL PREFRAGMENT-CHARGE-DIST. 
//    220	C               RECOMMENDED IS OPTCHA=1                                 
//    221	C     OPTCOL =0,1 COLLECTIVE ENHANCEMENT SWITCHED ON 1 OR OFF 0 IN DENSN
//    222	C     OPTAFAN=0,1 SWITCH FOR AF/AN = 1 IN DENSNIV 0 AF/AN>1 1 AF/AN=1   
//    223	C     AKAP =  REAL    ALWAYS EQUALS 10                                  
//    224	C     BET  =  REAL    REDUCED FRICTION COEFFICIENT / 10**(+21) S**(-1)  
//    225	C     HOMEGA = REAL   CURVATURE / MEV RECOMMENDED = 1. MEV              
//    226	C     KOEFF  = REAL   COEFFICIENT FOR FISSION BARRIER                   
//    227	C     OPTXFIS= INTEGER 0,1,2 FOR MYERS & SWIATECKI, DAHLINGER, ANDREYEV 
//    228	C              FISSILITY PARAMETER.                                     
//    229	C     EEFAC  = REAL EMPIRICAL FACTOR FOR THE EXCITATION ENERGY          
//    230	C                   RECOMMENDED 2.D0, STATISTICAL ABRASION MODELL 1.D0  
//    231	C     AV     = REAL KOEFFICIENTS FOR CALCULATION OF A(TILDE)            
//    232	C     AS     = REAL LEVEL DENSITY PARAMETER                             
//    233	C     AK     = REAL                                                     
//    234	C                                                                       
//    235	C This following inputs will be initialized in the main through the 
//    236	C         common /ABLAMAIN/  (A.B.)
//    237	

// switch-fission.1=on.0=off
	fiss->ifis = 1;

// shell+pairing.0-1-2-3
	fiss->optshp = 0;

// optemd =0,1  0 no emd, 1 incl. emd                                
	opt->optemd = 1;
// read(10,*,iostat=io) dum(10),optcha                               
	opt->optcha = 1;

// not.to.be.changed.(akap)
	fiss->akap = 10.;

// nuclear.viscosity.(beta)
	fiss->bet = 1.5;

// potential-curvature
	fiss->homega = 1.0;

// fission-barrier-coefficient
	fiss->koeff = 1.;

	//collective enhancement switched on 1 or off 0 in densn (qr=val or =1.)
	fiss->optcol = 0;

	// switch-for-low-energy-sys
	fiss->optles = 0;

	opt->eefac = 2.;

	ald->optafan = 0;

	ald->av = 0.073e0;
	ald->as = 0.095e0;
	ald->ak = 0.0e0;

	std::cout <<"ifis " << fiss->ifis << std::endl;
	std::cout <<"optshp " << fiss->optshp << std::endl;
	std::cout <<"optemd " << opt->optemd << std::endl;
	std::cout <<"optcha " << opt->optcha << std::endl;
	std::cout <<"akap " << fiss->akap << std::endl;
	std::cout <<"bet " << fiss->bet << std::endl;
	std::cout <<"homega " << fiss->homega << std::endl;
	std::cout <<"koeff " << fiss->koeff << std::endl;
	std::cout <<"optcol " << fiss->optcol << std::endl;
	std::cout <<"optles " << fiss->optles << std::endl;
	std::cout <<"eefac " << opt->eefac << std::endl;
	std::cout <<"optafan " << ald->optafan << std::endl;
	std::cout <<"av " << ald->av << std::endl;
	std::cout <<"as " << ald->as << std::endl;
	std::cout <<"ak " << ald->ak << std::endl;

	fiss->optxfis = 1;
	G4InclAblaHardcodedData *dataInterface = new G4InclAblaHardcodedData();
	dataInterface->readData();

	for(int z = 0; z < 98; z++) { //do 30  z = 0,98,1                                                 
		for(int n = 0; n < 154; n++) { //do 31  n = 0,153,1                                              
			ecld->ecfnz[n][z] = 0.e0;
			// Added. This extracts data from InclAblaData class (Hardcoded data in this case)
			ec2sub->ecnz[n][z] = dataInterface->getEcnz(n,z);
			ecld->ecgnz[n][z] = dataInterface->getEcnz(n,z);
			ecld->alpha[n][z] = dataInterface->getAlpha(n,z);
			ecld->vgsld[n][z] = dataInterface->getVgsld(n,z);
		} //31     continue                                                        
	} //30   continue                                                          

}

void G4Abla::qrot(G4double z, G4double a, G4double bet, G4double sig, G4double u, G4double *qr)
{
  // QROT INCLUDING DAMPING                                                
  // INPUT: Z,A,BET,SIG,U                                                  
  // OUTPUT: QR - COLLECTIVE ENHANCEMENT FACTOR                            
  // 
  // SEE  JUNGHANS ET AL., NUCL. PHYS. A 629 (1998) 635                    
  // 
  // 
  // FR(U) EXPONENTIAL FUNCTION TO DEFINE DAMPING                        
  // UCR   CRITICAL ENERGY FOR DAMPING                                   
  // DCR   WIDTH OF DAMPING                                              
  // BET   BETA-DEFORMATION !                                            
  // SIG   PERPENDICULAR SPIN CUTOFF FACTOR                              
  //   U   ENERGY                                                        
  //  QR   COEFFICIENT OF COLLECTIVE ENHANCEMENT                         
  //   A   MASS NUMBER                                                   
  //   Z   CHARGE NUMBER                                                 

  G4double ucr,dcr,ponq,dn,n,dz;

  dcr   = 10.0;

  ucr = 40.0;

  if(((fabs(bet)-1.15) < 0) || ((fabs(bet)-1.15) == 0)) {
    goto qrot10;
  }

  if((fabs(bet)-1.15) > 0) {
    goto qrot11;
  }

 qrot10:
  n = a - z;
  dz = fabs(z - 82.0);
  if (n > 104) {
    dn = fabs(n-126.e0);
  }
  else {
    dn = fabs(n - 82.0);
  }

  bet = 0.022 + 0.003*dn + 0.005*dz;

  sig = 25.0*pow(bet,2) * sig;

 qrot11:   
  ponq = (u - ucr)/dcr;

  if (ponq > 700.0) {
    ponq = 700.0;
  }
  if (sig < 1.0) {
    sig = 1.0;
  }
  (*qr) = 1.0/(1.0 + exp(ponq)) * (sig - 1.0) + 1.0;

  if ((*qr) < 1.0) {
    (*qr) = 1.0;
  }

 qrot12:   
  return;
}

void G4Abla::mglw(G4double a, G4double z, G4int refopt, G4double *el)
{
  // MODEL DE LA GOUTTE LIQUIDE DE C. F. WEIZSACKER.
  // USUALLY AN OBSOLETE OPTION

  G4int a1,z1;
  G4double xv,xs,xc,xa;                                   

  a1 = idnint(a);
  z1 = idnint(z);

  if ((a <= 0.01) || (z < 0.01)) {
    (*el) = 1.0e38;
  }
  else {
    xv = -15.56*a;
    xs = 17.23*pow(a,(2.0/3.0));

    if (a > 1.0) {
      xc = 0.7*z*(z-1.0)*pow((a-1.0),(-1.e0/3.e0));
    }
    else {
      xc = 0.0;
    }
  }

  xa = 23.6*(pow((a-2.0*z),2)/a);
  (*el) = xv+xs+xc+xa;
  return;	
}

void G4Abla::mglms(G4double a, G4double z, G4int refopt4, G4double *el)
{
  // USING FUNCTION EFLMAC(IA,IZ,0)                                    
  // 
  // REFOPT4 = 0 : WITHOUT MICROSCOPIC CORRECTIONS                     
  // REFOPT4 = 1 : WITH SHELL CORRECTION                               
  // REFOPT4 = 2 : WITH PAIRING CORRECTION                             
  // REFOPT4 = 3 : WITH SHELL- AND PAIRING CORRECTION                  

  //   1839	C-----------------------------------------------------------------------
  //   1840	C     A1       LOCAL    MASS NUMBER (INTEGER VARIABLE OF A)             
  //   1841	C     Z1       LOCAL    NUCLEAR CHARGE (INTEGER VARIABLE OF Z)          
  //   1842	C     REFOPT4           OPTION, SPECIFYING THE MASS FORMULA (SEE ABOVE) 
  //   1843	C     A                 MASS NUMBER                                     
  //   1844	C     Z                 NUCLEAR CHARGE                                  
  //   1845	C     DEL               PAIRING CORRECTION                              
  //   1846	C     EL                BINDING ENERGY                                  
  //   1847	C     ECNZ( , )         TABLE OF SHELL CORRECTIONS                      
  //   1848	C-----------------------------------------------------------------------
  //   1849	C                                                                       
  G4int a1 = idnint(a);
  G4int z1 = idnint(z);

  if ( (a1 < 0) || (z1 < 0) || ((a1-z1) < 0) )  { //then 
    // modif pour récupérer une masse p et n correcte:
    (*el) = 0.0;
    //    goto mglms50;
  }
  else {
    // binding energy incl. pairing contr. is calculated from                
    // function eflmac                                                       
    (*el) = eflmac(a1,z1,0,refopt4);

    if (refopt4 > 0) {
      if (refopt4 != 2) {
	(*el) = (*el) + ec2sub->ecnz[a1-z1][z1];
      }
    }
  }
  // mglms50: 
  return;
}

G4double G4Abla::spdef(G4int a, G4int z, G4int optxfis)
{

  // INPUT:  A,Z,OPTXFIS MASS AND CHARGE OF A NUCLEUS,                     
  // OPTION FOR FISSILITY                                          
  // OUTPUT: SPDEF                                                         

  // ALPHA2 SADDLE POINT DEF. COHEN&SWIATECKI ANN.PHYS. 22 (1963) 406      
  // RANGING FROM FISSILITY X=0.30 TO X=1.00 IN STEPS OF 0.02              

  G4int index;
  G4double x,fissilityValue,v,dx;
  G4double spdefResult;

  const G4int alpha2Size = 36;
  G4double alpha2[alpha2Size] = {2.5464e0, 2.4944e0, 2.4410e0, 2.3915e0, 2.3482e0,
			       2.3014e0, 2.2479e0, 2.1982e0, 2.1432e0, 2.0807e0, 2.0142e0, 1.9419e0,
			       1.8714e0, 1.8010e0, 1.7272e0, 1.6473e0, 1.5601e0, 1.4526e0, 1.3164e0,
			       1.1391e0, 0.9662e0, 0.8295e0, 0.7231e0, 0.6360e0, 0.5615e0, 0.4953e0,
			       0.4354e0, 0.3799e0, 0.3274e0, 0.2779e0, 0.2298e0, 0.1827e0, 0.1373e0,
			       0.0901e0, 0.0430e0, 0.0000e0};

  dx = 0.02;
  x  = fissility(a,z,optxfis);

  if (x > 1.0) {
    x = 1.0;
  }

  if (x < 0.0) {
    x = 0.0;
  }

  v  = (x - 0.3)/dx + 1.0;
  index = idnint(v);

  if (index < 1) {
    return(alpha2[0]);
  }

  if (index == 36) { //then // :::FIXME:: Possible off-by-one bug...                                            
    return(alpha2[35]); // alpha2(36)->alpha2[35]
  }
  else {
    return(alpha2[index] + (alpha2[index+1] - alpha2[index]) / dx * ( x - (0.3e0 + dx*(index-1)))); //:::FIXME::: Possible off-by-one
  }                                                       

  return alpha2[0]; // The algorithm is not supposed to reach this point.
}

G4double G4Abla::fissility(int a,int z, int optxfis)
{
  // CALCULATION OF FISSILITY PARAMETER                                 
  // 
  // INPUT: A,Z INTEGER MASS & CHARGE OF NUCLEUS                        
  // OPTXFIS = 0 : MYERS, SWIATECKI                              
  //           1 : DAHLINGER                                     
  //           2 : ANDREYEV                                      

  G4double aa,zz,i;
  G4double fissilityResult;

  aa = double(a);
  zz = double(z);
  i  = double(a-2*z) / aa;

  // myers & swiatecki droplet modell                        
  if (optxfis == 0) { //then                                            
    fissilityResult = pow(zz,2) / aa /50.8830e0 / (1.0e0 - 1.7826e0 * pow(i,2));
  }

  if (optxfis == 1) {
    // dahlinger fit:                                          
    fissilityResult = pow(zz,2) / aa * pow((49.22e0*(1.e0 - 0.3803e0*pow(i,2) - 20.489e0*pow(i,4))),(-1));
  }

  if (optxfis == 2) {
    // dubna fit:                                              
    fissilityResult = pow(zz,2) / aa  /(48.e0*(1.e0 - 17.22e0*pow(i,4)));
  }

  return fissilityResult;
}

void G4Abla::evapora(G4double zprf, G4double aprf, G4double ee, G4double jprf, 
		     G4double *zf_par, G4double *af_par, G4double *mtota_par,
		     G4double *pleva_par, G4double *pxeva_par)
{
  G4double zf = (*zf_par);
  G4double af = (*af_par);
  G4double mtota = (*mtota_par);
  G4double pleva = (*pleva_par);
  G4double pxeva = (*pxeva_par);
  
  //    533	C                                                                       
  //    534	C     INPUT:                                                            
  //    535	C                                                                       
  //    536	C     ZPRF, APRF, EE(EE IS MODIFIED!), JPRF                             
  //    537	C                                                                       
  //    538	C     PROJECTILE AND TARGET PARAMETERS + CROSS SECTIONS                 
  //    539	C     COMMON /ABRAMAIN/ AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,       
  //    540	C                       R_0,R_P,R_T, IMAX,IRNDM,PI,                     
  //    541	C                       BFPRO,SNPRO,SPPRO,SHELL                         
  //    542	C                                                                       
  //    543	C     AP,ZP,AT,ZT   - PROJECTILE AND TARGET MASSES                      
  //    544	C     EAP,BETA      - BEAM ENERGY PER NUCLEON, V/C                      
  //    545	C     BMAXNUC       - MAX. IMPACT PARAMETER FOR NUCL. REAC.             
  //    546	C     CRTOT,CRNUC   - TOTAL AND NUCLEAR REACTION CROSS SECTION          
  //    547	C     R_0,R_P,R_T,  - RADIUS PARAMETER, PROJECTILE+ TARGET RADII        
  //    548	C     IMAX,IRNDM,PI - MAXIMUM NUMBER OF EVENTS, DUMMY, 3.141...         
  //    549	C     BFPRO         - FISSION BARRIER OF THE PROJECTILE                 
  //    550	C     SNPRO         - NEUTRON SEPARATION ENERGY OF THE PROJECTILE       
  //    551	C     SPPRO         - PROTON    "           "   "    "   "              
  //    552	C     SHELL         - GROUND STATE SHELL CORRECTION                     
  //    553	C                                                                       
  //    554	C---------------------------------------------------------------------  
  //    555	C     FISSION BARRIERS                                                  
  //    556	C     COMMON /FB/     EFA                                               
  //    557	C     EFA    - ARRAY OF FISSION BARRIERS                                
  //    558	C---------------------------------------------------------------------  
  //    559	C     OUTPUT:                                                           
  //    560	C              ZF, AF, MTOTA, PLEVA, PTEVA, FF, INTTYPE, INUM           
  //    561	C                                                                       
  //    562	C     ZF,AF - CHARGE AND MASS OF FINAL FRAGMENT AFTER EVAPORATION       
  //    563	C     MTOTA _ NUMBER OF EVAPORATED ALPHAS                               
  //    564	C     PLEVA,PXEVA,PYEVA - MOMENTUM RECOIL BY EVAPORATION               
  //    565	C     INTTYPE - TYPE OF REACTION 0/1 NUCLEAR OR ELECTROMAGNETIC         
  //    566	C     FF      - 0/1 NO FISSION / FISSION EVENT                          
  //    567	C     INUM    - EVENTNUMBER                                             
  //    568	C   ____________________________________________________________________
  //    569	C  /                                                                    
  //    570	C  /  CALCUL DE LA MASSE ET CHARGE FINALES D'UNE CHAINE D'EVAPORATION   
  //    571	C  /                                                                    
  //    572	C  /  PROCEDURE FOR CALCULATING THE FINAL MASS AND CHARGE VALUES OF A   
  //    573	C  /  SPECIFIC EVAPORATION CHAIN, STARTING POINT DEFINED BY (APRF, ZPRF,
  //    574	C  /  EE)                                                               
  //    575	C  /  On ajoute les 3 composantes de l'impulsion (PXEVA,PYEVA,PLEVA)
  //    576	C  /    (actuellement PTEVA n'est pas correct; mauvaise norme...)                                               
  //    577	C  /____________________________________________________________________
  //    578	C                                                                       
  //    612	C                                                                       
  //    613	C-----------------------------------------------------------------------
  //    614	C     IRNDM             DUMMY ARGUMENT FOR RANDOM-NUMBER FUNCTION       
  //    615	C     SORTIE   LOCAL    HELP VARIABLE TO END THE EVAPORATION CHAIN      
  //    616	C     ZF                NUCLEAR CHARGE OF THE FRAGMENT                  
  //    617	C     ZPRF              NUCLEAR CHARGE OF THE PREFRAGMENT               
  //    618	C     AF                MASS NUMBER OF THE FRAGMENT                     
  //    619	C     APRF              MASS NUMBER OF THE PREFRAGMENT                  
  //    620	C     EPSILN            ENERGY BURNED IN EACH EVAPORATION STEP          
  //    621	C     MALPHA   LOCAL    MASS CONTRIBUTION TO MTOTA IN EACH EVAPORATION  
  //    622	C                        STEP                                           
  //    623	C     EE                EXCITATION ENERGY (VARIABLE)                    
  //    624	C     PROBP             PROTON EMISSION PROBABILITY                     
  //    625	C     PROBN             NEUTRON EMISSION PROBABILITY                    
  //    626	C     PROBA             ALPHA-PARTICLE EMISSION PROBABILITY             
  //    627	C     PTOTL             TOTAL EMISSION PROBABILITY                      
  //    628	C     E                 LOWEST PARTICLE-THRESHOLD ENERGY                
  //    629	C     SN                NEUTRON SEPARATION ENERGY                       
  //    630	C     SBP               PROTON SEPARATION ENERGY PLUS EFFECTIVE COULOMB 
  //    631	C                        BARRIER                                        
  //    632	C     SBA               ALPHA-PARTICLE SEPARATION ENERGY PLUS EFFECTIVE 
  //    633	C                        COULOMB BARRIER                                
  //    634	C     BP                EFFECTIVE PROTON COULOMB BARRIER                
  //    635	C     BA                EFFECTIVE ALPHA COULOMB BARRIER                 
  //    636	C     MTOTA             TOTAL MASS OF THE EVAPORATED ALPHA PARTICLES    
  //    637	C     X                 UNIFORM RANDOM NUMBER FOR NUCLEAR CHARGE        
  //    638	C     AMOINS   LOCAL    MASS NUMBER OF EVAPORATED PARTICLE              
  //    639	C     ZMOINS   LOCAL    NUCLEAR CHARGE OF EVAPORATED PARTICLE           
  //    640	C     ECP               KINETIC ENERGY OF PROTON WITHOUT COULOMB        
  //    641	C                        REPULSION                                      
  //    642	C     ECN               KINETIC ENERGY OF NEUTRON                       
  //    643	C     ECA               KINETIC ENERGY OF ALPHA PARTICLE WITHOUT COULOMB
  //    644	C                        REPULSION                                      
  //    645	C     PLEVA             TRANSVERSAL RECOIL MOMENTUM OF EVAPORATION      
  //    646	C     PTEVA             LONGITUDINAL RECOIL MOMENTUM OF EVAPORATION     
  //    647	C     FF                FISSION FLAG                                    
  //    648	C     INTTYPE           INTERACTION TYPE FLAG                           
  //    649	C     RNDX              RECOIL MOMENTUM IN X-DIRECTION IN A SINGLE STEP 
  //    650	C     RNDY              RECOIL MOMENTUM IN Y-DIRECTION IN A SINGLE STEP 
  //    651	C     RNDZ              RECOIL MOMENTUM IN Z-DIRECTION IN A SINGLE STEP 
  //    652	C     RNDN              NORMALIZATION OF RECOIL MOMENTUM FOR EACH STEP  
  //    653	C-----------------------------------------------------------------------
  //    654	C                                                                       
  //    655	      SAVE                                                              
  // SAVE -> static
	
  static G4int sortie;                            
  static G4double epsiln,alpha,probp,probn,proba,ptotl,e;  
  static G4double sn,sbp,sba,x,amoins,zmoins,ecn,ecp,eca,bp,ba;         
  static G4double pteva,rndx,rndy,rndz,rndn,pyeva;                       
  static G4double ff;

  static G4int itest;
  static G4double probf;
  static G4int inttype;

  static G4int inum;
  static G4int k, j, il;

  static G4double ctet1,stet1,phi1;                               
  static G4double sbfis,rnd;
  static G4double sel, selmax;
  static G4double segs;
  static G4double ef;
  static G4int irndm;

  static G4double pc, malpha;

  zf = zprf;
  af = aprf;
  pleva = 0.0;
  pteva = 0.0;
  pxeva = 0.0;
  pyeva = 0.0;

  sortie = 0;
  ff = 0;

  itest = 0;
  if (itest == 1) {
    std::cout << "***************************" << std::endl;
  }

 evapora10:

  if (itest == 1) {
    std::cout <<"------zf,af,ee------" << idnint(zf) << "," << idnint(af) << "," << ee << std::endl;
  }

  // calculation of the probabilities for the different decay channels     
  // plus separation energies and kinetic energies of the particles        
  direct(zf,af,ee,jprf,&probp,&probn,&proba,&probf,&ptotl,
	 &sn,&sbp,&sba,&ecn,&ecp,&eca,&bp,&ba,inttype,inum,itest); //:::FIXME::: Call

  k = idnint(zf);
  j = idnint(af-zf);

  // now ef is calculated from efa that depends on the subroutine
  // barfit which takes into account the modification on the ang. mom.
  // jb mvr 6-aug-1999
  // note *** shell correction! (ecgnz)  jb mvr 20-7-1999
  il  = idnint(jprf);
  barfit(k,k+j,il,&sbfis,&segs,&selmax);

  if ((fiss->optshp == 1) || (fiss->optshp == 3)) { //then                     
    fb->efa[k][j] = G4double(sbfis) -  ecld->ecgnz[j][k];
  }
  else {
    fb->efa[k][j] = G4double(sbfis);
  } //end if 
  ef = fb->efa[k][j];

  // here the final steps of the evaporation are calculated                
  if ((sortie == 1) || (ptotl == 0.e0)) {
    e = dmin1(sn,sbp,sba);
    if (e > 1.0e30) {
      std::cout <<"erreur a la sortie evapora,e>1.e30,af=" << af <<" zf=" << zf << std::endl;
    }
    if (zf <= 6.0) {
      goto evapora100;
    }
    if (e < 0.0) {
      if (sn == e) {
	af = af - 1.e0;
      }
      else if (sbp == e) {
	af = af - 1.0;
	zf = zf - 1.0;
      }
      else if (sba == e) {
	af = af - 4.0;
	zf = zf - 2.0;
      }
      if (af < 2.5) {
	goto evapora100;
      }
      goto evapora10;
    }
    goto evapora100;
  }
 evapora30:
  irndm = irndm + 1;

  // here the normal evaporation cascade starts                            

  // random number for the evaporation                                     
  //  x = double(Rndm(irndm))*ptotl;
  x = double(haz(1))*ptotl;
  
  if (x < proba) {
    // alpha evaporation                                                     
    if (itest == 1) {
      std::cout <<"< alpha evaporation >" << std::endl;
    }
    amoins = 4.0;
    zmoins = 2.0;
    epsiln = sba + eca;
    pc = sqrt(pow((1.0 + (eca+ba)/3.72834e3),2) - 1.0) * 3.72834e3;
    malpha = 4.0;

    // volant:
    volant->iv = volant->iv + 1;
    volant->acv[volant->iv] = 4.;
    volant->zpcv[volant->iv] = 2.;
    volant->pcv[volant->iv] = pc;
  }
  else if (x < proba+probp) {
    // proton evaporation                                                    
    if (itest == 1) {
      std::cout <<"< proton evaporation >" << std::endl;
      amoins = 1.0;
      zmoins = 1.0;
      epsiln = sbp + ecp;
      pc = sqrt(pow((1.0 + (ecp + bp)/9.3827e2),2) - 1.0) * 9.3827e2;
      malpha = 0.0;
      // volant:
      volant->iv = volant->iv + 1;
      volant->acv[volant->iv] = 1.;
      volant->zpcv[volant->iv] = 1.;
      volant->pcv[volant->iv] = pc;
    }
  }
  else if (x < proba+probp+probn) {
    // neutron evaporation                                                   
    if (itest == 1) {
      std::cout <<"< neutron evaporation >" << std::endl;
    }
    amoins = 1.0;
    zmoins = 0.0;
    epsiln = sn + ecn;
    pc = sqrt(pow((1.0 + (ecn)/9.3956e2),2) - 1.0) * 9.3956e2;
    malpha = 0.0;

    // volant:
    volant->iv = volant->iv + 1;
    volant->acv[volant->iv] = 1.;
    volant->zpcv[volant->iv] = 0.;
    volant->pcv[volant->iv] = pc;
  }
  else {
    // fission                                                               
    // in case of fission-events the fragment nucleus is the mother nucleus  
    // before fission occurs with excitation energy above the fis.- barrier. 
    // fission fragment mass distribution is calulated in subroutine fisdis  
    if (itest == 1) {
      std::cout <<"< fission >" << std::endl;
    }
    amoins = 0.0;
    zmoins = 0.0;
    epsiln = ef;

    malpha = 0.0;
    pc = 0.0;
    ff = 1;
  }

  if (itest == 1) {
    std::cout <<"sn,sbp,sba,ef" << sn << "," << sbp << "," << sba <<"," << ef << std::endl;
    std::cout <<"probn,probp,proba,probf,ptotl " <<","<< probn <<","<< probp <<","<< proba <<","<< probf <<","<< ptotl << std::endl;
  }

  // calculation of the daughter nucleus                                   
  af = af - amoins;
  zf = zf - zmoins;
  ee = ee - epsiln;
  if (ee <= 0.01) {
    ee = 0.01;
  }
  mtota = mtota + malpha;

  if(ff == 0) {
    standardRandom(&rnd,&(hazard->igraine[8]));
    ctet1 = 2.0*rnd - 1.0;
    standardRandom(&rnd,&(hazard->igraine[4]));
    phi1 = rnd*2.0*3.141592654;
    stet1 = sqrt(1.0 - pow(ctet1,2));
    volant->xcv[volant->iv] = stet1*cos(phi1);
    volant->ycv[volant->iv] = stet1*sin(phi1);
    volant->zcv[volant->iv] = ctet1;
    pxeva = pxeva - pc * volant->xcv[volant->iv];
    pyeva = pyeva - pc * volant->ycv[volant->iv];
    pleva = pleva - pc * ctet1;
  }

  // condition for end of evaporation                                   
  if ((af < 2.5) || (ff == 1)) {
    goto evapora100;
  }
  goto evapora10;

 evapora100:
  (*zf_par) = zf;
  (*af_par) = af;
  (*mtota_par) = mtota;
  (*pleva_par) = pleva;
  (*pxeva_par) = pxeva;
                                          
  return;
}

void G4Abla::direct(G4double zprf,G4double a, G4double ee, G4double jprf, 
		    G4double *probp_par, G4double *probn_par, G4double *proba_par, 
		    G4double *probf_par, G4double *ptotl_par, G4double *sn_par,
		    G4double *sbp_par, G4double *sba_par, G4double *ecn_par, 
		    G4double *ecp_par,G4double *eca_par, G4double *bp_par,
		    G4double *ba_par, G4int inttype, G4int inum, G4int itest)
{
  G4double probp = (*probp_par);
  G4double probn = (*probn_par);
  G4double proba = (*proba_par);
  G4double probf = (*probf_par);
  G4double ptool = (*ptotl_par);
  G4double sn = (*sn_par);
  G4double sbp = (*sbp_par);
  G4double sba = (*sba_par);
  G4double ecn = (*ecn_par);
  G4double ecp = (*ecp_par);
  G4double eca = (*eca_par);
  G4double bp = (*bp_par);
  G4double ba = (*ba_par);

  // CALCULATION OF PARTICLE-EMISSION PROBABILITIES & FISSION     / 
  // BASED ON THE SIMPLIFIED FORMULAS FOR THE DECAY WIDTH BY      / 
  // MORETTO, ROCHESTER MEETING TO AVOID COMPUTING TIME           / 
  // INTENSIVE INTEGRATION OF THE LEVEL DENSITIES                 / 
  // USES EFFECTIVE COULOMB BARRIERS AND AN AVERAGE KINETIC ENERGY/ 
  // OF THE EVAPORATED PARTICLES                                  / 
  // COLLECTIVE ENHANCMENT OF THE LEVEL DENSITY IS INCLUDED       / 
  // DYNAMICAL HINDRANCE OF FISSION IS INCLUDED BY A STEP FUNCTION/ 
  // APPROXIMATION. SEE A.R. JUNGHANS DIPLOMA THESIS              / 
  // SHELL AND PAIRING STRUCTURES IN THE LEVEL DENSITY IS INCLUDED/ 

  // INPUT:                                                            
  // ZPRF,A,EE  CHARGE, MASS, EXCITATION ENERGY OF COMPOUND     
  // NUCLEUS                                         
  // JPRF       ROOT-MEAN-SQUARED ANGULAR MOMENTUM                           

  // DEFORMATIONS AND G.S. SHELL EFFECTS                               
  // COMMON /ECLD/   ECGNZ,ECFNZ,VGSLD,ALPHA                           

  // ECGNZ - GROUND STATE SHELL CORR. FRLDM FOR A SPHERICAL G.S.       
  // ECFNZ - SHELL CORRECTION FOR THE SADDLE POINT (NOW: == 0)         
  // VGSLD - DIFFERENCE BETWEEN DEFORMED G.S. AND LDM VALUE            
  // ALPHA - ALPHA GROUND STATE DEFORMATION (THIS IS NOT BETA2!)       
  // BETA2 = SQRT((4PI)/5) * ALPHA                             

  //OPTIONS AND PARAMETERS FOR FISSION CHANNEL                        
  //COMMON /FISS/    AKAP,BET,HOMEGA,KOEFF,IFIS,                       
  //                 OPTSHP,OPTXFIS,OPTLES,OPTCOL               
  //
  // AKAP   - HBAR**2/(2* MN * R_0**2) = 10 MEV, R_0 = 1.4 FM          
  // BET    - REDUCED NUCLEAR FRICTION COEFFICIENT IN (10**21 S**-1)   
  // HOMEGA - CURVATURE OF THE FISSION BARRIER = 1 MEV                 
  // KOEFF  - COEFFICIENT FOR THE LD FISSION BARRIER == 1.0            
  // IFIS   - 0/1 FISSION CHANNEL OFF/ON                               
  // OPTSHP - INTEGER SWITCH FOR SHELL CORRECTION IN MASSES/ENERGY     
  //          = 0 NO MICROSCOPIC CORRECTIONS IN MASSES AND ENERGY        
  //          = 1 SHELL ,  NO PAIRING                                    
  //          = 2 PAIRING, NO SHELL                                      
  //          = 3 SHELL AND PAIRING                                      
  // OPTCOL - 0/1 COLLECTIVE ENHANCEMENT SWITCHED ON/OFF               
  // OPTXFIS- 0,1,2 FOR MYERS & SWIATECKI, DAHLINGER, ANDREYEV         
  //                FISSILITY PARAMETER.                                     
  // OPTLES - CONSTANT TEMPERATURE LEVEL DENSITY FOR A,Z > TH-224      
  // OPTCOL - 0/1 COLLECTIVE ENHANCEMENT OFF/ON                        

  // LEVEL DENSITY PARAMETERS                                          
  // COMMON /ALD/    AV,AS,AK,OPTAFAN                                  
  //                 AV,AS,AK - VOLUME,SURFACE,CURVATURE DEPENDENCE OF THE             
  //                            LEVEL DENSITY PARAMETER                                
  // OPTAFAN - 0/1  AF/AN >=1 OR AF/AN ==1                             
  //           RECOMMENDED IS OPTAFAN = 0                              

  // FISSION BARRIERS                                                  
  // COMMON /FB/     EFA                                               
  // EFA    - ARRAY OF FISSION BARRIERS                                


  // OUTPUT: PROBN,PROBP,PROBA,PROBF,PTOTL:                            
  // - EMISSION PROBABILITIES FOR N EUTRON, P ROTON,  A LPHA     
  // PARTICLES, F ISSION AND NORMALISATION                     
  // SN,SBP,SBA: SEPARATION ENERGIES N P A                     
  // INCLUDING EFFECTIVE BARRIERS                              
  // ECN,ECP,ECA,BP,BA                                         
  // - AVERAGE KINETIC ENERGIES (2*T) AND EFFECTIVE BARRIERS     

  G4double bk;
  G4double afp;
  G4double alpha;
  G4double at;
  G4double bet;
  G4double bs;
  G4double bshell;
  G4double call;
  G4double cf;
  G4double dconst;
  G4double defbet;
  G4double denomi;
  G4double densa;
  G4double densf;
  G4double densg;
  G4double densn;
  G4double densp;
  G4double edyn;
  G4double eer;
  G4double ef;
  G4double ft;
  G4double ga;
  G4double gf;
  G4double gn;
  G4double gngf;
  G4double gp;
  G4double gsum;
  G4double hbar;
  G4double homega;
  G4double iflag;
  G4int il;
  G4int ilast;
  G4int imaxwell;
  G4int in;
  G4int iz;
  G4int j;
  G4int k;
  G4int kkk;
  G4double ma1z;
  G4double ma1z1;
  G4double ma4z2;
  G4double maz;
  G4double nprf;
  G4double nt;
  G4double optcol;
  G4double optles;
  G4double optxfis;
  G4double parc;
  G4double pi;
  G4double pt;
  G4double ptotl;
  G4double ra;
  G4double rat;
  G4double refmod;
  G4double rf;
  G4double rn;
  G4double rnd;
  G4double rnt;
  G4double rp;
  G4double rpt;
  G4double sa;
  G4double sbf;
  G4double sbfis;
  G4double segs;
  G4double selmax;
  G4double sp;
  G4double tauc;
  G4double tconst;
  G4double temp;
  G4double ts1;
  G4double tsum;
  G4double wf;
  G4double wfex;
  G4double xx;
  G4double y;

  imaxwell = 1;

  // limiting of excitation energy where fission occurs                    
  // Note, this is not the dynamical hindrance (see end of routine)      
  edyn = 1000.0;

  // no limit if statistical model is calculated.                         
  if (bet <= 1.0e-16) {
    edyn = 10000.0;
  }

  // just a change of name until the end of this subroutine                
  eer = ee;
  if (inum == 1) {
    ilast = 1;
  }

  // calculation of masses                                           
  // refmod = 1 ==> myers,swiatecki model                              
  // refmod = 0 ==> weizsaecker model                                  
  refmod = 1;

  if (refmod == 1) {
    mglms(a,zprf,fiss->optshp,&maz);
    mglms(a-1.0,zprf,fiss->optshp,&ma1z);
    mglms(a-1.0,zprf-1.0,fiss->optshp,&ma1z1);
    mglms(a-4.0,zprf-2.0,fiss->optshp,&ma4z2);
  }
  else {
    mglw(a,zprf,fiss->optshp,&maz);
    mglw(a-1.0,zprf,fiss->optshp,&ma1z);
    mglw(a-1.0,zprf-1.0,fiss->optshp,&ma1z1);
    mglw(a-4.0,zprf-2.0,fiss->optshp,&ma4z2);
  }

  // separation energies and effective barriers                     
  sn = ma1z - maz;
  sp = ma1z1 - maz;
  sa = ma4z2 - maz-28.29688;
  if (zprf < 1.e0) {
    sbp = 1.0e75;
    goto direct30;
  }

  // parameterisation gaimard:
  // bp = 1.44*(zprf-1.d0)/(1.22*pow((a - 1.0),(1.0/3.0))+5.6)     
  // parameterisation khs (12-99)
  bp = 1.44*(zprf - 1.0)/(2.1*pow((a-1.0),(1.0/3.0)) + 0.0);

  sbp = sp + bp;
  if (a-4.0 <= 0.0) {
    sba = 1.0e+75;
    goto direct30;
  }

  // new effective barrier for alpha evaporation d=6.1: khs          
  // ba = 2.88d0*(zprf-2.d0)/(1.22d0*(a-4.d0)**(1.d0/3.d0)+6.1d0)
  // parametrisation khs (12-99)
  ba = 2.88*(zprf - 2.0)/(2.2*pow((a - 4.0),(1.0/3.0)) + 0.0);

  sba = sa + ba;
 direct30:

  // calculation of surface and curvature integrals needed to      
  // to calculate the level density parameter (in densniv)         
  if (fiss->ifis > 0) {
    k = idnint(zprf);
    j = idnint(a - zprf);

    // now ef is calculated from efa that depends on the subroutine
    // barfit which takes into account the modification on the ang. mom.
    // jb mvr 6-aug-1999
    // note *** shell correction! (ecgnz)  jb mvr 20-7-1999
    il = idnint(jprf);
    barfit(k,k+j,il,&sbfis,&segs,&selmax);
    if ((fiss->optshp == 1) || (fiss->optshp == 3)) {
      fb->efa[k][j] = G4double(sbfis) -  ecld->ecgnz[j][k];
    } 
    else {
      fb->efa[k][j] = G4double(sbfis);
    }
    ef = fb->efa[k][j];

    // to avoid negative values for impossible nuclei                        
    // the fission barrier is set to zero if smaller than zero.              
    if (fb->efa[k][j] < 0.0) {
      fb->efa[k][j] = 0.0;
    }

    // factor with jprf should be 0.0025d0 - 0.01d0 for                     
    // approximate influence of ang. momentum on bfis  a.j. 22.07.96        
    // 0.0 means no angular momentum                                       

    if (ef < 0.0) {
      ef = 0.0;
    }
    xx = fissility((k+j),k,fiss->optxfis);
    y = 1.00 - xx;
    if (y < 0.0) {
      y = 0.0;
    }
    if (y > 1.0) {
      y = 1.0;
    }
    bs = bipol(1,y);
    bk = bipol(2,y);
  }
  else {
    ef = 1.0e40;
    bs = 1.0;
    bk = 1.0;
  }
  sbf = ee - ef;

  afp = idnint(a);
  iz = idnint(zprf);
  in = afp - iz;
  bshell = ecld->ecfnz[in][iz];

  // ld saddle point deformation                                          
  // here: beta2 = sqrt(5/(4pi)) * alpha2                                  

  // for the ground state def. 1.5d0 should be used                        
  // because this was just the factor to produce the                       
  // alpha-deformation table 1.5d0 should be used                          
  // a.r.j. 6.8.97                                                         
  defbet = 1.58533e0 * spdef(idnint(a),idnint(zprf),fiss->optxfis);

  // level density and temperature at the saddle point                     
  densniv(a,zprf,ee,ef,&densf,bshell,bs,bk,&temp,fiss->optshp,optcol,defbet);
  ft = temp;
  if (iz >= 2) {
    bshell = ecld->ecgnz[in][iz-1] - ecld->vgsld[in][iz-1];
    defbet = 1.5 * (ecld->alpha[in][iz-1]);

    // level density and temperature in the proton daughter                  
    densniv(a-1.0,zprf-1.0e0,ee,sbp,&densp, bshell,1.e0,1.e0,&temp,fiss->optshp,optcol,defbet);

    pt = temp;
    if (imaxwell == 1) {
      // valentina - random kinetic energy in a maxwelliam distribution
      // modif juin/2002 a.b. c.v. for light targets; limit on the energy
      // from the maxwell distribution.
      rpt = pt;
      ecp = 2.0 * pt;
      if(rpt <= 1.0e-3) {
	goto direct2914;
      }
      iflag = 0;
    direct1914:
      ecp = fmaxhaz(kkk,rpt);
      iflag = iflag + 1;
      if(iflag >= 10) {
	standardRandom(&rnd,&(hazard->igraine[5]));
	ecp=sqrt(rnd)*(eer-sbp);
	goto direct2914;
      }
      if((ecp+sbp) > eer) {
	goto direct1914;
      }
    }
    else {
      ecp = 2.0 * pt;
    }

  direct2914:   
    std::cout <<""<<std::endl;
  }
  else {
    densp = 0.0;
    ecp = 0.0;
    pt = 0.0;
  }

  if (in >= 2) {
    bshell = ecld->ecgnz[in-1][iz] - ecld->vgsld[in-1][iz];
    defbet = 1.5e0 * (ecld->alpha[in-1][iz]);

    // level density and temperature in the neutron daughter                 
    densniv(a-1.0,zprf,ee,sn,&densn,bshell, 1.e0,1.e0,&temp,fiss->optshp,optcol,defbet);
    nt = temp;

    if (imaxwell == 1) {
      // valentina - random kinetic energy in a maxwelliam distribution
      // modif juin/2002 a.b. c.v. for light targets; limit on the energy
      // from the maxwell distribution.
      rnt = nt;
      ecn = 2.0 * nt;
      if(rnt <= 1.e-3) {
	goto direct2915;
      }

      iflag=0;
    direct1915:
      ecn = fmaxhaz(kkk,rnt);
      iflag=iflag+1;
      if(iflag >= 10) {
	standardRandom(&rnd,&(hazard->igraine[6]));
	ecn = sqrt(rnd)*(eer-sn);
	goto direct2915;
      }
//       if((ecn+sn) > eer) {
// 	goto direct1915;
//       }
//       else {
// 	ecn = 2.e0 * nt;
//       }
      if((ecn + sn) <= eer) {
	ecn = 2.0 * nt;
      }
    direct2915: 
      std::cout <<"" <<std::endl;
    }
  }
  else {
    densn = 0.0;
    ecn = 0.0;
    nt = 0.0;
  }

  if ((in >= 3) && (iz >= 3)) {
    bshell = ecld->ecgnz[in-2][iz-2] - ecld->vgsld[in-2][iz-2];
    defbet = 1.5 * (ecld->alpha[in-2][iz-2]);

    // level density and temperature in the alpha daughter                   
    densniv(a-4.0,zprf-2.0e0,ee,sba,&densa,bshell,1.e0,1.e0,&temp,fiss->optshp,optcol,defbet);

    // valentina - random kinetic energy in a maxwelliam distribution
    at = temp;
    if (imaxwell == 1) {
      // modif juin/2002 a.b. c.v. for light targets; limit on the energy
      // from the maxwell distribution.
      rat = at;
      eca= 2.e0 * at;
      if(rat <= 1.e-3) {
	goto direct2916;
      }
      iflag=0;
    direct1916:
      eca = fmaxhaz(kkk,rat);
      iflag=iflag+1;
      if(iflag >= 10) {
	standardRandom(&rnd,&(hazard->igraine[7]));
	eca=sqrt(rnd)*(eer-sba);
	goto direct2916;
      }
      if((eca+sba) > eer) {
	goto direct1916;
      }
      else {
	eca = 2.0 * at;
      }
    direct2916:
      std::cout <<"" << std::endl;
    }
    else {
      densa = 0.0;
      eca = 0.0;
      at = 0.0;
    }
  } // PK

  // special treatment for unbound nuclei                                                
  if (sn < 0.0) {
    probn = 1.0;
    probp = 0.0;
    proba = 0.0;
    probf = 0.0;
    goto direct70;
  }
  if (sbp < 0.0) {
    probp = 1.0;
    probn = 0.0;
    proba = 0.0;
    probf = 0.0;
    goto direct70;
  }

  if ((a < 50.e0) || (ee > edyn)) { // no fission if e*> edyn or mass < 50
    densf = 0.e0;
  }

  bshell = ecld->ecgnz[in][iz] - ecld->vgsld[in][iz];
  defbet = 1.5e0 * (ecld->alpha[in][iz]);

  // compound nucleus level density                                        
  densniv(a,zprf,ee,0.0e0,&densg,bshell,1.e0,1.e0,&temp,fiss->optshp,optcol,defbet);

  if ( densg > 0.e0) {
    // calculation of the partial decay width                                
    // used for both the time scale and the evaporation decay width
    gp = (pow(a,(2./3.))/fiss->akap)*densp/densg/pi*pow(pt,2);
    gn = (pow(a,(2./3.))/fiss->akap)*densn/densg/pi*pow(nt,2);
    ga = (pow(a,(2./3.))/fiss->akap)*densa/densg/pi*2.0e0*pow(at,2);
    gf= densf/densg/pi/2.0e0*ft;

    if (itest == 1) {
      std::cout <<"gn,gp,ga,gf " << gn <<","<< gp <<","<< ga <<","<< gf << std::endl;
    }
  }
  else {
    std::cout <<"direct: densg <= 0.e0 " << a <<","<< zprf <<","<< ee << std::endl;
  }

  gsum = ga + gp + gn;
  if (gsum > 0.0) {
    ts1  = hbar / gsum;
  }
  else {
    ts1  = 1.0e99;
  }

  if (inum > ilast) {  // new event means reset the time scale
    tsum = 0;
  }

  // calculate the relative probabilities for all decay channels        
  if (densf == 0.0) {
    if (densp == 0.0) {
      if (densn == 0.0) {
	if (densa == 0.0) {
	  // no reaction is possible                                               
	  probf = 0.0;
	  probp = 0.0;
	  probn = 0.0;
	  proba = 0.0;
	  goto direct70;
	}

	// alpha evaporation is the only open channel                            
	rf = 0.0;
	rp = 0.0;
	rn = 0.0;
	ra = 1.0;
	goto direct50;
      }

      // alpha emission and neutron emission                                   
      rf = 0.0;
      rp = 0.0;
      rn = 1.0;
      ra = densa*2.0/densn*pow((at/nt),2);
      goto direct50;
    }
    // alpha, proton and neutron emission                                    
    rf = 0.0;
    rp = 1.0;
    rn = densn/densp*pow((nt/pt),2);
    ra = densa*2.0/densp*pow((at/pt),2);
    goto direct50;
  }

  // here fission has taken place                                          
  rf = 1.0;
  // cramers and weidenmueller factors for the dynamical hindrances of     
  // fission                                                               
  if (bet <= 1.0e-16) {
    cf = 1.0;
    wf = 1.0;
  }
  else if (sbf > 0.0e0) {
    cf = cram(bet,homega);
    // if fission barrier ef=0.d0 then fission is the only possible      
    // channel. to avoid log(0) in function tau                          
    // a.j. 7/28/93                                                      
    if (ef <= 0.0) {
      rp = 0.0;
      rn = 0.0;
      ra = 0.0;
      goto direct50;
    }
    else {
      // transient time tau()                                                  
      tauc = tau(bet,homega,ef,ft);
    }
    wfex = (tauc - tsum)/ts1;

    if (wfex < 0.0) {
      wf = 1.0;
    }
    else {
      wf = exp( -wfex);
    }
  }
  else {
    cf=1.0;
    wf=1.0;
  }

  if (itest == 1) {
    std::cout <<"tsum,wf,cf " << tsum <<","<< wf <<","<< cf << std::endl;
  }

  tsum = tsum + ts1;

  // change by g.k. and a.h. 5.9.95                                       
  tconst = 0.7;
  dconst = 12.0/sqrt(a);
  nprf = a - zprf;

  if (fiss->optshp >= 2) { //then                                           
    parite(nprf,&parc);
    dconst = dconst*parc;
  }
  else {
    dconst= 0.0;
  }
  if ((ee <= 17.e0) && (optles == 1) && (iz >= 90) && (in >= 134)) { //then                              
    // constant changed to 5.0 accord to moretto & vandenbosch a.j. 19.3.96  
    gngf=pow(a,(2.e0/3.e0))*tconst/10.e0*exp((ef-sn+dconst)/tconst);

    // if the excitation energy is so low that densn=0 ==> gn = 0           
    // fission remains the only channel.                                    
    // a. j. 10.1.94                                                        
    if (gn == 0.0) {
      rn = 0.0;
      rp = 0.0;
      ra = 0.0;
    }
    else {
      rn=gngf;
      rp=gngf*gp/gn;
      ra=gngf*ga/gn;
    }
  } else {
    rn = gn/(gf*cf);
    rp = gp/(gf*cf);
    ra = ga/(gf*cf);
  }
 direct50:
  // relative decay probabilities                                          

  denomi = rp+rn+ra+rf;

  // decay probabilities after transient time
  probf = rf/denomi;
  probp = rp/denomi;
  probn = rn/denomi;
  proba = ra/denomi;

  // new treatment of grange-weidenmueller factor, 5.1.2000, khs !!!

  // decay probabilites with transient time included
  probf = probf * wf;
  probp = probp * (wf + (1.e0-wf)/(1.e0-probf));
  probn = probn * (wf + (1.e0-wf)/(1.e0-probf));
  proba = proba * (wf + (1.e0-wf)/(1.e0-probf));

 direct70:
  ptotl = probp+probn+proba+probf;

  ee = eer;
  ilast = inum;

  // Return values:
  (*probp_par) = probp;
  (*probn_par) = probn;
  (*proba_par) = proba;
  (*probf_par) = probf;
  (*ptotl_par) = ptotl;
  (*sn_par) = sn;
  (*sbp_par) = sbp;
  (*sba_par) = sba;
  (*ecn_par) = ecn;
  (*ecp_par) = ecp;
  (*eca_par) = eca;
  (*bp_par) = bp;
  (*ba_par) = ba;
}

void G4Abla::densniv(G4double a, G4double z, G4double ee, G4double esous, G4double *dens, G4double bshell, G4double bs, G4double bk, 
		     G4double *temp, G4int optshp, G4int optcol, G4double defbet)
{
  //   1498	C                                                                       
  //   1499	C     INPUT:                                                            
  //   1500	C             A,EE,ESOUS,OPTSHP,BS,BK,BSHELL,DEFBET                     
  //   1501	C                                                                       
  //   1502	C     LEVEL DENSITY PARAMETERS                                          
  //   1503	C     COMMON /ALD/    AV,AS,AK,OPTAFAN                                  
  //   1504	C     AV,AS,AK - VOLUME,SURFACE,CURVATURE DEPENDENCE OF THE             
  //   1505	C                LEVEL DENSITY PARAMETER                                
  //   1506	C     OPTAFAN - 0/1  AF/AN >=1 OR AF/AN ==1                             
  //   1507	C               RECOMMENDED IS OPTAFAN = 0                              
  //   1508	C---------------------------------------------------------------------  
  //   1509	C     OUTPUT: DENS,TEMP                                                 
  //   1510	C                                                                       
  //   1511	C   ____________________________________________________________________
  //   1512	C  /                                                                    
  //   1513	C  /  PROCEDURE FOR CALCULATING THE STATE DENSITY OF A COMPOUND NUCLEUS 
  //   1514	C  /____________________________________________________________________
  //   1515	C                                                                       
  //   1516	      INTEGER AFP,IZ,OPTSHP,OPTCOL,J,OPTAFAN                            
  //   1517	      REAL*8 A,EE,ESOUS,DENS,E,Y0,Y1,Y2,Y01,Y11,Y21,PA,BS,BK,TEMP       
  //   1518	C=====INSERTED BY KUDYAEV===============================================
  //   1519	      COMMON /ALD/ AV,AS,AK,OPTAFAN                                     
  //   1520	      REAL*8 ECR,ER,DELTAU,Z,DELTPP,PARA,PARZ,FE,HE,ECOR,ECOR1,Pi6      
  //   1521	      REAL*8 BSHELL,DELTA0,AV,AK,AS,PONNIV,PONFE,DEFBET,QR,SIG,FP       
  //   1522	C=======================================================================
  //   1523	C                                                                       
  //   1524	C                                                                       
  //   1525	C-----------------------------------------------------------------------
  //   1526	C     A                 MASS NUMBER OF THE DAUGHTER NUCLEUS             
  //   1527	C     EE                EXCITATION ENERGY OF THE MOTHER NUCLEUS         
  //   1528	C     ESOUS             SEPARATION ENERGY PLUS EFFECTIVE COULOMB BARRIER
  //   1529	C     DENS              STATE DENSITY OF DAUGHTER NUCLEUS AT EE-ESOUS-EC
  //   1530	C     BSHELL            SHELL CORRECTION                                
  //   1531	C     TEMP              NUCLEAR TEMPERATURE                             
  //   1532	C     E        LOCAL    EXCITATION ENERGY OF THE DAUGHTER NUCLEUS       
  //   1533	C     E1       LOCAL    HELP VARIABLE                                   
  //   1534	C     Y0,Y1,Y2,Y01,Y11,Y21                                              
  //   1535	C              LOCAL    HELP VARIABLES                                  
  //   1536	C     PA       LOCAL    STATE-DENSITY PARAMETER                         
  //   1537	C     EC                KINETIC ENERGY OF EMITTED PARTICLE WITHOUT      
  //   1538	C                        COULOMB REPULSION                              
  //   1539	C     IDEN              FAKTOR FOR SUBSTRACTING KINETIC ENERGY IDEN*TEMP
  //   1540	C     DELTA0            PAIRING GAP 12 FOR GROUND STATE                 
  //   1541	C                       14 FOR SADDLE POINT                             
  //   1542	C     EITERA            HELP VARIABLE FOR TEMPERATURE ITERATION         
  //   1543	C-----------------------------------------------------------------------
  //   1544	C                                                                       
  //   1545	C                                                                       
  G4double afp;
  G4double delta0;
  G4double deltau;
  G4double deltpp;
  G4double e;
  G4double ecor;
  G4double ecor1;
  G4double ecr;
  G4double er;
  G4double fe;
  G4double fp;
  G4double he;
  G4double iz;
  G4double pa;
  G4double para;
  G4double parz;
  G4double ponfe;
  G4double ponniv;
  G4double qr;
  G4double sig;
  G4double y01;
  G4double y11;
  G4double y2;
  G4double y21;
  G4double y1;
  G4double y0;

  G4double pi6 = pow(3.1415926535,2) / 6.0;
  ecr=10.0;
  er=28.0;
  afp=idnint(a);
  iz=idnint(z);

  // level density parameter                                               
  if((ald->optafan == 1)) {
    pa = (ald->av)*a + (ald->as)*pow(a,(2.e0/3.e0)) + (ald->ak)*pow(a,(1.e0/3.e0));
  }
  else {
    pa = (ald->av)*a + (ald->as)*bs*pow(a,(2.e0/3.e0)) + (ald->ak)*bk*pow(a,(1.e0/3.e0));
  }

  fp = 0.01377937231e0 * pow(a,(5.e0/3.e0)) * (1.e0 + defbet/3.e0);

  // pairing corrections                                                   
  if (bs > 1.0) {
    delta0 = 14;
  }
  else {
    delta0 = 12;
  }

  if (esous > 1.0e30) {
    (*dens) = 0.0;
    (*temp) = 0.0;
    goto densniv100;                                                       
  }
  e = ee - esous;

  if (e < 0.0) {
    (*dens) = 0.0;
    (*temp) = 0.0;
    goto densniv100;
  }

  // shell corrections                                                     
  if (optshp > 0) {
    deltau = bshell;
    if (optshp == 2) {
      deltau = 0.0;
    }
    if (optshp >= 2) {
      // pairing energy shift with condensation energy a.r.j. 10.03.97        
      deltpp = -0.25e0* (delta0/pow(sqrt(a),2)) * pa /pi6 + 2.e0*delta0/sqrt(a);

      parite(a,&para);
      if (para < 0.0) {
	e = e - delta0/sqrt(a);
      } else {                                                         
	parite(z,&parz);
	if (parz > 0.e0) {
	  e = e - 2.0*delta0/sqrt(a);
	} else {
	  e = e;
	}
      }
    } else {                                                          
      deltpp = 0.0;
    }
  } else {
    deltau = 0.0;
    deltpp = 0.0;
  }
  if (e < 0.0) {
    e = 0.0;
    (*temp) = 0.0;
  }

  // washing out is made stronger ! g.k. 3.7.96                           
  ponfe = -2.5*pa*e*pow(a,(-4.0/3.0));

  if (ponfe < -700.0)  {
    ponfe = -700.0;
  }
  fe = 1.0 - exp(ponfe);
  if (e < ecr) {
    // priv. comm. k.-h. schmidt                                         
    he = 1.0 - pow((1.0 - e/ecr),2);
  }
  else {
    he = 1.0;
  }

  // Excitation energy corrected for pairing and shell effects             
  // washing out with excitation energy is included.                        
  ecor = e + deltau*fe + deltpp*he;

  if (ecor <= 0.1) {
    ecor = 0.1;
  }

  // statt 170.d0 a.r.j. 8.11.97                                           

  // iterative procedure according to grossjean and feldmeier              
  // to avoid the singularity e = 0                                        
  if (ee < 5.0) {
    y1 = sqrt(pa*ecor);
    for(int j = 0; j < 5; j++) {
      y2 = pa*ecor*(1.e0-exp(-y1));
      y1 = sqrt(y2);
    }
    
    y0 = pa/y1;
    (*temp)=1.0/y0;
    (*dens) = exp(y0*ecor)/ (pow((pow(ecor,3)*y0),0.5)*pow((1.0-0.5*y0*ecor*exp(-y1)),0.5))* exp(y1)*(1.0-exp(-y1))*0.1477045;
    if (ecor < 1.0) {
      ecor1=1.0;
      y11 = sqrt(pa*ecor1);

      for(int j = 0; j < 7; j++) {
	y21 = pa*ecor1*(1.0-exp(-y11));
	y11 = sqrt(y21);
      }

      y01 = pa/y11;
      (*dens) = (*dens)*pow((y01/y0),1.5);
      (*temp) = (*temp)*pow((y01/y0),1.5);
    }
  }
  else {
    ponniv = 2.0*sqrt(pa*ecor);
    if (ponniv > 700.0) {
      ponniv = 700.0;
    }

    // fermi gas state density                                               
    (*dens) = pow(pa,(-0.25e0))*pow(ecor,(-1.25e0))*exp(ponniv) * 0.1477045e0;
    (*temp) = sqrt(ecor/pa);
  }
 densniv100:

  // spin cutoff parameter                                                 
  sig = fp * (*temp);

  // collective enhancement                                                
  if (optcol == 1) {
    qrot(z,a,defbet,sig,ecor,&qr);
  }
  else {
    qr   = 1.0;
  }

  (*dens) = (*dens) * qr;
}


G4double G4Abla::bfms67(G4double zms, G4double ams)
{
  // This subroutine calculates the fission barriers                                                                  
  // of the liquid-drop model of Myers and Swiatecki (1967).                                                                 
  // Analytic parameterization of Dahlinger 1982 
  // replaces tables. Barrier heights from Myers and Swiatecki !!!                                                                 

  G4double nms,ims,ksims,xms, ums;

  nms = ams - zms;
  ims = (nms-zms)/ams;
  ksims= 50.15e0 * (1.- 1.78 * pow(ims,2));
  xms = pow(zms,2) / (ams * ksims);
  ums = 0.368e0-5.057e0*xms+8.93e0*pow(xms,2)-8.71*pow(xms,3);
  return(0.7322e0*pow(zms,2)/pow(ams,(0.333333e0))*pow(10.e0,ums));
}

void G4Abla::lpoly(G4double x, G4int n, G4double pl[])
{
  // THIS SUBROUTINE CALCULATES THE ORDINARY LEGENDRE POLYNOMIALS OF   
  // ORDER 0 TO N-1 OF ARGUMENT X AND STORES THEM IN THE VECTOR PL.    
  // THEY ARE CALCULATED BY RECURSION RELATION FROM THE FIRST TWO      
  // POLYNOMIALS.                                                      
  // WRITTEN BY A.J.SIERK  LANL  T-9  FEBRUARY, 1984                   

  // NOTE: PL AND X MUST BE DOUBLE PRECISION ON 32-BIT COMPUTERS!      

  pl[0] = 1.0;
  pl[1] = x;

  for(int i = 2; i < n; i++) {
    pl[i] = ((2*i - 3)*x*pl[i-1] - (i - 2)*pl[i-2])/(i-1);
  }
}

G4double G4Abla::eflmac(G4int ia, G4int iz, G4int flag, G4int optshp)
{
  // CHANGED TO CALCULATE TOTAL BINDING ENERGY INSTEAD OF MASS EXCESS.     
  // SWITCH FOR PAIRING INCLUDED AS WELL.                                  
  // BINDING = EFLMAC(IA,IZ,0,OPTSHP)                                      
  // FORTRAN TRANSCRIPT OF /U/GREWE/LANG/EEX/FRLDM.C                       
  // A.J. 15.07.96                                                         

  // this function will calculate the liquid-drop nuclear mass for spheri
  // configuration according to the preprint NUCLEAR GROUND-STATE        
  // MASSES and DEFORMATIONS by P. M"oller et al. from August 16, 1993 p.
  // All constants are taken from this publication for consistency.      

  // Parameters:                                                         
  // a:    nuclear mass number                                         
  // z:    nuclear charge                                              
  // flag:     0       - return mass excess                            
  //       otherwise   - return pairing (= -1/2 dpn + 1/2 (Dp + Dn))   

  G4double eflmacResult;

  G4int in;
  G4double z,n,a,av,as,a0,c1,c4,b1,b3,f,ca,w,dp,dn,dpn,efl,pi;
  G4double rmac,bs,h,r0,kf,ks,kv,rp,ay,aden,x0,y0,mh,mn,esq,ael,i;
  pi = 3.141592653589793238e0;

  // fundamental constants
  // hydrogen-atom mass excess
  mh  = 7.289034;

  // neutron mass excess
  mn  = 8.071431;

  // electronic charge squared
  esq = 1.4399764;

  // constants from considerations other than nucl. masses
  // electronic binding
  ael = 1.433e-5;

  // proton rms radius
  rp  = 0.8;

  // nuclear radius constant
  r0  = 1.16;

  // range of yukawa-plus-expon. potential
  ay  = 0.68;

  // range of yukawa function used to generate                          
  // nuclear charge distribution
  aden= 0.70;

  // constants from considering odd-even mass differences
  // average pairing gap
  rmac= 4.80;

  // neutron-proton interaction
  h   = 6.6;

  // wigner constant
  w   = 30.0;

  // adjusted parameters
  // volume energy
  av  = 16.00126;

  // volume asymmetry
  kv  =  1.92240;

  // surface energy
  as  = 21.18466;

  // surface asymmetry
  ks  =  2.345;
  // a^0 constant
  a0  =  2.615;

  // charge asymmetry
  ca  =  0.10289;

  // we will account for deformation by using the microscopic           
  // corrections tabulated from p. 68ff */                               
  bs = 1.0;

  z   = double(iz);
  a   = double(ia);
  in  = ia - iz;                                                       
  n   = double(in);
  dn  = rmac*bs/pow(n,(1.0/3.0));
  dp  = rmac*bs/pow(z,(1.0/3.0));
  dpn = h/bs/pow(a,(2.0/3.0));

  c1  = 3.0/5.0*esq/r0;
  c4  = 5.0/4.0*pow((3.0/(2.0*pi)),(2.0/3.0)) * c1;

  kf  = pow((9.0*pi*z/(4.0*a)),(1.0/3.0))/r0;
  f = -1.0/8.0*rp*rp*esq/pow(r0,3) * (145.0/48.0 - 327.0/2880.0*pow(kf,2) * pow(rp,2) + 1527.0/1209600.0*pow(kf,4) * pow(rp,4));
  i   = (n-z)/a;

  x0  = r0 * pow(a,(1.0/3.0)) / ay;
  y0  = r0 * pow(a,(1.0/3.0)) / aden;

  b1  = 1.0 - 3.0/(pow(x0,2)) + (1.0 + x0) * (2.0 + 3.0/x0 + 3.0/pow(x0,2)) * exp(-2.0*x0);

  b3  = 1.0 - 5.0/pow(y0,2) * (1.0 - 15.0/(8.0*y0) + 21.0/(8.0 * pow(y0,3))
			       - 3.0/4.0 * (1.0 + 9.0/(2.0*y0) + 7.0/pow(y0,2)
					    + 7.0/(2.0 * pow(y0,3))) * exp(-2.0*y0));

  // now calulation of total binding energy a.j. 16.7.96                   

  efl = -1* av*(1.0 - kv*i*i)*a + as*(1.0 - ks*i*i)*b1 * pow(a,(2.0/3.0)) + a0
    + c1*z*z*b3/pow(a,(1.0/3.0)) - c4*pow(z,(4.0/3.0))/pow(a,(1.e0/3.e0))
    + f*pow(z,2)/a -ca*(n-z) - ael * pow(z,(2.39e0));

  if ((in == iz) && (mod(in,2) == 1) && (mod(iz,2) == 1)) {
    // n and z odd and equal
    efl = efl + w*(utilabs(i)+1.e0/a);
  }
  else {
    efl= efl + w* utilabs(i);
  }

  // pairing is made optional                                              
  if (optshp >= 2) {
    // average pairing
    if ((mod(in,2) == 1) && (mod(iz,2) == 1)) {
      efl = efl - dpn;
    }
    if (mod(in,2) == 1) {
      efl = efl + dn;
    }
    if (mod(iz,2) == 1) {
      efl    = efl + dp;
    }
    // end if for pairing term                                               
  }

  if (flag != 0) {
    eflmacResult =  (0.5*(dn + dp) - 0.5*dpn);
  }
  else {
    eflmacResult = efl;
  }

  return eflmacResult;
}

void G4Abla::appariem(G4double a, G4double z, G4double *del)
{
  // CALCUL DE LA CORRECTION, DUE A L'APPARIEMENT, DE L'ENERGIE DE     
  // LIAISON D'UN NOYAU                                                
  // PROCEDURE FOR CALCULATING THE PAIRING CORRECTION TO THE BINDING   
  // ENERGY OF A SPECIFIC NUCLEUS                                      

  double para,parz;
  // A                 MASS NUMBER                                     
  // Z                 NUCLEAR CHARGE                                  
  // PARA              HELP VARIABLE FOR PARITY OF A                   
  // PARZ              HELP VARIABLE FOR PARITY OF Z                   
  // DEL               PAIRING CORRECTION                              

  parite(a, &para);

  if (para < 0.0) {
    (*del) = 0.0;
  }
  else {
    parite(z, &parz);
    if (parz > 0.0) {
      (*del) = -12.0/sqrt(a);
    }
    else {
      (*del) = 12.0/sqrt(a);
    }
  }
}

void G4Abla::parite(G4double n, G4double *par)
{
  // CALCUL DE LA PARITE DU NOMBRE N                                   
  //
  // PROCEDURE FOR CALCULATING THE PARITY OF THE NUMBER N.             
  // RETURNS -1 IF N IS ODD AND +1 IF N IS EVEN                        

  G4double n1, n2, n3;

  // N                 NUMBER TO BE TESTED                             
  // N1,N2             HELP VARIABLES                                  
  // PAR               HELP VARIABLE FOR PARITY OF N                   

  n3 = double(idnint(n));
  n1 = n3/2.0;
  n2 = n1 - dint(n1);

  if (n2 > 0.0) {
    (*par) = -1.0;
  }
  else {
    (*par) = 1.0;
  }
}

G4double G4Abla::tau(G4double bet, G4double homega, G4double ef, G4double t)
{
  // INPUT : BET, HOMEGA, EF, T                                          
  // OUTPUT: TAU - RISE TIME IN WHICH THE FISSION WIDTH HAS REACHED      
  //               90 PERCENT OF ITS FINAL VALUE                               
  // 
  // BETA   - NUCLEAR VISCOSITY                                          
  // HOMEGA - CURVATURE OF POTENTIAL                                     
  // EF     - FISSION BARRIER                                            
  // T      - NUCLEAR TEMPERATURE                                        

  G4double tauResult;

  G4double tlim,tau1,tau2;
  tlim = 8.e0 * ef;
  if (t > tlim) {
    t = tlim;
  }

  // modified bj and khs 6.1.2000:
  if (bet/(sqrt(2.0)*10.0*(homega/6.582122)) <= 1.0) {
    tauResult = log(10.0*ef/t)/(bet*1.0e+21);
  }
  else {
    tauResult = log(10.0*ef/t)/ (2.0*pow((10.0*homega/6.582122),2))*(bet*1.0e-21);
  } //end if                                                            

  return tauResult;
}

G4double G4Abla::cram(G4double bet, G4double homega)
{
  // INPUT : BET, HOMEGA  NUCLEAR VISCOSITY + CURVATURE OF POTENTIAL      
  // OUTPUT: KRAMERS FAKTOR  - REDUCTION OF THE FISSION PROBABILITY       
  //                           INDEPENDENT OF EXCITATION ENERGY                             

  G4double rel = bet/(20.0*homega/6.582122);
  G4double cramResult = sqrt(1.0 + pow(rel,2)) - rel;
  // limitation introduced   6.1.2000  by  khs

  if (cramResult > 1.0) {
    cramResult = 1.0;
  }

  return cramResult;
}

G4double G4Abla::bipol(int iflag, G4double y)
{
  // CALCULATION OF THE SURFACE BS OR CURVATURE BK OF A NUCLEUS        
  // RELATIVE TO THE SPHERICAL CONFIGURATION                           
  // BASED ON  MYERS, DROPLET MODEL FOR ARBITRARY SHAPES               

  // INPUT: IFLAG - 0/1 BK/BS CALCULATION                              
  //         Y    - (1 - X) COMPLEMENT OF THE FISSILITY                

  // LINEAR INTERPOLATION OF BS BK TABLE                               

  int i;

  G4double bipolResult;

  const int bsbkSize = 52;

  G4double bk[bsbkSize] = {1.00000,1.00087,1.00352,1.00799,1.01433,1.02265,1.03306,
			   1.04576,1.06099,1.07910,1.10056,1.12603,1.15651,1.19348,
			   1.23915,1.29590,1.35951,1.41013,1.44103,1.46026,1.47339,
			   1.48308,1.49068,1.49692,1.50226,1.50694,1.51114,1.51502,
			   1.51864,1.52208,1.52539,1.52861,1.53177,1.53490,1.53803,
			   1.54117,1.54473,1.54762,1.55096,1.55440,1.55798,1.56173,
			   1.56567,1.56980,1.57413,1.57860,1.58301,1.58688,1.58688,
			   1.58688,1.58740,1.58740};

  G4double bs[bsbkSize] = {1.00000,1.00086,1.00338,1.00750,1.01319,
			   1.02044,1.02927,1.03974,
			   1.05195,1.06604,1.08224,1.10085,1.12229,1.14717,1.17623,1.20963,
			   1.24296,1.26532,1.27619,1.28126,1.28362,1.28458,1.28477,1.28450,
			   1.28394,1.28320,1.28235,1.28141,1.28042,1.27941,1.27837,1.27732,
			   1.27627,1.27522,1.27418,1.27314,1.27210,1.27108,1.27006,1.26906,
			   1.26806,1.26707,1.26610,1.26514,1.26418,1.26325,1.26233,1.26147,
			   1.26147,1.26147,1.25992,1.25992};

  i = idint(y/2.0e-02) + 1;

  if (iflag == 1) {
    bipolResult = bs[i] + (bs[i+1] - bs[i])/2.0e-02 * (y - 2.0e-02*(i - 1));
  }
  else {
    bipolResult = bk[i] + (bk[i+1] - bk[i])/2.e-02 * (y - 2.0e-02*(i - 1));
  }

  return bipolResult;
}

void G4Abla::barfit(G4int iz, G4int ia, G4int il, G4double *sbfis, G4double *segs, G4double *selmax)
{
  //   2223	C     VERSION FOR 32BIT COMPUTER                                        
  //   2224	C     THIS SUBROUTINE RETURNS THE BARRIER HEIGHT BFIS, THE              
  //   2225	C     GROUND-STATE ENERGY SEGS, IN MEV, AND THE ANGULAR MOMENTUM        
  //   2226	C     AT WHICH THE FISSION BARRIER DISAPPEARS, LMAX, IN UNITS OF        
  //   2227	C     H-BAR, WHEN CALLED WITH INTEGER AGUMENTS IZ, THE ATOMIC           
  //   2228	C     NUMBER, IA, THE ATOMIC MASS NUMBER, AND IL, THE ANGULAR           
  //   2229	C     MOMENTUM IN UNITS OF H-BAR. (PLANCK'S CONSTANT DIVIDED BY         
  //   2230	C     2*PI).                                                            
  //   2231	C                                                                       
  //   2232	C        THE FISSION BARRIER FO IL = 0 IS CALCULATED FROM A 7TH         
  //   2233	C     ORDER FIT IN TWO VARIABLES TO 638 CALCULATED FISSION              
  //   2234	C     BARRIERS FOR Z VALUES FROM 20 TO 110. THESE 638 BARRIERS ARE      
  //   2235	C     FIT WITH AN RMS DEVIATION OF 0.10 MEV BY THIS 49-PARAMETER        
  //   2236	C     FUNCTION.                                                         
  //   2237	C     IF BARFIT IS CALLED WITH (IZ,IA) VALUES OUTSIDE THE RANGE OF      
  //   2238	C     THE BARRIER HEIGHT IS SET TO 0.0, AND A MESSAGE IS PRINTED        
  //   2239	C     ON THE DEFAULT OUTPUT FILE.                                       
  //   2240	C                                                                       
  //   2241	C        FOR IL VALUES NOT EQUAL TO ZERO, THE VALUES OF L AT WHICH      
  //   2242	C     THE BARRIER IS 80% AND 20% OF THE L=0 VALUE ARE RESPECTIVELY      
  //   2243	C     FIT TO 20-PARAMETER FUNCTIONS OF Z AND A, OVER A MORE             
  //   2244	C     RESTRICTED RANGE OF A VALUES, THAN IS THE CASE FOR L = 0.         
  //   2245	C     THE VALUE OF L WHERE THE BARRIER DISAPPEARS, LMAX IS FIT TO       
  //   2246	C     A 24-PARAMETER FUNCTION OF Z AND A, WITH THE SAME RANGE OF        
  //   2247	C     Z AND A VALUES AS L-80 AND L-20.                                  
  //   2248	C        ONCE AGAIN, IF AN (IZ,IA) PAIR IS OUTSIDE OF THE RANGE OF      
  //   2249	C     VALIDITY OF THE FIT, THE BARRIER VALUE IS SET TO 0.0 AND A        
  //   2250	C     MESSAGE IS PRINTED. THESE THREE VALUES (BFIS(L=0),L-80, AND       
  //   2251	C     L-20) AND THE CONSTRINTS OF BFIS = 0 AND D(BFIS)/DL = 0 AT        
  //   2252	C     L = LMAX AND L=0 LEAD TO A FIFTH-ORDER FIT TO BFIS(L) FOR         
  //   2253	C     L>L-20. THE FIRST THREE CONSTRAINTS LEAD TO A THIRD-ORDER FIT     
  //   2254	C     FOR THE REGION L < L-20.                                          
  //   2255	C                                                                       
  //   2256	C        THE GROUND STATE ENERGIES ARE CALCULATED FROM A                
  //   2257	C     120-PARAMETER FIT IN Z, A, AND L TO 214 GROUND-STATE ENERGIES     
  //   2258	C     FOR 36 DIFFERENT Z AND A VALUES.                                  
  //   2259	C     (THE RANGE OF Z AND A IS THE SAME AS FOR L-80, L-20, AND          
  //   2260	C     L-MAX)                                                            
  //   2261	C                                                                       
  //   2262	C        THE CALCULATED BARRIERS FROM WHICH THE FITS WERE MADE WERE     
  //   2263	C     CALCULATED IN 1983-1984 BY A. J. SIERK OF LOS ALAMOS              
  //   2264	C     NATIONAL LABORATORY GROUP T-9, USING YUKAWA-PLUS-EXPONENTIAL      
  //   2265	C     G4DOUBLE FOLDED NUCLEAR ENERGY, EXACT COULOMB DIFFUSENESS           
  //   2266	C     CORRECTIONS, AND DIFFUSE-MATTER MOMENTS OF INERTIA.               
  //   2267	C     THE PARAMETERS OF THE MODEL R-0 = 1.16 FM, AS 21.13 MEV,          
  //   2268	C     KAPPA-S = 2.3, A = 0.68 FM.                                       
  //   2269	C     THE DIFFUSENESS OF THE MATTER AND CHARGE DISTRIBUTIONS USED       
  //   2270	C     CORRESPONDS TO A SURFACE DIFFUSENESS PARAMETER (DEFINED BY        
  //   2271	C     MYERS) OF 0.99 FM. THE CALCULATED BARRIERS FOR L = 0 ARE          
  //   2272	C     ACCURATE TO A LITTLE LESS THAN 0.1 MEV; THE OUTPUT FROM           
  //   2273	C     THIS SUBROUTINE IS A LITTLE LESS ACCURATE. WORST ERRORS MAY BE    
  //   2274	C     AS LARGE AS 0.5 MEV; CHARACTERISTIC UNCERTAINY IS IN THE RANGE    
  //   2275	C     OF 0.1-0.2 MEV. THE RMS DEVIATION OF THE GROUND-STATE FIT         
  //   2276	C     FROM THE 214 INPUT VALUES IS 0.20 MEV. THE MAXIMUM ERROR          
  //   2277	C     OCCURS FOR LIGHT NUCLEI IN THE REGION WHERE THE GROUND STATE      
  //   2278	C     IS PROLATE, AND MAY BE GREATER THAN 1.0 MEV FOR VERY NEUTRON      
  //   2279	C     DEFICIENT NUCLEI, WITH L NEAR LMAX. FOR MOST NUCLEI LIKELY TO     
  //   2280	C     BE ENCOUNTERED IN REAL EXPERIMENTS, THE MAXIMUM ERROR IS          
  //   2281	C     CLOSER TO 0.5 MEV, AGAIN FOR LIGHT NUCLEI AND L NEAR LMAX.        
  //   2282	C                                                                       
  //   2283	C     WRITTEN BY A. J. SIERK, LANL T-9                                  
  //   2284	C     VERSION 1.0 FEBRUARY, 1984                                        
  //   2285	C                                                                       
  //   2286	C     THE FOLLOWING IS NECESSARY FOR 32-BIT MACHINES LIKE DEC VAX,      
  //   2287	C     IBM, ETC                                                          

  G4double pa[7],pz[7],pl[10];
  G4double a,z,amin,amax,amin2,amax2,aa,zz,bfis;
  G4double bfis0,ell,el,egs,el80,el20,elmax,sel80,sel20,x,y,q,qa,qb;
  G4double aj,ak,a1,a2;

  G4int i,j,k,m;
  G4int l;

  G4double emncof[4][5] = {{-9.01100e+2,-1.40818e+3, 2.77000e+3,-7.06695e+2, 8.89867e+2}, 
			   {1.35355e+4,-2.03847e+4, 1.09384e+4,-4.86297e+3,-6.18603e+2},
			   {-3.26367e+3, 1.62447e+3, 1.36856e+3, 1.31731e+3, 1.53372e+2},
			   {7.48863e+3,-1.21581e+4, 5.50281e+3,-1.33630e+3, 5.05367e-2}};

  G4double elmcof[4][5] = {{1.84542e+3,-5.64002e+3, 5.66730e+3,-3.15150e+3, 9.54160e+2},
			   {-2.24577e+3, 8.56133e+3,-9.67348e+3, 5.81744e+3,-1.86997e+3},
			   {2.79772e+3,-8.73073e+3, 9.19706e+3,-4.91900e+3, 1.37283e+3},
			   {-3.01866e+1, 1.41161e+3,-2.85919e+3, 2.13016e+3,-6.49072e+2}};

  G4double emxcof[4][6] = {{19.43596e4,-2.241997e5,2.223237e5,-1.324408e5,4.68922e4,-8.83568e3},
			   {-1.655827e5,4.062365e5,-4.236128e5,2.66837e5,-9.93242e4,1.90644e4},
			   {1.705447e5,-4.032e5,3.970312e5,-2.313704e5,7.81147e4,-1.322775e4},
			   {-9.274555e4,2.278093e5,-2.422225e5,1.55431e5,-5.78742e4,9.97505e3}};

  G4double elzcof[7][7] = {{5.11819909e+5,-1.30303186e+6, 1.90119870e+6,-1.20628242e+6, 5.68208488e+5, 5.48346483e+4,-2.45883052e+4},
			   {-1.13269453e+6, 2.97764590e+6,-4.54326326e+6, 3.00464870e+6, -1.44989274e+6,-1.02026610e+5, 6.27959815e+4},
			   {1.37543304e+6,-3.65808988e+6, 5.47798999e+6,-3.78109283e+6, 1.84131765e+6, 1.53669695e+4,-6.96817834e+4},
			   {-8.56559835e+5, 2.48872266e+6,-4.07349128e+6, 3.12835899e+6, -1.62394090e+6, 1.19797378e+5, 4.25737058e+4},
			   {3.28723311e+5,-1.09892175e+6, 2.03997269e+6,-1.77185718e+6, 9.96051545e+5,-1.53305699e+5,-1.12982954e+4},
			   {4.15850238e+4, 7.29653408e+4,-4.93776346e+5, 6.01254680e+5, -4.01308292e+5, 9.65968391e+4,-3.49596027e+3},
			   {-1.82751044e+5, 3.91386300e+5,-3.03639248e+5, 1.15782417e+5, -4.24399280e+3,-6.11477247e+3, 3.66982647e+2}};

  const G4int sizex = 5;
  const G4int sizey = 6;
  const G4int sizez = 4;

  G4double egscof[sizey][sizey][sizez];

  G4double egs1[sizey][sizex] = {{1.927813e5, 7.666859e5, 6.628436e5, 1.586504e5,-7.786476e3},
				 {-4.499687e5,-1.784644e6,-1.546968e6,-4.020658e5,-3.929522e3},
				 {4.667741e5, 1.849838e6, 1.641313e6, 5.229787e5, 5.928137e4},
				 {-3.017927e5,-1.206483e6,-1.124685e6,-4.478641e5,-8.682323e4},
				 {1.226517e5, 5.015667e5, 5.032605e5, 2.404477e5, 5.603301e4},
				 {-1.752824e4,-7.411621e4,-7.989019e4,-4.175486e4,-1.024194e4}};

  G4double egs2[sizey][sizex] = {{-6.459162e5,-2.903581e6,-3.048551e6,-1.004411e6,-6.558220e4},
				 {1.469853e6, 6.564615e6, 6.843078e6, 2.280839e6, 1.802023e5},
				 {-1.435116e6,-6.322470e6,-6.531834e6,-2.298744e6,-2.639612e5},
				 {8.665296e5, 3.769159e6, 3.899685e6, 1.520520e6, 2.498728e5},      
				 {-3.302885e5,-1.429313e6,-1.512075e6,-6.744828e5,-1.398771e5},
				 {4.958167e4, 2.178202e5, 2.400617e5, 1.167815e5, 2.663901e4}};

  G4double egs3[sizey][sizex] = {{3.117030e5, 1.195474e6, 9.036289e5, 6.876190e4,-6.814556e4},
				 {-7.394913e5,-2.826468e6,-2.152757e6,-2.459553e5, 1.101414e5},
				 {7.918994e5, 3.030439e6, 2.412611e6, 5.228065e5, 8.542465e3},
				 {-5.421004e5,-2.102672e6,-1.813959e6,-6.251700e5,-1.184348e5},
				 {2.370771e5, 9.459043e5, 9.026235e5, 4.116799e5, 1.001348e5},
				 {-4.227664e4,-1.738756e5,-1.795906e5,-9.292141e4,-2.397528e4}};

  G4double egs4[sizey][sizex] = {{-1.072763e5,-5.973532e5,-6.151814e5, 7.371898e4, 1.255490e5},
				 {2.298769e5, 1.265001e6, 1.252798e6,-2.306276e5,-2.845824e5},
				 {-2.093664e5,-1.100874e6,-1.009313e6, 2.705945e5, 2.506562e5},
				 {1.274613e5, 6.190307e5, 5.262822e5,-1.336039e5,-1.115865e5},
				 {-5.715764e4,-2.560989e5,-2.228781e5,-3.222789e3, 1.575670e4},
				 {1.189447e4, 5.161815e4, 4.870290e4, 1.266808e4, 2.069603e3}};

  for(G4int i = 0; i < sizey; i++) {
    for(G4int j = 0; j < sizex; j++) {
      egscof[i][j][0] = egs1[i][j];
      egscof[i][j][1] = egs2[i][j];
      egscof[i][j][2] = egs3[i][j];
      egscof[i][j][3] = egs4[i][j];
    }
  }

  // the program starts here                                           
  if (iz < 19  ||  iz > 111) {
    goto barfit900;
  }

  if(iz > 102   &&  il > 0) {
    goto barfit902;
  }

  z=double(iz);
  a=double(ia);
  el=double(il);
  amin= 1.2e0*z + 0.01e0*z*z;
  amax= 5.8e0*z - 0.024e0*z*z;

  if(a  <  amin  ||  a  >  amax) {
    goto barfit910;
  }

  // angul.mom.zero barrier                 
  aa=2.5e-3*a;
  zz=1.0e-2*z;
  ell=1.0e-2*el;
  bfis0=0.e0;
  lpoly(zz,7,pz);
  lpoly(aa,7,pa);

  for(i = 0; i < 7; i++) { //do 10 i=1,7                                                       
    for(j = 0; j < 7; j++) { //do 10 j=1,7                                                       
      bfis0=bfis0+elzcof[j][i]*pz[j]*pa[i];
    }
  }

  bfis=bfis0;
  (*sbfis)=bfis;
  egs=0.0;
  (*segs)=egs;

  // values of l at which the barrier        
  // is 20%(el20) and 80%(el80) of l=0 value    
  amin2 = 1.4e0*z + 0.009e0*z*z;
  amax2 = 20.e0 + 3.0e0*z;

  if((a < amin2-5.e0  ||  a > amax2+10.e0) &&  il > 0) {
    goto barfit920;
  }

  lpoly(zz,5,pz);
  lpoly(aa,4,pa);
  el80=0.0;
  el20=0.0;
  elmax=0.0;

  for(i = 0; i < 4; i++) {
    for(j = 0; j < 5; j++) {
      el80 = el80 + elmcof[j][i]*pz[j]*pa[i];
      el20 = el20 + emncof[j][i]*pz[j]*pa[i];
    }
  }

  sel80 = el80;
  sel20 = el20;

  // value of l (elmax) where barrier disapp.
  lpoly(zz,6,pz);
  lpoly(ell,9,pl);

  for(i = 0; i < 4; i++) { //do 30 i= 1,4                                                      
    for(j = 0; j < 6; j++) { //do 30 j=1,6                                                       
      elmax = elmax + emxcof[j][i]*pz[j]*pa[i];
    }
  }

  (*selmax)=elmax;

  // value of barrier at ang.mom.  l          
  if(il < 1){
    return;                                                
  }

  x = sel20/(*selmax);
  y = sel80/(*selmax);

  if(el <= sel20) {
    // low l              
    q=0.2e0/(pow(sel20,2)*pow(sel80,2)*(sel20-sel80));
    qa=q*(4.e0*pow(sel80,3) - pow(sel20,3));
    qb=-q*(4.e0*pow(sel80,2) - pow(sel20,2));
    bfis=bfis*(1.e0 + qa*pow(el,2) + qb*pow(el,3));
  }
  else {
    // high l             
    aj=(-20.e0*pow(x,5) + 25.e0*pow(x,4) - 4.e0)*pow((y-1.e0),2)*y*y;
    ak=(-20.e0*pow(y,5) + 25.e0*pow(y,4) - 1.e0) * pow((x-1.e0),2)*x*x;
    q= 0.2e0/(pow((y-x)*((1.e0-x)*(1.e0-y)*x*y),2));
    qa=q*(aj*y - ak*x);
    qb=-q*(aj*(2.e0*y+1.e0) - ak*(2.e0*x+1.e0));
    z=el/(*selmax);
    a1=4.e0*pow(z,5) - 5.e0*pow(z,4) + 1.e0;
    a2=qa*(2.e0*z+1.e0);
    bfis=bfis*(a1 + (z-1.e0)*(a2 + qb*z)*z*z*(z-1.e0));
  }
  if(bfis <= 0.0) {
    bfis=0.0;
  }

  if(el > (*selmax)) {
    bfis=0.0;
  }
  (*sbfis)=bfis;

  // now calculate rotating ground state energy                        
  if(el > (*selmax)) {
    return;                                           
  }

  for(k = 0; k < 4; k++) {
    for(l = 0; l < 6; l++) {
      for(m = 0; m < 5; m++) {
	egs = egs + egscof[m][l][k]*pz[l]*pa[k]*pl[2*m-1];
      }
    }
  }

  (*segs)=egs;
  if((*segs) < 0.e0) {
    (*segs)=0.0e0;
  }

  return;                                                            

 barfit900:  //continue                                                          
  (*sbfis)=0.0;
  // for z<19 sbfis set to 1.0e3                                            
  if (iz < 19)  {
    (*sbfis) = 1.0e3;
  }
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                            

 barfit902:
  (*sbfis)=0.0;
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                            

 barfit910:
  (*sbfis)=0.0;
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                           

 barfit920:
  (*sbfis)=0.0;
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                            
}

G4double G4Abla::expohaz(G4int k, G4double T)
{
  // TIRAGE ALEATOIRE DANS UNE EXPONENTIELLLE : Y=EXP(-X/T)

  return (-1*T*log(haz(k)));
}

G4double G4Abla::fd(G4double E)
{
  // DISTRIBUTION DE MAXWELL

  return (E*exp(-E));
}

G4double G4Abla::f(G4double E)
{
  // FONCTION INTEGRALE DE FD(E)

  return (1 - (E + 1) * exp(-E));
}

G4double G4Abla::fmaxhaz(G4double k, G4double T)
{
  // tirage aleatoire dans une maxwellienne
  // t : temperature
  //
  // declaration des variables
  //

  static G4double fmaxhazResult;

  //   2565	        dimension p(100)
  const int pSize = 100;
  static G4double p[pSize];
  //   2566	       dimension iy(19)
  const G4int iySize = 19;
  static G4double iy[iySize];
  // ial generateur pour le cascade (et les iy pour eviter les correlations)

  static G4int itest = 0;
  // programme principal

  if (itest == 1) {
    goto fmaxhaz120;
  }
  // calcul des p(i) par approximation de newton
  p[pSize-1] = 8.0;
  static G4double x = 0.1;
  static G4double x1;
  static G4double y;
  static int i;
  for(i = 0; i < 99; i++) { //do i=1,99
  fmaxhaz20:
    x1 = x - (f(x) - i/100.0)/fd(x);
    x = x1;
    if (fabs(f(x) - G4double(i)/100.0) < 1e-5) {
      goto fmaxhaz100;
    }
    goto fmaxhaz20;
  fmaxhaz100:
    p[i] = x;
  } //end do

  itest = 1;
  //   2584	c------ tirage aleatoire et calcul du x correspondant 
  //   2585	c       par regression lineaire
  //   2586	c120     y=haz(k)
 fmaxhaz120:
  //call ribm(y,iy(18))
  standardRandom(&y, &(hazard->igraine[17]));
  i = nint(y*100);

  //   2590	c ici on evite froidement les depassements de tableaux....(a.b. 3/9/99)        
  if(i == 0) {
    goto fmaxhaz120;
  }

  if (i == 1) {
    x = p[i]*y*100;
  }
  else {
    x = (p[i] - p[i-1])*(y*100 - i) + p[i];
  }

  fmaxhazResult=x*T;

  return fmaxhazResult;
}

G4double G4Abla::pace2(G4double a, G4double z)
{
  // PACE2
  // Cette fonction retourne le defaut de masse du noyau A,Z en MeV
  // Révisée pour a, z flottants 25/4/2002	                       =

  const G4int u = 931500;

  G4double pace2;

  G4int ii = idint(a+0.5);
  G4int jj = idint(z+0.5);

  if(ii <= 0 || jj < 0) {
    pace2=0.;
    return pace2;
  }

  if(jj > 300) {
    pace2=0.0;
  }
  else {
    pace2=pace->dm[ii][jj];
  }
  pace2=pace2/1000.;

  if(pace->dm[ii][jj] == 0.) {
    if(ii < 12) {
      pace2=-500.;
    }
    else {
      guet(&a, &z, &pace2);
      pace2=pace2-ii*931.5;
      pace2=pace2/1000.;
    }
  }

  return pace2;
}

void G4Abla::guet(G4double *x_par, G4double *z_par, G4double *find_par)
{
  // TABLE DE MASSES ET FORMULE DE MASSE TIRE DU PAPIER DE BRACK-GUET
  // Gives the theoritical value for mass excess...
  // Révisée pour x, z flottants 25/4/2002

  //real*8 x,z
  //	dimension q(0:50,0:70)
  G4double x = (*x_par);
  G4double z = (*z_par);
  G4double find = (*find_par);

  const G4int qrows = 50;
  const G4int qcols = 70;
  G4double q[qrows][qcols];

  G4int ix=G4int(floor(x+0.5));
  G4int iz=G4int(floor(z+0.5));
  G4double zz = iz;
  G4double xx = ix;
  find = 0.0;
  G4double avol = 15.776;
  G4double asur = -17.22;
  G4double ac = -10.24;
  G4double azer = 8.0;
  G4double xjj = -30.03;
  G4double qq = -35.4;
  G4double c1 = -0.737;
  G4double c2 = 1.28;

  if(ix <= 7) {
    q[0][1]=939.50;
    q[1][1]=938.21;
    q[1][2]=1876.1;
    q[1][3]=2809.39;
    q[2][4]=3728.34;
    q[2][3]=2809.4;
    q[2][5]=4668.8;
    q[2][6]=5606.5;
    q[3][5]=4669.1;
    q[3][6]=5602.9;
    q[3][7]=6535.27;
    q[4][6]=5607.3;
    q[4][7]=6536.1;
    q[5][7]=6548.3;
    find=q[iz][ix];
  }
  else {
    G4double xneu=xx-zz;
    G4double si=(xneu-zz)/xx;
    G4double x13=pow(xx,.333);
    G4double ee1=c1*zz*zz/x13;
    G4double ee2=c2*zz*zz/xx;
    G4double aux=1.+(9.*xjj/4./qq/x13);
    G4double ee3=xjj*xx*si*si/aux;
    G4double ee4=avol*xx+asur*(pow(xx,.666))+ac*x13+azer;
    G4double tota = ee1 + ee2 + ee3 + ee4;
    find = 939.55*xneu+938.77*zz - tota;
  }

  (*x_par) = x;
  (*z_par) = z;
  (*find_par) = find;
}


// Fission code

void G4Abla::even_odd(G4double r_origin,G4double r_even_odd,G4int &i_out)     
{
  // Procedure to calculate I_OUT from R_IN in a way that
  // on the average a flat distribution in R_IN results in a
  // fluctuating distribution in I_OUT with an even-odd effect as
  // given by R_EVEN_ODD

//     /* ------------------------------------------------------------ */
//     /* EXAMPLES :                                                   */
//     /* ------------------------------------------------------------ */
//     /*    If R_EVEN_ODD = 0 :                                       */
//     /*           CEIL(R_IN)  ----                                   */
//     /*                                                              */
//     /*              R_IN ->                                         */
//     /*            (somewhere in between CEIL(R_IN) and FLOOR(R_IN)) */                                            */
//     /*                                                              */
//     /*           FLOOR(R_IN) ----       --> I_OUT                   */
//     /* ------------------------------------------------------------ */
//     /*    If R_EVEN_ODD > 0 :                                       */
//     /*      The interval for the above treatment is                 */
//     /*         larger for FLOOR(R_IN) = even and                    */
//     /*         smaller for FLOOR(R_IN) = odd                        */
//     /*    For R_EVEN_ODD < 0 : just opposite treatment              */
//     /* ------------------------------------------------------------ */

//     /* ------------------------------------------------------------ */
//     /* On input:   R_ORIGIN    nuclear charge (real number)         */
//     /*             R_EVEN_ODD  requested even-odd effect            */
//     /* Intermediate quantity: R_IN = R_ORIGIN + 0.5                 */
//     /* On output:  I_OUT       nuclear charge (integer)             */
//     /* ------------------------------------------------------------ */

//      G4double R_ORIGIN,R_IN,R_EVEN_ODD,R_REST,R_HELP;
  G4double r_in,r_rest,r_help;
  G4double r_floor;
  G4double r_middle;
  //      G4int I_OUT,N_FLOOR;
  G4int n_floor;

  r_in = r_origin + 0.5;
  r_floor = (float)((int)(r_in));
  if (r_even_odd < 0.001) {
    i_out = (int)(r_floor);
  } 
  else {
    r_rest = r_in - r_floor;
    r_middle = r_floor + 0.5;
    n_floor = (int)(r_floor);
    if (n_floor%2 == 0) {
      // even before modif.
      r_help = r_middle + (r_rest - 0.5) * (1.0 - r_even_odd);
    } 
    else {
      // odd before modification
      r_help = r_middle + (r_rest - 0.5) * (1.0 + r_even_odd);
    }
    i_out = (int)(r_help);
  }
}

G4double G4Abla::umass(G4double z,G4double n,G4double beta)
{
  // liquid-drop mass, Myers & Swiatecki, Lysekil, 1967
  // pure liquid drop, without pairing and shell effects

  // On input:    Z     nuclear charge of nucleus
  //              N     number of neutrons in nucleus
  //              beta  deformation of nucleus
  // On output:   binding energy of nucleus

  G4double a,umass;
  G4double alpha;
  G4double xcom,xvs,xe;
  const G4double pi = 3.1416;
     
  a = n + z;
  alpha = ( sqrt(5.0/(4.0*pi)) ) * beta;
  xcom = 1.0 - 1.7826 * ((a - 2.0*z)/a)*((a - 2.0*z)/a);

  // factor for asymmetry dependence of surface and volume term
  xvs = - xcom * ( 15.4941 * a - 
		   17.9439 * pow(a,0.66667) * (1.0+0.4*alpha*alpha) );
  // sum of volume and surface energy
  xe = z*z * (0.7053/(pow(a,0.33333)) * (1.0-0.2*alpha*alpha) - 1.1529/a);
  umass = xvs + xe;
      
  return umass;
}

G4double G4Abla::ecoul(G4double z1,G4double n1,G4double beta1,G4double z2,G4double n2,G4double beta2,G4double d)
{
// Coulomb potential between two nuclei
// surfaces are in a distance of d
// in a tip to tip configuration

// approximate formulation
// On input: Z1      nuclear charge of first nucleus
//           N1      number of neutrons in first nucleus
//           beta1   deformation of first nucleus
//           Z2      nuclear charge of second nucleus
//           N2      number of neutrons in second nucleus
//           beta2   deformation of second nucleus
//           d       distance of surfaces of the nuclei

//      G4double Z1,N1,beta1,Z2,N2,beta2,d,ecoul;
      G4double ecoul;
      G4double dtot;
      const G4double r0 = 1.16;

      dtot = r0 * ( pow((z1+n1),0.33333) * (1.0+(2.0/3.0)*beta1)
                  + pow((z2+n2),0.33333) * (1.0+(2.0/3.0)*beta2) ) + d;
      ecoul = z1 * z2 * 1.44 / dtot;

      return ecoul;
}

void G4Abla::fissionDistri(G4double a,G4double z,G4double e,
			   G4double &a1,G4double &z1,G4double &e1,G4double &v1,
			   G4double &a2,G4double &z2,G4double &e2,G4double &v2)
{
//  On input: A, Z, E (mass, atomic number and exc. energy of compound nucleus
//                     before fission)
//  On output: Ai, Zi, Ei (mass, atomic number and exc. energy of fragment 1 and 2
//                     after fission)

//  Additionally calculated but not put in the parameter list:
//  Kinetic energy of prefragments EkinR1, EkinR2

//  Translation of SIMFIS18.PLI (KHS, 2.1.2001)

// This program calculates isotopic distributions of fission fragments
// with a semiempirical model                      
// Copy from SIMFIS3, KHS, 8. February 1995        
// Modifications made by Jose Benlliure and KHS in August 1996
// Energy counted from lowest barrier (J. Benlliure, KHS 1997)
// Some bugs corrected (J. Benlliure, KHS 1997)         
// Version used for thesis S. Steinhaueser (August 1997)
// (Curvature of LD potential increased by factor of 2!)

// Weiter veraendert mit der Absicht, eine Version zu erhalten, die
// derjenigen entspricht, die von J. Benlliure et al.
// in Nucl. Phys. A 628 (1998) 458 verwendet wurde,
// allerdings ohne volle Neutronenabdampfung.

// The excitation energy was calculate now for each fission channel
// separately. The dissipation from saddle to scission was taken from
// systematics, the deformation energy at scission considers the shell
// effects in a simplified way, and the fluctuation is included. 
// KHS, April 1999 

// The width in N/Z was carefully adapted to values given by Lang et al.

// The width and eventually a shift in N/Z (polarization) follows the
// following rules:                                        

// The line N/Z following UCD has an angle of atan(Zcn/Ncn)
// to the horizontal axis on a chart of nuclides.
// (For 238U the angle is 32.2 deg.)

// The following relations hold: (from Armbruster)
// 
// sigma(N) (A=const) = sigma(Z) (A=const)
// sigma(A) (N=const) = sigma(Z) (N=const)
// sigma(A) (Z=const) = sigma(N) (Z=const)
// 
// From this we get:
// sigma(Z) (N=const) * N = sigma(N) (Z=const) * Z
// sigma(A) (Z=const) = sigma(Z) (A=const) * A/Z
// sigma(N) (Z=const) = sigma(Z) (A=const) * A/Z
// Z*sigma(N) (Z=const) = N*sigma(Z) (N=const) = A*sigma(Z) (A=const)

// Excitation energy now calculated above the lowest potential point
// Inclusion of a distribution of excitation energies             

// Several modifications, starting from SIMFIS12: KHS November 2000 
// This version seems to work quite well for 238U.                  
// The transition from symmetric to asymmetric fission around 226Th 
// is reasonably well reproduced, although St. I is too strong and St. II
// is too weak. St. I and St. II are also weakly seen for 208Pb.

// Extensions for an event generator of fission events (21.11.2000,KHS)

// Defalt parameters (IPARS) rather carefully adjusted to
// pre-neutron mass distributions of Vives et al. (238U + n)
// Die Parameter Fgamma1 und Fgamma2 sind kleiner als die resultierenden
// Breiten der Massenverteilungen!!!  
// Fgamma1 und Fgamma2 wurden angepa�, so da�
// Sigma-A(ST-I) = 3.3, Sigma-A(St-II) = 5.8 (nach Vives) 

// Parameters of the model carefully adjusted by KHS (2.2.2001) to
// 238U + 208Pb, 1000 A MeV, Timo Enqvist et al.         


      G4double     n;
      G4double     nlight1,nlight2;
      G4double     aheavy1,alight1,aheavy2,alight2;
      G4double     eheavy1,elight1,eheavy2,elight2;
      G4double     zheavy1_shell,zheavy2_shell;
      G4double     zlight1,zlight2;
      G4double     masscurv;
      G4double     sasymm1,sasymm2,ssymm,ysum,yasymm;
      G4double     ssymm_mode1,ssymm_mode2;
      G4double     cz_asymm1_saddle,cz_asymm2_saddle;
// Curvature at saddle, modified by ld-potential
      G4double     wzasymm1_saddle, wzasymm2_saddle, wzsymm_saddle;
      G4double     wzasymm1_scission, wzasymm2_scission, wzsymm_scission;
      G4double     wzasymm1,wzasymm2,wzsymm;
      G4double     nlight1_eff, nlight2_eff;
      G4int  imode;
      G4double     rmode;
      G4double     z1mean, z2mean, z1width, za1width;
//      G4double     Z1,Z2,N1R,N2R,A1R,A2R,N1,N2,A1,A2;
      G4double     n1r,n2r,a1r,a2r,n1,n2;

      G4double     zsymm,nsymm,asymm;
      G4double     n1mean, n2mean, n1width;
      G4double     dueff;
// effective shell effect at lowest barrier
      G4double     eld;
// Excitation energy with respect to ld barrier
      G4double     re1,re2,re3;
      G4double     eps1,eps2;
      G4double     n1ucd,n2ucd,z1ucd,z2ucd;
      G4double     beta,beta1,beta2;

      G4double     dn1_pol;
// shift of most probable neutron number for given Z,
// according to polarization
      G4int  i_help;
      const G4double pi = 3.141593;

//   /* Parameters of the semiempirical fission model */
      G4double a_levdens;
//           /* level-density parameter */
      G4double a_levdens_light1,a_levdens_light2;
      G4double a_levdens_heavy1,a_levdens_heavy2;
      const G4double r_null = 1.16;
//          /* radius parameter */
      G4double epsilon_1_saddle,epsilon0_1_saddle;
      G4double epsilon_2_saddle,epsilon0_2_saddle,epsilon_symm_saddle;
      G4double epsilon_1_scission,epsilon0_1_scission;
      G4double epsilon_2_scission,epsilon0_2_scission;
      G4double epsilon_symm_scission;
//                                   /* modified energy */
      G4double e_eff1_saddle,e_eff2_saddle;
      G4double epot0_mode1_saddle,epot0_mode2_saddle,epot0_symm_saddle;
      G4double epot_mode1_saddle,epot_mode2_saddle,epot_symm_saddle;
      G4double e_defo,e_defo1,e_defo2,e_scission,e_asym;
      G4double e1exc,e2exc;
      G4double e1exc_sigma,e2exc_sigma;
      G4double e1final,e2final;

      const G4double r0 = 1.16;
      G4double tker;
      G4double ekin1,ekin2;
//      G4double EkinR1,EkinR2,E1,E2,V1,V2;
      G4double ekinr1,ekinr2;
      G4int icz,k;

//   Input parameters:
//OMMENT(Nuclear charge number);
//      G4double Z;
//OMMENT(Nuclear mass number);
//      G4double A;
//OMMENT(Excitation energy above fission barrier);
//      G4double E;

//   Model parameters:
//OMMENT(position of heavy peak valley 1);
      const G4double nheavy1 = 83.0;
//OMMENT(position of heavy peak valley 2);
      const G4double nheavy2 = 90.0;
//OMMENT(Shell effect for valley 1);
      const G4double delta_u1_shell = -2.65;
//        Parameter (Delta_U1_shell = -2)
//OMMENT(Shell effect for valley 2);
      const G4double delta_u2_shell = -3.8;
//        Parameter (Delta_U2_shell = -3.2)
//OMMENT(I: used shell effect);
      G4double delta_u1;
//omment(I: used shell effect);
      G4double delta_u2;
//OMMENT(Curvature of asymmetric valley 1);
      const G4double cz_asymm1_shell = 0.7;
//OMMENT(Curvature of asymmetric valley 2);
      const G4double cz_asymm2_shell = 0.15;
//OMMENT(Factor for width of distr. valley 1);
      const G4double fwidth_asymm1 = 0.63;
//OMMENT(Factor for width of distr. valley 2);
      const G4double fwidth_asymm2 = 0.97;
//       Parameter (CZ_asymm2_scission = 0.12)
//OMMENT(Parameter x: a = A/x);
      const G4double xlevdens = 12.0;
//OMMENT(Factor to gamma_heavy1);
      const G4double fgamma1 = 2.0;
//OMMENT(I: fading of shells (general));
      G4double gamma;
//OMMENT(I: fading of shell 1);
      G4double gamma_heavy1;
//OMMENT(I: fading of shell 2);
      G4double gamma_heavy2;
//OMMENT(Calculate A =1  or Aprime =0);
      const G4int i_eva = 0;
//      real*8 I_NUM_EVTS) INIT(1000000) COMMENT(Number of events to calculate);
//OMMENT(Zero-point energy at saddle);
      const G4double e_zero_point = 0.5;
//OMMENT(I: friction from saddle to scission);
      G4double e_saddle_scission;
//OMMENT(Friction factor);
      const G4double friction_factor = 1.0;
//OMMENT(I: Internal counter for different modes); INIT(0,0,0)
//      Integer*4 I_MODE(3)
//OMMENT(I: Yield of symmetric mode);
      G4double ysymm;
//OMMENT(I: Yield of asymmetric mode 1);
      G4double yasymm1;
//OMMENT(I: Yield of asymmetric mode 2);
      G4double yasymm2;
//OMMENT(I: Effective position of valley 1);
      G4double nheavy1_eff;
//OMMENT(I: position of heavy peak valley 1);
      G4double zheavy1;
//omment(I: Effective position of valley 2);
      G4double nheavy2_eff;
//OMMENT(I: position of heavy peak valley 2);
      G4double zheavy2;
//omment(I: Excitation energy above saddle 1);
      G4double eexc1_saddle;
//omment(I: Excitation energy above saddle 2);
      G4double eexc2_saddle;
//omment(I: Excitation energy above lowest saddle);
      G4double eexc_max;
//omment(I: Effective mass mode 1);
      G4double aheavy1_mean;
//omment(I: Effective mass mode 2);
      G4double aheavy2_mean;
//omment(I: Width of symmetric mode);
      G4double wasymm_saddle;
//OMMENT(I: Width of asymmetric mode 1);
      G4double waheavy1_saddle;
//OMMENT(I: Width of asymmetric mode 2);
      G4double waheavy2_saddle;
//omment(I: Width of symmetric mode);
      G4double wasymm;
//OMMENT(I: Width of asymmetric mode 1);
      G4double waheavy1;
//OMMENT(I: Width of asymmetric mode 2);
      G4double waheavy2;
//OMMENT(I: Even-odd effect in Z);
      G4double r_e_o,r_e_o_exp;
//OMMENT(I: Curveture of symmetric valley);
      G4double cz_symm;
//OMMENT(I: Curvature of mass distribution for fixed Z);
      G4double cn;
//OMMENT(I: Curvature of Z distribution for fixed A);
      G4double cz;
//OMMENT(Minimum neutron width for constant Z);
      const G4double sigzmin = 1.16;
//OMMENT(Surface distance of scission configuration);
      const G4double d = 2.0;

//   /* Charge polarisation from Wagemanns p. 397: */
//OMMENT(Charge polarisation standard I);
      const G4double cpol1 = 0.65;
//OMMENT(Charge polarisation standard II);
      const G4double cpol2 = 0.55;
//OMMENT(=1: Polarisation simult. in N and Z);
      const G4int nzpol = 1;
//OMMENT(=1: test output, =0: no test output);
      const G4int itest = 0;
      
//      G4double UMASS, ECOUL, reps1, reps2, rn1_pol;
      G4double reps1, reps2, rn1_pol;
//      Float_t HAZ,GAUSSHAZ;
      G4int kkk;

//     I_MODE = 0;

      if (itest == 1) {
         std::cout << " cn mass " << a << std::endl;
	 std::cout << " cn charge " << z << std::endl;
	 std::cout << " cn energy " << e << std::endl;
      }

//     /* average Z of asymmetric and symmetric components: */
      n = a - z;  /* neutron number of the fissioning nucleus */

      k = 0;
      icz = 0;
      if ( (pow(z,2)/a < 25.0) || (n < nheavy2) || (e > 500.0) ) {
          icz = -1;
//          GOTO 1002;
	    goto milledeux;
      }

      nlight1 = n - nheavy1;
      nlight2 = n - nheavy2;
	
//    /* Polarisation assumed for standard I and standard II:
//      Z - Zucd = cpol (for A = const);
//      from this we get (see Armbruster)
//      Z - Zucd =  Acn/Ncn * cpol (for N = const)                        */

      zheavy1_shell = ((nheavy1/n) * z) - ((a/n) * cpol1);
      zheavy2_shell = ((nheavy2/n) * z) - ((a/n) * cpol2);

      e_saddle_scission = 
             (-24.0 + 0.02227 * (pow(z,2))/(pow(a,0.33333)) ) * friction_factor;
    
//      /* Energy dissipated from saddle to scission                        */
//      /* F. Rejmund et al., Nucl. Phys. A 678 (2000) 215, fig. 4 b        */
//      E_saddle_scission = DMAX1(0.,E_saddle_scission);
      if (e_saddle_scission > 0.) {
      	e_saddle_scission = e_saddle_scission;
      }
      else {
      	e_saddle_scission = 0.;
      }
//     /* Semiempirical fission model: */

//    /* Fit to experimental result on curvature of potential at saddle */
//           /* reference:                                              */
//    /* IF Z**2/A < 33.15E0 THEN
//       MassCurv = 30.5438538E0 - 4.00212049E0 * Z**2/A
//                               + 0.11983384E0 * Z**4 / (A**2) ;
//     ELSE
//       MassCurv = 10.E0 ** (7.16993332E0 - 0.26602401E0 * Z**2/A
//                               + 0.00283802E0 * Z**4 / (A**2)) ;  */
//  /* New parametrization of T. Enqvist according to Mulgin et al. 1998 */
      if ( (pow(z,2))/a < 34.0) {
       masscurv =  pow( 10.0,(-1.093364 + 0.082933 * (pow(z,2)/a)
                   - 0.0002602 * (pow(z,4)/pow(a,2))) );
      } else {
       masscurv = pow( 10.0,(3.053536 - 0.056477 * (pow(z,2)/a)
                  + 0.0002454 * (pow(z,4)/pow(a,2))) );
      }

      cz_symm = (8.0/pow(z,2)) * masscurv;

      if (itest == 1) {
         std::cout << "cz_symmetry= " << cz_symm << std::endl;
      }

      if (cz_symm < 0) {
          icz = -1;
//          GOTO 1002;
	  goto milledeux;
      }

//  /* proton number in symmetric fission (centre) */
      zsymm  = z/2.0;
      nsymm  = n/2.0;
      asymm = nsymm + zsymm;

      zheavy1 = (cz_symm*zsymm + cz_asymm1_shell*zheavy1_shell)/(cz_symm + cz_asymm1_shell);
      zheavy2 = (cz_symm*zsymm + cz_asymm2_shell*zheavy2_shell)/(cz_symm + cz_asymm2_shell);
//            /* position of valley due to influence of liquid-drop potential */
      nheavy1_eff = (zheavy1 + (a/n * cpol1))*(n/z);
      nheavy2_eff = (zheavy2 + (a/n * cpol2))*(n/z);
      nlight1_eff = n - nheavy1_eff;
      nlight2_eff = n - nheavy2_eff;
//  /* proton number of light fragments (centre) */
      zlight1 = z - zheavy1;
//  /* proton number of light fragments (centre) */
      zlight2 = z - zheavy2;
      aheavy1 = nheavy1_eff + zheavy1;
      aheavy2 = nheavy2_eff + zheavy2;
      aheavy1_mean = aheavy1;
      aheavy2_mean = aheavy2;
      alight1 = nlight1_eff + zlight1;
      alight2 = nlight2_eff + zlight2;

      a_levdens = a / xlevdens;
      a_levdens_heavy1 = aheavy1 / xlevdens;
      a_levdens_heavy2 = aheavy2 / xlevdens;
      a_levdens_light1 = alight1 / xlevdens;
      a_levdens_light2 = alight2 / xlevdens;
      gamma = a_levdens / (0.4 * (pow(a,1.3333)) );
      gamma_heavy1 = ( a_levdens_heavy1 / (0.4 * (pow(aheavy1,1.3333)) ) ) * fgamma1;
      gamma_heavy2 = a_levdens_heavy2 / (0.4 * (pow(aheavy2,1.3333)) );

      cz_asymm1_saddle = cz_asymm1_shell + cz_symm;
      cz_asymm2_saddle = cz_asymm2_shell + cz_symm;
	
// Up to here: Ok! Checked CS 10/10/05      	   

      cn = umass(zsymm,(nsymm+1.),0.0) + umass(zsymm,(nsymm-1.),0.0)
            + 1.44 * (pow(zsymm,2))/
                     ( (pow(r_null,2)) * 
		       ( pow((asymm+1.0),0.33333) + pow((asymm-1.0),0.33333) ) *
		       ( pow((asymm+1.0),0.33333) + pow((asymm-1.0),0.33333) ) )
            - 2.0 * umass(zsymm,nsymm,0.0)
            - 1.44 * (pow(zsymm,2))/
	             ( ( 2.0 * r_null * (pow(asymm,0.33333)) ) * 
		       ( 2.0 * r_null * (pow(asymm,0.33333)) ) );
	
// /* shell effect in valley of mode 1 */
      delta_u1 = delta_u1_shell + (pow((zheavy1_shell-zheavy1),2))*cz_asymm1_shell;
// /* shell effect in valley of mode 2 */
      delta_u2 = delta_u2_shell + (pow((zheavy2_shell-zheavy2),2))*cz_asymm2_shell;

//     /* liquid drop energies
//        at the centres of the different shell effects
//        with respect to liquid drop at symmetry: */
      epot0_mode1_saddle = (pow((zheavy1-zsymm),2)) * cz_symm;
      epot0_mode2_saddle = (pow((zheavy2-zsymm),2)) * cz_symm;
      epot0_symm_saddle = 0.0;
      
      if (itest == 1) {
	std::cout << "check zheavy1 = " << zheavy1  << std::endl;
	std::cout << "check zheavy2 = " << zheavy2  << std::endl;
	std::cout << "check zsymm = " << zsymm  << std::endl;
	std::cout << "check czsymm = " << cz_symm  << std::endl;
        std::cout << "check epot0_mode1_saddle = " << epot0_mode1_saddle  << std::endl;
	std::cout << "check epot0_mode2_saddle = " << epot0_mode2_saddle  << std::endl;
	std::cout << "check epot0_symm_saddle = " << epot0_symm_saddle  << std::endl;
	std::cout << "delta_u1 = " << delta_u1 << std::endl;
	std::cout << "delta_u2 = " << delta_u2 << std::endl;
      }
      
//     /* energies including shell effects
//        at the centres of the different shell effects
//        with respect to liquid drop at symmetry: */
      epot_mode1_saddle = epot0_mode1_saddle + delta_u1;
      epot_mode2_saddle = epot0_mode2_saddle + delta_u2;
      epot_symm_saddle = epot0_symm_saddle;
      if (itest == 1) {
        std::cout << "check epot_mode1_saddle = " << epot_mode1_saddle  << std::endl;
	std::cout << "check epot_mode2_saddle = " << epot_mode2_saddle  << std::endl;
	std::cout << "check epot_symm_saddle = " << epot_symm_saddle  << std::endl;
      }

//     /* Minimum of potential with respect to ld potential at symmetry */
      dueff = min(epot_mode1_saddle,epot_mode2_saddle);
      dueff = min(dueff,epot_symm_saddle);
      dueff = dueff - epot_symm_saddle;

      eld = e + dueff + e_zero_point;
      
      if (itest == 1) {
      	std::cout << "check dueff = " << dueff  << std::endl;
      	std::cout << "check e = " << e  << std::endl;
      	std::cout << "check e_zero_point = " << e_zero_point  << std::endl;
      	std::cout << "check eld = " << eld  << std::endl;
      }
// Up to here: Ok! Checked CS 10/10/05 
     
//          /* E = energy above lowest effective barrier */
//          /* Eld = energy above liquid-drop barrier */

//     /* Due to this treatment the energy E on input means the excitation  */
//     /* energy above the lowest saddle.                                   */

//  /* These energies are not used */
      eheavy1 = e * aheavy1 / a;
      eheavy2 = e * aheavy2 / a;
      elight1 = e * alight1 / a;
      elight2 = e * alight2 / a;

      epsilon0_1_saddle = eld - e_zero_point - epot0_mode1_saddle;
//            /* excitation energy at saddle mode 1 without shell effect */
      epsilon0_2_saddle = eld - e_zero_point - epot0_mode2_saddle;
//            /* excitation energy at saddle mode 2 without shell effect */

      epsilon_1_saddle = eld - e_zero_point - epot_mode1_saddle;
//            /* excitation energy at saddle mode 1 with shell effect */
      epsilon_2_saddle = eld - e_zero_point - epot_mode2_saddle;
//            /* excitation energy at saddle mode 2 with shell effect */
      epsilon_symm_saddle = eld - e_zero_point - epot_symm_saddle;

//  /* global parameters */
      eexc1_saddle = epsilon_1_saddle;
      eexc2_saddle = epsilon_2_saddle;
      eexc_max = max(eexc1_saddle,eexc2_saddle);
      eexc_max = max(eexc_max,eld);
       
//         /* EEXC_MAX is energy above the lowest saddle */


      epsilon0_1_scission = eld + e_saddle_scission - epot0_mode1_saddle;
//                    /* excitation energy without shell effect */
      epsilon0_2_scission = eld + e_saddle_scission - epot0_mode2_saddle;
//                    /* excitation energy without shell effect */

      epsilon_1_scission = eld + e_saddle_scission - epot_mode1_saddle;
//                    /* excitation energy at scission */
      epsilon_2_scission = eld+ e_saddle_scission - epot_mode2_saddle;
//                    /* excitation energy at scission */
      epsilon_symm_scission = eld + e_saddle_scission - epot_symm_saddle;
//           /* excitation energy of symmetric fragment at scission */

//     /* Calculate widhts at the saddle:                */

      e_eff1_saddle = epsilon0_1_saddle - delta_u1 * (exp((-epsilon_1_saddle*gamma)));
   
      if (e_eff1_saddle > 0.0) {
         wzasymm1_saddle = sqrt( (0.5 * 
         (sqrt(1.0/a_levdens*e_eff1_saddle)) /
         (cz_asymm1_shell * exp((-epsilon_1_saddle*gamma)) + cz_symm) ) );
      } 
      else {
         wzasymm1_saddle = 1.0;
      }

      e_eff2_saddle = epsilon0_2_saddle - delta_u2 * (exp((-epsilon_2_saddle*gamma)));
      if (e_eff2_saddle > 0.0) {
       wzasymm2_saddle = sqrt( (0.5 * 
         (sqrt(1.0/a_levdens*e_eff2_saddle)) /
         (cz_asymm2_shell * exp((-epsilon_2_saddle*gamma)) + cz_symm) ) );
      } 
      else {
         wzasymm2_saddle = 1.0;
      }

      if (eld > e_zero_point) {
       if ( (eld + epsilon_symm_saddle) < 0.0)  {
       	std::cout << "<e> eld + epsilon_symm_saddle < 0" << std::endl;
       }
       wzsymm_saddle = sqrt( (0.5 * 
         (sqrt(1.0/a_levdens*(eld+epsilon_symm_saddle))) / cz_symm ) );
      } else {
       wzsymm_saddle = 1.0;
      }

      if (itest == 1) {
      	std::cout << "wz1(saddle) = " << wzasymm1_saddle << std::endl;
	std::cout << "wz2(saddle) = " << wzasymm2_saddle << std::endl;
	std::cout << "wzsymm(saddle) = " << wzsymm_saddle << std::endl;
      }
     
//     /* Calculate widhts at the scission point: */
//     /* fits of ref. Beizin 1991 (Plots brought to GSI by Sergei Zhdanov) */

      wzsymm_scission = wzsymm_saddle;

      if (e_saddle_scission == 0.0) {

        wzasymm1_scission = wzasymm1_saddle;
        wzasymm2_scission = wzasymm2_saddle;

      } 
      else {

        if (nheavy1_eff > 75.0) {
           wzasymm1_scission = (sqrt(21.0)) * z/a;
           wzasymm2_scission = (sqrt (max( (70.0-28.0)/3.0*(z*z/a-35.0)+28.,0.0 )) ) * z/a;
        } 
	else {
           wzasymm1_scission = wzasymm1_saddle;
           wzasymm2_scission = wzasymm2_saddle;
        }

      }

      wzasymm1_scission = max(wzasymm1_scission,wzasymm1_saddle);
      wzasymm2_scission = max(wzasymm2_scission,wzasymm2_saddle);

      wzasymm1 = wzasymm1_scission * fwidth_asymm1;
      wzasymm2 = wzasymm2_scission * fwidth_asymm2;
      wzsymm = wzsymm_scission;

/*      if (ITEST == 1) {
      	std::cout << "WZ1(scission) = " << WZasymm1_scission << std::endl;
	std::cout << "WZ2(scission) = " << WZasymm2_scission << std::endl;
	std::cout << "WZsymm(scission) = " << WZsymm_scission << std::endl;
      }
      if (ITEST == 1) {
      	std::cout << "WZ1(scission) final= " << WZasymm1 << std::endl;
	std::cout << "WZ2(scission) final= " << WZasymm2 << std::endl;
	std::cout << "WZsymm(scission) final= " << WZsymm << std::endl;
      } */
      
      wasymm = wzsymm * a/z;
      waheavy1 = wzasymm1 * a/z;
      waheavy2 = wzasymm2 * a/z;

      wasymm_saddle = wzsymm_saddle * a/z;
      waheavy1_saddle = wzasymm1_saddle * a/z;
      waheavy2_saddle = wzasymm2_saddle * a/z;

      if (itest == 1) {
      	std::cout << "wasymm = " << wzsymm << std::endl;
	std::cout << "waheavy1 = " << waheavy1 << std::endl;
	std::cout << "waheavy2 = " << waheavy2 << std::endl;
      }
// Up to here: Ok! Checked CS 11/10/05
            
      if ( (epsilon0_1_saddle - delta_u1*exp((-epsilon_1_saddle*gamma_heavy1))) < 0.0) {
       sasymm1 = -10.0;
      } 
      else {
       sasymm1 = 2.0 * sqrt( a_levdens * (epsilon0_1_saddle - 
                 delta_u1*(exp((-epsilon_1_saddle*gamma_heavy1))) ) );
      }

      if ( (epsilon0_2_saddle - delta_u2*exp((-epsilon_2_saddle*gamma_heavy2))) < 0.0) {
       sasymm2 = -10.0;
      } 
      else {
       sasymm2 = 2.0 * sqrt( a_levdens * (epsilon0_2_saddle - 
                 delta_u2*(exp((-epsilon_2_saddle*gamma_heavy2))) ) );
      }
              
      if (epsilon_symm_saddle > 0.0) {
       ssymm = 2.0 * sqrt( a_levdens*(epsilon_symm_saddle) );
      } 
      else {
       ssymm = -10.0;
      }
      
      if (ssymm > -10.0) {
       ysymm = 1.0;

       if (epsilon0_1_saddle < 0.0) {
//  /* low energy */
         yasymm1 = exp((sasymm1-ssymm)) * wzasymm1_saddle / wzsymm_saddle * 2.0;
//           /* factor of 2 for symmetry classes */
       } 
       else {
//        /* high energy */
         ssymm_mode1 = 2.0 * sqrt( a_levdens*(epsilon0_1_saddle) );
         yasymm1 = ( exp((sasymm1-ssymm)) - exp((ssymm_mode1 - ssymm)) )  
	           * wzasymm1_saddle / wzsymm_saddle * 2.0;
       }

       if (epsilon0_2_saddle < 0.0) {
//  /* low energy */
         yasymm2 = exp((sasymm2-ssymm)) * wzasymm2_saddle / wzsymm_saddle * 2.0;
//           /* factor of 2 for symmetry classes */
       } 
       else {
//        /* high energy */
         ssymm_mode2 = 2.0 * sqrt( a_levdens*(epsilon0_2_saddle) );
         yasymm2 = ( exp((sasymm2-ssymm)) - exp((ssymm_mode2 - ssymm)) )  
	           * wzasymm2_saddle / wzsymm_saddle * 2.0;
       }       
//                            /* difference in the exponent in order */
//                            /* to avoid numerical overflow         */

      } 
      else {
       if ( (sasymm1 > -10.0) && (sasymm2 > -10.0) ) {
          ysymm = 0.0;
          yasymm1 = exp(sasymm1) * wzasymm1_saddle * 2.0;
          yasymm2 = exp(sasymm2) * wzasymm2_saddle * 2.0;
       }
      }
        
//  /* normalize */
      ysum = ysymm + yasymm1 + yasymm2;
      if (ysum > 0.0) {
       ysymm = ysymm / ysum;
       yasymm1 = yasymm1 / ysum;
       yasymm2 = yasymm2 / ysum;
       yasymm = yasymm1 + yasymm2;
      } 
      else {
       ysymm = 0.0;
       yasymm1 = 0.0;
       yasymm2 = 0.0;
//        /* search minimum threshold and attribute all events to this mode */
       if ( (epsilon_symm_saddle < epsilon_1_saddle) && (epsilon_symm_saddle < epsilon_2_saddle) ) {
        ysymm = 1.0;
       } 
       else {
          if (epsilon_1_saddle < epsilon_2_saddle) {
            yasymm1 = 1.0;
          } 
	  else {
            yasymm2 = 1.0;
          }
       }
      }

      if (itest == 1) {
      	std::cout << "ysymm normalized= " << ysymm  << std::endl;
	std::cout << "yasymm1 normalized= " << yasymm1  << std::endl;
	std::cout << "yasymm2 normalized= " << yasymm2  << std::endl;
      }
// Up to here: Ok! Ckecked CS 11/10/05      
      
//      /* even-odd effect */
//      /* simple parametrization KHS, Nov. 2000. From Rejmund et al. */
      if ((int)(z) % 2 == 0) {
      r_e_o_exp = -0.017 * (e_saddle_scission + eld) * (e_saddle_scission + eld);
         if ( r_e_o_exp < -307.0) {
	 	r_e_o_exp = -307.0;
         	r_e_o = pow(10.0,r_e_o_exp);
	}
	else {
		r_e_o = pow(10.0,r_e_o_exp);
	}
      } 
      else {
         r_e_o = 0.0;
      }
    
//      $LOOP;    /* event loop */
//     I_COUNT = I_COUNT + 1;

//     /* random decision: symmetric or asymmetric */
//     /* IMODE = 1 means asymmetric fission, mode 1,
//        IMODE = 2 means asymmetric fission, mode 2,
//        IMODE = 3 means symmetric  */
//      RMODE = dble(HAZ(k));
//      rmode = rnd.rndm();  
      rmode = haz(k);
// Cast for test CS 11/10/05
//      RMODE = 0.54;    
      if (rmode < yasymm1) {
         imode = 1;
      }
      if ( (rmode > yasymm1) && (rmode < (yasymm1+yasymm2)) ) {
      	 imode = 2;
      }
      if ( (rmode > yasymm1) && (rmode > (yasymm1+yasymm2)) ) {
	 imode = 3;
      }
  
//     /* determine parameters of the Z distribution */
	if (imode == 1) {
	 z1mean = zheavy1;
         z1width = wzasymm1;
	}
	if (imode == 2) {
	 z1mean = zheavy2;
         z1width = wzasymm2;
	}
	if (imode == 3) {
	 z1mean = zsymm;
         z1width = wzsymm;
	}

      if (itest == 1) {
      	std::cout << "nbre aleatoire tire " << rmode << std::endl;
	std::cout << "fission mode " << imode << std::endl;
	std::cout << "z1mean= " << z1mean << std::endl;
	std::cout << "z1width= " << z1width << std::endl;
      }
		      
//     /* random decision: Z1 and Z2 at scission: */
      z1 = 1.0;
      z2 = 1.0;
      while  ( (z1<5.0) || (z2<5.0) ) {
//         Z1 = dble(GAUSSHAZ(K,sngl(Z1mean),sngl(Z1width)));
//	 z1 = rnd.gaus(z1mean,z1width);
	z1 = gausshaz(k, z1mean, z1width);
	z2 = z - z1;
      }
      if (itest == 1) {
      	std::cout << "ff charge sample " << std::endl;
	std::cout << "z1 =  " << z1 << std::endl;
	std::cout << "z2 = " << z2 << std::endl;
      }

//     CALL EVEN_ODD(Z1,R_E_O,I_HELP);
//         /* Integer proton number with even-odd effect */
//     Z1 = REAL(I_HELP)
//      /* Z1 = INT(Z1+0.5E0); */
      z2 = z - z1;

//     /* average N of both fragments: */
      if (imode == 1) {
         n1mean = (z1 + cpol1 * a/n) * n/z;
      }
      if (imode == 2) { 
         n1mean = (z1 + cpol2 * a/n) * n/z;
      }
/*       CASE(99)   ! only for testing;
         N1UCD = Z1 * N/Z;
         N2UCD = Z2 * N/Z;
         re1 = UMASS(Z1,N1UCD,0.6) +;
     &         UMASS(Z2,N2UCD,0.6) +;
     &         ECOUL(Z1,N1UCD,0.6,Z2,N2UCD,0.6,d);
         re2 = UMASS(Z1,N1UCD+1.,0.6) +;
     &         UMASS(Z2,N2UCD-1.,0.6) +;
     &         ECOUL(Z1,N1UCD+1.,0.6,Z2,N2UCD-1.,0.6,d);
         re3 = UMASS(Z1,N1UCD+2.,0.6) +;
     &         UMASS(Z2,N2UCD-2.,0.6) +;
     &         ECOUL(Z1,N1UCD+2.,0.6,Z2,N2UCD-2.,0.6,d);
         eps2 = (re1-2.0*re2+re3) / 2.0;
         eps1 = re2 - re1 - eps2;
         DN1_POL = - eps1 / (2.0 * eps2);
         N1mean = N1UCD + DN1_POL; */
       if (imode == 3) {
         n1ucd = z1 * n/z;
         n2ucd = z2 * n/z;
         re1 = umass(z1,n1ucd,0.6) + umass(z2,n2ucd,0.6) + ecoul(z1,n1ucd,0.6,z2,n2ucd,0.6,d);
         re2 = umass(z1,n1ucd+1.,0.6) + umass(z2,n2ucd-1.,0.6) + ecoul(z1,n1ucd+1.,0.6,z2,n2ucd-1.,0.6,d);
         re3 = umass(z1,n1ucd+2.,0.6) + umass(z2,n2ucd-2.,0.6) + ecoul(z1,n1ucd+2.,0.6,z2,n2ucd-2.,0.6,d);
         eps2 = (re1-2.0*re2+re3) / 2.0;
         eps1 = re2 - re1 - eps2;
         dn1_pol = - eps1 / (2.0 * eps2);
         n1mean = n1ucd + dn1_pol;
       }
// all fission modes features have been checked CS 11/10/05
      n2mean = n - n1mean;
      z2mean = z - z1mean;
      
//   /* Excitation energies */
//     /* formulated in energies in close consistency with the fission model */

//     /* E_defo = UMASS(Z*0.5E0,N*0.5E0,0.6E0) -
//                 UMASS(Z*0.5E0,N*0.5E0,0);   */
//       /* calculates the deformation energy of the liquid drop for
//          deformation beta = 0.6 which is most probable at scission */

//    /* N1R and N2R provisionaly taken without fluctuations in
//       polarisation:                                              */
       n1r = n1mean;
       n2r = n2mean;
       a1r = n1r + z1;
       a2r = n2r + z2;

         if (imode == 1) { /* N = 82 */;
//! /* Eexc at scission */
           e_scission = max(epsilon_1_scission,1.0);
           if (n1mean > (n * 0.5) ) {
//! /* 1. fragment is spherical */
             beta1 = 0.0;
             beta2 = 0.6;
             e1exc = epsilon_1_scission * a1r / a;
             e_defo = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
             e2exc = epsilon_1_scission * a2r / a + e_defo;
	   }
           else {
//!                       /* 2. fragment is spherical */
             beta1 = 0.6;
             beta2 = 0.0;
             e_defo = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
             e1exc = epsilon_1_scission * a1r / a + e_defo;
             e2exc = epsilon_1_scission * a2r / a;
           }
	 }
	   
         if (imode == 2) { 
//! /* N appr. 86 */
           e_scission = max(epsilon_2_scission,1.0);
           if (n1mean >  (n * 0.5) ) { 
//! /* 2. fragment is spherical */
             beta1 = (n1r - nheavy2) * 0.034 + 0.3;
             e_defo = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
             e1exc = epsilon_2_scission * a1r / a + e_defo;
             beta2 = 0.6 - beta1;
             e_defo = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
             e2exc = epsilon_2_scission * a2r / a + e_defo;
           }
	   else {
//!                      /* 1. fragment is spherical */
             beta2 = (n2r - nheavy2) * 0.034 + 0.3;
             e_defo = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
             e2exc = epsilon_2_scission * a2r / a + e_defo;
             beta1 = 0.6 - beta2;
             e_defo = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
             e1exc = epsilon_2_scission * a1r / a + e_defo;
           }
	  }
         
	  if (imode == 3) { 
// ! /* Symmetric fission channel */

//             /* the fit function for beta is the deformation for
//                optimum energy at the scission point, d = 2 */
//             /* beta  : deformation of symmetric fragments */
//             /* beta1 : deformation of first fragment */
//             /* beta2 : deformation of second fragment */
             beta =  0.177963 + 0.0153241 * zsymm - 0.000162037 * zsymm*zsymm;
             beta1 = 0.177963 + 0.0153241 * z1 - 0.000162037 * z1*z1;
//            beta1 = 0.6
             e_defo1 = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
             beta2 = 0.177963 + 0.0153241 * z2 - 0.000162037 * z2*z2;
//            beta2 = 0.6
             e_defo2 = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
             e_asym = umass(z1 , n1r, beta1) + umass(z2, n2r ,beta2)
     		      + ecoul(z1,n1r,beta1,z2,n2r,beta2,2.0)
     		      - 2.0 * umass(zsymm,nsymm,beta)
                      - ecoul(zsymm,nsymm,beta,zsymm,nsymm,beta,2.0);
//            E_asym = CZ_symm * (Z1 - Zsymm)**2
             e_scission = max((epsilon_symm_scission - e_asym),1.0);
//         /*  $LIST(Z1,N1R,Z2,N2R,E_asym,E_scission); */
             e1exc = e_scission * a1r / a + e_defo1;
             e2exc = e_scission * a2r / a + e_defo2;
          }
// Energies checked for all the modes CS 11/10/05
	  
//    /* random decision: N1R and N2R at scission, before evaporation: */
//    /* CN =       UMASS(Zsymm , Nsymm + 1.E0,0) +
//                UMASS(Zsymm, Nsymm - 1.E0,0)
//                 + 1.44E0 * (Zsymm)**2 /
//                 (r_null**2 * ((Asymm+1)**1/3 + (Asymm-1)**1/3)**2 )
//                 - 2.E0 * UMASS(Zsymm,Nsymm,0)
//                 - 1.44E0 * (Zsymm)**2 / (r_null * 2.E0 * (Asymm)**1/3)**2; */


//    /* N1width = sqrt(0.5E0 * sqrt(1.E0/A_levdens*(Eld+E_saddle_scission)) / CN); */
//    /* 8. 9. 1998: KHS (see also consideration in the first comment block)
//       sigma_N(Z=const) = A/Z  * sigma_Z(A=const)
//       sigma_Z(A=const) = 0.4 to 0.5  (from Lang paper Nucl Phys. A345 (1980) 34)
//       sigma_N(Z=const) = 0.45 * A/Z  (= 1.16 for 238U)
//       therefore: SIGZMIN = 1.16                                              */

         if ( (imode == 1) || (imode == 2) ) {
           cn=(umass(z1,n1mean+1.,beta1) + umass(z1,n1mean-1.,beta1)
             + umass(z2,n2mean+1.,beta2) + umass(z2,n2mean-1.,beta2)
             + ecoul(z1,n1mean+1.,beta1,z2,n2mean-1.,beta2,2.0)
             + ecoul(z1,n1mean-1.,beta1,z2,n2mean+1.,beta2,2.0)
             - 2.0 * ecoul(z1,n1mean,beta1,z2,n2mean,beta2,2.0)
             - 2.0 * umass(z1, n1mean, beta1)
             - 2.0 * umass(z2, n2mean, beta2) ) * 0.5;
//          /* Coulomb energy neglected for the moment! */
//           IF (E_scission.lt.0.) Then
//             write(6,*)'<E> E_scission < 0, MODE 1,2'
//           ENDIF
//           IF (CN.lt.0.) Then
//             write(6,*)'CN < 0, MODE 1,2'
//           ENDIF
          n1width=sqrt( (0.5 * (sqrt(1.0/a_levdens*(e_scission)))/cn) );
          n1width=max(n1width, sigzmin);

//          /* random decision: N1R and N2R at scission, before evaporation: */
           n1r = 1.0;
           n2r = 1.0;
           while  ( (n1r<5.0) || (n2r<5.0) ) {
//             n1r = dble(gausshaz(k,sngl(n1mean),sngl(n1width)));
//	   	n1r = rnd.gaus(n1mean,n1width);
	     n1r = gausshaz(k, n1mean, n1width);
	     n2r = n - n1r;
           }	   
//           N1R = GAUSSHAZ(K,N1mean,N1width)
	   if (itest == 1) {
	   	std::cout << "after neutron sample " << n1r << std::endl;
	   }	   
           n1r = (float)( (int)((n1r+0.5)) );
           n2r = n - n1r;
 
           even_odd(z1,r_e_o,i_help);
//         /*  proton number with even-odd effect */
           z1 = (float)(i_help);
           z2 = z - z1;

           a1r = z1 + n1r;
           a2r = z2 + n2r;
	}
	
         if (imode == 3) {
//!  /* When(3) */
          if (nzpol > 0.0) {
//           /* We treat a simultaneous split in Z and N to determine polarisation */
           cz = ( umass(z1-1., n1mean+1.,beta1)
               +  umass(z2+1., n2mean-1.,beta1) 
               +  umass(z1+1., n1mean-1.,beta2)
               +  umass(z2 - 1., n2mean + 1.,beta2)
               +  ecoul(z1-1.,n1mean+1.,beta1,z2+1.,n2mean-1.,beta2,2.0)
               +  ecoul(z1+1.,n1mean-1.,beta1,z2-1.,n2mean+1.,beta2,2.0)
               -  2.0 * ecoul(z1,n1mean,beta1,z2,n2mean,beta2,2.0)
               -  2.0 * umass(z1, n1mean,beta1)
               -  2.0 * umass(z2, n2mean,beta2) ) * 0.5;
//           IF (E_scission.lt.0.) Then
//             write(6,*) '<E> E_scission < 0, MODE 1,2'
//           ENDIF
//           IF (CZ.lt.0.) Then
//             write(6,*) 'CZ < 0, MODE 1,2'
//           ENDIF
         za1width=sqrt( (0.5 * sqrt(1.0/a_levdens*(e_scission)) / cz) );
         za1width=sqrt( (max((za1width*za1width-(1.0/12.0)),0.1)) );
//                        /* Check the value of 0.1 ! */
//                        /* Shephard correction */
           a1r = z1 + n1mean;
           a1r = (float)((int)((a1r+0.5)));
           a2r = a - a1r;
//           /* A1R and A2R are integer numbers now */
//        /* $LIST(A1R,A2R,ZA1WIDTH); */

           n1ucd = n/a * a1r;
           n2ucd = n/a * a2r;
           z1ucd = z/a * a1r;
           z2ucd = z/a * a2r;

           re1 = umass(z1ucd-1.,n1ucd+1.,beta1) + umass(z2ucd+1.,n2ucd-1.,beta2)
                 + ecoul(z1ucd-1.,n1ucd+1.,beta1,z2ucd+1.,n2ucd-1.,beta2,d);
           re2 = umass(z1ucd,n1ucd,beta1) + umass(z2ucd,n2ucd,beta2)
                 + ecoul(z1ucd,n1ucd,beta1,z2ucd,n2ucd,beta2,d);
           re3 = umass(z1ucd+1.,n1ucd-1.,beta1) + umass(z2ucd-1.,n2ucd+1.,beta2) +
                 + ecoul(z1ucd+1.,n1ucd-1.,beta1,z2ucd-1.,n2ucd+1.,beta2,d);
		 
           eps2 = (re1-2.0*re2+re3) / 2.0;
           eps1 = (re3 - re1)/2.0;
           dn1_pol = - eps1 / (2.0 * eps2);
           z1 = z1ucd + dn1_pol;
	   if (itest == 1) {
	   	std::cout << "before proton sample " << z1 << std::endl;
	   }  
//           Z1 = dble(GAUSSHAZ(k,sngl(Z1),sngl(ZA1width)));
//	   z1 = rnd.gaus(z1,za1width);
	   z1 = gausshaz(k, z1, za1width);
	   if (itest == 1) {
	   	std::cout << "after proton sample " << z1 << std::endl;
	   }	   
           even_odd(z1,r_e_o,i_help);
//         /* proton number with even-odd effect */
           z1 = (float)(i_help);
           z2 = (float)((int)( (z - z1 + 0.5)) );

           n1r = a1r - z1;
           n2r = n - n1r;
          } 
	  else {
//           /* First division of protons, then adjustment of neutrons */
           cn = ( umass(z1, n1mean+1.,beta1) + umass(z1, n1mean-1., beta1)
                  + umass(z2, n2mean+1.,beta2) + umass(z2, n2mean-1., beta2)
                  + ecoul(z1,n1mean+1.,beta1,z2,n2mean-1.,beta2,2.0)
                  + ecoul(z1,n1mean-1.,beta1,z2,n2mean+1.,beta2,2.0)
                  - 2.0 * ecoul(z1,n1mean,beta1,z2,n2mean,beta2,2.0)
                  - 2.0 * umass(z1, n1mean, 0.6)
                  - 2.0 * umass(z2, n2mean, 0.6) ) * 0.5;
//          /* Coulomb energy neglected for the moment! */
//           IF (E_scission.lt.0.) Then
//             write(6,*) '<E> E_scission < 0, MODE 1,2'
//           Endif
//           IF (CN.lt.0.) Then
//             write(6,*) 'CN < 0, MODE 1,2'
//           Endif
          n1width=sqrt( (0.5 * sqrt(1.0/a_levdens*(e_scission)) / cn) );
          n1width=max(n1width, sigzmin);

//          /* random decision: N1R and N2R at scission, before evaporation: */
//           N1R = dble(GAUSSHAZ(k,sngl(N1mean),sngl(N1width)));
//	   n1r = rnd.gaus(n1mean,n1width);
	  n1r = gausshaz(k, n1mean, n1width);
           n1r = (float)( (int)((n1r+0.5)) );
           n2r = n - n1r;

           even_odd(z1,r_e_o,i_help);
//         /* Integer proton number with even-odd effect */
           z1 = (float)(i_help);
           z2 = z - z1;

           a1r = z1 + n1r;
           a2r = z2 + n2r;
          
	  }
	}

      if (itest == 1) {
      	std::cout << "remid imode = " << imode << std::endl;
	std::cout << "n1width =  " << n1width << std::endl;
	std::cout << "n1r = " << n1r << std::endl;
	std::cout << "a1r = " << a1r << std::endl;
	std::cout << "n2r = " << n2r << std::endl;
	std::cout << "a2r = " << a2r << std::endl;
      }
// Up to here: checked CS 11/10/05	
	
//      /* Extracted from Lang et al. Nucl. Phys. A 345 (1980) 34 */
       e1exc_sigma = 5.5;
       e2exc_sigma = 5.5;

neufcentquatrevingtsept:
//       E1final = dble(Gausshaz(k,sngl(E1exc),sngl(E1exc_sigma)));
//       E2final = dble(Gausshaz(k,sngl(E2exc),sngl(E2exc_sigma)));
//        e1final = rnd.gaus(e1exc,e1exc_sigma);
//        e2final = rnd.gaus(e2exc,e2exc_sigma);
       e1final = gausshaz(k, e1exc, e1exc_sigma);
       e2final = gausshaz(k, e2exc, e2exc_sigma);
       if ( (e1final < 0.0) || (e2final < 0.0) ) goto neufcentquatrevingtsept;
       if (itest == 1) {
	 std::cout << "sampled exc 1 " << e1final << std::endl;
	 std::cout << "sampled exc 2 " << e2final << std::endl;
       }

//      /* OUTPUT QUANTITIES OF THE EVENT GENERATOR:        */

//      /* Quantities before neutron evaporation            */

//      /* Neutron number of prefragments: N1R and N2R      */
//      /* Atomic number of fragments: Z1 and Z2            */
//      /* Kinetic energy of fragments: EkinR1, EkinR2      *7

//      /* Quantities after neutron evaporation:            */

//      /* Neutron number of fragments: N1 and N2           */
//      /* Mass number of fragments: A1 and A2              */
//      /* Atomic number of fragments: Z1 and Z2            */
//      /* Number of evaporated neutrons: N1R-N1 and N2R-N2 */
//      /* Kinetic energy of fragments: EkinR1*A1/A1R and
//                                      EkinR2*A2/A2R       */

       n1 = n1r;
       n2 = n2r;
       a1 = n1 + z1;
       a2 = n2 + z2;
       e1 = e1final;
       e2 = e2final;

//       /* Pre-neutron-emission total kinetic energy: */
       tker = (z1 * z2 * 1.44) /
              ( r0 * pow(a1,0.33333) * (1.0 + 2.0/3.0 * beta1) +
                r0 * pow(a2,0.33333) * (1.0 + 2.0/3.0 * beta2) + 2.0 );
//       /* Pre-neutron-emission kinetic energy of 1. fragment: */
       ekinr1 = tker * a2 / a;
//       /* Pre-neutron-emission kinetic energy of 2. fragment: */
       ekinr2 = tker * a1 / a;

       v1 = sqrt( (ekinr1/a1) ) * 1.3887;
       v2 = sqrt( (ekinr2/a2) ) * 1.3887;

       if (itest == 1) {
	 std::cout << "ekinr1 " << ekinr1 << std::endl;
	 std::cout << "ekinr2 " << ekinr2 << std::endl;
       }

milledeux:       
//**************************
//*** only symmetric fission
//**************************
// Symmetric fission: Ok! Checked CS 10/10/05
      if ( (icz == -1) || (a1 < 0.0) || (a2 < 0.0) ) {
//           IF (z.eq.92) THEN
//              write(6,*)'symmetric fission'
//              write(6,*)'Z,A,E,A1,A2,icz,Atot',Z,A,E,A1,A2,icz,Atot
//           END IF

      if (itest == 1) {
	std::cout << "milledeux: liquid-drop option "  << std::endl;
      }

      n = a-z;
//  proton number in symmetric fission (centre) *
     	zsymm  = z / 2.0;
     	nsymm  = n / 2.0;
     	asymm = nsymm + zsymm;

     	a_levdens = a / xlevdens;

      masscurv = 2.0;
     	cz_symm = 8.0 / pow(z,2) * masscurv;

     	wzsymm = sqrt( (0.5 * sqrt(1.0/a_levdens*e) / cz_symm) ) ;

      if (itest == 1) {
        std::cout << " symmetric high energy fission " << std::endl;
	std::cout << "wzsymm " << wzsymm << std::endl;
      }

        z1mean = zsymm;
        z1width = wzsymm;

// random decision: Z1 and Z2 at scission: */
        z1 = 1.0;
        z2 = 1.0;
        while  ( (z1 < 5.0) || (z2 < 5.0) ) {
//           z1 = dble(gausshaz(kkk,sngl(z1mean),sngl(z1width)));
//	   z1 = rnd.gaus(z1mean,z1width);
	  z1 = gausshaz(kkk, z1mean, z1width);
     	   z2 = z - z1;
        }

      if (itest == 1) {
	std::cout << " z1 " << z1 << std::endl;
	std::cout << " z2 " << z2 << std::endl;
      }
      if (itest == 1) {
	std::cout << " zsymm " << zsymm << std::endl;
	std::cout << " nsymm " << nsymm << std::endl;
	std::cout << " asymm " << asymm << std::endl;
      }
//    	CN =  UMASS(Zsymm , Nsymm + 1.E0) + UMASS(Zsymm, Nsymm - 1.E0)
//    #            + 1.44E0 * (Zsymm)**2 /
//    #            (r_null**2 * ((Asymm+1)**(1./3.) +
//    #            (Asymm-1)**(1./3.))**2 )
//    #            - 2.E0 * UMASS(Zsymm,Nsymm)
//    #            - 1.44E0 * (Zsymm)**2 /
//    #            (r_null * 2.E0 * (Asymm)**(1./3.))**2

        n1ucd = z1 * n/z;
        n2ucd = z2 * n/z;
        re1 = umass(z1,n1ucd,0.6) + umass(z2,n2ucd,0.6) +
              ecoul(z1,n1ucd,0.6,z2,n2ucd,0.6,2.0);
        re2 = umass(z1,n1ucd+1.,0.6) + umass(z2,n2ucd-1.,0.6) +
              ecoul(z1,n1ucd+1.,0.6,z2,n2ucd-1.,0.6,2.0);
        re3 = umass(z1,n1ucd+2.,0.6) + umass(z2,n2ucd-2.,0.6) +
              ecoul(z1,n1ucd+2.,0.6,z2,n2ucd-2.,0.6,2.0);
        reps2 = (re1-2.0*re2+re3)/2.0;
        reps1 = re2 - re1 -reps2;
        rn1_pol = -reps1/(2.0*reps2);
        n1mean = n1ucd + rn1_pol;
        n2mean = n - n1mean;

      if (itest == 1) {
	std::cout << " n1mean " << n1mean << std::endl;
	std::cout << " n2mean " << n2mean << std::endl;
      }

     	cn = (umass(z1,n1mean+1.,0.0) + umass(z1,n1mean-1.,0.0) +
              + umass(z2,n2mean+1.,0.0) + umass(z2,n2mean-1.,0.0)
              - 2.0 * umass(z1,n1mean,0.0) +
              - 2.0 * umass(z2,n2mean,0.0) ) * 0.5;
//      This is an approximation! Coulomb energy is neglected.

     	n1width = sqrt( (0.5 * sqrt(1.0/a_levdens*e) / cn) );

      if (itest == 1) {
	std::cout << " cn " << cn << std::endl;
	std::cout << " n1width " << n1width << std::endl;
      }
	
// random decision: N1R and N2R at scission, before evaporation: */
//       N1R = dfloat(NINT(GAUSSHAZ(KKK,sngl(N1mean),sngl(N1width))));
//      n1r = (float)( (int)(rnd.gaus(n1mean,n1width)) );
      n1r = (float)( (int)(gausshaz(k, n1mean,n1width)) );
      n2r = n - n1r;
// Mass of first and second fragment */
     	a1 = z1 + n1r;
     	a2 = z2 + n2r;

        e1 = e*a1/(a1+a2);
        e2 = e - e*a1/(a1+a2);
      if (itest == 1) {
	std::cout << " n1r " << n1r << std::endl;
	std::cout << " n2r " << n2r << std::endl;
      }

        }

      if (itest == 1) {
	 std::cout << " a1 " << a1 << std::endl;
	 std::cout << " z1 " << z1 << std::endl;
	 std::cout << " a2 " << a2 << std::endl;
	 std::cout << " z2 " << z2 << std::endl;
	 std::cout << " e1 " << e1 << std::endl;
	 std::cout << " e2 " << e << std::endl;
      }

//       /* Pre-neutron-emission total kinetic energy: */
       tker = (z1 * z2 * 1.44) /
              ( r0 * pow(a1,0.33333) * (1.0 + 2.0/3.0 * beta1) +
                r0 * pow(a2,0.33333) * (1.0 + 2.0/3.0 * beta2) + 2.0 );
//       /* Pre-neutron-emission kinetic energy of 1. fragment: */
       ekin1 = tker * a2 / a;
//       /* Pre-neutron-emission kinetic energy of 2. fragment: */
       ekin2 = tker * a1 / a;

       v1 = sqrt( (ekin1/a1) ) * 1.3887;
       v2 = sqrt( (ekin2/a2) ) * 1.3887;

       if (itest == 1) {
       	 std::cout << " kinetic energies " << std::endl;
	 std::cout << " ekin1 " << ekin1 << std::endl;
	 std::cout << " ekin2 " << ekin2 << std::endl;
       }
}

// Methods related to the internal ABLA random number generator. In
// the future the random number generation must be factored into its
// own class

void G4Abla::standardRandom(G4double *rndm, G4int *seed)
{
  // IBM standard random number generator imported from INCL.

  (*seed) = (*seed) * 65539;

  if((*seed) < 0) {
    (*seed) = (*seed) + 2147483647+1;
  }

  (*rndm) = (*seed) * 0.4656613e-9;
}


G4double G4Abla::haz(G4int k)
{
  const G4int pSize = 110;
  static G4double p[pSize];
  static G4int ix,i;
  static G4double x,y,a,haz;
  //    973       c***      si k=<-1 on initialise                                        
  //    974       c***      si k=-1 c'est reproductible                                   
  //    975       c***      si k<>-1 ce n'est pas reproductible

  //    977              save   --> static!!!                        

  if (k <= -1) { //then                                             
    if(k == -1) { //then                                            
      ix=0;
    }
    else {
      x=0.0;
      y=secnds(x);
      ix=int(y*100+43543000);
      if(mod(ix,2) == 0) {
	ix=ix+1;
      }
    }

    //    x=ranf(ix);
    
    // Here we are using random number generator copied from INCL code
    // instead of the CERNLIB one! This causes difficulties for
    // automatic testing since the random number generators, and thus
    // the behavior of the routines in C++ and FORTRAN versions is no
    // longer exactly the same!
    standardRandom(&x, &ix);
    
    for(G4int i = 0; i < pSize; i++) { //do i=1,110                                                 
      //p[i]=ranf(ix);
      standardRandom(&(p[i]), &ix);
    } //end do                                                     

    //    a=ranf(ix);
    standardRandom(&a, &ix);
    k=0;
  } //end if                                                        

  i=nint(100*a)+1;
  haz=p[i];
  //  a=ranf(ix);
  standardRandom(&a, &ix);
  p[i]=a;

  return haz;
}


G4double G4Abla::gausshaz(int k, double xmoy, double sig)
{
  // Gaussian random numbers:

  //   1005       C*** TIRAGE ALEATOIRE DANS UNE GAUSSIENNE DE LARGEUR SIG ET MOYENNE XMOY
  static G4int  iset = 0;
  static G4double v1,v2,r,fac,gset,gausshaz;

  if(iset == 0) { //then                                              
    do {
      v1 = 2.0*haz(k) - 1.0;
      v2 = 2.0*haz(k) - 1.0;
      r = pow(v1,2) + pow(v2,2);
    } while(r >= 1);

    fac = sqrt(-2.*log(r)/r);
    gset = v1*fac;
    gausshaz = v2*fac*sig+xmoy;
    iset = 1;
  }
  else {
    gausshaz=gset*sig+xmoy;
    iset=0;
  }

  return gausshaz;                                                         
}


// Utilities

G4double G4Abla::min(G4double a, G4double b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4Abla::min(G4int a, G4int b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4double G4Abla::max(G4double a, G4double b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4Abla::max(G4int a, G4int b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4Abla::nint(G4double number)
{
  G4double intpart;
  G4double fractpart;
  fractpart = modf(number, &intpart);
  if(number == 0) {
    return 0;
  }
  if(number > 0) {
    if(fractpart < 0.5) {
      return int(floor(number));
    }
    else {
      return int(ceil(number));
    }
  }
  if(number < 0) {
    if(fractpart < -0.5) {
      return int(floor(number));
    }
    else {
      return int(ceil(number));
    }
  }
}

G4int G4Abla::secnds(G4double x)
{
  time_t mytime;
  tm *mylocaltime;

  time(&mytime);
  mylocaltime = localtime(&mytime);

  return(mylocaltime->tm_hour*60*60 + mylocaltime->tm_min*60 + mylocaltime->tm_sec);
}

G4int G4Abla::mod(G4int a, G4int b)
{
  if(b != 0) {
    return (a - (a/b)*b);
  }
  else {
    return 0;
  } 
}

G4double G4Abla::dint(G4double a)
{
  G4double value = 0.0;

  if(a < 0.0) {
    value = double(ceil(a));
  }
  else {
    value = double(floor(a));
  }

  return value;
}

G4int G4Abla::idint(G4double a)
{
  G4int value = 0;

  if(a < 0) {
    value = int(ceil(a));
  }
  else {
    value = int(floor(a));
  }

  return value;
}

G4int G4Abla::idnint(G4double a)
{
  G4int value = 0;

  G4int valueCeil = int(ceil(a));
  G4int valueFloor = int(floor(a));

  if(abs(value - valueCeil) < abs(value - valueFloor)) {
    return valueCeil;
  }
  else {
    return valueFloor;
  }
}

G4double G4Abla::dmin1(G4double a, G4double b, G4double c)
{
  if(a < b && a < c) {
    return a;
  }
  if(b < a && b < c) {
    return b;
  }
  if(c < a && c < b) {
    return c;
  }
}

G4double G4Abla::utilabs(G4double a)
{
  if(a > 0) {
    return a;
  }
  if(a < 0) {
    return (-1*a);
  }
  if(a == 0) {
    return a;
  }
}
