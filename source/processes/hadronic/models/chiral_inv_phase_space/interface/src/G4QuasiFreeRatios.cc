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
//
// $Id: G4QuasiFreeRatios.cc,v 1.1 2007-03-20 17:49:39 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Physics class: G4QuasiFreeRatios for N+A elastic cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Oct-06
// 
//================================================================================

//#define debug
//#define isodebug
//#define pdebug
//#define ppdebug
//#define tdebug
//#define sdebug

#include "G4QuasiFreeRatios.hh"

// Initialization of the static parameters
const G4double G4QuasiFreeRatios::tolerance=.001; // relTolerance in dMom to get old retio
const G4int G4QuasiFreeRatios::nPoints=128;//#of points in AMDB tables (>anyPar)(D)
const G4int G4QuasiFreeRatios::nLast=nPoints-1;// the Last element in the table   (D)
G4double  G4QuasiFreeRatios::lPMin=-8.; // Min tabulated logarithmic Momentum     (D)
G4double  G4QuasiFreeRatios::lPMax= 8.; // Max tabulated logarithmic Momentum     (D)
G4double  G4QuasiFreeRatios::dlnP=(lPMax-lPMin)/nLast;// Log step in the table    (D)
std::pair<G4double,G4double> G4QuasiFreeRatios::lastRAT(0.,0.);// Last QE,QF rat  (L)
G4double  G4QuasiFreeRatios::lastLP=-10.;// Last log(mom_of_the_incident_hadron)  (L)
G4double  G4QuasiFreeRatios::theQFP=0.; // The Last p-dependent QF/IN function    (L)
G4double  G4QuasiFreeRatios::theETR=0.; // The Last  p-dependent El/Tot Ratio     (L)
G4int     G4QuasiFreeRatios::lastTZ=0;  // Last atomic number of the target				
G4int     G4QuasiFreeRatios::lastTN=0;  // Last number of neutrons of the target
G4double  G4QuasiFreeRatios::lastPIN=0.;// Last initialized max momentum
std::pair<G4double,G4double>* G4QuasiFreeRatios::lastRST=0;// Last ratio(EQ/QF,QF/IN) Table
G4double* G4QuasiFreeRatios::lastPAR=0; // Parameters for functional calculation					
G4int     G4QuasiFreeRatios::lastPDG=0; // The last PDG code of the projectile
G4int     G4QuasiFreeRatios::lastN=0;   // The last N of calculated nucleus
G4int     G4QuasiFreeRatios::lastZ=0;   // The last Z of calculated nucleus
G4double  G4QuasiFreeRatios::lastP=0.;  // Last used in cross section Momentum
std::pair<G4double,G4double> G4QuasiFreeRatios::lastQIR(0.,0.); // Last ratio(EQ/QF,QF/IN)
G4int     G4QuasiFreeRatios::lastI=0;   // The last position in the DAMDB

// Returns Pointer to the G4VQCrossSection class
G4QuasiFreeRatios* G4QuasiFreeRatios::GetPointer()
{
  static G4QuasiFreeRatios theRatios;   // *** Static body of the QEl Cross Section ***
  return &theRatios;
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
 std::pair<G4double,G4double> G4QuasiFreeRatios::GetRatios(G4double pMom, G4int tgZ,
                                                           G4int tgN, G4int pPDG)
{
  static std::vector <G4int>    colPDG;// Vector of the projectile PDG code
  static std::vector <G4int>    colN;  // Vector of N for calculated nuclei (isotops)
  static std::vector <G4int>    colZ;  // Vector of Z for calculated nuclei (isotops)
  static std::vector <G4double> colP;  // Vector of last momenta for the reaction
  static std::vector <std::pair<G4double,G4double> > colQIR; // Vector of last ratios
  // ***---*** End of the mandatory Static Definitions of the Associative Memory ***---***
#ifdef pdebug
  G4cout<<"G4QElCS::GetCS:>>> f="<<fR<<", p="<<pMom<<", Z="<<tgZ<<"("<<lastZ<<") ,N="<<tgN
        <<"("<<lastN<<"),PDG="<<pPDG<<"("<<lastPDG<<") ,Sz="<<colN.size()<<G4endl;
#endif
  if(!pPDG)
  {
#ifdef pdebug
    G4cout<<"G4QElCS::GetR: *** Found pPDG="<<pPDG<<" ====> CS=0"<<G4endl;
#endif
    return std::make_pair(0.,0.);      // projectile PDG=0 is a mistake (?!) @@
  }
  G4bool in=false;                     // By default the isotope must be found in the AMDB
  if(tgN!=lastN || tgZ!=lastZ || pPDG!=lastPDG)// The nucleus was not the last used isotope
  {
    in = false;                        // By default the isotope haven't be found in AMDB  
    lastP   = 0.;                      // New momentum history (nothing to compare with)
    lastPDG = pPDG;                    // The last PDG of the projectile
    lastN   = tgN;                     // The last N of the calculated nucleus
    lastZ   = tgZ;                     // The last Z of the calculated nucleus
    lastI   = colN.size();             // Size of the Associative Memory DB in the heap
    if(lastI) for(G4int i=0; i<lastI; i++) // Loop over proj/tgZ/tgN lines of DB
	   {                                  // The nucleus with projPDG is found in AMDB
      if(colPDG[i]==pPDG && colN[i]==tgN && colZ[i]==tgZ)
						{
        lastI=i;
#ifdef pdebug
        G4cout<<"G4QElCS::GetCS:*Found* P="<<pMom<<",i="<<i<<G4endl;
#endif
        lastP  =colP [i];               // Last Momentum  (A-dependent)
        lastQIR =colQIR[i];             // Last Ratios (A-dependent)
        if(std::fabs(lastP/pMom-1.)<tolerance)
        {
#ifdef pdebug
          G4cout<<"G4QElCS::GetCS:P="<<pMom<<",R="<<lastQIR<<G4endl;
#endif
          CalculateRatios(-1,i,lastPDG,lastZ,lastN,pMom); // Update param's only
          return lastQIR;      // Use theLastRatios
        }
        in = true;                      // This is the case when the isotop is found in DB
        // Momentum pMom is in IU ! @@ Units
#ifdef pdebug
        G4cout<<"G4QElCS::GR:UpdateDB P="<<pMom<<",f="<<fR<<",I="<<lastI<<",i="<<i<<G4endl;
#endif
        lastQIR=CalculateRatios(-1,i,lastPDG,lastZ,lastN,pMom); // read & update
#ifdef pdebug
        G4cout<<"G4QElCS::GetR: *****> New (inDB) Calculated R="<<lastQIR<<G4endl;
#endif
        break;                           // Go out of the LOOP with found lastI
      }
#ifdef pdebug
      G4cout<<"---G4QElCrossSec::GetCrosSec:pPDG="<<pPDG<<",i="<<i<<",N="<<colN[i]
            <<",Z["<<i<<"]="<<colZ[i]<<",cPDG="<<colPDG[i]<<G4endl;
#endif
	   }
	   if(!in)                            // This nucleus has not been calculated previously
	   {
#ifdef pdebug
      G4cout<<"G4QElCS::GetCrosSec:CalcNew P="<<pMom<<",f="<<fR<<",lastI="<<lastI<<G4endl;
#endif
      //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
      lastQIR=CalculateRatios(0,lastI,lastPDG,lastZ,lastN,pMom);//calculate&create
#ifdef pdebug
      G4cout<<"G4QElCS::GetCrosSec: New R="<<lastQIR<<",lZ="<<lastN<<",lN="<<lastZ<<G4endl;
#endif
      colN.push_back(tgN);
      colZ.push_back(tgZ);
      colPDG.push_back(pPDG);
      colP.push_back(pMom);
      colQIR.push_back(lastQIR);
#ifdef pdebug
      G4cout<<"G4QElCS::GetCS:1st,P="<<pMom<<"(MeV),R="<<lastQIR<<"(mb)"<<G4endl;
#endif
      return lastQIR;
	   } // End of creation of the new set of parameters
    else
				{
#ifdef pdebug
      G4cout<<"G4QElCS::GetCS: Update lastI="<<lastI<<G4endl;
#endif
      colP[lastI]=pMom;
      colPDG[lastI]=pPDG;
      colQIR[lastI]=lastQIR;
    }
  } // End of parameters udate
  else if(std::fabs(lastP/pMom-1.)<tolerance)
  {
#ifdef pdebug
    G4cout<<"G4QElCS::GetCS:OldCur P="<<pMom<<"="<<pMom<<", R="<<lastQIR<<G4endl;
#endif
    return lastQIR;     // Use theLastRatios
  }
  else
  {
#ifdef pdebug
    G4cout<<"G4QElCS::GetCS:UpdatCur P="<<pMom<<",f="<<fR<<",I="<<lastI<<",j="<<j<<G4endl;
#endif
    lastQIR=CalculateRatios(1,lastI,lastPDG,lastZ,lastN,pMom); // Only UpdateDB
    lastP=pMom;
  }
#ifdef pdebug
  G4cout<<"G4QElCS::GetCrSec:End,P="<<pMom<<"(MeV),R="<<lastQIR<<G4endl;
#endif
  return lastQIR;
}

// Calculation of total elastic cross section (p in IU, CS in mb) @@ Units (?)
// F=0 - create AMDB, F=-1 - read&update AMDB, F=1 - update AMDB (sinchro with higher AMDB)
std::pair<G4double,G4double> G4QuasiFreeRatios::CalculateRatios(G4int F,G4int I, G4int PDG,
                                                        G4int tgZ, G4int tgN, G4double pIU)
{
  // *** Begin of Associative Memory DB for acceleration of the cross section calculations
  static std::vector <G4double>  PIN;   // Vector of max initialized log(P) in the table
  static std::vector <std::pair<G4double,G4double>*> RST; // Vector of ratio pai tables
  static std::vector <G4double*> SST;   // Vector of the first squared slope
  static std::vector <G4double*> S1T;   // Vector of the first mantissa
  static std::vector <G4double*> B1T;   // Vector of the first slope
  static std::vector <G4double*> S2T;   // Vector of the secon mantissa
  static std::vector <G4double*> B2T;   // Vector of the second slope
  static std::vector <G4double*> S3T;   // Vector of the third mantissa
  static std::vector <G4double*> B3T;   // Vector of the third slope
  static std::vector <G4double*> S4T;   // Vector of the 4-th mantissa (gloria)
  static std::vector <G4double*> B4T;   // Vector of the 4-th slope    (gloria)
  // *** End of Static Definitions (Associative Memory Data Base) ***
  G4double pMom=pIU/GeV;                // All calculations are in GeV
#ifdef pdebug
  G4cout<<"G4QuasiFreeRatios::CalcCS:->F="<<F<<",p="<<pIU<<G4endl;
#endif
  lastLP=std::log(pMom);                // Make a logarithm of the momentum for calculation
  if(F)                                 // This isotope was found in AMDB =>RETRIEVE/UPDATE
		{
    if(F<0)                             // the AMDB must be loded
    {
      lastPIN = PIN[I];                 // Max log(P) initialised for this table set
      lastRST = RST[I];                 // Pointer to the total sross-section table
#ifdef pdebug
      G4cout<<"G4QElasticCS::CalcCS: DB is updated for I="<<I<<",*,PIN4="<<PIN[4]<<G4endl;
#endif
    }
#ifdef pdebug
    G4cout<<"G4QuasiFreeRatios::CalcCS:*read*, LP="<<lastLP<<",PIN="<<lastPIN<<G4endl;
#endif
    if(lastLP>lastPIN && lastLP<lPMax)
    {
      lastPIN = GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);// Can update upper logP-Limit in tabs
#ifdef pdebug
						G4cout<<"G4QElCS::CalcCS:*updated(I)*,LP="<<lastLP<<"<IN["<<I<<"]="<<lastPIN<<G4endl;
#endif
      PIN[I]=lastPIN;                   // Remember the new P-Limit of the tables
    }
	 }
	 else                                  // This isotope wasn't initialized => CREATE
	 {
    lastRST = new std::pair<G4double,G4double>[nPoints];// Allocate memory for TabulatedRat
    lastRST[0].first=0.;                // A flag that the Table was not initiated
#ifdef pdebug
    G4cout<<"G4QuasiFreeRatios::CalcCS:*ini*,lastLP="<<lastLP<<",min="<<lPMin<<G4endl;
#endif
    lastPIN = GetPTables(lastLP,lPMin,PDG,tgZ,tgN); // Returns the new P-limit for tables
#ifdef pdebug
    G4cout<<"G4QElCS::CalcCS:*i*Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<",LP"<<lastPIN<<G4endl;
#endif
    PIN.push_back(lastPIN);             // Fill parameters of CS function to AMDB
    RST.push_back(lastRST);             // Fill Tabulated CS function to AMDB				
	 } // End of creation/update of the new set of parameters and tables
  // ============= NOW Update (if necessary) and Calculate the Cross Section ===========
#ifdef pdebug
  G4cout<<"G4QElCS::CalcCS:?update?,LP="<<lastLP<<",IN="<<lastPIN<<",ML="<<lPMax<<G4endl;
#endif
  if(lastLP>lastPIN && lastLP<lPMax)
  {
    lastPIN = GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);
#ifdef pdebug
    G4cout<<"G4QElCS::CalcCS: *updated(O)*, LP="<<lastLP<<" < IN="<<lastPIN<<G4endl;
#endif
  }
#ifdef pdebug
  G4cout<<"G4QElastCS::CalcCS: lastLP="<<lastLP<<",lPM="<<lPMin<<",lPIN="<<lastPIN<<G4endl;
#endif
#ifdef pdebug
  G4cout<<"G4QElasticCrosSec::CalcCS: p="<<lastLP<<G4endl;
#endif
  if(lastLP>lPMin && lastLP<=lastPIN)   // Linear fit is made using precalculated tables
  {
    if(std::fabs(lastLP-lastPIN)<.0001) // Just take the highest tabulated value
    {
      G4double shift=(lastLP-lPMin)/dlnP+.000001; // Log distance from lPMin
      G4int    blast=static_cast<int>(shift); // this is a bin number of the lower edge (0)
      if(blast<0 || blast>=nLast) G4cout<<"G4QEleastCS::CCS:b="<<blast<<","<<nLast<<G4endl;
      lastRAT  = lastRST[blast];
#ifdef pdebug
      G4cout<<"G4QuasiFreeRatios::CalculateCS:(E) S1="<<theS1<<", B1="<<theB1<<G4endl;
#endif
				}
    else
    {
      G4double shift=(lastLP-lPMin)/dlnP;        // a shift from the beginning of the table
      G4int    blast=static_cast<int>(shift);    // the lower bin number
      if(blast<0)   blast=0;
      if(blast>=nLast) blast=nLast-1;            // low edge of the last bin
      shift-=blast;                              // step inside the unit bin
      G4int lastL=blast+1;                       // the upper bin number
      G4double QEL=lastRST[blast].first;         // the basic value of Q-el/Q-free rat
      lastRAT.first= QEL+shift*(lastRST[lastL].first-QEL); // calculated Q-el/Q-free rat
      G4double QFL=lastRST[blast].second;        // the basic value of Q-free/inelastic rat
      lastRAT.second= QFL+shift*(lastRST[lastL].second-QFL); // calculated Q-free/inel rat
#ifdef pdebug
      G4cout<<"G4QElCS::CalcCrossSection: Rat="<<lastRAT<<", P="<<pMom<<", Z="<<tgZ<<", N="
            <<tgN<<", PDG="<<PDG<<G4endl;
#endif
#ifdef pdebug
      G4cout<<"G4QuasiFreeRatios::CalculateCS:(I) S1="<<theS1<<", B1="<<theB1<<G4endl;
#endif
    }
  }
  else lastRAT=GetTabValues(lastLP, PDG, tgZ, tgN); // Direct calculation beyond the table
#ifdef pdebug
  G4cout<<"G4QuasiFreeRatios::CalculateCS: END result="<<lastRAT<<G4endl;
#endif
  return lastRAT;
}

// It has parameter sets for all tZ/tN/PDG, using them the tables can be created/updated
G4double G4QuasiFreeRatios::GetPTables(G4double LP,G4double ILP, G4int PDG, G4int tgZ,
                                                                                 G4int tgN)
{
  // ---------> Each parameter set can have not more than nPoints=128 parameters
  //static const G4double lmi=3.5;       // min of (lnP-lmi)^2 parabola (common,@@ global?)
  static const G4int nnn=26;           // #of parameters for nn  (<nPoints=128)
  static const G4int npi=37;           // #of parameters for pin (<nPoints=128)
  static const G4int nkm=18;           // #of parameters for kmn (<nPoints=128)
  static const G4int nkp=16;           // #of parameters for kpn (<nPoints=128)
  static const G4int n_h=13;           // #of parameters for HyperonN (<nPoints=128)
  static const G4int n_a=10;           // #of parameters for antibaryonN (<nPoints=128)
  // p2=p*p;p3=p2*p;lp=log(p);dl=lp-3.5
  //   np     pp
		//  0/ 1 | 14/15: parabola coefficient (e/t)
  //  2/ 3 | 16/17: constant term        (e/t)
  //  4/ 5 | 18/19: Hiperbola            (e/t)
  //  6/ 7 | 20/21: LEcutOfHE            (e/t)
  //  8/ 9 | 22/23: ReversedLETop        (e/t)
  // 10/11 | 24/25: LESlope              (e/t)
  // 12/13 |      : LESuperSlope         (e/t)
  //                       -0- -1-| -2-  -3-|-4--5-|-6- -7-|  -8-    -9- |-10- -11-|12-13|
  static G4double nn[nnn]={.557,.3,6.72,38.2,30.,0.,.54,.54,.00012,.00012,.051,.051,.1,.1,
                           .557,.3,6.72,38.2,32.6,52.7,1.,2.72,.00012,.00012,.2,.2};
  //                       -14--15|-16- -17-|-18- -19-|-20-21-| -22-  -23-  |24-25|
  // --- nn elastic/total
		// E=(#0*dl*dl+#2+#4/p)/(1+#6/p2/p)+1./(#8+p2*(#10+p2*#12);
		// T=(#1*dl*dl+#3+#5/p2)/(1+#7/p2/p2)+1./(#9+p2*(#11+p2*#13));
  // ======================================================================================
  //  hp-N   ap-N (@@)
		//  0/ 1 |  0/ 1: parabola coefficient (e/t)
  //  2/ 3 |  2/ 3: constant term        (e/t)
  //  4/ 5 |  4/ 5: Hiperbola            (e/t)
  //  6/ 7 |  6/ 7: SQRT_Suppress | power(e/t)
  //  8/ 9 |  8/ 9: LEcutOfHE | norm     (e/t)
  // 10/11 |      : Norm0/Slope          (e=t)
  // 12    |      : Suppression          (e=t)
  //                        -0- -1-| -2- -3- |-4- -5- |-6--7-|-8-9-|-10--11-|12|
  static G4double hpn[n_h]={.557,.3,6.72,38.2,99.,900.,2.,27.,2.,3.,.002,.12,1.};
  //                        -0- -1-| -2- -3- |-4- -5-|-6-  -7-|-8-9-|
  static G4double ann[n_a]={.557,.3,6.72,38.2,80.,80.,1.25,.35,1.,.3};

  // --- hyperon-N elastic/total: sp=SQRT(p)
		// E=(#0*dl*dl+#2)/(1+#4/sp+#6/p2/p2)+#8/(dd^2+#11^2)+#12/(dl2^2+#11^2)+[#12/p/sp];
		// T=(#1*dl*dl+#3)/(1+#5/sp+#7/p2/p2)+#9/(dd^2+#11^2)+#13/(dl2^2+#11^2)+[#12/p/sp];
  // --- antibaryon-N elastic/total: sp=SQRT(p)
		// E=#0*dl*dl+#2+#4/(p^#6+#8)
		// T=#1*dl*dl+#3+(#5/p^#7+#8)/p^#7
  // ======================================================================================
  // pim-p pip-p
		//  0/ 1 | 20/21: parabola coefficient (e/t)
  //  2/ 3 | 22/23: constant term        (e/t)
  //  4/ 5 | 24/25: Hiperbola            (e/t)
  //  6/ 7 | 26/27: LEcutOfHE            (e/t)
  //  8/ 9 | 28/29: Norm 1st resonance   (e/t)
  // 10/11 | 30/31: Pos/Wid 1st resonance(e=t)
  // 12/13 | 32/33: Norm 2nd resonance   (e/t)
  // 14/15 | 34/35: Pos/Wid 2nd resonance(e=t)
  // 16/17 | 36   : Norm 3d | Sup 1st    (e/t)
  // 18/19 |      : Pos/Wid 3d resonance (e/t)
  //                        -0- -1-|-2-  -3-|-4--5-|-6-7-| -8- -9-| -10--11-|12-13|-14-15-|
  //                        -16-17-|-18--19-|
  static G4double pin[npi]={.557,.3,2.4,22.3,7.,12.,.7,.4,1.53,4.7,-1.27,.26,.6,1.,-.36,.2,
																												.05,.06,.017,.05,
                            .557,.3,2.4,22.3,6.,5.,3.,1.,13.,13.,-1.27,.26,.7,.8,.32,.24,1.
																											};
  //                         -20-21|-22--23-|24-25|26-27|-28-29-| -30--31-|32-33|-34-35-|36
  // --- pin elastic/total: sp=SQRT(p), dd=lp-#10, dl2=lp-#14, dl3=lp-#18 
		// E=(#0*dl*dl+#2+#4/sp)/(1+#6/p2/p2)+#8/(dd^2*[1+#36*dd^2]+#11^2)+#12/(dl2^2+#11^2)+[3];
		// T=(#1*dl*dl+#3+#5/sp)/(1+#7/p2/p2)+#9/(dd^2*[1+#36*dd^2]+#11^2)+#13/(dl2^2+#11^2)+[3];
  // ======================================================================================
  //  km-N   kp-N (@@)
		//  0/ 1 |  0/ 1: parabola coefficient (e/t)
  //  2/ 3 |  2/ 3: constant term        (e/t)
  //  4/ 5 |  4/ 5: SQRT_Suppression     (e/t)
  //  6/ 7 |  6/ 7: LEcutOfHE            (e/t)
  //  8/ 9 |  8/ 9: Norm 1st resonance   (e/t)
  // 10/11 | 10/11: Pos/Wid 1st resonance(e=t)
  // 12/13 | 12/13: Norm 2nd resonance   (e/t)
  // 14/15 | 14/15: Pos/Wid 2nd resonance(e=t)
  // 16/17 |      : LAMBDA_Hyperbola     (e/t)
  //                        -0- -1-| -2- -3- |-4- -5- | -6- -7-| -8- -9- |-10--11-|-12-13|
  //                        -14-15-|-16-17-|
  static G4double kmn[nkm]={.557,.3,2.23,19.5,-.7,-.21,.075,.52,.004,.006,.39,.0125,.14,.3,
                            1.,.125,5.2,14.};
  //                        -0- -1-| -2- -3- |-4- -5-|-6--7-|-8--9-|10-11-|12-13|-14-15-|
  static G4double kpn[nkp]={.557,.3,2.23,19.5,-.7,.46,.1,1.6,2.,2.6,1.,.61,.7,.7,.38,.26};

  // --- pin elastic/total: sp=SQRT(p), dd=lp-#10, dl2=lp-#14
		// E=(#0*dl*dl+#2)/(1+#4/sp+#6/p2/p2)+#8/(dd^2+#11^2)+#12/(dl2^2+#11^2)+[#12/p/sp];
		// T=(#1*dl*dl+#3)/(1+#5/sp+#7/p2/p2)+#9/(dd^2+#11^2)+#13/(dl2^2+#11^2)+[#12/p/sp];
  // ======================================================================================
  G4bool kfl=true;                             // Flag of K0/aK0 oscillation
  G4bool kf=false;
  if(PDG==130||PDG==310)
  {
    kf=true;
    if(G4UniformRand()>.5) kfl=false;
  }
  if     (PDG==2212 || PDG==2112) lastPAR=nn;  // np/pp parameters
  else if(PDG== 211 || PDG==-112) lastPAR=pin; // pipp/pimp parameters
		else if(PDG== 321 || PDG== 311 || (kf && kfl)) lastPAR=kpn; // KN parameters
  else if(PDG==-321 || PDG==-311 || (kf &&!kfl)) lastPAR=kmn; // antiKN parameters
  else if(PDG> 3000 && PDG< 3335) lastPAR=hpn; // @@ for all hyperons - take Lambda
  else if(PDG<-2000 && PDG>-3335) lastPAR=ann; // @@ for all anti-baryons - anti-proton/an
  else
  {
    G4cout<<"*Error*G4QuasiFreeRatios::GetPTables: PDG="<<PDG
          <<", while it is defined only for p,n,hyperons,anti-baryons,pi,K/antiK"<<G4endl;
    throw G4QException("G4QuasiFreeRatios::GetPTables: projectile is not implemented");
  }
    if(lastRST[0].first==0) // A unique flag to avoid the repeatable initialization
    {
      // and initialize the zero element of the table
      G4double lp=lPMin;                                      // ln(momentum)
      lastRST[0]=GetTabValues(lp, PDG, tgZ, tgN);             // Calculate AMDB tables
#ifdef pdebug
      G4cout<<"G4QuasiFreeRatios::GetPTables:ip=0(init), lp="<<lp<<",QI="<<theQFP
            <<",B3="<<theB3<<",S4="<<theS4<<",B4="<<theB4<<G4endl;
#endif
    }
    if(LP>ILP)
				{
      G4int ini = static_cast<int>((ILP-lPMin+.000001)/dlnP)+1; // already inited till this
      if(ini<0) ini=0;
      if(ini<nPoints)
      {
        G4int fin = static_cast<int>((LP-lPMin)/dlnP)+1; // final bin of initialization
        if(fin>nPoints) fin=nLast;                // Limit of the tabular initialization
        if(fin>=ini)
        {
          G4double lp=0.;
          for(G4int ip=ini; ip<=fin; ip++)        // Calculate tabular CS,S1,B1,S2,B2,S3,B3
										{
            lp=lPMin+ip*dlnP;                     // ln(momentum)
            lastRST[ip]=GetTabValues(lp, PDG, tgZ, tgN); // Calculate AMDB tables (ret CS)
#ifdef pdebug
            G4cout<<"G4QuasiFreeRatios::GetPTables:ip="<<ip<<",lp="<<lp<<",S1="<<theS1
                  <<",B1="<<theB1<<",S2="<<theS2<<",B2="<<theB2<<",S3="
                  <<theS3<<",B3="<<theB3<<",S4="<<theS4<<",B4="<<theB4<<G4endl;
#endif
          }
          return lp;
        }
        else G4cout<<"*Warning*G4QuasiFreeRatios::GetPTables: PDG="<<PDG<<", Z="<<tgZ
                   <<", N="<<tgN<<", i="<<ini<<" > fin="<<fin<<", LP="<<LP<<" > ILP="<<ILP
                   <<" nothing is done!"<<G4endl;
      }
      else G4cout<<"*Warning*G4QuasiFreeRatios::GetPTables: PDG="<<PDG<<", Z="<<tgZ
                 <<", N="<<tgN<<", i="<<ini<<">= max="<<nPoints<<", LP="<<LP<<" > ILP="
                 <<ILP<<", lPMax="<<lPMax<<" nothing is done!"<<G4endl;
    }
#ifdef pdebug
    else G4cout<<"*Warning*G4QuasiFreeRatios::GetPTables: PDG="<<PDG<<", Z="<<tgZ
               <<", N="<<tgN<<", LP="<<LP<<" <= ILP="<<ILP<<" nothing is done!"<<G4endl;
#endif
  return ILP;
}

// lastLP is used, so calculating tables, one need to remember and then recover lastLP
std::pair<G4double,G4double> G4QuasiFreeRatios::GetTabValues(G4double lp, G4int PDG,
                                                                   G4int tgZ,G4int tgN)
{
  //static const G4double GeVSQ=gigaelectronvolt*gigaelectronvolt;
  if(tgZ<1 || tgZ>92)
  {
    G4cout<<"*Warning*G4QElasticCS::GetTabValue: (1-92) No isotopes for Z="<<tgZ<<G4endl;
    return std::make_pair(0.,0.);
  }
#ifdef pdebug
  G4cout<<"G4QElasticCS::GetTabVal: lp="<<lp<<",Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<G4endl;
#endif
  G4double p=std::exp(lp);              // momentum
  G4double sp=std::sqrt(p);             // sqrt(p)
  if(PDG==2112 && tgZ==1 && tgN==0)       // np
  {
    theQFP=lastPAR[1]*p;
    theETR=lastPAR[2]*sp;
#ifdef tdebug
    G4cout<<"G4QElasticCS::GetTableValues:(np) S1="<<theS1<<",B1="<<theB1
          <<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS1<<",B3="<<theB1<<G4endl;
#endif
    // Returns the total elastic pp cross-section (to avoid spoiling lastRAT)
    G4double RQEQF=0.;
    G4double RQFIN=0.;
    return std::make_pair(RQEQF,RQFIN);

  }
  else if(PDG==2212 && tgZ==1 && tgN==0) // pp
  {
  }
  else if(PDG==2212 || PDG==2112) // n/p+A (isotope/projectile invariant)
  {
    // Returns the ratios (to avoid spoiling lastRAT)
#ifdef tdebug
    G4cout<<"G4QElCS::GetTabV: PDG="<<PDG<<",P="<<p<<",N="<<tgN<<",Z="<<tgZ<<G4endl;
#endif
    G4double RQEQF=0.;
    G4double RQFIN=0.;
    return std::make_pair(RQEQF,RQFIN);
  }
  else
  {
    G4cout<<"*Error*G4QuasiFreeRatios::GetTabValues: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for PDG=2212, Z=1, N=0"<<G4endl;
    throw G4QException("G4QuasiFreeRatios::GetTabValues: only pp is implemented");
  }
  return std::make_pair(0.,0.);
} // End of GetTableValues

// Calculatio QasiFree/Inelastic Ratio as a function of total hN cross-section and A
G4double G4QuasiFreeRatios::QF2IN_Ratio(G4double s, G4double A)
{
  static const G4double C=1.246;
		G4double s2=s*s;
  G4double s4=s2*s2;
		G4double ss=std::sqrt(std::sqrt(s));
  G4double P=7.48e-5*s2/(1.+8.77e12/s4/s4/s2);
  G4double E=.2644+.016/(1.+std::exp((29.54-s)/2.49));
  G4double F=ss*.1526*std::exp(-s2*ss*.0000859);
	 return C*std::exp(-E*pow(A-1.,F))/std::pow(A,P);
} // End of QF2IN_Ratio

// Calculates Mean Elastic and Total Cross-Sections (mb) for PDG+(Z,N) at P=p[GeV/c]
std::pair<G4double,G4double> G4QuasiFreeRatios::GetElTot(G4double p, G4int PDG,
                                                         G4int Z, G4int N, G4bool F)
{
  G4double A=0;
  if(F) A=Z+N;
  return std::make_pair(Z+PDG,N+p);
} // End of GetHpHnTotal
