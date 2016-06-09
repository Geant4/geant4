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
//
// $Id$
//
//
// G4 Physics class: G4QNuNuNuclearCrossSection for (nu,nu)A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-Oct-03
// 
// ****************************************************************************************
// ********** This CLASS is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************
// -----------------------------------------------------------
// Short description: Neutral Current nu -> nu nuclear XS
// -----------------------------------------------------------

//#define debug
//#define edebug
//#define pdebug
//#define ppdebug
//#define tdebug
//#define sdebug

#include "G4QNuNuNuclearCrossSection.hh"
#include "G4SystemOfUnits.hh"

// Initialization of the
G4bool    G4QNuNuNuclearCrossSection::onlyCS=true;// Flag to calculate only CS (not QE)
G4double  G4QNuNuNuclearCrossSection::lastSig=0.;// Last calculated total cross section
G4double  G4QNuNuNuclearCrossSection::lastQEL=0.;// Last calculated quasi-el. cross section
G4int     G4QNuNuNuclearCrossSection::lastL=0;   // Last used in cross section TheLastBin
G4double  G4QNuNuNuclearCrossSection::lastE=0.;  // Last used in cross section TheEnergy
G4double* G4QNuNuNuclearCrossSection::lastEN=0;  // Pointer to the Energy Scale of TX & QE
G4double* G4QNuNuNuclearCrossSection::lastTX=0;  // Pointer to the LastArray of TX function
G4double* G4QNuNuNuclearCrossSection::lastQE=0;  // Pointer to the LastArray of QE function
G4int     G4QNuNuNuclearCrossSection::lastPDG=0; // The last PDG code of the projectile
G4int     G4QNuNuNuclearCrossSection::lastN=0;   // The last N of calculated nucleus
G4int     G4QNuNuNuclearCrossSection::lastZ=0;   // The last Z of calculated nucleus
G4double  G4QNuNuNuclearCrossSection::lastP=0.;  // Last used in cross section Momentum
G4double  G4QNuNuNuclearCrossSection::lastTH=0.; // Last threshold momentum
G4double  G4QNuNuNuclearCrossSection::lastCS=0.; // Last value of the Cross Section
G4int     G4QNuNuNuclearCrossSection::lastI=0;   // The last position in the DAMDB
std::vector<G4double*>* G4QNuNuNuclearCrossSection::TX = new std::vector<G4double*>;
std::vector<G4double*>* G4QNuNuNuclearCrossSection::QE = new std::vector<G4double*>;

// Returns Pointer to the G4VQCrossSection class
G4VQCrossSection* G4QNuNuNuclearCrossSection::GetPointer()
{
  static G4QNuNuNuclearCrossSection theCrossSection; //**Static body of the Cross Section**
  return &theCrossSection;
}

G4QNuNuNuclearCrossSection::~G4QNuNuNuclearCrossSection()
{
  G4int lens=TX->size();
  for(G4int i=0; i<lens; ++i) delete[] (*TX)[i];
  delete TX;
  G4int hens=QE->size();
  for(G4int i=0; i<hens; ++i) delete[] (*QE)[i];
  delete QE;
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4QNuNuNuclearCrossSection::GetCrossSection(G4bool fCS, G4double pMom,
                                                     G4int tgZ, G4int tgN, G4int pPDG)
{
  static G4int j;                      // A#0f records found in DB for this projectile
  static std::vector <G4int>    colPDG;// Vector of the projectile PDG code
  static std::vector <G4int>    colN;  // Vector of N for calculated nuclei (isotops)
  static std::vector <G4int>    colZ;  // Vector of Z for calculated nuclei (isotops)
  static std::vector <G4double> colP;  // Vector of last momenta for the reaction
  static std::vector <G4double> colTH; // Vector of energy thresholds for the reaction
  static std::vector <G4double> colCS; // Vector of last cross sections for the reaction
  // ***---*** End of the mandatory Static Definitions of the Associative Memory ***---***
  G4double pEn=pMom;
#ifdef debug
  G4cout<<"G4QNMNCS::GetCS:>> f="<<fCS<<", p="<<pMom<<", Z="<<tgZ<<"("<<lastZ<<") ,N="<<tgN
        <<"("<<lastN<<"),PDG="<<pPDG<<"("<<lastPDG<<"), T="<<pEn<<"("<<lastTH<<")"<<",Sz="
        <<colN.size()<<G4endl;
  //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  if(pPDG!=14)
  {
#ifdef pdebug
    G4cout<<"G4QNMNCS::GetCS: *** Found pPDG="<<pPDG<<" =--=> CS=0"<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return 0.;                         // projectile PDG=0 is a mistake (?!) @@
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
    j  = 0;                            // A#0f records found in DB for this projectile
    if(lastI) for(G4int i=0; i<lastI; i++) if(colPDG[i]==pPDG) // The partType is found
    {                                  // The nucleus with projPDG is found in AMDB
      if(colN[i]==tgN && colZ[i]==tgZ)
      {
        lastI=i;
        lastTH =colTH[i];                // Last THreshold (A-dependent)
#ifdef pdebug
        G4cout<<"G4QNMNCS::GetCS:*Found*P="<<pMom<<",Threshold="<<lastTH<<",j="<<j<<G4endl;
        //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(pEn<=lastTH)
        {
#ifdef pdebug
          G4cout<<"G4QNMNCS::GetCS:Found T="<<pEn<<" < Threshold="<<lastTH<<",X=0"<<G4endl;
          //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
          return 0.;                     // Energy is below the Threshold value
        }
        lastP  =colP [i];                // Last Momentum  (A-dependent)
        lastCS =colCS[i];                // Last CrossSect (A-dependent)
        if(std::fabs(lastP/pMom-1.)<tolerance)
        {
#ifdef pdebug
          G4cout<<"G4QNMNCS::GetCS:P="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
          //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
          return lastCS*millibarn;     // Use theLastCS
        }
        in = true;                       // This is the case when the isotop is found in DB
        // Momentum pMom is in IU ! @@ Units
#ifdef pdebug
        G4cout<<"G4QNMNCS::G:UpdaDB P="<<pMom<<",f="<<fCS<<",lI="<<lastI<<",j="<<j<<G4endl;
#endif
        lastCS=CalculateCrossSection(fCS,-1,j,lastPDG,lastZ,lastN,pMom); // read & update
#ifdef pdebug
        G4cout<<"G4QNMNCS::GetCrosSec: *****> New (inDB) Calculated CS="<<lastCS<<G4endl;
        //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(lastCS<=0. && pEn>lastTH)    // Correct the threshold
        {
#ifdef pdebug
          G4cout<<"G4QNMNCS::GetCS: New T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
        break;                           // Go out of the LOOP
      }
#ifdef pdebug
      G4cout<<"---G4QNMNCrossSec::GetCrosSec:pPDG="<<pPDG<<",j="<<j<<",N="<<colN[i]
            <<",Z["<<i<<"]="<<colZ[i]<<",cPDG="<<colPDG[i]<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      j++;                             // Increment a#0f records found in DB for this pPDG
    }
    if(!in)                            // This nucleus has not been calculated previously
    {
#ifdef pdebug
      G4cout<<"G4QNMNCS::GetCrSec: CalcNew P="<<pMom<<",f="<<fCS<<",lastI="<<lastI<<G4endl;
#endif
      //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
      lastCS=CalculateCrossSection(fCS,0,j,lastPDG,lastZ,lastN,pMom); //calculate & create
      if(lastCS<=0.)
      {
        lastTH = ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last
#ifdef pdebug
        G4cout<<"G4QNMNCrossSection::GetCrossSect:NewThresh="<<lastTH<<",T="<<pEn<<G4endl;
#endif
        if(pEn>lastTH)
        {
#ifdef pdebug
          G4cout<<"G4QNMNCS::GetCS: First T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
      }
#ifdef pdebug
      G4cout<<"G4QNMNCS::GetCrosSec:New CS="<<lastCS<<",lZ="<<lastN<<",lN="<<lastZ<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      colN.push_back(tgN);
      colZ.push_back(tgZ);
      colPDG.push_back(pPDG);
      colP.push_back(pMom);
      colTH.push_back(lastTH);
      colCS.push_back(lastCS);
#ifdef pdebug
      G4cout<<"G4QNMNCS::GetCS:1st,P="<<pMom<<"(MeV),X="<<lastCS*millibarn<<"(mb)"<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      return lastCS*millibarn;
    } // End of creation of the new set of parameters
    else
    {
#ifdef pdebug
      G4cout<<"G4QNMNCS::GetCS: Update lastI="<<lastI<<",j="<<j<<G4endl;
#endif
      colP[lastI]=pMom;
      colPDG[lastI]=pPDG;
      colCS[lastI]=lastCS;
    }
  } // End of parameters udate
  else if(pEn<=lastTH)
  {
#ifdef pdebug
    G4cout<<"G4QNMNCS::GetCS: Current T="<<pEn<<" < Threshold="<<lastTH<<", CS=0"<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return 0.;                         // Momentum is below the Threshold Value -> CS=0
  }
  else if(std::fabs(lastP/pMom-1.)<tolerance)
  {
#ifdef pdebug
    G4cout<<"G4QNMNCS::GetCS:OldCur P="<<pMom<<"="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return lastCS*millibarn;     // Use theLastCS
  }
  else
  {
#ifdef pdebug
    G4cout<<"G4QNMNCS::GetCS:UpdaCur P="<<pMom<<",f="<<fCS<<",I="<<lastI<<",j="<<j<<G4endl;
#endif
    lastCS=CalculateCrossSection(fCS,1,j,lastPDG,lastZ,lastN,pMom); // Only UpdateDB
    lastP=pMom;
  }
#ifdef pdebug
  G4cout<<"G4QNMNCS::GetCrSec:End,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
  //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  return lastCS*millibarn;
}

// Gives the threshold energy = the same for all nuclei (@@ can be reduced for hevy nuclei)
G4double G4QNuNuNuclearCrossSection::ThresholdEnergy(G4int, G4int, G4int) {return 0.;}

// The main member function giving the gamma-A cross section (E_kin in MeV, CS in mb)
G4double G4QNuNuNuclearCrossSection::CalculateCrossSection(G4bool CS, G4int F, G4int I,
                                        G4int, G4int targZ, G4int targN, G4double Momentum)
{
  static const G4double mb38=1.E-11; // Conversion 10^-38 cm^2 to mb=10^-27 cm^2
  static const G4int nE=65;   // !! If change this, change it in GetFunctions() (*.hh) !!
  static const G4int mL=nE-1;
  static const G4double EMi=0.;      // Universal threshold of the reaction in GeV
  static const G4double EMa=300.;    // Maximum tabulated Energy of nu_mu in GeV 
  // *** Begin of the Associative memory for acceleration of the cross section calculations
  static std::vector <G4double> colH;//?? Vector of HighEnergyCoefficients (functional)
  static G4bool first=true;          // Flag of initialization of the energy axis
  // *** End of Static Definitions (Associative Memory) ***
  //const G4double Energy = aPart->GetKineticEnergy()/MeV; // Energy of the Muon
  //G4double TotEnergy2=Momentum;
  onlyCS=CS;                         // Flag to calculate only CS (not TX & QE)
  lastE=Momentum/GeV;                // Kinetic energy of the muon neutrino (in GeV!)
  if (lastE<=EMi)                    // Energy is below the minimum energy in the table
  {
    lastE=0.;
    lastSig=0.;
    return 0.;
  }
  G4int Z=targZ;                     // New Z, which can change the sign
  if(F<=0)                           // This isotope was not the last used isotop
  {
    if(F<0)                          // This isotope was found in DAMDB =-------=> RETRIEVE
    {
      lastTX =(*TX)[I];              // Pointer to the prepared TX function (same isotope)
      lastQE =(*QE)[I];              // Pointer to the prepared QE function (same isotope)
   }
   else                              // This isotope wasn't calculated previously => CREATE
   {
      if(first)
      {
        lastEN = new G4double[nE];   // This must be done only once!
        Z=-Z;                        // To explain GetFunctions that E-axis must be filled
        first=false;                 // To make it only once
      }
      lastTX = new G4double[nE];     // Allocate memory for the new TX function
      lastQE = new G4double[nE];     // Allocate memory for the new QE function
      G4int res=GetFunctions(Z,targN,lastTX,lastQE,lastEN);//@@analize(0=first,-1=bad,1=OK)
      if(res<0) G4cerr<<"*W*G4NuMuNuclearCS::CalcCrossSect:Bad Function Retrieve"<<G4endl;
      // *** The synchronization check ***
      G4int sync=TX->size();
      if(sync!=I) G4cerr<<"***G4NuMuNuclearCS::CalcCrossSect:Sync.="<<sync<<"#"<<I<<G4endl;
      TX->push_back(lastTX);
      QE->push_back(lastQE);
    } // End of creation of the new set of parameters
  } // End of parameters udate
  // =-------------= NOW Calculate the Cross Section =---------------------=
  if (lastE<=EMi)                   // Check that the neutrinoEnergy is higher than ThreshE
  {
    lastE=0.;
    lastSig=0.;
    return 0.;
  }
  if(lastE<EMa) // Linear fit is made explicitly to fix the last bin for the randomization
  {
    G4int chk=1;
    G4int ran=mL/2;
    G4int sep=ran;  // as a result = an index of the left edge of the interval
    while(ran>=2)
    {
      G4int newran=ran/2;
      if(lastE<=lastEN[sep]) sep-=newran;
      else                   sep+=newran;
      ran=newran;
      chk=chk+chk; 
    }
    if(chk+chk!=mL) G4cerr<<"*Warn*G4NuMuNuclearCS::CalcCS:Table! mL="<<mL<<G4endl;
    G4double lowE=lastEN[sep];
    G4double highE=lastEN[sep+1];
    G4double lowTX=lastTX[sep];
    if(lastE<lowE||sep>=mL||lastE>highE)
      G4cerr<<"*Warn*G4NuMuNuclearCS::CalcCS:Bin! "<<lowE<<" < "<<lastE<<" < "<<highE
            <<", sep="<<sep<<", mL="<<mL<<G4endl;
    lastSig=lastE*(lastE-lowE)*(lastTX[sep+1]-lowTX)/(highE-lowE)+lowTX; // Recover *E
    if(!onlyCS)                       // Skip the differential cross-section parameters
    {
      G4double lowQE=lastQE[sep];
      lastQEL=(lastE-lowE)*(lastQE[sep+1]-lowQE)/(highE-lowE)+lowQE;
#ifdef pdebug
      G4cout<<"G4NuMuNuclearCS::CalcCS: T="<<lastSig<<",Q="<<lastQEL<<",E="<<lastE<<G4endl;
#endif
    }
  }
  else
  {
    lastSig=lastTX[mL]; // @@ No extrapolation, just a const, while it looks shrinking...
    lastQEL=lastQE[mL];
  }
  if(lastQEL<0.) lastQEL = 0.;
  if(lastSig<0.) lastSig = 0.;
  // The cross-sections are expected to be in mb
  lastSig*=mb38;
  if(!onlyCS) lastQEL*=mb38;
  return lastSig;
}

// Calculate the cros-section functions
// ****************************************************************************************
// *** This tables are the same for all lepto-nuclear reactions, only mass is different ***
// ***@@ IT'S REASONABLE TO MAKE ADDiTIONAL VIRTUAL CLASS FOR LEPTO-NUCLEAR WITH THIS@@ ***
// ****************************************************************************************
G4int G4QNuNuNuclearCrossSection::GetFunctions(G4int z, G4int n,
                                                     G4double* t, G4double* q, G4double* e)
{
  static const G4int nE=65; // !! If change this, change it in GetCrossSection() (*.cc) !!
  static const G4double nuEn[nE]={0.,
  1.00463e-5,1.05336e-5,1.10692e-5,1.16592e-5,1.23109e-5,1.30323e-5,1.38331e-5,1.47245e-5,
  1.57194e-5,1.68335e-5,1.80848e-5,1.94948e-5,2.10894e-5,2.28991e-5,2.49608e-5,2.73189e-5,
  3.00273e-5,3.31516e-5,3.67722e-5,4.09881e-5,4.59217e-5,5.17255e-5,5.85908e-5,6.67583e-5,
  7.65338e-5,8.83078e-5,.000102583,.000120011,.000141441,.000167995,.000201160,.000242926,
  .000295985,.000364008,.000452051,.000567152,.000719210,.000922307,.001196710,.001571930,
  .002091530,.002820590,.003857810,.005354930,.007548840,.010815300,.015760100,.023376900,
  .035325600,.054430800,.085595700,.137508000,.225898000,.379892000,.654712000,1.15767000,
  2.10277000,3.92843000,7.55861000,14.9991000,30.7412000,65.1734000,143.155000,326.326000};
  static const G4double TOTX[nE]={0.,
  3.18319e-5,3.33759e-5,3.50729e-5,3.69425e-5,3.90071e-5,4.12928e-5,4.38300e-5,4.66541e-5,
  4.98065e-5,5.33360e-5,5.73004e-5,6.17678e-5,6.68196e-5,7.25529e-5,7.90845e-5,8.65551e-5,
  9.51355e-5,.000105033,.000116503,.000129858,.000145486,.000163871,.000185616,.000211484,
  .000242444,.000279730,.000324932,.000380110,.000447953,.000531999,.000636946,.000769076,
  .000936873,.001151900,.001430050,.001793410,.002272980,.002912690,.003775490,.004952530,
  .006577150,.008846400,.012054100,.017404800,.027443300,.040285600,.058251500,.084054600,
  .121309000,.173273000,.240200000,.312865000,.367478000,.379071000,.348461000,.301143000,
  .259802000,.231449000,.214007000,.203827000,.198265000,.195186000,.192821000,.190555000};
  static const G4double QELX[nE]={0.,
  .319793e-9,.351567e-9,.388228e-9,.430721e-9,.480212e-9,.538141e-9,.606306e-9,.686956e-9,
  .782930e-9,.897832e-9,1.03626e-9,1.20415e-9,1.40918e-9,1.66139e-9,1.97401e-9,2.36458e-9,
  2.85666e-9,3.48202e-9,4.28408e-9,5.32264e-9,6.68098e-9,8.47630e-9,1.08754e-8,1.41183e-8,
  1.85552e-8,2.47023e-8,3.33325e-8,4.56172e-8,6.33590e-8,8.93735e-8,1.28128e-7,1.86829e-7,
  2.77301e-7,4.19300e-7,6.46455e-7,1.01714e-6,1.63475e-6,2.68640e-6,4.51816e-6,7.78503e-6,
  1.37563e-5,2.49521e-5,4.65025e-5,9.32015e-5,.000207165,.000435700,.000918046,.001964940,
  .004285300,.009431390,.020560100,.043021600,.083012600,.142280000,.194239000,.222213000,
  .231903000,.234030000,.234691000,.235204000,.235575000,.235793000,.235902000,.235946000};
  // --------------------------------
  G4int first=0;
  if(z<0.)
  {
    first=1;
    z=-z;
  }
  if(z<1 || z>92)             // neutron and plutonium are forbidden
  {
    G4cout<<"***G4QNuNuNuclearCrossSection::GetFunctions:Z="<<z<<".No CS returned"<<G4endl;
    return -1;
  }
  for(G4int k=0; k<nE; k++)
  {
    G4double a=n+z;
    G4double na=n+a;
    G4double dn=n+n;
    G4double da=a+a;
    G4double ta=da+a;
    if(first) e[k]=nuEn[k];       // Energy of neutrino E (first bin k=0 can be modified)
    t[k]=TOTX[k]*nuEn[k]*(na+na)/ta+QELX[k]*(dn+dn-da)/ta; // TotalCrossSection
    q[k]=QELX[k]*dn/a;                                     // QuasiElasticCrossSection
  }
  return first;
}

// Randomize Q2 from neutrino to the scattered muon when the scattering is quasi-elastic
G4double G4QNuNuNuclearCrossSection::GetQEL_ExchangeQ2()
{
  static const double MN=.931494043;   // Nucleon mass (inside nucleus, atomicMassUnit,GeV)
  static const G4double power=-3.5;    // direct power for the magic variable
  static const G4double pconv=1./power;// conversion power for the magic variable
  static const G4int nQ2=101;          // #Of point in the Q2l table (in GeV^2)
  static const G4int lQ2=nQ2-1;        // index of the last in the Q2l table
  static const G4int bQ2=lQ2-1;        // index of the before last in the Q2 ltable
  // Reversed table
  static const G4double Xl[nQ2]={1.87905e-10,
 .005231, .010602, .016192, .022038, .028146, .034513, .041130, .047986, .055071, .062374,
 .069883, .077587, .085475, .093539, .101766, .110150, .118680, .127348, .136147, .145069,
 .154107, .163255, .172506, .181855, .191296, .200825, .210435, .220124, .229886, .239718,
 .249617, .259578, .269598, .279675, .289805, .299986, .310215, .320490, .330808, .341169,
 .351568, .362006, .372479, .382987, .393527, .404099, .414700, .425330, .435987, .446670,
 .457379, .468111, .478866, .489643, .500441, .511260, .522097, .532954, .543828, .554720,
 .565628, .576553, .587492, .598447, .609416, .620398, .631394, .642403, .653424, .664457,
 .675502, .686557, .697624, .708701, .719788, .730886, .741992, .753108, .764233, .775366,
 .786508, .797658, .808816, .819982, .831155, .842336, .853524, .864718, .875920, .887128,
 .898342, .909563, .920790, .932023, .943261, .954506, .965755, .977011, .988271, .999539};
  // Direct table
  static const G4double Xmax=Xl[lQ2];
  static const G4double Xmin=Xl[0];
  static const G4double dX=(Xmax-Xmin)/lQ2;  // step in X(Q2, GeV^2)
  static const G4double inl[nQ2]={0,
 1.88843, 3.65455, 5.29282, 6.82878, 8.28390, 9.67403, 11.0109, 12.3034, 13.5583, 14.7811,
 15.9760, 17.1466, 18.2958, 19.4260, 20.5392, 21.6372, 22.7215, 23.7933, 24.8538, 25.9039,
 26.9446, 27.9766, 29.0006, 30.0171, 31.0268, 32.0301, 33.0274, 34.0192, 35.0058, 35.9876,
 36.9649, 37.9379, 38.9069, 39.8721, 40.8337, 41.7920, 42.7471, 43.6992, 44.6484, 45.5950,
 46.5390, 47.4805, 48.4197, 49.3567, 50.2916, 51.2245, 52.1554, 53.0846, 54.0120, 54.9377,
 55.8617, 56.7843, 57.7054, 58.6250, 59.5433, 60.4603, 61.3761, 62.2906, 63.2040, 64.1162,
 65.0274, 65.9375, 66.8467, 67.7548, 68.6621, 69.5684, 70.4738, 71.3784, 72.2822, 73.1852,
 74.0875, 74.9889, 75.8897, 76.7898, 77.6892, 78.5879, 79.4860, 80.3835, 81.2804, 82.1767,
 83.0724, 83.9676, 84.8622, 85.7563, 86.6499, 87.5430, 88.4356, 89.3277, 90.2194, 91.1106,
 92.0013, 92.8917, 93.7816, 94.6711, 95.5602, 96.4489, 97.3372, 98.2252, 99.1128, 100.000};
  G4double Enu=lastE;                 // Get energy of the last calculated cross-section
  G4double dEnu=Enu+Enu;              // doubled energy of nu/anu
  G4double Enu2=Enu*Enu;              // squared energy of nu/anu
  G4double ME=Enu*MN;                 // M*E
  G4double dME=ME+ME;                 // 2*M*E
  G4double dEMN=(dEnu+MN)*ME;
  G4double sqE=Enu*ME;
  G4double E2M=MN*Enu2;
  G4double ymax=(E2M+sqE)/dEMN;
  G4double Q2mi=0; // Q2_min(E_nu)
  G4double Q2ma=dME*ymax;                                                  // Q2_max(E_nu)
  G4double Xma=std::pow((1.+Q2mi),power);  // X_max(E_nu)
  G4double Xmi=std::pow((1.+Q2ma),power);  // X_min(E_nu)
  // Find the integral values integ(Xmi) & integ(Xma) using the direct table
  G4double rXi=(Xmi-Xmin)/dX;
  G4int    iXi=static_cast<int>(rXi);
  if(iXi<0) iXi=0;
  if(iXi>bQ2) iXi=bQ2;
  G4double dXi=rXi-iXi;
  G4double bnti=inl[iXi];
  G4double inti=bnti+dXi*(inl[iXi+1]-bnti);
  //
  G4double rXa=(Xma-Xmin)/dX;
  G4int    iXa=static_cast<int>(rXa);
  if(iXa<0) iXa=0;
  if(iXa>bQ2) iXa=bQ2;
  G4double dXa=rXa-iXa;
  G4double bnta=inl[iXa];
  G4double inta=bnta+dXa*(inl[iXa+1]-bnta);
  // *** Find X using the reversed table ***
  G4double intx=inti+(inta-inti)*G4UniformRand();
  G4int    intc=static_cast<int>(intx);
  if(intc<0) intc=0;
  if(intc>bQ2) intc=bQ2;         // If it is more than max, then the BAD extrapolation
  G4double dint=intx-intc;
  G4double mX=Xl[intc];
  G4double X=mX+dint*(Xl[intc+1]-mX);
  G4double Q2=std::pow(X,pconv)-1.;
  return Q2*GeV*GeV;
}

// Randomize Q2 from neutrino to the scattered muon when the scattering is not quasiElastic
G4double G4QNuNuNuclearCrossSection::GetNQE_ExchangeQ2()
{
  static const double mpi=.13957018;    // charged pi meson mass in GeV
  static const double MN=.931494043;    // Nucleon mass (inside nucleus,atomicMassUnit,GeV)
  static const double dMN=MN+MN;        // 2*M_N in GeV
  static const double mcV=(dMN+mpi)*mpi;// constant of W>M+mc cut for Quasi-Elastic
  static const G4int power=7;           // direct power for the magic variable
  static const G4double pconv=1./power; // conversion power for the magic variable
  static const G4int nX=21;             // #Of point in the Xl table (in GeV^2)
  static const G4int lX=nX-1;           // index of the last in the Xl table
  static const G4int bX=lX-1;           // @@ index of the before last in the Xl table
  static const G4int nE=20;             // #Of point in the El table (in GeV^2)
  static const G4int bE=nE-1;           // index of the last in the El table
  static const G4int pE=bE-1;           // index of the before last in the El table
  // Reversed table
  static const G4double X0[nX]={6.14081e-05,
 .413394, .644455, .843199, 1.02623, 1.20032, 1.36916, 1.53516, 1.70008, 1.86539, 2.03244,
 2.20256, 2.37723, 2.55818, 2.74762, 2.94857, 3.16550, 3.40582, 3.68379, 4.03589, 4.77419};
  static const G4double X1[nX]={.00125268,
 .861178, 1.34230, 1.75605, 2.13704, 2.49936, 2.85072, 3.19611, 3.53921, 3.88308, 4.23049,
 4.58423, 4.94735, 5.32342, 5.71700, 6.13428, 6.58447, 7.08267, 7.65782, 8.38299, 9.77330};
  static const G4double X2[nX]={.015694,
 1.97690, 3.07976, 4.02770, 4.90021, 5.72963, 6.53363, 7.32363, 8.10805, 8.89384, 9.68728,
 10.4947, 11.3228, 12.1797, 13.0753, 14.0234, 15.0439, 16.1692, 17.4599, 19.0626, 21.7276};
  static const G4double X3[nX]={.0866877,
 4.03498, 6.27651, 8.20056, 9.96931, 11.6487, 13.2747, 14.8704, 16.4526, 18.0351, 19.6302,
 21.2501, 22.9075, 24.6174, 26.3979, 28.2730, 30.2770, 32.4631, 34.9243, 37.8590, 41.9115};
  static const G4double X4[nX]={.160483,
 5.73111, 8.88884, 11.5893, 14.0636, 16.4054, 18.6651, 20.8749, 23.0578, 25.2318, 27.4127,
 29.6152, 31.8540, 34.1452, 36.5074, 38.9635, 41.5435, 44.2892, 47.2638, 50.5732, 54.4265};
  static const G4double X5[nX]={.0999307,
 5.25720, 8.11389, 10.5375, 12.7425, 14.8152, 16.8015, 18.7296, 20.6194, 22.4855, 24.3398,
 26.1924, 28.0527, 29.9295, 31.8320, 33.7699, 35.7541, 37.7975, 39.9158, 42.1290, 44.4649};
  static const G4double X6[nX]={.0276367,
 3.53378, 5.41553, 6.99413, 8.41629, 9.74057, 10.9978, 12.2066, 13.3796, 14.5257, 15.6519,
 16.7636, 17.8651, 18.9603, 20.0527, 21.1453, 22.2411, 23.3430, 24.4538, 25.5765, 26.7148};
  static const G4double X7[nX]={.00472383,
 2.08253, 3.16946, 4.07178, 4.87742, 5.62140, 6.32202, 6.99034, 7.63368, 8.25720, 8.86473,
 9.45921, 10.0430, 10.6179, 11.1856, 11.7475, 12.3046, 12.8581, 13.4089, 13.9577, 14.5057};
  static const G4double X8[nX]={.000630783,
 1.22723, 1.85845, 2.37862, 2.84022, 3.26412, 3.66122, 4.03811, 4.39910, 4.74725, 5.08480,
 5.41346, 5.73457, 6.04921, 6.35828, 6.66250, 6.96250, 7.25884, 7.55197, 7.84232, 8.13037};
  static const G4double X9[nX]={7.49179e-05,
 .772574, 1.16623, 1.48914, 1.77460, 2.03586, 2.27983, 2.51069, 2.73118, 2.94322, 3.14823,
 3.34728, 3.54123, 3.73075, 3.91638, 4.09860, 4.27779, 4.45428, 4.62835, 4.80025, 4.97028};
  static const G4double XA[nX]={8.43437e-06,
 .530035, .798454, 1.01797, 1.21156, 1.38836, 1.55313, 1.70876, 1.85712, 1.99956, 2.13704,
 2.27031, 2.39994, 2.52640, 2.65007, 2.77127, 2.89026, 3.00726, 3.12248, 3.23607, 3.34823};
  static const G4double XB[nX]={9.27028e-07,
 .395058, .594211, .756726, .899794, 1.03025, 1.15167, 1.26619, 1.37523, 1.47979, 1.58059,
 1.67819, 1.77302, 1.86543, 1.95571, 2.04408, 2.13074, 2.21587, 2.29960, 2.38206, 2.46341};
  static const G4double XC[nX]={1.00807e-07,
 .316195, .474948, .604251, .717911, .821417, .917635, 1.00829, 1.09452, 1.17712, 1.25668,
 1.33364, 1.40835, 1.48108, 1.55207, 1.62150, 1.68954, 1.75631, 1.82193, 1.88650, 1.95014};
  static const G4double XD[nX]={1.09102e-08,
 .268227, .402318, .511324, .606997, .694011, .774803, .850843, .923097, .992243, 1.05878,
 1.12309, 1.18546, 1.24613, 1.30530, 1.36313, 1.41974, 1.47526, 1.52978, 1.58338, 1.63617};
  static const G4double XE[nX]={1.17831e-09,
 .238351, .356890, .453036, .537277, .613780, .684719, .751405, .814699, .875208, .933374,
 .989535, 1.04396, 1.09685, 1.14838, 1.19870, 1.24792, 1.29615, 1.34347, 1.38996, 1.43571};
  static const G4double XF[nX]={1.27141e-10,
 .219778, .328346, .416158, .492931, .562525, .626955, .687434, .744761, .799494, .852046,
 .902729, .951786, .999414, 1.04577, 1.09099, 1.13518, 1.17844, 1.22084, 1.26246, 1.30338};
  static const G4double XG[nX]={1.3713e-11,
 .208748, .310948, .393310, .465121, .530069, .590078, .646306, .699515, .750239, .798870,
 .845707, .890982, .934882, .977559, 1.01914, 1.05973, 1.09941, 1.13827, 1.17637, 1.21379};
  static const G4double XH[nX]={1.47877e-12,
 .203089, .301345, .380162, .448646, .510409, .567335, .620557, .670820, .718647, .764421,
 .808434, .850914, .892042, .931967, .970812, 1.00868, 1.04566, 1.08182, 1.11724, 1.15197};
  static const G4double XI[nX]={1.59454e-13,
 .201466, .297453, .374007, .440245, .499779, .554489, .605506, .653573, .699213, .742806,
 .784643, .824952, .863912, .901672, .938353, .974060, 1.00888, 1.04288, 1.07614, 1.10872};
  static const G4double XJ[nX]={1.71931e-14,
 .202988, .297870, .373025, .437731, .495658, .548713, .598041, .644395, .688302, .730147,
 .770224, .808762, .845943, .881916, .916805, .950713, .983728, 1.01592, 1.04737, 1.07813};
  // Direct table
  static const G4double Xmin[nE]={X0[0],X1[0],X2[0],X3[0],X4[0],X5[0],X6[0],X7[0],X8[0],
                        X9[0],XA[0],XB[0],XC[0],XD[0],XE[0],XF[0],XG[0],XH[0],XI[0],XJ[0]};
  static const G4double dX[nE]={
    (X0[lX]-X0[0])/lX, (X1[lX]-X1[0])/lX, (X2[lX]-X2[0])/lX, (X3[lX]-X3[0])/lX,
    (X4[lX]-X4[0])/lX, (X5[lX]-X5[0])/lX, (X6[lX]-X6[0])/lX, (X7[lX]-X7[0])/lX,
    (X8[lX]-X8[0])/lX, (X9[lX]-X9[0])/lX, (XA[lX]-XA[0])/lX, (XB[lX]-XB[0])/lX,
    (XC[lX]-XC[0])/lX, (XD[lX]-XD[0])/lX, (XE[lX]-XE[0])/lX, (XF[lX]-XF[0])/lX,
    (XG[lX]-XG[0])/lX, (XH[lX]-XH[0])/lX, (XI[lX]-XI[0])/lX, (XJ[lX]-XJ[0])/lX};
  static const G4double* Xl[nE]=
                             {X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,XA,XB,XC,XD,XE,XF,XG,XH,XI,XJ};
  static const G4double I0[nX]={0,
 .411893, 1.25559, 2.34836, 3.60264, 4.96046, 6.37874, 7.82342, 9.26643, 10.6840, 12.0555,
 13.3628, 14.5898, 15.7219, 16.7458, 17.6495, 18.4217, 19.0523, 19.5314, 19.8501, 20.0000};
  static const G4double I1[nX]={0,
 .401573, 1.22364, 2.28998, 3.51592, 4.84533, 6.23651, 7.65645, 9.07796, 10.4780, 11.8365,
 13.1360, 14.3608, 15.4967, 16.5309, 17.4516, 18.2481, 18.9102, 19.4286, 19.7946, 20.0000};
  static const G4double I2[nX]={0,
 .387599, 1.17339, 2.19424, 3.37090, 4.65066, 5.99429, 7.37071, 8.75427, 10.1232, 11.4586,
 12.7440, 13.9644, 15.1065, 16.1582, 17.1083, 17.9465, 18.6634, 19.2501, 19.6982, 20.0000};
  static const G4double I3[nX]={0,
 .366444, 1.09391, 2.04109, 3.13769, 4.33668, 5.60291, 6.90843, 8.23014, 9.54840, 10.8461,
 12.1083, 13.3216, 14.4737, 15.5536, 16.5512, 17.4573, 18.2630, 18.9603, 19.5417, 20.0000};
  static const G4double I4[nX]={0,
 .321962, .959681, 1.79769, 2.77753, 3.85979, 5.01487, 6.21916, 7.45307, 8.69991, 9.94515,
 11.1759, 12.3808, 13.5493, 14.6720, 15.7402, 16.7458, 17.6813, 18.5398, 19.3148, 20.0000};
  static const G4double I5[nX]={0,
 .257215, .786302, 1.49611, 2.34049, 3.28823, 4.31581, 5.40439, 6.53832, 7.70422, 8.89040,
 10.0865, 11.2833, 12.4723, 13.6459, 14.7969, 15.9189, 17.0058, 18.0517, 19.0515, 20.0000};
  static const G4double I6[nX]={0,
 .201608, .638914, 1.24035, 1.97000, 2.80354, 3.72260, 4.71247, 5.76086, 6.85724, 7.99243,
 9.15826, 10.3474, 11.5532, 12.7695, 13.9907, 15.2117, 16.4275, 17.6337, 18.8258, 20.0000};
  static const G4double I7[nX]={0,
 .168110, .547208, 1.07889, 1.73403, 2.49292, 3.34065, 4.26525, 5.25674, 6.30654, 7.40717,
 8.55196, 9.73492, 10.9506, 12.1940, 13.4606, 14.7460, 16.0462, 17.3576, 18.6767, 20.0000};
  static const G4double I8[nX]={0,
 .150652, .497557, .990048, 1.60296, 2.31924, 3.12602, 4.01295, 4.97139, 5.99395, 7.07415,
 8.20621, 9.38495, 10.6057, 11.8641, 13.1561, 14.4781, 15.8267, 17.1985, 18.5906, 20.0000};
  static const G4double I9[nX]={0,
 .141449, .470633, .941304, 1.53053, 2.22280, 3.00639, 3.87189, 4.81146, 5.81837, 6.88672,
 8.01128, 9.18734, 10.4106, 11.6772, 12.9835, 14.3261, 15.7019, 17.1080, 18.5415, 20.0000};
  static const G4double IA[nX]={0,
 .136048, .454593, .912075, 1.48693, 2.16457, 2.93400, 3.78639, 4.71437, 5.71163, 6.77265,
 7.89252, 9.06683, 10.2916, 11.5631, 12.8780, 14.2331, .625500, 17.0525, 18.5115, 20.0000};
  static const G4double IB[nX]={0,
 .132316, .443455, .891741, 1.45656, 2.12399, 2.88352, 3.72674, 4.64660, 5.63711, 6.69298,
 7.80955, 8.98262, 10.2084, 11.4833, 12.8042, 14.1681, 15.5721, 17.0137, 18.4905, 20.0000};
  static const G4double IC[nX]={0,
 .129197, .434161, .874795, 1.43128, 2.09024, 2.84158, 3.67721, 4.59038, 5.57531, 6.62696,
 7.74084, 8.91291, 10.1395, 11.4173, 12.7432, 14.1143, 15.5280, 16.9817, 18.4731, 20.0000};
  static const G4double ID[nX]={0,
 .126079, .424911, .857980, 1.40626, 2.05689, 2.80020, 3.62840, 4.53504, 5.51456, 6.56212,
 7.67342, 8.84458, 10.0721, 11.3527, 12.6836, 14.0618, 15.4849, 16.9504, 18.4562, 20.0000};
  static const G4double IE[nX]={0,
 .122530, .414424, .838964, 1.37801, 2.01931, 2.75363, 3.57356, 4.47293, 5.44644, 6.48949,
 7.59795, 8.76815, 9.99673, 11.2806, 12.6170, 14.0032, 15.4369, 16.9156, 18.4374, 20.0000};
  static const G4double IF[nX]={0,
 .118199, .401651, .815838, 1.34370, 1.97370, 2.69716, 3.50710, 4.39771, 5.36401, 6.40164,
 7.50673, 8.67581, 9.90572, 11.1936, 12.5367, 13.9326, 15.3790, 16.8737, 18.4146, 20.0000};
  static const G4double IG[nX]={0,
 .112809, .385761, .787075, 1.30103, 1.91700, 2.62697, 3.42451, 4.30424, 5.26158, 6.29249,
 7.39341, 8.56112, 9.79269, 11.0855, 12.4369, 13.8449, 15.3071, 16.8216, 18.3865, 20.0000};
  static const G4double IH[nX]={0,
 .106206, .366267, .751753, 1.24859, 1.84728, 2.54062, 3.32285, .189160, 5.13543, 6.15804,
 7.25377, 8.41975, 9.65334, 10.9521, 12.3139, 13.7367, 15.2184, 16.7573, 18.3517, 20.0000};
  static const G4double II[nX]={0,
 .098419, .343194, .709850, 1.18628, 1.76430, 2.43772, 3.20159, 4.05176, 4.98467, 5.99722,
 7.08663, 8.25043, 9.48633, 10.7923, 12.1663, 13.6067, 15.1118, 16.6800, 18.3099, 20.0000};
  static const G4double IJ[nX]={0,
 .089681, .317135, .662319, 1.11536, 1.66960, 2.32002, 3.06260, 3.89397, 4.81126, 5.81196,
 6.89382, 8.05483, 9.29317, 10.6072, 11.9952, 13.4560, 14.9881, 16.5902, 18.2612, 20.0000};
  static const G4double* Il[nE]=
                             {I0,I1,I2,I3,I4,I5,I6,I7,I8,I9,IA,IB,IC,ID,IE,IF,IG,IH,II,IJ};
  static const G4double lE[nE]={
-1.98842,-1.58049,-1.17256,-.764638,-.356711, .051215, .459141, .867068, 1.27499, 1.68292,
 2.09085, 2.49877, 2.90670, 3.31463, 3.72255, 4.13048, 4.53840, 4.94633, 5.35426, 5.76218};
  static const G4double lEmi=lE[0];
  static const G4double lEma=lE[nE-1];
  static const G4double dlE=(lEma-lEmi)/bE;
  //***************************************************************************************
  G4double Enu=lastE;                 // Get energy of the last calculated cross-section
  G4double lEn=std::log(Enu);         // log(E) for interpolation
  G4double rE=(lEn-lEmi)/dlE;         // Position of the energy
  G4int fE=static_cast<int>(rE);      // Left bin for interpolation
  if(fE<0) fE=0;
  if(fE>pE)fE=pE;
  G4int    sE=fE+1;                   // Right bin for interpolation
  G4double dE=rE-fE;                  // relative log shift from the left bin
  G4double dEnu=Enu+Enu;              // doubled energy of nu/anu
  G4double Enu2=Enu*Enu;              // squared energy of nu/anu
  G4double Emu=Enu;                   // Free Energy of neutrino/anti-neutrino
  G4double ME=Enu*MN;                 // M*E
  G4double dME=ME+ME;                 // 2*M*E
  G4double dEMN=(dEnu+MN)*ME;
  G4double sqE=Enu*ME;
  G4double E2M=MN*Enu2;
  G4double ymax=(E2M+sqE)/dEMN;
  G4double Q2mi=0.;                   // Q2_min(E_nu)
  G4double Q2ma=dME*ymax;             // Q2_max(E_nu)
  G4double Q2nq=Emu*dMN-mcV;
  if(Q2ma>Q2nq) Q2ma=Q2nq;            // Correction for Non Quasi Elastic
  // --- now r_min=Q2mi/Q2ma and r_max=1.; when r is randomized -> Q2=r*Q2ma ---
  G4double Rmi=Q2mi/Q2ma;
  G4double shift=1.+.9673/(1.+.323/Enu/Enu)/std::pow(Enu,.78); //@@ different for anti-nu
  // --- E-interpolation must be done in a log scale ---
  G4double Xmi=std::pow((shift-Rmi),power);// X_min(E_nu)
  G4double Xma=std::pow((shift-1.),power); // X_max(E_nu)
  // Find the integral values integ(Xmi) & integ(Xma) using the direct table
  G4double idX=dX[fE]+dE*(dX[sE]-dX[fE]); // interpolated X step
  G4double iXmi=Xmin[fE]+dE*(Xmin[sE]-Xmin[fE]); // interpolated X minimum
  G4double rXi=(Xmi-iXmi)/idX;
  G4int    iXi=static_cast<int>(rXi);
  if(iXi<0) iXi=0;
  if(iXi>bX) iXi=bX;
  G4double dXi=rXi-iXi;
  G4double bntil=Il[fE][iXi];
  G4double intil=bntil+dXi*(Il[fE][iXi+1]-bntil);
  G4double bntir=Il[sE][iXi];
  G4double intir=bntir+dXi*(Il[sE][iXi+1]-bntir);
  G4double inti=intil+dE*(intir-intil);// interpolated begin of the integral
  //
  G4double rXa=(Xma-iXmi)/idX;
  G4int    iXa=static_cast<int>(rXa);
  if(iXa<0) iXa=0;
  if(iXa>bX) iXa=bX;
  G4double dXa=rXa-iXa;
  G4double bntal=Il[fE][iXa];
  G4double intal=bntal+dXa*(Il[fE][iXa+1]-bntal);
  G4double bntar=Il[sE][iXa];
  G4double intar=bntar+dXa*(Il[sE][iXa+1]-bntar);
  G4double inta=intal+dE*(intar-intal);// interpolated end of the integral
  //
  // *** Find X using the reversed table ***
  G4double intx=inti+(inta-inti)*G4UniformRand(); 
  G4int    intc=static_cast<int>(intx);
  if(intc<0) intc=0;
  if(intc>bX) intc=bX;
  G4double dint=intx-intc;
  G4double mXl=Xl[fE][intc];
  G4double Xlb=mXl+dint*(Xl[fE][intc+1]-mXl);
  G4double mXr=Xl[sE][intc];
  G4double Xrb=mXr+dint*(Xl[sE][intc+1]-mXr);
  G4double X=Xlb+dE*(Xrb-Xlb);        // interpolated X value
  G4double R=shift-std::pow(X,pconv);
  G4double Q2=R*Q2ma;
  return Q2*GeV*GeV;
}

// It returns a fraction of the direct interaction of the neutrino with quark-partons
G4double G4QNuNuNuclearCrossSection::GetDirectPart(G4double Q2)
{
  G4double f=Q2/4.62;
  G4double ff=f*f;
  G4double r=ff*ff;
  G4double s_value=std::pow((1.+.6/Q2),(-1.-(1.+r)/(12.5+r/.3)));
  //@@ It is the same for nu/anu, but for nu it is a bit less, and for anu a bit more (par)
  return 1.-s_value*(1.-s_value/2);
}

// #of quark-partons in the nonperturbative phase space is the same for neut and anti-neut
G4double G4QNuNuNuclearCrossSection::GetNPartons(G4double Q2)
{
  return 3.+.3581*std::log(1.+Q2/.04); // a#of partons in the nonperturbative phase space
}

// This class can provide only virtual exchange gamma (a substitute for Z0 boson)
G4int G4QNuNuNuclearCrossSection::GetExchangePDGCode() {return 22;}
