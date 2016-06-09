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
// G4 Physics class: G4QANuMuNuclearCrossSection for (anu_mu,mu+)A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-Oct-03
// 
// ****************************************************************************************
// ********** This CLASS is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************
// ----------------------------------------------------------------------------------------
// Short desctription: antinu_mu -> mu cross-section
// ----------------------------------------------------------------------------------------

//#define debug
//#define edebug
//#define pdebug
//#define ppdebug
//#define tdebug
//#define sdebug

#include "G4QANuMuNuclearCrossSection.hh"
#include "G4SystemOfUnits.hh"

// Initialization of the
G4bool    G4QANuMuNuclearCrossSection::onlyCS=true;//Flag to calculate only CS (not QE)
G4double  G4QANuMuNuclearCrossSection::lastSig=0.;//Last calculated total cross section
G4double  G4QANuMuNuclearCrossSection::lastQEL=0.;//Last calculated quasi-el. cross section
G4int     G4QANuMuNuclearCrossSection::lastL=0;   //Last used in cross section TheLastBin
G4double  G4QANuMuNuclearCrossSection::lastE=0.;  //Last used in cross section TheEnergy
G4double* G4QANuMuNuclearCrossSection::lastEN=0;  //Pointer to the Energy Scale of TX & QE
G4double* G4QANuMuNuclearCrossSection::lastTX=0;  //Pointer to the LastArray of TX function
G4double* G4QANuMuNuclearCrossSection::lastQE=0;  //Pointer to the LastArray of QE function
G4int     G4QANuMuNuclearCrossSection::lastPDG=0; // The last PDG code of the projectile
G4int     G4QANuMuNuclearCrossSection::lastN=0;   // The last N of calculated nucleus
G4int     G4QANuMuNuclearCrossSection::lastZ=0;   // The last Z of calculated nucleus
G4double  G4QANuMuNuclearCrossSection::lastP=0.;  // Last used in cross section Momentum
G4double  G4QANuMuNuclearCrossSection::lastTH=0.; // Last threshold momentum
G4double  G4QANuMuNuclearCrossSection::lastCS=0.; // Last value of the Cross Section
G4int     G4QANuMuNuclearCrossSection::lastI=0;   // The last position in the DAMDB
std::vector<G4double*>* G4QANuMuNuclearCrossSection::TX = new std::vector<G4double*>;
std::vector<G4double*>* G4QANuMuNuclearCrossSection::QE = new std::vector<G4double*>;

// Returns Pointer to the G4VQCrossSection class
G4VQCrossSection* G4QANuMuNuclearCrossSection::GetPointer()
{
  static G4QANuMuNuclearCrossSection theCrossSection;//**Static body of the Cross Section**
  return &theCrossSection;
}

G4QANuMuNuclearCrossSection::~G4QANuMuNuclearCrossSection()
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
G4double G4QANuMuNuclearCrossSection::GetCrossSection(G4bool fCS, G4double pMom,
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
  G4cout<<"G4QAMNCS::GetCS:>> f="<<fCS<<", p="<<pMom<<", Z="<<tgZ<<"("<<lastZ<<") ,N="<<tgN
        <<"("<<lastN<<"),PDG="<<pPDG<<"("<<lastPDG<<"), T="<<pEn<<"("<<lastTH<<")"<<",Sz="
        <<colN.size()<<G4endl;
  //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  if(pPDG!=-14)
  {
#ifdef debug
    G4cout<<"G4QAMNCS::GetCS: *** Found pPDG="<<pPDG<<" =--=> CS=0"<<G4endl;
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
        G4cout<<"G4QAMNCS::GetCS:*Found*P="<<pMom<<",Threshold="<<lastTH<<",j="<<j<<G4endl;
        //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(pEn<=lastTH)
        {
#ifdef pdebug
          G4cout<<"G4QAMNCS::GetCS:Found T="<<pEn<<" < Threshold="<<lastTH<<",X=0"<<G4endl;
          //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
          return 0.;                     // Energy is below the Threshold value
        }
        lastP  =colP [i];                // Last Momentum  (A-dependent)
        lastCS =colCS[i];                // Last CrossSect (A-dependent)
        if(std::fabs(lastP/pMom-1.)<tolerance)
        {
#ifdef pdebug
          G4cout<<"G4QAMNCS::GetCS:P="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
          //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
          return lastCS*millibarn;     // Use theLastCS
        }
        in = true;                       // This is the case when the isotop is found in DB
        // Momentum pMom is in IU ! @@ Units
#ifdef pdebug
        G4cout<<"G4QAMNCS::G:UpdaDB P="<<pMom<<",f="<<fCS<<",lI="<<lastI<<",j="<<j<<G4endl;
#endif
        lastCS=CalculateCrossSection(fCS,-1,j,lastPDG,lastZ,lastN,pMom); // read & update
#ifdef pdebug
        G4cout<<"G4QAMNCS::GetCrosSec: *****> New (inDB) Calculated CS="<<lastCS<<G4endl;
        //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(lastCS<=0. && pEn>lastTH)    // Correct the threshold
        {
#ifdef pdebug
          G4cout<<"G4QAMNCS::GetCS: New T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
        break;                           // Go out of the LOOP
      }
#ifdef pdebug
      G4cout<<"---G4QAMNCrossSec::GetCrosSec:pPDG="<<pPDG<<",j="<<j<<",N="<<colN[i]
            <<",Z["<<i<<"]="<<colZ[i]<<",cPDG="<<colPDG[i]<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      j++;                             // Increment a#0f records found in DB for this pPDG
    }
    if(!in)                            // This nucleus has not been calculated previously
    {
#ifdef pdebug
      G4cout<<"G4QAMNCS::GetCrosSec:CalcNew P="<<pMom<<",f="<<fCS<<",lstI="<<lastI<<G4endl;
#endif
      //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
      lastCS=CalculateCrossSection(fCS,0,j,lastPDG,lastZ,lastN,pMom); //calculate & create
      if(lastCS<=0.)
      {
        lastTH = ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last
#ifdef pdebug
        G4cout<<"G4QAMNCrossSection::GetCrossSect: NewThresh="<<lastTH<<",T="<<pEn<<G4endl;
#endif
        if(pEn>lastTH)
        {
#ifdef pdebug
          G4cout<<"G4QAMNCS::GetCS: First T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
      }
#ifdef pdebug
      G4cout<<"G4QAMNCS::GetCrosSec:New CS="<<lastCS<<",lZ="<<lastN<<",lN="<<lastZ<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      colN.push_back(tgN);
      colZ.push_back(tgZ);
      colPDG.push_back(pPDG);
      colP.push_back(pMom);
      colTH.push_back(lastTH);
      colCS.push_back(lastCS);
#ifdef pdebug
      G4cout<<"G4QAMNCS::GetCS:1st,P="<<pMom<<"(MeV),X="<<lastCS*millibarn<<"(mb)"<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      return lastCS*millibarn;
    } // End of creation of the new set of parameters
    else
    {
#ifdef pdebug
      G4cout<<"G4QAMNCS::GetCS: Update lastI="<<lastI<<",j="<<j<<G4endl;
#endif
      colP[lastI]=pMom;
      colPDG[lastI]=pPDG;
      colCS[lastI]=lastCS;
    }
  } // End of parameters udate
  else if(pEn<=lastTH)
  {
#ifdef pdebug
    G4cout<<"G4QAMNCS::GetCS: Current T="<<pEn<<" < Threshold="<<lastTH<<", CS=0"<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return 0.;                         // Momentum is below the Threshold Value -> CS=0
  }
  else if(std::fabs(lastP/pMom-1.)<tolerance)
  {
#ifdef pdebug
    G4cout<<"G4QAMNCS::GetCS:OldCur P="<<pMom<<"="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return lastCS*millibarn;     // Use theLastCS
  }
  else
  {
#ifdef pdebug
    G4cout<<"G4QAMNCS::GetCS:UpdaCur P="<<pMom<<",f="<<fCS<<",I="<<lastI<<",j="<<j<<G4endl;
#endif
    lastCS=CalculateCrossSection(fCS,1,j,lastPDG,lastZ,lastN,pMom); // Only UpdateDB
    lastP=pMom;
  }
#ifdef pdebug
  G4cout<<"G4QAMNCS::GetCrSec:End,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
  //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  return lastCS*millibarn;
}

// Gives the threshold energy = the same for all nuclei (@@ can be reduced for hevy nuclei)
G4double G4QANuMuNuclearCrossSection::ThresholdEnergy(G4int Z, G4int N, G4int)
{
  //static const G4double mNeut = G4NucleiProperties::GetNuclearMass(1,0)/GeV;
  //static const G4double mProt = G4NucleiProperties::GetNuclearMass(1,1)/GeV;
  //static const G4double mDeut = G4NucleiProperties::GetNuclearMass(2,1)/GeV/2.;
  static const G4double mN=.931494043;// Nucleon mass (inside nucleus, AtomicMassUnit, GeV)
  static const G4double dmN=mN+mN;    // Doubled nucleon mass (2*AtomicMassUnit, GeV)
  static const G4double mmu=.105658369; // Mass of a muon in GeV
  static const G4double mmu2=mmu*mmu;   // Squared mass of a muon in GeV^2
  static const G4double thresh=mmu+mmu2/dmN; // Universal threshold in GeV
  // ---------
  //static const G4double infEn = 9.e27;
  G4double dN=0.;
  if(Z>0||N>0) dN=thresh*GeV; // @@ if upgraded, change it in a total cross section
  //@@ "dN=mmu+mmu2/G4NucleiProperties::GetNuclearMass(<G4double>(Z+N),<G4double>(Z)/GeV"
  return dN;
}

// The main member function giving the gamma-A cross section (E_kin in MeV, CS in mb)
G4double G4QANuMuNuclearCrossSection::CalculateCrossSection(G4bool CS, G4int F, G4int I,
                                       G4int , G4int targZ, G4int targN, G4double Momentum)
{
  static const G4double mb38=1.E-11; // Conversion 10^-38 cm^2 to mb=10^-27 cm^2
  static const G4int nE=65;   // !! If change this, change it in GetFunctions() (*.hh) !!
  static const G4int mL=nE-1;
  static const G4double mN=.931494043;// Nucleon mass (inside nucleus, AtomicMassUnit, GeV)
  static const G4double dmN=mN+mN;   // Doubled nucleon mass (2*AtomicMassUnit, GeV)
  static const G4double mmu=.105658369; // Mass of a muon in GeV
  static const G4double mmu2=mmu*mmu;// Squared mass of a muon in GeV^2
  static const G4double EMi=mmu+mmu2/dmN; // Universal threshold of the reaction in GeV
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
  // =---------------= NOW Calculate the Cross Section =-------------------=
  if (lastE<=EMi)                    // Check that antiNuEnergy is higher than ThreshE
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
G4int G4QANuMuNuclearCrossSection::GetFunctions (G4int z, G4int n,
                                                     G4double* t, G4double* q, G4double* e)
{
  static const G4double mN=.931494043;// Nucleon mass (inside nucleus, AtomicMassUnit, GeV)
  static const G4double dmN=mN+mN;    // Doubled nucleon mass (2*AtomicMassUnit, GeV)
  static const G4double mmu=.105658369; // Mass of a muon in GeV
  static const G4double mmu2=mmu*mmu;   // Squared mass of a muon in GeV^2
  static const G4double thresh=mmu+mmu2/dmN; // Universal threshold in GeV
  static const G4int nE=65; // !! If change this, change it in GetCrossSection() (*.cc) !!
  static const G4double nuEn[nE]={thresh,
    .112039,.116079,.120416,.125076,.130090,.135494,.141324,.147626,.154445,.161838,
    .169864,.178594,.188105,.198485,.209836,.222272,.235923,.250941,.267497,.285789,
    .306045,.328530,.353552,.381466,.412689,.447710,.487101,.531538,.581820,.638893,
    .703886,.778147,.863293,.961275,1.07445,1.20567,1.35843,1.53701,1.74667,1.99390,
    2.28679,2.63542,3.05245,3.55386,4.15990,4.89644,5.79665,6.90336,8.27224,9.97606,
    12.1106,14.8029,18.2223,22.5968,28.2351,35.5587,45.1481,57.8086,74.6682,97.3201,
    128.036,170.085,228.220,309.420};
  static const G4double TOTX[nE]={0.,
    .077498,.247583,.329691,.386384,.429087,.462699,.489899,.512316,.530996,.546614,
    .559616,.570292,.578840,.585395,.590053,.593083,.594197,.593614,.591396,.587611,
    .582335,.575653,.567667,.558490,.548417,.537270,.525352,.512825,.499857,.486620,
    .473283,.460014,.446970,.434294,.422116,.410656,.399782,.389665,.380349,.371860,
    .364207,.357387,.351388,.346192,.341778,.338122,.335198,.332980,.331439,.330544,
    .330263,.330558,.331391,.332718,.334494,.336667,.339182,.341697,.344470,.348125,
    .351322,.354481,.357507,.359239};
  static const G4double QELX[nE]={0.,
    .008683,.028739,.039700,.048327,.055820,.062693,.069235,.075631,.082010,.088463,
    .095059,.101851,.108883,.116192,.123814,.131826,.140185,.148962,.158197,.167933,
    .178221,.189119,.200700,.213045,.226326,.240454,.255277,.270612,.286388,.302608,
    .319318,.336582,.354468,.373031,.392427,.412445,.433146,.454448,.476222,.498289,
    .520430,.542558,.564130,.585003,.604928,.623680,.641266,.657255,.671704,.684586,
    .696111,.706028,.714553,.721951,.728085,.733182,.737348,.740958,.743716,.746059,
    .747806,.749129,.750331,.751100};

  // --------------------------------
  G4int first=0;
  if(z<0.)
  {
    first=1;
    z=-z;
  }
  if(z<1 || z>92)             // neutron & plutonium are forbidden
  {
    G4cout<<"**G4QANuMuNuclearCrossSection::GetFunctions:Z="<<z<<".No CS returned"<<G4endl;
    return -1;
  }
  for(G4int k=0; k<nE; k++)
  {
    G4double a=n+z;
    G4double za=z+a;
    G4double dz=z+z;
    G4double da=a+a;
    G4double ta=da+a;
    if(first) e[k]=nuEn[k];       // Energy of neutrino E (first bin k=0 can be modified)
    t[k]=TOTX[k]*nuEn[k]*(za+za)/ta+QELX[k]*(dz+dz-da)/ta; // TotalCrossSection
    q[k]=QELX[k]*dz/a;                                     // QuasiElasticCrossSection
  }
  return first;
}

// Randomize Q2 from neutrino to the scattered muon when the scattering is quasi-elastic
G4double G4QANuMuNuclearCrossSection::GetQEL_ExchangeQ2()
{
  static const G4double mmu=.105658369;// Mass of muon in GeV
  static const G4double mmu2=mmu*mmu;  // Squared Mass of muon in GeV^2
  static const double hmmu2=mmu2/2;    // .5*m_mu^2 in GeV^2
  static const double MN=.931494043;   // Nucleon mass (inside nucleus, atomicMassUnit,GeV)
  static const double MN2=MN*MN;       // M_N^2 in GeV^2
  static const G4double power=-3.5;    // direct power for the magic variable
  static const G4double pconv=1./power;// conversion power for the magic variable
  static const G4int nQ2=101;          // #Of point in the Q2l table (in GeV^2)
  static const G4int lQ2=nQ2-1;        // index of the last in the Q2l table
  static const G4int bQ2=lQ2-1;        // index of the before last in the Q2 ltable
  // Reversed table
  static const G4double Xl[nQ2]={5.20224e-16,
 .006125,.0137008,.0218166,.0302652,.0389497,.0478144,.0568228,.0659497,.0751768,.0844898,
 .093878, .103332, .112844, .122410, .132023, .141680, .151376, .161109, .170875, .180672,
 .190499, .200352, .210230, .220131, .230055, .239999, .249963, .259945, .269944, .279960,
 .289992, .300039, .310099, .320173, .330260, .340359, .350470, .360592, .370724, .380867,
 .391019, .401181, .411352, .421531, .431719, .441915, .452118, .462329, .472547, .482771,
 .493003, .503240, .513484, .523734, .533989, .544250, .554517, .564788, .575065, .585346,
 .595632, .605923, .616218, .626517, .636820, .647127, .657438, .667753, .678072, .688394,
 .698719, .709048, .719380, .729715, .740053, .750394, .760738, .771085, .781434, .791786,
 .802140, .812497, .822857, .833219, .843582, .853949, .864317, .874687, .885060, .895434,
 .905810, .916188, .926568, .936950, .947333, .957719, .968105, .978493, .988883, .999275};
   // Direct table
  static const G4double Xmax=Xl[lQ2];
  static const G4double Xmin=Xl[0];
  static const G4double dX=(Xmax-Xmin)/lQ2;  // step in X(Q2, GeV^2)
  static const G4double inl[nQ2]={0,
 1.52225, 2.77846, 3.96651, 5.11612, 6.23990, 7.34467, 8.43466, 9.51272, 10.5809, 11.6406,
 12.6932, 13.7394, 14.7801, 15.8158, 16.8471, 17.8743, 18.8979, 19.9181, 20.9353, 21.9496,
 22.9614, 23.9707, 24.9777, 25.9826, 26.9855, 27.9866, 28.9860, 29.9837, 30.9798, 31.9745,
 32.9678, 33.9598, 34.9505, 35.9400, 36.9284, 37.9158, 38.9021, 39.8874, 40.8718, 41.8553,
 42.8379, 43.8197, 44.8007, 45.7810, 46.7605, 47.7393, 48.7174, 49.6950, 50.6718, 51.6481,
 52.6238, 53.5990, 54.5736, 55.5476, 56.5212, 57.4943, 58.4670, 59.4391, 60.4109, 61.3822,
 62.3531, 63.3236, 64.2937, 65.2635, 66.2329, 67.2019, 68.1707, 69.1390, 70.1071, 71.0748,
 72.0423, 73.0095, 73.9763, 74.9429, 75.9093, 76.8754, 77.8412, 78.8068, 79.7721, 80.7373,
 81.7022, 82.6668, 83.6313, 84.5956, 85.5596, 86.5235, 87.4872, 88.4507, 89.4140, 90.3771,
 91.3401, 92.3029, 93.2656, 94.2281, 95.1904, 96.1526, 97.1147, 98.0766, 99.0384, 100.000};
  G4double Enu=lastE;                 // Get energy of the last calculated cross-section
  G4double dEnu=Enu+Enu;              // doubled energy of nu/anu
  G4double Enu2=Enu*Enu;              // squared energy of nu/anu
  G4double ME=Enu*MN;                 // M*E
  G4double dME=ME+ME;                 // 2*M*E
  G4double dEMN=(dEnu+MN)*ME;
  G4double MEm=ME-hmmu2;
  G4double sqE=Enu*std::sqrt(MEm*MEm-mmu2*MN2);
  G4double E2M=MN*Enu2-(Enu+MN)*hmmu2;
  G4double ymax=(E2M+sqE)/dEMN;
  G4double ymin=(E2M-sqE)/dEMN;
  G4double rmin=1.-ymin;
  G4double rhm2E=hmmu2/Enu2;
  G4double Q2mi=(Enu2+Enu2)*(rmin-rhm2E-std::sqrt(rmin*rmin-rhm2E-rhm2E)); // Q2_min(E_nu)
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
G4double G4QANuMuNuclearCrossSection::GetNQE_ExchangeQ2()
{
  static const double mpi=.13957018;    // charged pi meson mass in GeV
  static const G4double mmu=.105658369; // Mass of muon in GeV
  static const G4double mmu2=mmu*mmu;   // Squared Mass of muon in GeV^2
  static const double hmmu2=mmu2/2;     // .5*m_mu^2 in GeV^2
  static const double MN=.931494043;    // Nucleon mass (inside nucleus,atomicMassUnit,GeV)
  static const double MN2=MN*MN;        // M_N^2 in GeV^2
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
  static const G4double X0[nX]={5.21412e-05,
 .437860, .681908, .891529, 1.08434, 1.26751, 1.44494, 1.61915, 1.79198, 1.96493, 2.13937,
 2.31664, 2.49816, 2.68559, 2.88097, 3.08705, 3.30774, 3.54917, 3.82233, 4.15131, 4.62182};
  static const G4double X1[nX]={.00102591,
 1.00443, 1.55828, 2.03126, 2.46406, 2.87311, 3.26723, 3.65199, 4.03134, 4.40835, 4.78561,
 5.16549, 5.55031, 5.94252, 6.34484, 6.76049, 7.19349, 7.64917, 8.13502, 8.66246, 9.25086};
  static const G4double X2[nX]={.0120304,
 2.59903, 3.98637, 5.15131, 6.20159, 7.18024, 8.10986, 9.00426, 9.87265, 10.7217, 11.5564,
 12.3808, 13.1983, 14.0116, 14.8234, 15.6359, 16.4515, 17.2723, 18.1006, 18.9386, 19.7892};
  static const G4double X3[nX]={.060124,
 5.73857, 8.62595, 10.9849, 13.0644, 14.9636, 16.7340, 18.4066, 20.0019, 21.5342, 23.0142,
 24.4497, 25.8471, 27.2114, 28.5467, 29.8564, 31.1434, 32.4102, 33.6589, 34.8912, 36.1095};
  static const G4double X4[nX]={.0992363,
 8.23746, 12.1036, 15.1740, 17.8231, 20.1992, 22.3792, 24.4092, 26.3198, 28.1320, 29.8615,
 31.5200, 33.1169, 34.6594, 36.1536, 37.6044, 39.0160, 40.3920, 41.7353, 43.0485, 44.3354};
  static const G4double X5[nX]={.0561127,
 7.33661, 10.5694, 13.0778, 15.2061, 17.0893, 18.7973, 20.3717, 21.8400, 23.2211, 24.5291,
 25.7745, 26.9655, 28.1087, 29.2094, 30.2721, 31.3003, 32.2972, 33.2656, 34.2076, 35.1265};
  static const G4double X6[nX]={.0145859,
 4.81774, 6.83565, 8.37399, 9.66291, 10.7920, 11.8075, 12.7366, 13.5975, 14.4025, 15.1608,
 15.8791, 16.5628, 17.2162, 17.8427, 18.4451, 19.0259, 19.5869, 20.1300, 20.6566, 21.1706};
  static const G4double X7[nX]={.00241155,
 2.87095, 4.02492, 4.89243, 5.61207, 6.23747, 6.79613, 7.30433, 7.77270, 8.20858, 8.61732,
 9.00296, 9.36863, 9.71682, 10.0495, 10.3684, 10.6749, 10.9701, 11.2550, 11.5306, 11.7982};
  static const G4double X8[nX]={.000316863,
 1.76189, 2.44632, 2.95477, 3.37292, 3.73378, 4.05420, 4.34415, 4.61009, 4.85651, 5.08666,
 5.30299, 5.50738, 5.70134, 5.88609, 6.06262, 6.23178, 6.39425, 6.55065, 6.70149, 6.84742};
  static const G4double X9[nX]={3.73544e-05,
 1.17106, 1.61289, 1.93763, 2.20259, 2.42976, 2.63034, 2.81094, 2.97582, 3.12796, 3.26949,
 3.40202, 3.52680, 3.64482, 3.75687, 3.86360, 3.96557, 4.06323, 4.15697, 4.24713, 4.33413};
  static const G4double XA[nX]={4.19131e-06,
 .849573, 1.16208, 1.38955, 1.57379, 1.73079, 1.86867, 1.99221, 2.10451, 2.20770, 2.30332,
 2.39252, 2.47622, 2.55511, 2.62977, 2.70066, 2.76818, 2.83265, 2.89437, 2.95355, 3.01051};
  static const G4double XB[nX]={4.59981e-07,
 .666131, .905836, 1.07880, 1.21796, 1.33587, 1.43890, 1.53080, 1.61399, 1.69011, 1.76040,
 1.82573, 1.88682, 1.94421, 1.99834, 2.04959, 2.09824, 2.14457, 2.18878, 2.23107, 2.27162};
  static const G4double XC[nX]={4.99861e-08,
 .556280, .752730, .893387, 1.00587, 1.10070, 1.18317, 1.25643, 1.32247, 1.38269, 1.43809,
 1.48941, 1.53724, 1.58203, 1.62416, 1.66391, 1.70155, 1.73728, 1.77128, 1.80371, 1.83473};
  static const G4double XD[nX]={5.40832e-09,
 .488069, .657650, .778236, .874148, .954621, 1.02432, 1.08599, 1.14138, 1.19172, 1.23787,
 1.28049, 1.32008, 1.35705, 1.39172, 1.42434, 1.45514, 1.48429, 1.51197, 1.53829, 1.56339};
  static const G4double XE[nX]={5.84029e-10,
 .445057, .597434, .705099, .790298, .861468, .922865, .976982, 1.02542, 1.06930, 1.10939,
 1.14630, 1.18050, 1.21233, 1.24208, 1.27001, 1.29630, 1.32113, 1.34462, 1.36691, 1.38812};
  static const G4double XF[nX]={6.30137e-11,
 .418735, .560003, .659168, .737230, .802138, .857898, .906854, .950515, .989915, 1.02580,
 1.05873, 1.08913, 1.11734, 1.14364, 1.16824, 1.19133, 1.21306, 1.23358, 1.25298, 1.27139};
  static const G4double XG[nX]={6.79627e-12,
 .405286, .539651, .633227, .706417, .766929, .818642, .863824, .903931, .939963, .972639,
 1.00250, 1.02995, 1.05532, 1.07887, 1.10082, 1.12134, 1.14058, 1.15867, 1.17572, 1.19183};
  static const G4double XH[nX]={7.32882e-13,
 .404391, .535199, .625259, .695036, .752243, .800752, .842823, .879906, .912994, .942802,
 .969862, .994583, 1.01729, 1.03823, 1.05763, 1.07566, 1.09246, 1.10816, 1.12286, 1.13667};
  static const G4double XI[nX]={7.90251e-14,
 .418084, .548382, .636489, .703728, .758106, .803630, .842633, .876608, .906576, .933269,
 .957233, .978886, .998556, 1.01651, 1.03295, 1.04807, 1.06201, 1.07489, 1.08683, 1.09792};
  static const G4double XJ[nX]={8.52083e-15,
 .447299, .579635, .666780, .731788, .783268, .825512, .861013, .891356, .917626, .940597,
 .960842, .978802, .994820, 1.00917, 1.02208, 1.03373, 1.04427, 1.05383, 1.06253, 1.07046};
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
 .354631, 1.08972, 2.05138, 3.16564, 4.38343, 5.66828, 6.99127, 8.32858, 9.65998, 10.9680,
 12.2371, 13.4536, 14.6050, 15.6802, 16.6686, 17.5609, 18.3482, 19.0221, 19.5752, 20.0000};
  static const G4double I1[nX]={0,
 .281625, .877354, 1.67084, 2.60566, 3.64420, 4.75838, 5.92589, 7.12829, 8.34989, 9.57708,
 10.7978, 12.0014, 13.1781, 14.3190, 15.4162, 16.4620, 17.4496, 18.3724, 19.2245, 20.0000};
  static const G4double I2[nX]={0,
 .201909, .642991, 1.24946, 1.98463, 2.82370, 3.74802, 4.74263, 5.79509, 6.89474, 8.03228,
 9.19947, 10.3889, 11.5938, 12.8082, 14.0262, 15.2427, 16.4527, 17.6518, 18.8356, 20.0000};
  static const G4double I3[nX]={0,
 .140937, .461189, .920216, 1.49706, 2.17728, 2.94985, 3.80580, 4.73758, 5.73867, 6.80331,
 7.92637, 9.10316, 10.3294, 11.6013, 12.9150, 14.2672, 15.6548, 17.0746, 18.5239, 20.0000};
  static const G4double I4[nX]={0,
 .099161, .337358, .694560, 1.16037, 1.72761, 2.39078, 3.14540, 3.98768, 4.91433, 5.92245,
 7.00942, 8.17287, 9.41060, 10.7206, 12.1010, 13.5500, 15.0659, 16.6472, 18.2924, 20.0000};
  static const G4double I5[nX]={0,
 .071131, .255084, .543312, .932025, 1.41892, 2.00243, 2.68144, 3.45512, 4.32283, 5.28411,
 6.33859, 7.48602, 8.72621, 10.0590, 11.4844, 13.0023, 14.6128, 16.3158, 18.1115, 20.0000};
  static const G4double I6[nX]={0,
 .053692, .202354, .443946, .778765, 1.20774, 1.73208, 2.35319, 3.07256, 3.89177, 4.81249,
 5.83641, 6.96528, 8.20092, 9.54516, 10.9999, 12.5670, 14.2486, 16.0466, 17.9630, 20.0000};
  static const G4double I7[nX]={0,
 .043065, .168099, .376879, .672273, 1.05738, 1.53543, 2.10973, 2.78364, 3.56065, 4.44429,
 5.43819, 6.54610, 7.77186, 9.11940, 10.5928, 12.1963, 13.9342, 15.8110, 17.8313, 20.0000};
  static const G4double I8[nX]={0,
 .036051, .143997, .327877, .592202, .941572, 1.38068, 1.91433, 2.54746, 3.28517, 4.13277,
 5.09574, 6.17984, 7.39106, 8.73568, 10.2203, 11.8519, 13.6377, 15.5854, 17.7033, 20.0000};
  static const G4double I9[nX]={0,
 .030977, .125727, .289605, .528146, .846967, 1.25183, 1.74871, 2.34384, 3.04376, 3.85535,
 4.78594, 5.84329, 7.03567, 8.37194, 9.86163, 11.5150, 13.3430, 15.3576, 17.5719, 20.0000};
  static const G4double IA[nX]={0,
 .027129, .111420, .258935, .475812, .768320, 1.14297, 1.60661, 2.16648, 2.83034, 3.60650,
 4.50394, 5.53238, 6.70244, 8.02569, 9.51488, 11.1841, 13.0488, 15.1264, 17.4362, 20.0000};
  static const G4double IB[nX]={0,
 .024170, .100153, .234345, .433198, .703363, 1.05184, 1.48607, 2.01409, 2.64459, 3.38708,
 4.25198, 5.25084, 6.39647, 7.70319, 9.18708, 10.8663, 12.7617, 14.8968, 17.2990, 20.0000};
  static const G4double IC[nX]={0,
 .021877, .091263, .214670, .398677, .650133, .976322, 1.38510, 1.88504, 2.48555, 3.19709,
 4.03129, 5.00127, 6.12184, 7.40989, 8.88482, 10.5690, 12.4888, 14.6748, 17.1638, 20.0000};
  static const G4double ID[nX]={0,
 .020062, .084127, .198702, .370384, .606100, .913288, 1.30006, 1.77535, 2.34912, 3.03253,
 3.83822, 4.78063, 5.87634, 7.14459, 8.60791, 10.2929, 12.2315, 14.4621, 17.0320, 20.0000};
  static const G4double IE[nX]={0,
 .018547, .078104, .185102, .346090, .567998, .858331, 1.22535, 1.67824, 2.22735, 2.88443,
 3.66294, 4.57845, 5.64911, 6.89637, 8.34578, 10.0282, 11.9812, 14.2519, 16.8993, 20.0000};
  static const G4double IF[nX]={0,
 .017143, .072466, .172271, .323007, .531545, .805393, 1.15288, 1.58338, 2.10754, 2.73758,
 3.48769, 4.37450, 5.41770, 6.64092, 8.07288, 9.74894, 11.7135, 14.0232, 16.7522, 20.0000};
  static const G4double IG[nX]={0,
 .015618, .066285, .158094, .297316, .490692, .745653, 1.07053, 1.47479, 1.96931, 2.56677,
 3.28205, 4.13289, 5.14068, 6.33158, 7.73808, 9.40133, 11.3745, 13.7279, 16.5577, 20.0000};
  static const G4double IH[nX]={0,
 .013702, .058434, .139923, .264115, .437466, .667179, .961433, 1.32965, 1.78283, 2.33399,
 2.99871, 3.79596, 4.74916, 5.88771, 7.24937, 8.88367, 10.8576, 13.2646, 16.2417, 20.0000};
  static const G4double II[nX]={0,
 .011264, .048311, .116235, .220381, .366634, .561656, .813132, 1.13008, 1.52322, 2.00554,
 2.59296, 3.30542, 4.16834, 5.21490, 6.48964, 8.05434, 9.99835, 12.4580, 15.6567, 20.0000};
  static const G4double IJ[nX]={0,
 .008628, .037206, .089928, .171242, .286114, .440251, .640343, .894382, 1.21208, 1.60544,
 2.08962, 2.68414, 3.41486, 4.31700, 5.44048, 6.85936, 8.69067, 11.1358, 14.5885, 20.0000};
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
  G4double Emu=Enu-mmu;               // Free Energy of neutrino/anti-neutrino
  G4double ME=Enu*MN;                 // M*E
  G4double dME=ME+ME;                 // 2*M*E
  G4double dEMN=(dEnu+MN)*ME;
  G4double MEm=ME-hmmu2;
  G4double sqE=Enu*std::sqrt(MEm*MEm-mmu2*MN2);
  G4double E2M=MN*Enu2-(Enu+MN)*hmmu2;
  G4double ymax=(E2M+sqE)/dEMN;
  G4double ymin=(E2M-sqE)/dEMN;
  G4double rmin=1.-ymin;
  G4double rhm2E=hmmu2/Enu2;
  G4double Q2mi=(Enu2+Enu2)*(rmin-rhm2E-std::sqrt(rmin*rmin-rhm2E-rhm2E)); // Q2_min(E_nu)
  G4double Q2ma=dME*ymax;                                                  // Q2_max(E_nu)
  G4double Q2nq=Emu*dMN-mcV;
  if(Q2ma>Q2nq) Q2ma=Q2nq;            // Correction for Non Quasi Elastic
  // --- now r_min=Q2mi/Q2ma and r_max=1.; when r is randomized -> Q2=r*Q2ma ---
  G4double Rmi=Q2mi/Q2ma;
  G4double shift=.875/(1.+.2977/Enu/Enu)/std::pow(Enu,.78);
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
G4double G4QANuMuNuclearCrossSection::GetDirectPart(G4double Q2)
{
  G4double f=Q2/4.62;
  G4double ff=f*f;
  G4double r=ff*ff;
  G4double s_value=std::pow((1.+.6/Q2),(-1.-(1.+r)/(12.5+r/.3)));
  //@@ It is the same for nu/anu, but for nu it is a bit less, and for anu a bit more (par)
  return 1.-s_value*(1.-s_value/2);
}

// #of quark-partons in the nonperturbative phase space is the same for neut and anti-neut
G4double G4QANuMuNuclearCrossSection::GetNPartons(G4double Q2)
{
  return 3.+.3581*std::log(1.+Q2/.04); // a#of partons in the nonperturbative phase space
}

// This class can provide only virtual exchange pi- (a substitute for W- boson)
G4int G4QANuMuNuclearCrossSection::GetExchangePDGCode() {return -211;}
