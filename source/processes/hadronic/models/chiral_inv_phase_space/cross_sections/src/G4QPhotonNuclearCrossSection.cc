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
// The lust update: M.V. Kossov, CERN/ITEP(Moscow) 17-June-02
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Physics class: G4QPhotonNuclearCrossSection for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-03
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Feb-04
// ****************************************************************************************
// ********** This CLASS is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************
// Short description: This is an original CHIPS process for photo-nuclear
// interactions, which does not include "fast and dirty" corrections for
// reactions near threshold, with respect to the GHAD application of CHIPS.
// ------------------------------------------------------------------------

//#define debug
//#define pdebug
//#define debug3
//#define debugn
//#define debugs

#include "G4QPhotonNuclearCrossSection.hh"
#include "G4SystemOfUnits.hh"

// Initialization of the static variables
G4bool    G4QPhotonNuclearCrossSection::onlyCS=true;// Flag to calculate only CS
G4double  G4QPhotonNuclearCrossSection::lastSig=0.;// Last value of the Cross Section
G4double* G4QPhotonNuclearCrossSection::lastGDR=0; // Pointer to the lastArray of GDR CS
G4double* G4QPhotonNuclearCrossSection::lastHEN=0; // Pointer to the last array of HEn CS
G4double  G4QPhotonNuclearCrossSection::lastE=0.;  // LastUsed in CrossSections TheEnergy
G4double  G4QPhotonNuclearCrossSection::lastSP=0.; // Last value of ShadowingPomeron(A-dep)
G4int     G4QPhotonNuclearCrossSection::lastPDG=0; // The last PDG code of the projectile
G4int     G4QPhotonNuclearCrossSection::lastN=0;   // The last N of calculated nucleus
G4int     G4QPhotonNuclearCrossSection::lastZ=0;   // The last Z of calculated nucleus
G4double  G4QPhotonNuclearCrossSection::lastP=0.;  // Last used in cross section Momentum
G4double  G4QPhotonNuclearCrossSection::lastTH=0.; // Last threshold momentum
G4double  G4QPhotonNuclearCrossSection::lastCS=0.; // Last value of the Cross Section
G4int     G4QPhotonNuclearCrossSection::lastI=0;   // The last position in the DAMDB
std::vector<G4double*>* G4QPhotonNuclearCrossSection::GDR = new std::vector<G4double*>;
std::vector<G4double*>* G4QPhotonNuclearCrossSection::HEN = new std::vector<G4double*>;

// Returns Pointer to the G4VQCrossSection class
G4VQCrossSection* G4QPhotonNuclearCrossSection::GetPointer()
{
  static G4QPhotonNuclearCrossSection theCrossSection; //**Static body of Cross Section**
  return &theCrossSection;
}

G4QPhotonNuclearCrossSection::~G4QPhotonNuclearCrossSection()
{
  G4int lens=GDR->size();
  for(G4int i=0; i<lens; ++i) delete[] (*GDR)[i];
  delete GDR;
  G4int hens=HEN->size();
  for(G4int i=0; i<hens; ++i) delete[] (*HEN)[i];
  delete HEN;
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4QPhotonNuclearCrossSection::GetCrossSection(G4bool fCS, G4double pMom,
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
  G4cout<<"G4QPhCS::GetCS:>>> f="<<fCS<<", p="<<pMom<<", Z="<<tgZ<<"("<<lastZ<<") ,N="<<tgN
        <<"("<<lastN<<"),PDG="<<pPDG<<"("<<lastPDG<<"), T="<<pEn<<"("<<lastTH<<")"<<",Sz="
        <<colN.size()<<G4endl;
  //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  if(!pPDG)
  {
#ifdef debug
    G4cout<<"G4QPhCS::GetCS: *** Found pPDG="<<pPDG<<" =--=> CS=0"<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return 0.;                         // projectile PDG=0 is a mistake (?!) @@
  }
  if(pPDG==22 && tgZ==1 && !tgN && pMom<150.*MeV) return 0.; // Examle of pre-threshold (A)
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
        lastTH =colTH[i];              // Last THreshold (A-dependent)
#ifdef debug
        G4cout<<"G4QPhCS::GetCS:*Found* P="<<pMom<<",Threshold="<<lastTH<<",j="<<j<<G4endl;
        //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(pEn<=lastTH)
        {
#ifdef debug
          G4cout<<"G4QPhCS::GetCS:Found T="<<pEn<<" < Threshold="<<lastTH<<",CS=0"<<G4endl;
          //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
          return 0.;                   // Energy is below the Threshold value
        }
        lastP  =colP [i];              // Last Momentum  (A-dependent)
        lastCS =colCS[i];              // Last CrossSect (A-dependent)
        if(lastP == pMom)
        {
#ifdef debug
          G4cout<<"G4QPhCS::GetCS:P="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
#endif
          CalculateCrossSection(fCS,-1,j,lastPDG,lastZ,lastN,pMom); // Update param's only
          return lastCS*millibarn;     // Use theLastCS
        }
        in = true;                     // This is the case when the isotop is found in DB
        // Momentum pMom is in IU ! @@ Units
#ifdef pdebug
        G4cout<<"G4QPhCS::G:UpdatDB P="<<pMom<<",f="<<fCS<<",lI="<<lastI<<",j="<<j<<G4endl;
#endif
        lastCS=CalculateCrossSection(fCS,-1,j,lastPDG,lastZ,lastN,pMom); // read & update
#ifdef debug
        G4cout<<"G4QPhCS::GetCrosSec: *****> New (inDB) Calculated CS="<<lastCS<<G4endl;
        //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(lastCS<=0. && pEn>lastTH)   // Correct the threshold
        {
#ifdef debug
          G4cout<<"G4QPhCS::GetCS: New T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
        break;                         // Go out of the LOOP
      }
#ifdef debug
      G4cout<<"---G4QPhCrossSec::GetCrosSec:pPDG="<<pPDG<<",j="<<j<<",N="<<colN[i]
            <<",Z["<<i<<"]="<<colZ[i]<<",cPDG="<<colPDG[i]<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      j++;                             // Increment a#0f records found in DB for this pPDG
    }
    if(!in)                            // This nucleus has not been calculated previously
    {
#ifdef debug
      G4cout<<"G4QPhCS::GetCrosSec:CalcNew P="<<pMom<<",f="<<fCS<<",lastI="<<lastI<<G4endl;
#endif
      //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
      lastCS=CalculateCrossSection(fCS,0,j,lastPDG,lastZ,lastN,pMom); //calculate & create
      if(lastCS<=0.)
      {
        lastTH = ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last
#ifdef debug
        G4cout<<"G4QPhCrossSection::GetCrossSect: NewThresh="<<lastTH<<",T="<<pEn<<G4endl;
#endif
        if(pEn>lastTH)
        {
#ifdef debug
          G4cout<<"G4QPhCS::GetCS: First T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
      }
#ifdef debug
      G4cout<<"G4QPhCS::GetCrosSec: New CS="<<lastCS<<",lZ="<<lastN<<",lN="<<lastZ<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      colN.push_back(tgN);
      colZ.push_back(tgZ);
      colPDG.push_back(pPDG);
      colP.push_back(pMom);
      colTH.push_back(lastTH);
      colCS.push_back(lastCS);
#ifdef debug
      G4cout<<"G4QPhCS::GetCS:1st,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      return lastCS*millibarn;
    } // End of creation of the new set of parameters
    else
    {
#ifdef debug
      G4cout<<"G4QPrCS::GetCS: Update lastI="<<lastI<<",j="<<j<<G4endl;
#endif
      colP[lastI]=pMom;
      colPDG[lastI]=pPDG;
      colCS[lastI]=lastCS;
    }
  } // End of parameters udate
  else if(pEn<=lastTH)
  {
#ifdef debug
    G4cout<<"G4QPhCS::GetCS: Current T="<<pEn<<" < Threshold="<<lastTH<<", CS=0"<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return 0.;                         // Momentum is below the Threshold Value -> CS=0
  }
  else if(lastP == pMom)
  {
#ifdef debug
    G4cout<<"G4QPhCS::GetCS:OldCur P="<<pMom<<"="<<pMom<<", CS="<<lastCS*millibarn<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return lastCS*millibarn;     // Use theLastCS
  }
  else
  {
#ifdef debug
    G4cout<<"G4QPhCS::GetCS:UpdatCur P="<<pMom<<",f="<<fCS<<",I="<<lastI<<",j="<<j<<G4endl;
#endif
    lastCS=CalculateCrossSection(fCS,1,j,lastPDG,lastZ,lastN,pMom); // Only UpdateDB
    lastP=pMom;
  }
#ifdef debug
  G4cout<<"G4QPhCS::GetCroSec:End,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
  //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  return lastCS*millibarn;
}

// Gives the threshold energy for different nuclei (min of p- and n-threshold)
// *******************************************************************************
// *** This function is the same for all lepto- & photo-nuclear reactions, for ***
// *** (nu,l) reactions the mass value of the final state lepton must be added ***
// ***@@ IT IS REASONABLE TO MAKE ADDITIONAL VIRTUAL CLASS FOR LEPTO-NUCLEAR @@***
// *******************************************************************************
G4double G4QPhotonNuclearCrossSection::ThresholdEnergy(G4int Z, G4int N, G4int)
{
  // CHIPS - Direct GEANT
  //static const G4double mNeut = G4QPDGCode(2112).GetMass();
  //static const G4double mProt = G4QPDGCode(2212).GetMass();
  static const G4double mNeut = G4NucleiProperties::GetNuclearMass(1,0)/MeV;
  static const G4double mProt = G4NucleiProperties::GetNuclearMass(1,1)/MeV;
  static const G4double mAlph = G4NucleiProperties::GetNuclearMass(4,2)/MeV;
  // ---------
  static const G4double infEn = 9.e27;

  G4int A=Z+N;
  if(A<1) return infEn;
  else if(A==1) return 134.9766; // Pi0 threshold for the nucleon
  // CHIPS - Direct GEANT
  //G4double mT= G4QPDGCode(111).GetNuclMass(Z,N,0);
  G4double mT= 0.;
  if(Z&&G4NucleiProperties::IsInStableTable(A,Z))
                           mT = G4NucleiProperties::GetNuclearMass(A,Z)/MeV;
  else return 0.;       // If it is not in the Table of Stable Nuclei, then the Threshold=0
  G4double mP= infEn;
  if(A>1&&Z&&G4NucleiProperties::IsInStableTable(A-1,Z-1))
          mP = G4NucleiProperties::GetNuclearMass(A-1,Z-1)/MeV;// ResNucMass for a proton

  G4double mN= infEn;
  if(A>1&&G4NucleiProperties::IsInStableTable(A-1,Z))
          mN = G4NucleiProperties::GetNuclearMass(A-1,Z)/MeV;  // ResNucMass for a neutron

  G4double mA= infEn;
  if(A>4&&Z>1&&G4NucleiProperties::IsInStableTable(A-4,Z-2))
          mA=G4NucleiProperties::GetNuclearMass(A-4,Z-2)/MeV;  // ResNucMass for an alpha

  G4double dP= mP +mProt - mT;
  G4double dN= mN +mNeut - mT;
  G4double dA= mA +mAlph - mT;
#ifdef debug
  G4cout<<"G4QPhotoNucCS::ThreshEn: mP="<<mP<<",dP="<<dP<<",mN="<<mN<<",dN="<<dN<<",mA="
        <<mA<<",dA="<<dA<<",mT="<<mT<<",A="<<A<<",Z="<<Z<<G4endl;
#endif
  if(dP<dN)dN=dP;
  if(dA<dN)dN=dA;
  return dN;
}

// The main member function giving the gamma-A cross section (E in GeV, CS in mb)
G4double G4QPhotonNuclearCrossSection::CalculateCrossSection(G4bool CS, G4int F, G4int I,
                                          G4int, G4int targZ, G4int targN, G4double Energy)
{
#ifdef debug
  G4cout<<"G4QPhotonNucCrossSection::CalculateCrossSection: ***Called***"<<G4endl;
#endif
  static const G4double THmin=2.;    // minimum Energy Threshold
  static const G4double dE=1.;       // step for the GDR table
  static const G4int    nL=105;      // A#of GDResonance points in E (1 MeV from 2 to 106)
  static const G4double Emin=THmin+(nL-1)*dE; // minE for the HighE part
  static const G4double Emax=50000.;       // maxE for the HighE part
  static const G4int    nH=224;            // A#of HResonance points in lnE
  static const G4double milE=std::log(Emin);    // Low logarithm energy for the HighE part
  static const G4double malE=std::log(Emax);    // High logarithm energy (each 2.75 %)
  static const G4double dlE=(malE-milE)/(nH-1); // Step in log energy in the HighE part
  //
  //static const G4double shd=1.075-.0023*log(2.);  // HE PomShadowing(D)
  static const G4double shd=1.0734;  // HE PomShadowing(D)
  static const G4double shc=0.072;   // HE Shadowing constant
  static const G4double poc=0.0375;  // HE Pomeron coefficient
  static const G4double pos=16.5;    // HE Pomeron shift
  static const G4double reg=.11;     // HE Reggeon slope
  static const G4double shp=1.075;   // HE PomShadowing(P)
  //
  // Associative memory for acceleration
  static std::vector <G4double> spA; // shadowing coefficients (A-dependent)
  //
  onlyCS=CS;                         // Flag to calculate only CS (not Si/Bi)
#ifdef pdebug
  G4cout<<"G4QPhotonNucCS::CalcCS: P="<<Energy<<", F="<<F<<", I="<<I<<", Z="<<targZ
        <<", N="<<targN<<", onlyCS="<<CS<<",E="<<Energy<<",th="<<THmin<<G4endl;
  if(F==-27) return 0.;
#endif
  //if (Energy<THmin)
  //{
  //  lastE=0.;
  //  lastSig=0.;
#ifdef debug
  //  G4cout<<"---> G4QMuonNucCS::CalcCS: CS=0  as E="<<Energy<<" < "<<THmin<<G4endl;
#endif
  //  return 0.;                       // @@ This can be dangerouse for the heaviest nuc.!
  //}
  G4double sigma=0.;
  G4double A=targN+targZ;
  if(F<=0)                           // This isotope was not the last used isotop
  {
    if(F<0)                          // This isotope was found in DAMDB =-------=> RETRIEVE
    {
      lastGDR=(*GDR)[I];             // Pointer to prepared GDR cross sections
      lastHEN=(*HEN)[I];             // Pointer to prepared High Energy cross sections
      lastSP =spA[I];                // Shadowing coefficient for UHE
    }
    else                             // This isotope wasn't calculated previously => CREATE
    {
      G4double lnA=std::log(A);      // The nucleus is not found in DB. It is new.
      if(A==1.) lastSP=1.;           // The Reggeon shadowing (A=1)
      else      lastSP=A*(1.-shc*lnA); // The Reggeon shadowing
#ifdef debug
      G4cout<<">>>G4QPhotonNuclearCrossSect::CalcCS:lnA="<<lnA<<",lastSP="<<lastSP<<G4endl;
#endif
#ifdef debug3
      if(A==3) G4cout<<"G4QPhotonNuclearCrossSection::CalcCS: lastSP="<<lastSP<<G4endl;
#endif
      lastGDR = new G4double[nL];    // Allocate memory for the new GDR cross sections
      lastHEN = new G4double[nH];    // Allocate memory for the new HEN cross sections
      G4int er=GetFunctions(A,lastGDR,lastHEN);// set newZeroPosition and fill theFunctions
      if(er<1) G4cerr<<"***G4QPhotNucCrosSec::CalcCrossSection: A="<<A<<" failed"<<G4endl;
#ifdef pdebug
      G4cout<<">>G4QPhotNucCS::CalcCS:**GDR/HEN're made for A="<<A<<"** GFEr="<<er<<G4endl;
#endif
      // *** The synchronization check ***
      G4int sync=GDR->size();
      if(sync!=I) G4cerr<<"***G4QPhotoNuclearCS::CalcCS: PDG=22, S="<<sync<<"#"<<I<<G4endl;
      GDR->push_back(lastGDR);       // added GDR, found by AH 10/7/02
      HEN->push_back(lastHEN);       // added HEN, found by AH 10/7/02
      spA.push_back(lastSP);         // Pomeron Shadowing
    } // End of creation of the new set of parameters
  } // End of parameters udate
  // =-------------------= NOW the Magic Formula =--------------------------=
  if (Energy<lastTH) return 0.;             // It must be already checked in the interface
  else if (Energy<=Emin)                    // GDR region (approximated in E, not in lnE)
  {
#ifdef debug
    G4cout<<"G4QPhNCS::CalcCS:bGDR A="<<A<<", nL="<<nL<<",TH="<<THmin<<",dE="<<dE<<G4endl;
#endif
    if(A<=1. || dE <= THmin) sigma=0.;      // No GDR for A=1
    else      sigma=EquLinearFit(Energy,nL,THmin,dE,lastGDR);
#ifdef debugn
    if(sigma<0.)
      G4cout<<"G4QPhoNucCS::CalcCS:A="<<A<<",E="<<Energy<<",T="<<THmin<<",dE="<<dE<<G4endl;
#endif
  }
  else if (Energy<Emax)                     // High Energy region
  {
    G4double lE=std::log(Energy);
#ifdef pdebug
    G4cout<<"G4QPhotNucCS::CalcCS:lE="<<milE<<",n="<<nH<<",iE="<<milE<<",dE="<<dlE<<",HEN="
          <<lastHEN<<",A="<<A<<G4endl;
#endif
    if(lE > milE) sigma=EquLinearFit(lE,nH,milE,dlE,lastHEN);
    else          sigma=0.;
  }
  else                                      // UHE region (calculation, not frequent)
  {
    G4double lE=std::log(Energy);
    G4double sh=shd;
    if(A==1.)sh=shp;
    sigma=lastSP*(poc*(lE-pos)+sh*std::exp(-reg*lE));
  }
#ifdef debug
  G4cout<<"G4QPhotonNuclearCrossSection::CalcCS: sigma="<<sigma<<G4endl;
#endif
#ifdef debug
  if(Energy>45000.&&Energy<60000.)
    G4cout<<"G4QPhotoNucCS::GetCS: A="<<A<<", E="<<Energy<<",CS="<<sigma<<G4endl;
#endif
  if(sigma<0.) return 0.;
  return sigma;
}

// Linear fit for YN[N] tabulated (from X0 with fixed step DX) function to X point

// Calculate the functions for the log(A)
G4int G4QPhotonNuclearCrossSection::GetFunctions(G4double a, G4double* y, G4double* z)
{
  static const G4int nLA=49;       // A#of GDResonance basic nuclei
  static const G4double LA[nLA]={ 2.,4.,6.,7.,9.,12.,14.,15.,16.,19.,23.,24.,27.,28.,32.,
   34.,40.,54., 55.,56.,58.7,58.9,63.5,65.4,76.,82.,107.9,112.4,118.7,126.9,154.,156.,159.,
   165.,168.,174.,178.,180.,181.,184.,186.,197.,204.4,207.2,209.,232.,235.,238.,239.};
  static const G4int nL=105;       // A#of GDResonance points in E (each MeV from 2 to 106)
  static const G4int nHA=14;       // A#of HResonance basic nuclei
  static const G4double HA[nHA]={1.,2.,3.,4.,6.,7.,9.,12.,16.,27.,63.5,118.7,207.2,238.};
  static const G4int nH=224;       // A#of HResonance points in lnE (each 2.75 percents)
  // If the cross section approximation formula is changed - replace from file.
  static const G4double SL0[nL]={
    7.094260e-1,1.532987e+0,2.449381e+0,2.785790e+0,2.525673e+0,2.128172e+0,1.780549e+0,
    1.506934e+0,1.294560e+0,1.128048e+0,9.953850e-1,8.879274e-1,7.995356e-1,7.258111e-1,
    6.635555e-1,6.104038e-1,5.645786e-1,5.247229e-1,4.897864e-1,4.589445e-1,4.315429e-1,
    4.070560e-1,3.850576e-1,3.651990e-1,3.471920e-1,3.307971e-1,3.158133e-1,3.020711e-1,
    2.894266e-1,2.777569e-1,2.669563e-1,2.569336e-1,2.476099e-1,2.389161e-1,2.307920e-1,
    2.231848e-1,2.160475e-1,2.093390e-1,2.030225e-1,1.970653e-1,1.914383e-1,1.861152e-1,
    1.810725e-1,1.762891e-1,1.717459e-1,1.674254e-1,1.633120e-1,1.593914e-1,1.556505e-1,
    1.520775e-1,1.486616e-1,1.453926e-1,1.422615e-1,1.392599e-1,1.363800e-1,1.336147e-1,
    1.309573e-1,1.284017e-1,1.259423e-1,1.235738e-1,1.212914e-1,1.190904e-1,1.169666e-1,
    1.149161e-1,1.129353e-1,1.110206e-1,1.091688e-1,1.073770e-1,1.056423e-1,1.039619e-1,
    1.023336e-1,1.007548e-1,9.922335e-2,9.773724e-2,9.629446e-2,9.489316e-2,9.353161e-2,
    9.220814e-2,9.092120e-2,8.966931e-2,8.845106e-2,8.726514e-2,8.611027e-2,8.498527e-2,
    8.388900e-2,8.282039e-2,8.177841e-2,8.076208e-2,7.977047e-2,7.880271e-2,7.785794e-2,
    7.693536e-2,7.603421e-2,7.515376e-2,7.429330e-2,7.345216e-2,7.262971e-2,7.182534e-2,
    7.103847e-2,7.026852e-2,6.951498e-2,6.877732e-2,6.805505e-2,6.734772e-2,6.665486e-2};
  static const G4double SL1[nL]={
    2.017310e-4,9.866847e-4,3.081371e-3,7.486476e-3,1.550083e-2,2.873865e-2,4.915763e-2,
    7.909925e-2,1.213409e-1,1.791563e-1,2.563778e-1,3.574504e-1,4.874600e-1,6.521048e-1,
    8.575237e-1,1.109763e+0,1.413389e+0,1.768398e+0,2.164804e+0,2.576439e+0,2.960166e+0,
    3.267863e+0,3.467401e+0,3.555207e+0,3.550465e+0,3.480379e+0,3.369351e+0,3.235278e+0,
    3.090040e+0,2.941162e+0,2.793315e+0,2.649362e+0,2.511013e+0,2.379237e+0,2.254508e+0,
    2.136980e+0,2.026587e+0,1.923127e+0,1.826306e+0,1.735781e+0,1.651181e+0,1.572129e+0,
    1.498250e+0,1.429182e+0,1.364578e+0,1.304111e+0,1.247474e+0,1.194383e+0,1.144574e+0,
    1.097802e+0,1.053842e+0,1.012488e+0,9.735481e-1,9.368486e-1,9.022288e-1,8.695413e-1,
    8.386509e-1,8.094332e-1,7.817741e-1,7.555686e-1,7.307199e-1,7.071392e-1,6.847445e-1,
    6.634600e-1,6.432161e-1,6.239481e-1,6.055964e-1,5.881058e-1,5.714251e-1,5.555068e-1,
    5.403069e-1,5.257844e-1,5.119013e-1,4.986222e-1,4.859140e-1,4.737459e-1,4.620893e-1,
    4.509173e-1,4.402048e-1,4.299283e-1,4.200659e-1,4.105971e-1,4.015025e-1,3.927640e-1,
    3.843647e-1,3.762885e-1,3.685206e-1,3.610467e-1,3.538537e-1,3.469290e-1,3.402609e-1,
    3.338384e-1,3.276510e-1,3.216889e-1,3.159429e-1,3.104042e-1,3.050648e-1,2.999168e-1,
    2.949531e-1,2.901668e-1,2.855516e-1,2.811014e-1,2.768107e-1,2.726743e-1,2.686872e-1};
  static const G4double SL2[nL]={
    4.776434e-3,2.412116e-2,7.595870e-2,1.835144e-1,3.703569e-1,6.466818e-1,9.877908e-1,
    1.324697e+0,1.575559e+0,1.699764e+0,1.715038e+0,1.669943e+0,1.614318e+0,1.588675e+0,
    1.625472e+0,1.751419e+0,1.984310e+0,2.321229e+0,2.721619e+0,3.102285e+0,3.366107e+0,
    3.454042e+0,3.371971e+0,3.170282e+0,2.907055e+0,2.626333e+0,2.354348e+0,2.103834e+0,
    1.879323e+0,1.681008e+0,1.507070e+0,1.354935e+0,1.221899e+0,1.105415e+0,1.003205e+0,
    9.132844e-1,8.339468e-1,7.637380e-1,7.014215e-1,6.459465e-1,5.964196e-1,5.520797e-1,
    5.122769e-1,4.764549e-1,4.441358e-1,4.149083e-1,3.884171e-1,3.643543e-1,3.424525e-1,
    3.224787e-1,3.042291e-1,2.875255e-1,2.722111e-1,2.581479e-1,2.452141e-1,2.333017e-1,
    2.223151e-1,2.121693e-1,2.027886e-1,1.941052e-1,1.860586e-1,1.785946e-1,1.716644e-1,
    1.652242e-1,1.592345e-1,1.536595e-1,1.484670e-1,1.436278e-1,1.391153e-1,1.349053e-1,
    1.309760e-1,1.273073e-1,1.238809e-1,1.206803e-1,1.176900e-1,1.148962e-1,1.122861e-1,
    1.098477e-1,1.075703e-1,1.054441e-1,1.034596e-1,1.016087e-1,9.988340e-2,9.827659e-2,
    9.678165e-2,9.539245e-2,9.410337e-2,9.290919e-2,9.180512e-2,9.078673e-2,8.984997e-2,
    8.899108e-2,8.820664e-2,8.749353e-2,8.684888e-2,8.627010e-2,8.575488e-2,8.530112e-2,
    8.490697e-2,8.457084e-2,8.429132e-2,8.406729e-2,8.389779e-2,8.378214e-2,8.371985e-2};
  static const G4double SL3[nL]={
    1.375991e-4,6.420490e-4,2.009594e-3,5.073626e-3,1.137383e-2,2.408187e-2,5.091978e-2,
    1.151175e-1,2.955817e-1,8.132651e-1,1.635125e+0,1.931572e+0,2.185333e+0,2.701264e+0,
    3.269689e+0,3.632210e+0,3.708366e+0,3.594398e+0,3.418556e+0,3.260141e+0,3.149899e+0,
    3.091216e+0,3.075568e+0,3.090472e+0,3.123146e+0,3.162196e+0,3.198373e+0,3.224873e+0,
    3.237305e+0,3.233448e+0,3.212852e+0,3.176382e+0,3.125768e+0,3.063213e+0,2.991084e+0,
    2.911700e+0,2.827186e+0,2.739409e+0,2.649943e+0,2.560078e+0,2.470840e+0,2.383021e+0,
    2.297216e+0,2.213858e+0,2.133241e+0,2.055557e+0,1.980911e+0,1.909346e+0,1.840852e+0,
    1.775386e+0,1.712877e+0,1.653235e+0,1.596357e+0,1.542133e+0,1.490448e+0,1.441186e+0,
    1.394230e+0,1.349469e+0,1.306789e+0,1.266085e+0,1.227254e+0,1.190196e+0,1.154820e+0,
    1.121035e+0,1.088758e+0,1.057908e+0,1.028411e+0,1.000196e+0,9.731948e-1,9.473457e-1,
    9.225887e-1,8.988681e-1,8.761312e-1,8.543286e-1,8.334136e-1,8.133424e-1,7.940736e-1,
    7.755685e-1,7.577903e-1,7.407046e-1,7.242788e-1,7.084823e-1,6.932862e-1,6.786633e-1,
    6.645878e-1,6.510355e-1,6.379834e-1,6.254100e-1,6.132949e-1,6.016187e-1,5.903633e-1,
    5.795116e-1,5.690472e-1,5.589548e-1,5.492201e-1,5.398295e-1,5.307700e-1,5.220296e-1,
    5.135969e-1,5.054613e-1,4.976128e-1,4.900419e-1,4.827400e-1,4.756989e-1,4.689110e-1};
  static const G4double SL4[nL]={
    1.531367e-4,6.750684e-4,2.023434e-3,4.818832e-3,9.866691e-3,1.816857e-2,3.094217e-2,
    4.965477e-2,7.607934e-2,1.123974e-1,1.614108e-1,2.270208e-1,3.153403e-1,4.372460e-1,
    6.139880e-1,8.886525e-1,1.345605e+0,2.121366e+0,3.298049e+0,4.533310e+0,5.172459e+0,
    5.243522e+0,5.175754e+0,5.149633e+0,5.156364e+0,5.151144e+0,5.108382e+0,5.025027e+0,
    4.909480e+0,4.772279e+0,4.621981e+0,4.464473e+0,4.303590e+0,4.141874e+0,3.981115e+0,
    3.822656e+0,3.667551e+0,3.516631e+0,3.370536e+0,3.229738e+0,3.094556e+0,2.965180e+0,
    2.841688e+0,2.724066e+0,2.612228e+0,2.506035e+0,2.405305e+0,2.309830e+0,2.219381e+0,
    2.133721e+0,2.052608e+0,1.975802e+0,1.903066e+0,1.834170e+0,1.768894e+0,1.707024e+0,
    1.648361e+0,1.592714e+0,1.539903e+0,1.489759e+0,1.442123e+0,1.396846e+0,1.353788e+0,
    1.312819e+0,1.273817e+0,1.236668e+0,1.201266e+0,1.167510e+0,1.135308e+0,1.104573e+0,
    1.075223e+0,1.047183e+0,1.020382e+0,9.947538e-1,9.702356e-1,9.467696e-1,9.243013e-1,
    9.027797e-1,8.821569e-1,8.623879e-1,8.434307e-1,8.252457e-1,8.077958e-1,7.910459e-1,
    7.749634e-1,7.595173e-1,7.446784e-1,7.304196e-1,7.167151e-1,7.035406e-1,6.908733e-1,
    6.786917e-1,6.669756e-1,6.557060e-1,6.448649e-1,6.344356e-1,6.244022e-1,6.147497e-1,
    6.054644e-1,5.965332e-1,5.879438e-1,5.796850e-1,5.717461e-1,5.641173e-1,5.567897e-1};
  static const G4double SL5[nL]={
    1.905569e-4,7.771730e-4,2.250919e-3,5.273053e-3,1.071640e-2,1.969996e-2,3.365091e-2,
    5.440813e-2,8.439169e-2,1.268914e-1,1.866020e-1,2.707115e-1,3.912405e-1,5.701376e-1,
    8.501724e-1,1.317340e+0,2.143911e+0,3.657987e+0,6.387255e+0,1.074352e+1,1.571664e+1,
    1.840405e+1,1.776700e+1,1.557514e+1,1.329204e+1,1.138076e+1,9.874227e+0,8.700723e+0,
    7.781216e+0,7.050490e+0,6.458855e+0,5.969695e+0,5.556515e+0,5.200371e+0,4.887807e+0,
    4.609287e+0,4.358030e+0,4.129172e+0,3.919172e+0,3.725403e+0,3.545861e+0,3.378977e+0,
    3.223486e+0,3.078336e+0,2.942636e+0,2.815610e+0,2.696573e+0,2.584914e+0,2.480080e+0,
    2.381572e+0,2.288930e+0,2.201736e+0,2.119606e+0,2.042187e+0,1.969152e+0,1.900204e+0,
    1.835064e+0,1.773480e+0,1.715215e+0,1.660054e+0,1.607794e+0,1.558252e+0,1.511256e+0,
    1.466647e+0,1.424278e+0,1.384013e+0,1.345726e+0,1.309299e+0,1.274624e+0,1.241600e+0,
    1.210133e+0,1.180134e+0,1.151523e+0,1.124223e+0,1.098163e+0,1.073277e+0,1.049504e+0,
    1.026785e+0,1.005065e+0,9.842955e-1,9.644275e-1,9.454166e-1,9.272213e-1,9.098022e-1,
    8.931224e-1,8.771474e-1,8.618447e-1,8.471837e-1,8.331356e-1,8.196734e-1,8.067717e-1,
    7.944065e-1,7.825554e-1,7.711972e-1,7.603122e-1,7.498817e-1,7.398883e-1,7.303156e-1,
    7.211483e-1,7.123722e-1,7.039741e-1,6.959417e-1,6.882635e-1,6.809293e-1,6.739294e-1};
  static const G4double SL6[nL]={
    2.222448e-4,8.620556e-4,2.444896e-3,5.705453e-3,1.171159e-2,2.205349e-2,3.918281e-2,
    6.696997e-2,1.115720e-1,1.827533e-1,2.959155e-1,4.753435e-1,7.596938e-1,1.211738e+0,
    1.936099e+0,3.111254e+0,5.043478e+0,8.232698e+0,1.330416e+1,2.019140e+1,2.638709e+1,
    2.859878e+1,2.728600e+1,2.464338e+1,2.186072e+1,1.931943e+1,1.710886e+1,1.522482e+1,
    1.363142e+1,1.228407e+1,1.113950e+1,1.015995e+1,9.314220e+0,8.577403e+0,7.929898e+0,
    7.356396e+0,6.844921e+0,6.386042e+0,5.972257e+0,5.597518e+0,5.256885e+0,4.946267e+0,
    4.662228e+0,4.401855e+0,4.162646e+0,3.942438e+0,3.739345e+0,3.551711e+0,3.378077e+0,
    3.217151e+0,3.067782e+0,2.928947e+0,2.799729e+0,2.679306e+0,2.566940e+0,2.461965e+0,
    2.363783e+0,2.271852e+0,2.185682e+0,2.104827e+0,2.028885e+0,1.957488e+0,1.890304e+0,
    1.827026e+0,1.767378e+0,1.711104e+0,1.657972e+0,1.607769e+0,1.560299e+0,1.515382e+0,
    1.472853e+0,1.432558e+0,1.394357e+0,1.358121e+0,1.323729e+0,1.291070e+0,1.260042e+0,
    1.230549e+0,1.202502e+0,1.175819e+0,1.150425e+0,1.126246e+0,1.103218e+0,1.081279e+0,
    1.060370e+0,1.040438e+0,1.021433e+0,1.003308e+0,9.860183e-1,9.695234e-1,9.537847e-1,
    9.387662e-1,9.244342e-1,9.107573e-1,8.977058e-1,8.852523e-1,8.733710e-1,8.620378e-1,
    8.512302e-1,8.409275e-1,8.311102e-1,8.217603e-1,8.128613e-1,8.043977e-1,7.963557e-1};
  static const G4double SL7[nL]={
    2.400132e-4,9.082999e-4,2.545511e-3,5.912609e-3,1.214175e-2,2.297237e-2,4.117454e-2,
    7.124517e-2,1.204927e-1,2.006898e-1,3.306145e-1,5.401144e-1,8.769596e-1,1.418938e+0,
    2.295653e+0,3.727247e+0,6.087430e+0,9.967584e+0,1.601149e+1,2.371712e+1,2.968929e+1,
    3.091057e+1,2.878613e+1,2.564897e+1,2.255784e+1,1.981778e+1,1.747914e+1,1.551277e+1,
    1.386591e+1,1.248298e+1,1.131384e+1,1.031654e+1,9.457362e+0,8.709876e+0,8.053561e+0,
    7.472543e+0,6.954497e+0,6.489768e+0,6.070707e+0,5.691171e+0,5.346149e+0,5.031499e+0,
    4.743750e+0,4.479955e+0,4.237588e+0,4.014464e+0,3.808676e+0,3.618553e+0,3.442619e+0,
    3.279569e+0,3.128239e+0,2.987593e+0,2.856703e+0,2.734737e+0,2.620949e+0,2.514665e+0,
    2.415276e+0,2.322236e+0,2.235046e+0,2.153255e+0,2.076456e+0,2.004275e+0,1.936374e+0,
    1.872444e+0,1.812202e+0,1.755391e+0,1.701773e+0,1.651134e+0,1.603273e+0,1.558007e+0,
    1.515169e+0,1.474604e+0,1.436168e+0,1.399730e+0,1.365168e+0,1.332368e+0,1.301226e+0,
    1.271645e+0,1.243535e+0,1.216812e+0,1.191399e+0,1.167222e+0,1.144214e+0,1.122313e+0,
    1.101460e+0,1.081599e+0,1.062680e+0,1.044655e+0,1.027479e+0,1.011109e+0,9.955079e-1,
    9.806377e-1,9.664643e-1,9.529557e-1,9.400818e-1,9.278146e-1,9.161277e-1,9.049967e-1,
    8.943989e-1,8.843129e-1,8.747192e-1,8.655994e-1,8.569368e-1,8.487159e-1,8.409226e-1};
  static const G4double SL8[nL]={
    2.590923e-4,9.573672e-4,2.651275e-3,6.130118e-3,1.259782e-2,2.396311e-2,4.335926e-2,
    7.599430e-2,1.304582e-1,2.206539e-1,3.685911e-1,6.084138e-1,9.922345e-1,1.598590e+0,
    2.544422e+0,4.001018e+0,6.212589e+0,9.507804e+0,1.423676e+1,2.030771e+1,2.598385e+1,
    2.841920e+1,2.739643e+1,2.481830e+1,2.199791e+1,1.941582e+1,1.718791e+1,1.530626e+1,
    1.372572e+1,1.239454e+1,1.126539e+1,1.029864e+1,9.462554e+0,8.732318e+0,8.088729e+0,
    7.516956e+0,7.005488e+0,6.545306e+0,6.129254e+0,5.751565e+0,5.407521e+0,5.093205e+0,
    4.805315e+0,4.541037e+0,4.297946e+0,4.073932e+0,3.867148e+0,3.675966e+0,3.498945e+0,
    3.334802e+0,3.182394e+0,3.040697e+0,2.908794e+0,2.785859e+0,2.671150e+0,2.563996e+0,
    2.463791e+0,2.369986e+0,2.282085e+0,2.199636e+0,2.122228e+0,2.049489e+0,1.981076e+0,
    1.916681e+0,1.856017e+0,1.798827e+0,1.744870e+0,1.693929e+0,1.645803e+0,1.600307e+0,
    1.557271e+0,1.516539e+0,1.477966e+0,1.441418e+0,1.406772e+0,1.373914e+0,1.342737e+0,
    1.313142e+0,1.285040e+0,1.258344e+0,1.232976e+0,1.208862e+0,1.185934e+0,1.164128e+0,
    1.143384e+0,1.123646e+0,1.104863e+0,1.086985e+0,1.069967e+0,1.053766e+0,1.038344e+0,
    1.023662e+0,1.009685e+0,9.963805e-1,9.837187e-1,9.716705e-1,9.602093e-1,9.493103e-1,
    9.389503e-1,9.291078e-1,9.197629e-1,9.108970e-1,9.024933e-1,8.945360e-1,8.870112e-1};
  static const G4double SL9[nL]={
    3.243985e-4,1.122034e-3,3.000932e-3,6.850212e-3,1.414720e-2,2.751937e-2,5.204925e-2,
    9.887958e-2,1.966468e-1,4.282973e-1,1.041076e+0,2.706630e+0,6.509565e+0,1.085114e+1,
    1.162472e+1,1.124054e+1,1.202416e+1,1.402207e+1,1.659634e+1,1.891975e+1,2.032292e+1,
    2.059083e+1,1.993672e+1,1.873926e+1,1.732572e+1,1.590211e+1,1.457097e+1,1.336993e+1,
    1.230272e+1,1.135820e+1,1.052046e+1,9.773672e+0,9.103884e+0,8.499562e+0,7.951408e+0,
    7.451996e+0,6.995366e+0,6.576679e+0,6.191930e+0,5.837748e+0,5.511235e+0,5.209864e+0,
    4.931401e+0,4.673850e+0,4.435420e+0,4.214488e+0,4.009588e+0,3.819384e+0,3.642664e+0,
    3.478323e+0,3.325356e+0,3.182848e+0,3.049964e+0,2.925943e+0,2.810093e+0,2.701782e+0,
    2.600432e+0,2.505518e+0,2.416558e+0,2.333114e+0,2.254783e+0,2.181197e+0,2.112021e+0,
    2.046943e+0,1.985682e+0,1.927976e+0,1.873586e+0,1.822292e+0,1.773891e+0,1.728195e+0,
    1.685032e+0,1.644242e+0,1.605677e+0,1.569201e+0,1.534686e+0,1.502017e+0,1.471082e+0,
    1.441781e+0,1.414020e+0,1.387711e+0,1.362772e+0,1.339127e+0,1.316705e+0,1.295438e+0,
    1.275266e+0,1.256130e+0,1.237976e+0,1.220753e+0,1.204413e+0,1.188912e+0,1.174209e+0,
    1.160265e+0,1.147042e+0,1.134507e+0,1.122628e+0,1.111376e+0,1.100721e+0,1.090639e+0,
    1.081106e+0,1.072098e+0,1.063597e+0,1.055582e+0,1.048036e+0,1.040943e+0,1.034290e+0};
  static const G4double SL10[nL]={
    4.311217e-4,1.384716e-3,3.549518e-3,7.988549e-3,1.667330e-2,3.341344e-2,6.552895e-2,
    1.266167e-1,2.409191e-1,4.501490e-1,8.243911e-1,1.480280e+0,2.612343e+0,4.545249e+0,
    7.790746e+0,1.287033e+1,1.909053e+1,2.392952e+1,2.652790e+1,2.742592e+1,2.690891e+1,
    2.536786e+1,2.330746e+1,2.113183e+1,1.907371e+1,1.723144e+1,1.562538e+1,1.423904e+1,
    1.304262e+1,1.200459e+1,1.109667e+1,1.029534e+1,9.581841e+0,8.941546e+0,8.363124e+0,
    7.837784e+0,7.358628e+0,6.920153e+0,6.517878e+0,6.148072e+0,5.807568e+0,5.493625e+0,
    5.203837e+0,4.936070e+0,4.688412e+0,4.459144e+0,4.246710e+0,4.049702e+0,3.866841e+0,
    3.696964e+0,3.539013e+0,3.392026e+0,3.255128e+0,3.127519e+0,3.008473e+0,2.897326e+0,
    2.793474e+0,2.696364e+0,2.605493e+0,2.520398e+0,2.440659e+0,2.365890e+0,2.295737e+0,
    2.229875e+0,2.168007e+0,2.109859e+0,2.055179e+0,2.003736e+0,1.955316e+0,1.909723e+0,
    1.866773e+0,1.826299e+0,1.788145e+0,1.752167e+0,1.718232e+0,1.686214e+0,1.655999e+0,
    1.627479e+0,1.600555e+0,1.575133e+0,1.551128e+0,1.528457e+0,1.507047e+0,1.486825e+0,
    1.467726e+0,1.449689e+0,1.432656e+0,1.416572e+0,1.401389e+0,1.387057e+0,1.373533e+0,
    1.360776e+0,1.348747e+0,1.337409e+0,1.326730e+0,1.316677e+0,1.307222e+0,1.298337e+0,
    1.289997e+0,1.282179e+0,1.274863e+0,1.268027e+0,1.261656e+0,1.255732e+0,1.250242e+0};
  static const G4double SL11[nL]={
    4.614524e-4,1.458509e-3,3.702639e-3,8.309380e-3,1.740590e-2,3.519535e-2,6.986551e-2,
    1.367187e-1,2.630019e-1,4.950763e-1,9.087988e-1,1.624204e+0,2.825210e+0,4.782440e+0,
    7.867272e+0,1.250247e+1,1.878669e+1,2.530271e+1,2.928727e+1,3.015114e+1,2.903038e+1,
    2.689359e+1,2.438858e+1,2.190927e+1,1.964845e+1,1.767000e+1,1.597064e+1,1.451878e+1,
    1.327514e+1,1.220222e+1,1.126790e+1,1.044614e+1,9.716524e+0,9.063269e+0,8.474255e+0,
    7.940129e+0,7.453592e+0,7.008848e+0,6.601198e+0,6.226749e+0,5.882205e+0,5.564731e+0,
    5.271842e+0,5.001344e+0,4.751274e+0,4.519872e+0,4.305549e+0,4.106867e+0,3.922524e+0,
    3.751336e+0,3.592227e+0,3.444221e+0,3.306426e+0,3.178034e+0,3.058307e+0,2.946572e+0,
    2.842215e+0,2.744679e+0,2.653450e+0,2.568064e+0,2.488092e+0,2.413143e+0,2.342861e+0,
    2.276915e+0,2.215004e+0,2.156852e+0,2.102204e+0,2.050824e+0,2.002497e+0,1.957023e+0,
    1.914217e+0,1.873910e+0,1.835942e+0,1.800170e+0,1.766455e+0,1.734674e+0,1.704709e+0,
    1.676451e+0,1.649799e+0,1.624659e+0,1.600943e+0,1.578569e+0,1.557460e+0,1.537545e+0,
    1.518758e+0,1.501034e+0,1.484317e+0,1.468550e+0,1.453684e+0,1.439670e+0,1.426463e+0,
    1.414022e+0,1.402307e+0,1.391282e+0,1.380912e+0,1.371166e+0,1.362014e+0,1.353430e+0,
    1.345386e+0,1.337862e+0,1.330834e+0,1.324283e+0,1.318193e+0,1.312546e+0,1.307330e+0};
  static const G4double SL12[nL]={
    5.615148e-4,1.700309e-3,4.203181e-3,9.368359e-3,1.987519e-2,4.133574e-2,8.507565e-2,
    1.726852e-1,3.430025e-1,6.623201e-1,1.238631e+0,2.240098e+0,3.915001e+0,6.601693e+0,
    1.070034e+1,1.656745e+1,2.430795e+1,3.323297e+1,4.042222e+1,4.203499e+1,3.892326e+1,
    3.426056e+1,2.971854e+1,2.578645e+1,2.251925e+1,1.984150e+1,1.764928e+1,1.584399e+1,
    1.434242e+1,1.307807e+1,1.199938e+1,1.106710e+1,1.025167e+1,9.530877e+0,8.888043e+0,
    8.310522e+0,7.788614e+0,7.314751e+0,6.882912e+0,6.488204e+0,6.126573e+0,5.794594e+0,
    5.489325e+0,5.208212e+0,4.949006e+0,4.709718e+0,4.488572e+0,4.283980e+0,4.094514e+0,
    3.918886e+0,3.755935e+0,3.604610e+0,3.463959e+0,3.333120e+0,3.211310e+0,3.097817e+0,
    2.991993e+0,2.893250e+0,2.801049e+0,2.714902e+0,2.634361e+0,2.559016e+0,2.488493e+0,
    2.422448e+0,2.360566e+0,2.302559e+0,2.248160e+0,2.197125e+0,2.149229e+0,2.104262e+0,
    2.062033e+0,2.022365e+0,1.985093e+0,1.950064e+0,1.917138e+0,1.886183e+0,1.857078e+0,
    1.829709e+0,1.803970e+0,1.779763e+0,1.756997e+0,1.735586e+0,1.715451e+0,1.696516e+0,
    1.678712e+0,1.661974e+0,1.646241e+0,1.631456e+0,1.617565e+0,1.604520e+0,1.592272e+0,
    1.580780e+0,1.570002e+0,1.559901e+0,1.550441e+0,1.541590e+0,1.533317e+0,1.525595e+0,
    1.518397e+0,1.511701e+0,1.505485e+0,1.499729e+0,1.494416e+0,1.489531e+0,1.485061e+0};
  static const G4double SL13[nL]={
    5.979521e-4,1.787895e-3,4.384312e-3,9.755476e-3,2.079561e-2,4.366898e-2,9.094059e-2,
    1.867226e-1,3.746609e-1,7.299098e-1,1.376720e+0,2.513601e+0,4.446871e+0,7.627694e+0,
    1.267423e+1,2.032656e+1,3.102537e+1,4.279863e+1,4.924268e+1,4.764583e+1,4.223031e+1,
    3.635559e+1,3.114336e+1,2.680306e+1,2.327037e+1,2.041256e+1,1.809450e+1,1.619914e+1,
    1.463175e+1,1.331842e+1,1.220266e+1,1.124188e+1,1.040420e+1,9.665816e+0,9.008883e+0,
    8.419927e+0,7.888645e+0,7.407021e+0,6.968700e+0,6.568538e+0,6.202287e+0,5.866370e+0,
    5.557733e+0,5.273725e+0,5.012027e+0,4.770588e+0,4.547585e+0,4.341389e+0,4.150539e+0,
    3.973720e+0,3.809747e+0,3.657550e+0,3.516159e+0,3.384698e+0,3.262371e+0,3.148455e+0,
    3.042294e+0,2.943289e+0,2.850895e+0,2.764617e+0,2.683999e+0,2.608629e+0,2.538126e+0,
    2.472142e+0,2.410358e+0,2.352483e+0,2.298245e+0,2.247399e+0,2.199714e+0,2.154982e+0,
    2.113006e+0,2.073607e+0,2.036619e+0,2.001888e+0,1.969270e+0,1.938632e+0,1.909852e+0,
    1.882814e+0,1.857411e+0,1.833544e+0,1.811120e+0,1.790052e+0,1.770261e+0,1.751669e+0,
    1.734206e+0,1.717807e+0,1.702410e+0,1.687957e+0,1.674395e+0,1.661672e+0,1.649742e+0,
    1.638561e+0,1.628089e+0,1.618286e+0,1.609117e+0,1.600551e+0,1.592555e+0,1.585102e+0,
    1.578167e+0,1.571726e+0,1.565757e+0,1.560241e+0,1.555161e+0,1.550502e+0,1.546251e+0};
  static const G4double SL14[nL]={
    7.595609e-4,2.174487e-3,5.184472e-3,1.148979e-2,2.501660e-2,5.458957e-2,1.187206e-1,
    2.534357e-1,5.242273e-1,1.043266e+0,1.992371e+0,3.648981e+0,6.401444e+0,1.071384e+1,
    1.696937e+1,2.517085e+1,3.466133e+1,4.362657e+1,4.818786e+1,4.632665e+1,4.110517e+1,
    3.547368e+1,3.049881e+1,2.637019e+1,2.301312e+1,2.029373e+1,1.808104e+1,1.626396e+1,
    1.475365e+1,1.348118e+1,1.239409e+1,1.145289e+1,1.062807e+1,9.897589e+0,9.244965e+0,
    8.657732e+0,8.126345e+0,7.643355e+0,7.202834e+0,6.799954e+0,6.430699e+0,6.091668e+0,
    5.779927e+0,5.492917e+0,5.228379e+0,4.984303e+0,4.758895e+0,4.550538e+0,4.357779e+0,
    4.179304e+0,4.013924e+0,3.860561e+0,3.718237e+0,3.586064e+0,3.463234e+0,3.349014e+0,
    3.242733e+0,3.143782e+0,3.051603e+0,2.965688e+0,2.885571e+0,2.810827e+0,2.741066e+0,
    2.675929e+0,2.615087e+0,2.558239e+0,2.505107e+0,2.455433e+0,2.408981e+0,2.365534e+0,
    2.324889e+0,2.286859e+0,2.251272e+0,2.217967e+0,2.186795e+0,2.157618e+0,2.130308e+0,
    2.104744e+0,2.080815e+0,2.058418e+0,2.037456e+0,2.017838e+0,1.999480e+0,1.982303e+0,
    1.966234e+0,1.951205e+0,1.937150e+0,1.924011e+0,1.911731e+0,1.900259e+0,1.889545e+0,
    1.879545e+0,1.870216e+0,1.861520e+0,1.853420e+0,1.845884e+0,1.838880e+0,1.832382e+0,
    1.826362e+0,1.820798e+0,1.815670e+0,1.810960e+0,1.806650e+0,1.802729e+0,1.799183e+0};
  static const G4double SL15[nL]={
    8.500963e-4,2.390172e-3,5.632030e-3,1.247632e-2,2.747950e-2,6.109914e-2,1.355108e-1,
    2.941224e-1,6.161245e-1,1.237476e+0,2.378852e+0,4.376594e+0,7.697785e+0,1.288755e+1,
    2.037233e+1,3.017649e+1,4.195065e+1,5.485775e+1,6.453653e+1,6.432845e+1,5.643212e+1,
    4.707792e+1,3.899572e+1,3.258608e+1,2.760535e+1,2.373102e+1,2.068965e+1,1.827124e+1,
    1.631941e+1,1.471906e+1,1.338580e+1,1.225785e+1,1.128994e+1,1.044881e+1,9.709837e+0,
    9.054659e+0,8.469380e+0,7.943327e+0,7.468155e+0,7.037214e+0,6.645102e+0,6.287353e+0,
    5.960212e+0,5.660481e+0,5.385398e+0,5.132559e+0,4.899852e+0,4.685413e+0,4.487583e+0,
    4.304884e+0,4.135994e+0,3.979725e+0,3.835009e+0,3.700883e+0,3.576477e+0,3.461003e+0,
    3.353748e+0,3.254064e+0,3.161362e+0,3.075105e+0,2.994806e+0,2.920017e+0,2.850330e+0,
    2.785372e+0,2.724800e+0,2.668301e+0,2.615584e+0,2.566383e+0,2.520455e+0,2.477573e+0,
    2.437529e+0,2.400129e+0,2.365194e+0,2.332560e+0,2.302073e+0,2.273590e+0,2.246978e+0,
    2.222116e+0,2.198888e+0,2.177188e+0,2.156916e+0,2.137979e+0,2.120292e+0,2.103773e+0,
    2.088348e+0,2.073945e+0,2.060501e+0,2.047953e+0,2.036245e+0,2.025323e+0,2.015140e+0,
    2.005648e+0,1.996806e+0,1.988574e+0,1.980915e+0,1.973798e+0,1.967190e+0,1.961065e+0,
    1.955397e+0,1.950164e+0,1.945345e+0,1.940923e+0,1.936883e+0,1.933213e+0,1.929901e+0};
  static const G4double SL16[nL]={
    1.161977e-3,3.130797e-3,7.178175e-3,1.596595e-2,3.647036e-2,8.543942e-2,1.991615e-1,
    4.493705e-1,9.672255e-1,1.976461e+0,3.832134e+0,7.044334e+0,1.221939e+1,1.977226e+1,
    2.929890e+1,3.906811e+1,4.690664e+1,5.234861e+1,5.669474e+1,5.908286e+1,5.608983e+1,
    4.880825e+1,4.096475e+1,3.429691e+1,2.901943e+1,2.490513e+1,2.168321e+1,1.913006e+1,
    1.707654e+1,1.539807e+1,1.400363e+1,1.282689e+1,1.181939e+1,1.094565e+1,1.017948e+1,
    9.501359e+0,8.896578e+0,8.353847e+0,7.864355e+0,7.421093e+0,7.018379e+0,6.651519e+0,
    6.316575e+0,6.010193e+0,5.729483e+0,5.471930e+0,5.235327e+0,5.017723e+0,4.817386e+0,
    4.632774e+0,4.462505e+0,4.305337e+0,4.160155e+0,4.025952e+0,3.901818e+0,3.786931e+0,
    3.680544e+0,3.581977e+0,3.490615e+0,3.405895e+0,3.327304e+0,3.254375e+0,3.186678e+0,
    3.123820e+0,3.065443e+0,3.011214e+0,2.960830e+0,2.914011e+0,2.870498e+0,2.830053e+0,
    2.792457e+0,2.757506e+0,2.725010e+0,2.694797e+0,2.666703e+0,2.640580e+0,2.616286e+0,
    2.593694e+0,2.572682e+0,2.553138e+0,2.534958e+0,2.518045e+0,2.502310e+0,2.487667e+0,
    2.474040e+0,2.461354e+0,2.449544e+0,2.438546e+0,2.428302e+0,2.418759e+0,2.409866e+0,
    2.401579e+0,2.393854e+0,2.386654e+0,2.379943e+0,2.373689e+0,2.367864e+0,2.362442e+0,
    2.357400e+0,2.352718e+0,2.348380e+0,2.344371e+0,2.340680e+0,2.337299e+0,2.334221e+0};
  static const G4double SL17[nL]={
    2.137065e-3,5.442007e-3,1.210645e-2,2.774945e-2,6.888202e-2,1.771769e-1,4.450546e-1,
    1.057471e+0,2.354951e+0,4.918482e+0,9.652965e+0,1.776486e+1,3.037627e+1,4.763569e+1,
    6.860085e+1,9.419933e+1,1.267075e+2,1.511993e+2,1.442063e+2,1.180149e+2,9.193020e+1,
    7.155618e+1,5.654375e+1,4.555533e+1,3.744235e+1,3.137117e+1,2.675825e+1,2.319641e+1,
    2.039992e+1,1.816677e+1,1.635312e+1,1.485589e+1,1.360078e+1,1.253381e+1,1.161550e+1,
    1.081664e+1,1.011536e+1,9.495083e+0,8.943002e+0,8.449085e+0,8.005329e+0,7.605245e+0,
    7.243487e+0,6.915590e+0,6.617778e+0,6.346818e+0,6.099922e+0,5.874662e+0,5.668909e+0,
    5.480788e+0,5.308639e+0,5.150985e+0,5.006505e+0,4.874020e+0,4.752467e+0,4.640891e+0,
    4.538429e+0,4.444301e+0,4.357799e+0,4.278281e+0,4.205160e+0,4.137906e+0,4.076030e+0,
    4.019090e+0,3.966677e+0,3.918419e+0,3.873973e+0,3.833027e+0,3.795290e+0,3.760498e+0,
    3.728406e+0,3.698788e+0,3.671436e+0,3.646160e+0,3.622782e+0,3.601139e+0,3.581082e+0,
    3.562470e+0,3.545177e+0,3.529084e+0,3.514082e+0,3.500071e+0,3.486960e+0,3.474663e+0,
    3.463103e+0,3.452209e+0,3.441916e+0,3.432164e+0,3.422901e+0,3.414078e+0,3.405652e+0,
    3.397582e+0,3.389836e+0,3.382382e+0,3.375196e+0,3.368254e+0,3.361538e+0,3.355034e+0,
    3.348730e+0,3.342620e+0,3.336699e+0,3.330967e+0,3.325427e+0,3.320085e+0,3.314951e+0};
  static const G4double SL18[nL]={
    2.220534e-3,5.640053e-3,1.253572e-2,2.881392e-2,7.191580e-2,1.859408e-1,4.687157e-1,
    1.115760e+0,2.485562e+0,5.183559e+0,1.013008e+1,1.847496e+1,3.103145e+1,4.701870e+1,
    6.345164e+1,7.777111e+1,8.950804e+1,9.321427e+1,8.410731e+1,6.975786e+1,5.670984e+1,
    4.641759e+1,3.856198e+1,3.257293e+1,2.796698e+1,2.438084e+1,2.154901e+1,1.927832e+1,
    1.742802e+1,1.589540e+1,1.460538e+1,1.350313e+1,1.254846e+1,1.171188e+1,1.097157e+1,
    1.031123e+1,9.718498e+0,9.183826e+0,8.699693e+0,8.260038e+0,7.859873e+0,7.495011e+0,
    7.161876e+0,6.857372e+0,6.578785e+0,6.323715e+0,6.090025e+0,5.875801e+0,5.679326e+0,
    5.499048e+0,5.333567e+0,5.181614e+0,5.042039e+0,4.913795e+0,4.795932e+0,4.687583e+0,
    4.587960e+0,4.496341e+0,4.412068e+0,4.334539e+0,4.263201e+0,4.197551e+0,4.137125e+0,
    4.081496e+0,4.030275e+0,3.983101e+0,3.939643e+0,3.899598e+0,3.862684e+0,3.828641e+0,
    3.797233e+0,3.768237e+0,3.741451e+0,3.716686e+0,3.693770e+0,3.672542e+0,3.652854e+0,
    3.634571e+0,3.617565e+0,3.601721e+0,3.586931e+0,3.573099e+0,3.560132e+0,3.547947e+0,
    3.536470e+0,3.525629e+0,3.515361e+0,3.505610e+0,3.496321e+0,3.487449e+0,3.478950e+0,
    3.470787e+0,3.462928e+0,3.455342e+0,3.448006e+0,3.440898e+0,3.434002e+0,3.427303e+0,
    3.420792e+0,3.414463e+0,3.408314e+0,3.402345e+0,3.396560e+0,3.390968e+0,3.385579e+0};
  static const G4double SL19[nL]={
    2.305897e-3,5.842654e-3,1.297593e-2,2.991119e-2,7.506153e-2,1.950960e-1,4.938019e-1,
    1.179632e+0,2.638978e+0,5.539887e+0,1.095013e+1,2.037657e+1,3.550284e+1,5.759776e+1,
    8.715375e+1,1.188643e+2,1.303680e+2,1.150932e+2,9.265204e+1,7.339629e+1,5.867930e+1,
    4.768315e+1,3.944920e+1,3.322398e+1,2.845918e+1,2.476205e+1,2.185087e+1,1.952250e+1,
    1.762968e+1,1.606532e+1,1.475137e+1,1.363085e+1,1.266212e+1,1.181463e+1,1.106578e+1,
    1.039872e+1,9.800670e+0,9.261781e+0,8.774299e+0,8.331989e+0,7.929726e+0,7.563221e+0,
    7.228815e+0,6.923346e+0,6.644047e+0,6.388477e+0,6.154464e+0,5.940066e+0,5.743540e+0,
    5.563316e+0,5.397976e+0,5.246237e+0,5.106936e+0,4.979015e+0,4.861515e+0,4.753562e+0,
    4.654358e+0,4.563177e+0,4.479354e+0,4.402282e+0,4.331405e+0,4.266213e+0,4.206241e+0,
    4.151059e+0,4.100274e+0,4.053523e+0,4.010474e+0,3.970819e+0,3.934277e+0,3.900587e+0,
    3.869510e+0,3.840823e+0,3.814323e+0,3.789821e+0,3.767142e+0,3.746128e+0,3.726629e+0,
    3.708509e+0,3.691642e+0,3.675912e+0,3.661213e+0,3.647445e+0,3.634520e+0,3.622354e+0,
    3.610872e+0,3.600004e+0,3.589688e+0,3.579867e+0,3.570489e+0,3.561507e+0,3.552881e+0,
    3.544574e+0,3.536552e+0,3.528789e+0,3.521261e+0,3.513947e+0,3.506832e+0,3.499903e+0,
    3.493151e+0,3.486572e+0,3.480165e+0,3.473930e+0,3.467875e+0,3.462007e+0,3.456340e+0};
  static const G4double SL20[nL]={
    2.545914e-3,6.412659e-3,1.422001e-2,3.303967e-2,8.409149e-2,2.213646e-1,5.649122e-1,
    1.354715e+0,3.029540e+0,6.323258e+0,1.232016e+1,2.225805e+1,3.662567e+1,5.344971e+1,
    6.796031e+1,7.669870e+1,8.176394e+1,8.725461e+1,8.966246e+1,8.202204e+1,6.857177e+1,
    5.574957e+1,4.551980e+1,3.772808e+1,3.182031e+1,2.730095e+1,2.379636e+1,2.103601e+1,
    1.882555e+1,1.702523e+1,1.553424e+1,1.427965e+1,1.320850e+1,1.228212e+1,1.147206e+1,
    1.075719e+1,1.012157e+1,9.553019e+0,9.042024e+0,8.581021e+0,8.163877e+0,7.785521e+0,
    7.441691e+0,7.128756e+0,6.843578e+0,6.583419e+0,6.345874e+0,6.128812e+0,5.930337e+0,
    5.748752e+0,5.582537e+0,5.430325e+0,5.290879e+0,5.163086e+0,5.045934e+0,4.938508e+0,
    4.839973e+0,4.749573e+0,4.666616e+0,4.590473e+0,4.520568e+0,4.456376e+0,4.397414e+0,
    4.343241e+0,4.293452e+0,4.247676e+0,4.205572e+0,4.166825e+0,4.131148e+0,4.098275e+0,
    4.067961e+0,4.039982e+0,4.014130e+0,3.990215e+0,3.968060e+0,3.947505e+0,3.928399e+0,
    3.910607e+0,3.894002e+0,3.878468e+0,3.863899e+0,3.850198e+0,3.837275e+0,3.825050e+0,
    3.813448e+0,3.802402e+0,3.791850e+0,3.781739e+0,3.772017e+0,3.762641e+0,3.753573e+0,
    3.744778e+0,3.736225e+0,3.727891e+0,3.719753e+0,3.711795e+0,3.704003e+0,3.696369e+0,
    3.688887e+0,3.681555e+0,3.674374e+0,3.667351e+0,3.660494e+0,3.653816e+0,3.647333e+0};
  static const G4double SL21[nL]={
    2.564250e-3,6.456227e-3,1.431544e-2,3.328130e-2,8.479343e-2,2.234161e-1,5.704940e-1,
    1.368585e+0,3.061075e+0,6.389076e+0,1.244551e+1,2.247024e+1,3.692501e+1,5.375094e+1,
    6.804530e+1,7.622389e+1,8.041323e+1,8.592812e+1,9.222352e+1,8.976000e+1,7.746020e+1,
    6.319441e+1,5.115788e+1,4.188960e+1,3.488950e+1,2.958086e+1,2.550619e+1,2.233120e+1,
    1.981637e+1,1.779045e+1,1.613068e+1,1.474867e+1,1.358052e+1,1.257968e+1,1.171205e+1,
    1.095233e+1,1.028156e+1,9.685270e+0,9.152260e+0,8.673689e+0,8.242451e+0,7.852734e+0,
    7.499706e+0,7.179292e+0,6.888008e+0,6.622849e+0,6.381195e+0,6.160747e+0,5.959474e+0,
    5.775573e+0,5.607438e+0,5.453629e+0,5.312856e+0,5.183957e+0,5.065885e+0,4.957691e+0,
    4.858517e+0,4.767585e+0,4.684186e+0,4.607676e+0,4.537466e+0,4.473021e+0,4.413851e+0,
    4.359506e+0,4.309575e+0,4.263683e+0,4.221483e+0,4.182658e+0,4.146917e+0,4.113990e+0,
    4.083632e+0,4.055616e+0,4.029732e+0,4.005789e+0,3.983609e+0,3.963030e+0,3.943902e+0,
    3.926087e+0,3.909459e+0,3.893901e+0,3.879307e+0,3.865579e+0,3.852627e+0,3.840370e+0,
    3.828735e+0,3.817652e+0,3.807062e+0,3.796909e+0,3.787144e+0,3.777723e+0,3.768606e+0,
    3.759759e+0,3.751154e+0,3.742764e+0,3.734569e+0,3.726552e+0,3.718700e+0,3.711004e+0,
    3.703458e+0,3.696062e+0,3.688816e+0,3.681726e+0,3.674803e+0,3.668058e+0,3.661508e+0};
  static const G4double SL22[nL]={
    3.007427e-3,7.510236e-3,1.663782e-2,3.922952e-2,1.022486e-1,2.747505e-1,7.108351e-1,
    1.719543e+0,3.868385e+0,8.112125e+0,1.586986e+1,2.877567e+1,4.762099e+1,7.091149e+1,
    9.555037e+1,1.153230e+2,1.159518e+2,9.952423e+1,8.059227e+1,6.490265e+1,5.285116e+1,
    4.371769e+1,3.676918e+1,3.143378e+1,2.728996e+1,2.403042e+1,2.143067e+1,1.932609e+1,
    1.759555e+1,1.614984e+1,1.492332e+1,1.386771e+1,1.294753e+1,1.213663e+1,1.141562e+1,
    1.076997e+1,1.018862e+1,9.662970e+0,9.186212e+0,8.752798e+0,8.358120e+0,7.998266e+0,
    7.669857e+0,7.369937e+0,7.095886e+0,6.845367e+0,6.616282e+0,6.406739e+0,6.215024e+0,
    6.039585e+0,5.879011e+0,5.732019e+0,5.597439e+0,5.474206e+0,5.361347e+0,5.257976e+0,
    5.163281e+0,5.076520e+0,4.997014e+0,4.924143e+0,4.857338e+0,4.796075e+0,4.739877e+0,
    4.688303e+0,4.640952e+0,4.597451e+0,4.557461e+0,4.520669e+0,4.486786e+0,4.455548e+0,
    4.426712e+0,4.400054e+0,4.375369e+0,4.352466e+0,4.331174e+0,4.311331e+0,4.292793e+0,
    4.275425e+0,4.259104e+0,4.243719e+0,4.229166e+0,4.215352e+0,4.202192e+0,4.189610e+0,
    4.177536e+0,4.165907e+0,4.154667e+0,4.143766e+0,4.133160e+0,4.122809e+0,4.112681e+0,
    4.102746e+0,4.092980e+0,4.083363e+0,4.073880e+0,4.064520e+0,4.055276e+0,4.046143e+0,
    4.037123e+0,4.028220e+0,4.019441e+0,4.010799e+0,4.002309e+0,3.993991e+0,3.985867e+0};
  static const G4double SL23[nL]={
    3.202591e-3,7.975022e-3,1.767003e-2,4.191214e-2,1.102179e-1,2.983388e-1,7.754839e-1,
    1.881188e+0,4.239060e+0,8.896849e+0,1.740204e+1,3.149771e+1,5.193761e+1,7.708634e+1,
    1.037149e+2,1.232794e+2,1.201643e+2,1.011825e+2,8.136267e+1,6.536884e+1,5.318799e+1,
    4.398670e+1,3.699713e+1,3.163451e+1,2.747158e+1,2.419807e+1,2.158786e+1,1.947527e+1,
    1.773851e+1,1.628791e+1,1.505750e+1,1.399877e+1,1.307609e+1,1.226319e+1,1.154059e+1,
    1.089372e+1,1.031145e+1,9.785155e+0,9.307994e+0,8.874394e+0,8.479723e+0,8.120049e+0,
    7.791975e+0,7.492525e+0,7.219066e+0,6.969243e+0,6.740944e+0,6.532261e+0,6.341471e+0,
    6.167007e+0,6.007447e+0,5.861497e+0,5.727978e+0,5.605816e+0,5.494030e+0,5.391724e+0,
    5.298080e+0,5.212351e+0,5.133851e+0,5.061955e+0,4.996087e+0,4.935723e+0,4.880379e+0,
    4.829614e+0,4.783020e+0,4.740226e+0,4.700887e+0,4.664691e+0,4.631348e+0,4.600593e+0,
    4.572181e+0,4.545889e+0,4.521511e+0,4.498859e+0,4.477759e+0,4.458053e+0,4.439594e+0,
    4.422252e+0,4.405903e+0,4.390437e+0,4.375754e+0,4.361761e+0,4.348375e+0,4.335521e+0,
    4.323132e+0,4.311147e+0,4.299511e+0,4.288177e+0,4.277103e+0,4.266251e+0,4.255591e+0,
    4.245096e+0,4.234743e+0,4.224516e+0,4.214402e+0,4.204390e+0,4.194477e+0,4.184662e+0,
    4.174947e+0,4.165339e+0,4.155849e+0,4.146491e+0,4.137283e+0,4.128247e+0,4.119408e+0};
  static const G4double SL24[nL]={
    4.424391e-3,1.089365e-2,2.425454e-2,5.950323e-2,1.636625e-1,4.585667e-1,1.218434e+0,
    3.000059e+0,6.849073e+0,1.459387e+1,2.913751e+1,5.434530e+1,9.302757e+1,1.344961e+2,
    1.455690e+2,1.285470e+2,1.059779e+2,8.603243e+1,7.002571e+1,5.756765e+1,4.795028e+1,
    4.051569e+1,3.473537e+1,3.020519e+1,2.662094e+1,2.375420e+1,2.143316e+1,1.952845e+1,
    1.794276e+1,1.660309e+1,1.545498e+1,1.445799e+1,1.358216e+1,1.280530e+1,1.211087e+1,
    1.148642e+1,1.092239e+1,1.041132e+1,9.947182e+0,9.525039e+0,9.140708e+0,8.790584e+0,
    8.471504e+0,8.180659e+0,7.915525e+0,7.673825e+0,7.453491e+0,7.252639e+0,7.069553e+0,
    6.902666e+0,6.750545e+0,6.611879e+0,6.485472e+0,6.370227e+0,6.265143e+0,6.169303e+0,
    6.081868e+0,6.002071e+0,5.929208e+0,5.862637e+0,5.801770e+0,5.746067e+0,5.695037e+0,
    5.648227e+0,5.605225e+0,5.565653e+0,5.529167e+0,5.495449e+0,5.464212e+0,5.435193e+0,
    5.408150e+0,5.382864e+0,5.359136e+0,5.336785e+0,5.315644e+0,5.295565e+0,5.276412e+0,
    5.258062e+0,5.240405e+0,5.223341e+0,5.206783e+0,5.190650e+0,5.174872e+0,5.159387e+0,
    5.144141e+0,5.129087e+0,5.114183e+0,5.099396e+0,5.084697e+0,5.070063e+0,5.055476e+0,
    5.040923e+0,5.026397e+0,5.011893e+0,4.997412e+0,4.982959e+0,4.968544e+0,4.954178e+0,
    4.939879e+0,4.925668e+0,4.911569e+0,4.897612e+0,4.883829e+0,4.870257e+0,4.856937e+0};
  static const G4double SL25[nL]={
    5.218262e-3,1.279812e-2,2.863691e-2,7.159626e-2,2.012920e-1,5.725263e-1,1.533263e+0,
    3.785207e+0,8.620686e+0,1.819011e+1,3.550779e+1,6.346712e+1,1.028732e+2,1.498888e+2,
    1.795448e+2,1.653851e+2,1.332773e+2,1.044804e+2,8.246369e+1,6.611941e+1,5.397243e+1,
    4.486079e+1,3.794689e+1,3.263665e+1,2.850697e+1,2.525335e+1,2.265438e+1,2.054771e+1,
    1.881372e+1,1.736415e+1,1.613394e+1,1.507520e+1,1.415268e+1,1.334041e+1,1.261911e+1,
    1.197431e+1,1.139498e+1,1.087253e+1,1.040012e+1,9.972149e+0,9.583943e+0,9.231507e+0,
    8.911366e+0,8.620462e+0,8.356070e+0,8.115746e+0,7.897283e+0,7.698685e+0,7.518138e+0,
    7.353992e+0,7.204745e+0,7.069027e+0,6.945589e+0,6.833293e+0,6.731099e+0,6.638057e+0,
    6.553303e+0,6.476044e+0,6.405559e+0,6.341189e+0,6.282333e+0,6.228441e+0,6.179012e+0,
    6.133589e+0,6.091756e+0,6.053132e+0,6.017372e+0,5.984159e+0,5.953209e+0,5.924259e+0,
    5.897073e+0,5.871438e+0,5.847158e+0,5.824059e+0,5.801982e+0,5.780784e+0,5.760337e+0,
    5.740528e+0,5.721253e+0,5.702421e+0,5.683953e+0,5.665777e+0,5.647831e+0,5.630063e+0,
    5.612426e+0,5.594882e+0,5.577398e+0,5.559949e+0,5.542515e+0,5.525080e+0,5.507636e+0,
    5.490177e+0,5.472703e+0,5.455220e+0,5.437735e+0,5.420261e+0,5.402815e+0,5.385418e+0,
    5.368096e+0,5.350876e+0,5.333791e+0,5.316880e+0,5.300182e+0,5.283744e+0,5.267615e+0};
  static const G4double SL26[nL]={
    9.533418e-3,2.324917e-2,5.364098e-2,1.447139e-1,4.381268e-1,1.303754e+0,3.571583e+0,
    8.890991e+0,2.014541e+1,4.138069e+1,7.546159e+1,1.179996e+2,1.568622e+2,1.907924e+2,
    2.305942e+2,2.457159e+2,2.095925e+2,1.607399e+2,1.215717e+2,9.340695e+1,7.340074e+1,
    5.903474e+1,4.855591e+1,4.078885e+1,3.494131e+1,3.047098e+1,2.700048e+1,2.426339e+1,
    2.206960e+1,2.028240e+1,1.880303e+1,1.755986e+1,1.650082e+1,1.558780e+1,1.479275e+1,
    1.409475e+1,1.347801e+1,1.293034e+1,1.244217e+1,1.200582e+1,1.161497e+1,1.126435e+1,
    1.094945e+1,1.066638e+1,1.041174e+1,1.018252e+1,9.976076e+0,9.790007e+0,9.622186e+0,
    9.470694e+0,9.333803e+0,9.209956e+0,9.097750e+0,8.995915e+0,8.903308e+0,8.818895e+0,
    8.741745e+0,8.671018e+0,8.605956e+0,8.545878e+0,8.490170e+0,8.438281e+0,8.389716e+0,
    8.344033e+0,8.300834e+0,8.259766e+0,8.220512e+0,8.182794e+0,8.146361e+0,8.110993e+0,
    8.076498e+0,8.042705e+0,8.009465e+0,7.976650e+0,7.944147e+0,7.911862e+0,7.879713e+0,
    7.847632e+0,7.815563e+0,7.783460e+0,7.751287e+0,7.719017e+0,7.686631e+0,7.654117e+0,
    7.621469e+0,7.588688e+0,7.555781e+0,7.522758e+0,7.489637e+0,7.456437e+0,7.423184e+0,
    7.389908e+0,7.356640e+0,7.323419e+0,7.290286e+0,7.257286e+0,7.224468e+0,7.191884e+0,
    7.159594e+0,7.127657e+0,7.096141e+0,7.065118e+0,7.034663e+0,7.004858e+0,6.975791e+0};
  static const G4double SL27[nL]={
    1.043535e-2,2.545247e-2,5.908485e-2,1.613411e-1,4.935021e-1,1.477078e+0,4.058976e+0,
    1.012392e+1,2.297013e+1,4.720890e+1,8.609450e+1,1.355874e+2,1.876865e+2,2.480184e+2,
    2.968858e+2,2.729073e+2,2.110699e+2,1.571879e+2,1.183563e+2,9.114148e+1,7.189981e+1,
    5.807038e+1,4.795957e+1,4.044610e+1,3.477555e+1,3.043070e+1,2.705064e+1,2.437987e+1,
    2.223556e+1,2.048595e+1,1.903566e+1,1.781544e+1,1.677484e+1,1.587694e+1,1.509453e+1,
    1.440731e+1,1.379991e+1,1.326049e+1,1.277969e+1,1.235001e+1,1.196525e+1,1.162022e+1,
    1.131048e+1,1.103217e+1,1.078194e+1,1.055679e+1,1.035408e+1,1.017143e+1,1.000671e+1,
    9.858006e+0,9.723600e+0,9.601935e+0,9.491608e+0,9.391355e+0,9.300037e+0,9.216625e+0,
    9.140195e+0,9.069912e+0,9.005027e+0,8.944866e+0,8.888826e+0,8.836365e+0,8.786998e+0,
    8.740293e+0,8.695865e+0,8.653370e+0,8.612505e+0,8.573001e+0,8.534620e+0,8.497153e+0,
    8.460418e+0,8.424256e+0,8.388530e+0,8.353120e+0,8.317926e+0,8.282862e+0,8.247857e+0,
    8.212851e+0,8.177798e+0,8.142660e+0,8.107410e+0,8.072027e+0,8.036502e+0,8.000828e+0,
    7.965007e+0,7.929047e+0,7.892960e+0,7.856763e+0,7.820479e+0,7.784135e+0,7.747759e+0,
    7.711389e+0,7.675061e+0,7.638818e+0,7.602707e+0,7.566777e+0,7.531083e+0,7.495683e+0,
    7.460641e+0,7.426023e+0,7.391901e+0,7.358353e+0,7.325461e+0,7.293313e+0,7.262004e+0};
  static const G4double SL28[nL]={
    1.177612e-2,2.873855e-2,6.729897e-2,1.868027e-1,5.790453e-1,1.745770e+0,4.814971e+0,
    1.203191e+1,2.730811e+1,5.598968e+1,1.014774e+2,1.592489e+2,2.225285e+2,2.776163e+2,
    2.688536e+2,2.132293e+2,1.608280e+2,1.222362e+2,9.485335e+1,7.529219e+1,6.110006e+1,
    5.063436e+1,4.279795e+1,3.684661e+1,3.226548e+1,2.869156e+1,2.586467e+1,2.359614e+1,
    2.174808e+1,2.021937e+1,1.893580e+1,1.784291e+1,1.690071e+1,1.607977e+1,1.535827e+1,
    1.471984e+1,1.415200e+1,1.364501e+1,1.319113e+1,1.278403e+1,1.241840e+1,1.208972e+1,
    1.179407e+1,1.152799e+1,1.128840e+1,1.107255e+1,1.087797e+1,1.070243e+1,1.054390e+1,
    1.040057e+1,1.027077e+1,1.015302e+1,1.004596e+1,9.948371e+0,9.859147e+0,9.777296e+0,
    9.701922e+0,9.632220e+0,9.567466e+0,9.507013e+0,9.450281e+0,9.396751e+0,9.345963e+0,
    9.297506e+0,9.251017e+0,9.206173e+0,9.162691e+0,9.120322e+0,9.078848e+0,9.038079e+0,
    8.997851e+0,8.958023e+0,8.918474e+0,8.879101e+0,8.839819e+0,8.800557e+0,8.761258e+0,
    8.721877e+0,8.682378e+0,8.642738e+0,8.602939e+0,8.562975e+0,8.522843e+0,8.482548e+0,
    8.442103e+0,8.401522e+0,8.360828e+0,8.320046e+0,8.279205e+0,8.238341e+0,8.197491e+0,
    8.156696e+0,8.116003e+0,8.075461e+0,8.035123e+0,7.995046e+0,7.955292e+0,7.915926e+0,
    7.877018e+0,7.838644e+0,7.800883e+0,7.763820e+0,7.727547e+0,7.692162e+0,7.657767e+0};
  static const G4double SL29[nL]={
    1.365967e-2,3.337537e-2,7.906748e-2,2.239693e-1,7.052155e-1,2.143189e+0,5.929377e+0,
    1.480860e+1,3.341772e+1,6.741315e+1,1.178596e+2,1.735378e+2,2.243020e+2,2.786430e+2,
    3.008511e+2,2.541744e+2,1.921484e+2,1.437224e+2,1.094876e+2,8.542953e+1,6.828273e+1,
    5.584793e+1,4.667796e+1,3.980970e+1,3.458954e+1,3.056493e+1,2.741696e+1,2.491771e+1,
    2.290273e+1,2.125263e+1,1.988051e+1,1.872296e+1,1.773367e+1,1.687863e+1,1.613270e+1,
    1.547709e+1,1.489750e+1,1.438287e+1,1.392442e+1,1.351505e+1,1.314884e+1,1.282080e+1,
    1.252665e+1,1.226263e+1,1.202543e+1,1.181213e+1,1.162010e+1,1.144700e+1,1.129071e+1,
    1.114933e+1,1.102117e+1,1.090467e+1,1.079846e+1,1.070131e+1,1.061208e+1,1.052979e+1,
    1.045353e+1,1.038251e+1,1.031602e+1,1.025342e+1,1.019413e+1,1.013767e+1,1.008359e+1,
    1.003150e+1,9.981059e+0,9.931961e+0,9.883946e+0,9.836788e+0,9.790289e+0,9.744282e+0,
    9.698622e+0,9.653189e+0,9.607882e+0,9.562616e+0,9.517324e+0,9.471954e+0,9.426464e+0,
    9.380825e+0,9.335018e+0,9.289031e+0,9.242863e+0,9.196518e+0,9.150008e+0,9.103348e+0,
    9.056561e+0,9.009673e+0,8.962717e+0,8.915727e+0,8.868742e+0,8.821806e+0,8.774965e+0,
    8.728269e+0,8.681773e+0,8.635535e+0,8.589616e+0,8.544082e+0,8.499003e+0,8.454454e+0,
    8.410514e+0,8.367269e+0,8.324807e+0,8.283225e+0,8.242627e+0,8.203121e+0,8.164824e+0};
  static const G4double SL30[nL]={
    2.103117e-2,5.172350e-2,1.273601e-1,3.830439e-1,1.257967e+0,3.901184e+0,1.087470e+1,
    2.709930e+1,6.005843e+1,1.151063e+2,1.807622e+2,2.248182e+2,2.337093e+2,2.305565e+2,
    2.282780e+2,2.025860e+2,1.607669e+2,1.243617e+2,9.741054e+1,7.794111e+1,6.376177e+1,
    5.328718e+1,4.543843e+1,3.947900e+1,3.489731e+1,3.133117e+1,2.851997e+1,2.627388e+1,
    2.445368e+1,2.295694e+1,2.170826e+1,2.065218e+1,1.974793e+1,1.896539e+1,1.828223e+1,
    1.768163e+1,1.715075e+1,1.667954e+1,1.625998e+1,1.588548e+1,1.555055e+1,1.525048e+1,
    1.498116e+1,1.473901e+1,1.452086e+1,1.432386e+1,1.414548e+1,1.398346e+1,1.383576e+1,
    1.370056e+1,1.357624e+1,1.346135e+1,1.335459e+1,1.325483e+1,1.316103e+1,1.307231e+1,
    1.298788e+1,1.290703e+1,1.282916e+1,1.275374e+1,1.268031e+1,1.260848e+1,1.253790e+1,
    1.246828e+1,1.239939e+1,1.233100e+1,1.226295e+1,1.219511e+1,1.212735e+1,1.205959e+1,
    1.199175e+1,1.192379e+1,1.185567e+1,1.178738e+1,1.171889e+1,1.165022e+1,1.158138e+1,
    1.151237e+1,1.144324e+1,1.137400e+1,1.130470e+1,1.123538e+1,1.116608e+1,1.109684e+1,
    1.102773e+1,1.095880e+1,1.089011e+1,1.082171e+1,1.075367e+1,1.068607e+1,1.061897e+1,
    1.055244e+1,1.048657e+1,1.042144e+1,1.035714e+1,1.029376e+1,1.023139e+1,1.017015e+1,
    1.011014e+1,1.005149e+1,9.994313e+0,9.938752e+0,9.884949e+0,9.833059e+0,9.783245e+0};
  static const G4double SL31[nL]={
    2.164664e-2,5.326850e-2,1.315360e-1,3.972007e-1,1.307919e+0,4.061184e+0,1.132688e+1,
    2.822812e+1,6.251535e+1,1.195506e+2,1.871899e+2,2.337421e+2,2.500245e+2,2.569217e+2,
    2.417968e+2,1.964670e+2,1.512314e+2,1.169992e+2,9.234431e+1,7.451650e+1,6.143973e+1,
    5.170376e+1,4.435567e+1,3.874082e+1,3.440004e+1,3.100496e+1,2.831703e+1,2.616109e+1,
    2.440774e+1,2.296126e+1,2.175087e+1,2.072434e+1,1.984318e+1,1.907891e+1,1.841036e+1,
    1.782157e+1,1.730033e+1,1.683705e+1,1.642406e+1,1.605504e+1,1.572468e+1,1.542842e+1,
    1.516228e+1,1.492277e+1,1.470678e+1,1.451152e+1,1.433452e+1,1.417353e+1,1.402657e+1,
    1.389184e+1,1.376774e+1,1.365284e+1,1.354587e+1,1.344569e+1,1.335132e+1,1.326186e+1,
    1.317653e+1,1.309467e+1,1.301566e+1,1.293899e+1,1.286422e+1,1.279096e+1,1.271887e+1,
    1.264768e+1,1.257715e+1,1.250708e+1,1.243731e+1,1.236770e+1,1.229815e+1,1.222857e+1,
    1.215890e+1,1.208909e+1,1.201911e+1,1.194895e+1,1.187861e+1,1.180808e+1,1.173739e+1,
    1.166655e+1,1.159559e+1,1.152454e+1,1.145345e+1,1.138236e+1,1.131131e+1,1.124036e+1,
    1.116955e+1,1.109896e+1,1.102863e+1,1.095863e+1,1.088903e+1,1.081989e+1,1.075129e+1,
    1.068331e+1,1.061602e+1,1.054952e+1,1.048389e+1,1.041922e+1,1.035562e+1,1.029319e+1,
    1.023205e+1,1.017232e+1,1.011413e+1,1.005761e+1,1.000292e+1,9.950209e+0,9.899647e+0};
  static const G4double SL32[nL]={
    2.258863e-2,5.563670e-2,1.379665e-1,4.191065e-1,1.385411e+0,4.309706e+0,1.202998e+1,
    2.998579e+1,6.634936e+1,1.265293e+2,1.975930e+2,2.495457e+2,2.807339e+2,3.044005e+2,
    2.792187e+2,2.183968e+2,1.646189e+2,1.256518e+2,9.814320e+1,7.851348e+1,6.426612e+1,
    5.375330e+1,4.588048e+1,3.990576e+1,3.531490e+1,3.174407e+1,2.893143e+1,2.668630e+1,
    2.486877e+1,2.337590e+1,2.213188e+1,2.108096e+1,2.018214e+1,1.940513e+1,1.872746e+1,
    1.813222e+1,1.760647e+1,1.714012e+1,1.672508e+1,1.635475e+1,1.602359e+1,1.572684e+1,
    1.546043e+1,1.522073e+1,1.500457e+1,1.480910e+1,1.463180e+1,1.447040e+1,1.432289e+1,
    1.418745e+1,1.406248e+1,1.394655e+1,1.383837e+1,1.373683e+1,1.364092e+1,1.354977e+1,
    1.346262e+1,1.337878e+1,1.329768e+1,1.321880e+1,1.314171e+1,1.306603e+1,1.299144e+1,
    1.291767e+1,1.284450e+1,1.277173e+1,1.269921e+1,1.262681e+1,1.255444e+1,1.248202e+1,
    1.240949e+1,1.233681e+1,1.226397e+1,1.219095e+1,1.211774e+1,1.204437e+1,1.197085e+1,
    1.189720e+1,1.182346e+1,1.174967e+1,1.167586e+1,1.160209e+1,1.152840e+1,1.145484e+1,
    1.138148e+1,1.130838e+1,1.123558e+1,1.116317e+1,1.109121e+1,1.101977e+1,1.094893e+1,
    1.087876e+1,1.080935e+1,1.074079e+1,1.067317e+1,1.060658e+1,1.054114e+1,1.047695e+1,
    1.041413e+1,1.035280e+1,1.029310e+1,1.023517e+1,1.017916e+1,1.012525e+1,1.007359e+1};
  static const G4double SL33[nL]={
    2.454062e-2,6.055745e-2,1.514382e-1,4.653812e-1,1.549692e+0,4.836226e+0,1.351062e+1,
    3.362441e+1,7.394409e+1,1.387856e+2,2.098544e+2,2.505891e+2,2.554933e+2,2.445847e+2,
    2.124190e+2,1.687167e+2,1.317569e+2,1.042289e+2,8.404509e+1,6.911577e+1,5.793269e+1,
    4.945432e+1,4.295680e+1,3.792803e+1,3.399868e+1,3.089784e+1,2.842428e+1,2.642742e+1,
    2.479430e+1,2.344031e+1,2.230232e+1,2.133341e+1,2.049880e+1,1.977269e+1,1.913579e+1,
    1.857355e+1,1.807474e+1,1.763054e+1,1.723383e+1,1.687870e+1,1.656015e+1,1.627387e+1,
    1.601607e+1,1.578341e+1,1.557289e+1,1.538185e+1,1.520789e+1,1.504885e+1,1.490281e+1,
    1.476805e+1,1.464303e+1,1.452637e+1,1.441687e+1,1.431346e+1,1.421517e+1,1.412119e+1,
    1.403079e+1,1.394334e+1,1.385828e+1,1.377515e+1,1.369354e+1,1.361311e+1,1.353357e+1,
    1.345466e+1,1.337620e+1,1.329802e+1,1.321997e+1,1.314197e+1,1.306393e+1,1.298578e+1,
    1.290750e+1,1.282905e+1,1.275042e+1,1.267163e+1,1.259268e+1,1.251358e+1,1.243438e+1,
    1.235509e+1,1.227578e+1,1.219646e+1,1.211721e+1,1.203806e+1,1.195908e+1,1.188033e+1,
    1.180186e+1,1.172374e+1,1.164603e+1,1.156881e+1,1.149216e+1,1.141614e+1,1.134083e+1,
    1.126633e+1,1.119272e+1,1.112009e+1,1.104854e+1,1.097817e+1,1.090910e+1,1.084144e+1,
    1.077532e+1,1.071087e+1,1.064824e+1,1.058758e+1,1.052905e+1,1.047282e+1,1.041909e+1};
  static const G4double SL34[nL]={
    2.555084e-2,6.311099e-2,1.584856e-1,4.897911e-1,1.636776e+0,5.116560e+0,1.430535e+1,
    3.561279e+1,7.827264e+1,1.466546e+2,2.221104e+2,2.722029e+2,3.010875e+2,3.106231e+2,
    2.658424e+2,2.035783e+2,1.541383e+2,1.188706e+2,9.384637e+1,7.582817e+1,6.263287e+1,
    5.281838e+1,4.541822e+1,3.977013e+1,3.540993e+1,3.200563e+1,2.931600e+1,2.716381e+1,
    2.541809e+1,2.398190e+1,2.278352e+1,2.177003e+1,2.090240e+1,2.015174e+1,1.949657e+1,
    1.892070e+1,1.841172e+1,1.795992e+1,1.755750e+1,1.719807e+1,1.687624e+1,1.658741e+1,
    1.632757e+1,1.609322e+1,1.588122e+1,1.568882e+1,1.551353e+1,1.535316e+1,1.520573e+1,
    1.506950e+1,1.494290e+1,1.482456e+1,1.471325e+1,1.460789e+1,1.450755e+1,1.441139e+1,
    1.431870e+1,1.422884e+1,1.414128e+1,1.405556e+1,1.397128e+1,1.388810e+1,1.380574e+1,
    1.372397e+1,1.364260e+1,1.356147e+1,1.348045e+1,1.339945e+1,1.331839e+1,1.323723e+1,
    1.315592e+1,1.307446e+1,1.299284e+1,1.291106e+1,1.282914e+1,1.274710e+1,1.266499e+1,
    1.258282e+1,1.250066e+1,1.241854e+1,1.233652e+1,1.225466e+1,1.217300e+1,1.209162e+1,
    1.201058e+1,1.192993e+1,1.184976e+1,1.177014e+1,1.169113e+1,1.161282e+1,1.153529e+1,
    1.145863e+1,1.138292e+1,1.130827e+1,1.123477e+1,1.116254e+1,1.109168e+1,1.102231e+1,
    1.095457e+1,1.088860e+1,1.082454e+1,1.076255e+1,1.070280e+1,1.064547e+1,1.059075e+1};
  static const G4double SL35[nL]={
    2.764032e-2,6.840683e-2,1.732180e-1,5.412140e-1,1.820749e+0,5.707551e+0,1.596471e+1,
    3.965730e+1,8.650314e+1,1.589997e+2,2.314210e+2,2.624498e+2,2.501423e+2,2.222984e+2,
    1.892896e+2,1.529510e+2,1.218148e+2,9.800381e+1,8.015642e+1,6.671593e+1,5.649992e+1,
    4.866270e+1,4.259935e+1,3.787106e+1,3.415439e+1,3.120761e+1,2.884831e+1,2.693819e+1,
    2.537244e+1,2.407192e+1,2.297720e+1,2.204394e+1,2.123915e+1,2.053827e+1,1.992295e+1,
    1.937927e+1,1.889649e+1,1.846614e+1,1.808135e+1,1.773642e+1,1.742651e+1,1.714742e+1,
    1.689546e+1,1.666738e+1,1.646026e+1,1.627150e+1,1.609875e+1,1.593994e+1,1.579319e+1,
    1.565683e+1,1.552939e+1,1.540955e+1,1.529615e+1,1.518817e+1,1.508473e+1,1.498504e+1,
    1.488843e+1,1.479431e+1,1.470218e+1,1.461162e+1,1.452226e+1,1.443381e+1,1.434600e+1,
    1.425863e+1,1.417154e+1,1.408459e+1,1.399767e+1,1.391072e+1,1.382367e+1,1.373650e+1,
    1.364917e+1,1.356169e+1,1.347407e+1,1.338633e+1,1.329848e+1,1.321057e+1,1.312263e+1,
    1.303472e+1,1.294688e+1,1.285916e+1,1.277162e+1,1.268433e+1,1.259735e+1,1.251074e+1,
    1.242457e+1,1.233891e+1,1.225383e+1,1.216942e+1,1.208575e+1,1.200290e+1,1.192096e+1,
    1.184003e+1,1.176019e+1,1.168155e+1,1.160422e+1,1.152831e+1,1.145394e+1,1.138123e+1,
    1.131034e+1,1.124140e+1,1.117458e+1,1.111003e+1,1.104795e+1,1.098853e+1,1.093197e+1};
  static const G4double SL36[nL]={
    2.908478e-2,7.207879e-2,1.835222e-1,5.774992e-1,1.951299e+0,6.129692e+0,1.716666e+1,
    4.268370e+1,9.316524e+1,1.715413e+2,2.537232e+2,3.118305e+2,3.698348e+2,4.142092e+2,
    3.599490e+2,2.689394e+2,1.969985e+2,1.471686e+2,1.129308e+2,8.901259e+1,7.196919e+1,
    5.959877e+1,5.047252e+1,4.364148e+1,3.846005e+1,3.447932e+1,3.138128e+1,2.893751e+1,
    2.698235e+1,2.539504e+1,2.408733e+1,2.299466e+1,2.206971e+1,2.127766e+1,2.059274e+1,
    1.999559e+1,1.947149e+1,1.900901e+1,1.859907e+1,1.823428e+1,1.790856e+1,1.761673e+1,
    1.735437e+1,1.711766e+1,1.690324e+1,1.670817e+1,1.652986e+1,1.636602e+1,1.621462e+1,
    1.607388e+1,1.594222e+1,1.581828e+1,1.570083e+1,1.558882e+1,1.548134e+1,1.537759e+1,
    1.527688e+1,1.517863e+1,1.508233e+1,1.498755e+1,1.489394e+1,1.480120e+1,1.470908e+1,
    1.461738e+1,1.452595e+1,1.443465e+1,1.434339e+1,1.425209e+1,1.416072e+1,1.406924e+1,
    1.397765e+1,1.388593e+1,1.379410e+1,1.370220e+1,1.361024e+1,1.351827e+1,1.342633e+1,
    1.333447e+1,1.324275e+1,1.315122e+1,1.305994e+1,1.296897e+1,1.287838e+1,1.278824e+1,
    1.269862e+1,1.260959e+1,1.252123e+1,1.243361e+1,1.234682e+1,1.226095e+1,1.217608e+1,
    1.209230e+1,1.200973e+1,1.192845e+1,1.184859e+1,1.177025e+1,1.169358e+1,1.161869e+1,
    1.154574e+1,1.147488e+1,1.140627e+1,1.134009e+1,1.127653e+1,1.121578e+1,1.115808e+1};
  static const G4double SL37[nL]={
    2.982256e-2,7.395762e-2,1.888215e-1,5.962481e-1,2.018837e+0,6.347438e+0,1.778020e+1,
    4.418698e+1,9.625560e+1,1.763759e+2,2.587902e+2,3.153492e+2,3.717654e+2,4.152507e+2,
    3.605764e+2,2.693717e+2,1.973299e+2,1.474418e+2,1.131676e+2,8.922544e+1,7.216593e+1,
    5.978476e+1,5.065157e+1,4.381637e+1,3.863284e+1,3.465151e+1,3.155398e+1,2.911148e+1,
    2.715812e+1,2.557296e+1,2.426760e+1,2.317735e+1,2.225483e+1,2.146515e+1,2.078248e+1,
    2.018745e+1,1.966529e+1,1.920456e+1,1.879617e+1,1.843272e+1,1.810811e+1,1.781717e+1,
    1.755549e+1,1.731924e+1,1.710506e+1,1.691002e+1,1.673154e+1,1.656733e+1,1.641538e+1,
    1.627392e+1,1.614137e+1,1.601639e+1,1.589776e+1,1.578445e+1,1.567554e+1,1.557025e+1,
    1.546791e+1,1.536794e+1,1.526984e+1,1.517321e+1,1.507768e+1,1.498297e+1,1.488885e+1,
    1.479512e+1,1.470162e+1,1.460824e+1,1.451488e+1,1.442149e+1,1.432802e+1,1.423444e+1,
    1.414075e+1,1.404695e+1,1.395306e+1,1.385910e+1,1.376512e+1,1.367114e+1,1.357723e+1,
    1.348342e+1,1.338977e+1,1.329635e+1,1.320321e+1,1.311042e+1,1.301805e+1,1.292616e+1,
    1.283483e+1,1.274413e+1,1.265413e+1,1.256493e+1,1.247660e+1,1.238923e+1,1.230290e+1,
    1.221773e+1,1.213380e+1,1.205122e+1,1.197011e+1,1.189058e+1,1.181277e+1,1.173682e+1,
    1.166286e+1,1.159105e+1,1.152158e+1,1.145460e+1,1.139033e+1,1.132895e+1,1.127071e+1};
  static const G4double SL38[nL]={
    3.019534e-2,7.490782e-2,1.915089e-1,6.057947e-1,2.053522e+0,6.462026e+0,1.812420e+1,
    4.516020e+1,9.894848e+1,1.839628e+2,2.814112e+2,3.765563e+2,4.501148e+2,4.042364e+2,
    3.029151e+2,2.212396e+2,1.647273e+2,1.259548e+2,9.887598e+1,7.956617e+1,6.553083e+1,
    5.515881e+1,4.738384e+1,4.148198e+1,3.694970e+1,3.342931e+1,3.066239e+1,2.845999e+1,
    2.668305e+1,2.522895e+1,2.402196e+1,2.300629e+1,2.214081e+1,2.139515e+1,2.074672e+1,
    2.017852e+1,1.967755e+1,1.923363e+1,1.883866e+1,1.848598e+1,1.817002e+1,1.788607e+1,
    1.763002e+1,1.739830e+1,1.718776e+1,1.699564e+1,1.681946e+1,1.665704e+1,1.650645e+1,
    1.636599e+1,1.623415e+1,1.610960e+1,1.599119e+1,1.587790e+1,1.576885e+1,1.566328e+1,
    1.556054e+1,1.546006e+1,1.536137e+1,1.526406e+1,1.516780e+1,1.507230e+1,1.497734e+1,
    1.488273e+1,1.478833e+1,1.469401e+1,1.459971e+1,1.450534e+1,1.441089e+1,1.431633e+1,
    1.422164e+1,1.412685e+1,1.403197e+1,1.393703e+1,1.384207e+1,1.374712e+1,1.365224e+1,
    1.355748e+1,1.346290e+1,1.336855e+1,1.327450e+1,1.318081e+1,1.308756e+1,1.299481e+1,
    1.290264e+1,1.281111e+1,1.272032e+1,1.263033e+1,1.254124e+1,1.245313e+1,1.236609e+1,
    1.228022e+1,1.219562e+1,1.211239e+1,1.203067e+1,1.195055e+1,1.187219e+1,1.179570e+1,
    1.172125e+1,1.164898e+1,1.157908e+1,1.151171e+1,1.144709e+1,1.138541e+1,1.132690e+1};
  static const G4double SL39[nL]={
    3.132934e-2,7.780169e-2,1.997199e-1,6.350103e-1,2.159047e+0,6.802554e+0,1.908355e+1,
    4.750616e+1,1.037394e+2,1.913355e+2,2.887491e+2,3.778240e+2,4.291402e+2,3.695233e+2,
    2.761193e+2,2.038392e+2,1.536530e+2,1.188331e+2,9.423086e+1,7.650062e+1,6.349328e+1,
    5.380369e+1,4.649042e+1,4.090667e+1,3.659746e+1,3.323617e+1,3.058461e+1,2.846722e+1,
    2.675392e+1,2.534816e+1,2.417843e+1,2.319185e+1,2.234937e+1,2.162211e+1,2.098852e+1,
    2.043239e+1,1.994127e+1,1.950544e+1,1.911710e+1,1.876982e+1,1.845825e+1,1.817778e+1,
    1.792446e+1,1.769479e+1,1.748570e+1,1.729449e+1,1.711874e+1,1.695632e+1,1.680535e+1,
    1.666414e+1,1.653123e+1,1.640532e+1,1.628529e+1,1.617015e+1,1.605904e+1,1.595123e+1,
    1.584607e+1,1.574304e+1,1.564166e+1,1.554157e+1,1.544242e+1,1.534396e+1,1.524598e+1,
    1.514829e+1,1.505077e+1,1.495332e+1,1.485585e+1,1.475831e+1,1.466069e+1,1.456296e+1,
    1.446512e+1,1.436719e+1,1.426920e+1,1.417118e+1,1.407316e+1,1.397520e+1,1.387735e+1,
    1.377966e+1,1.368219e+1,1.358501e+1,1.348818e+1,1.339176e+1,1.329584e+1,1.320047e+1,
    1.310574e+1,1.301172e+1,1.291849e+1,1.282613e+1,1.273474e+1,1.264439e+1,1.255519e+1,
    1.246724e+1,1.238063e+1,1.229548e+1,1.221190e+1,1.213003e+1,1.204999e+1,1.197192e+1,
    1.189599e+1,1.182235e+1,1.175118e+1,1.168267e+1,1.161701e+1,1.155444e+1,1.149517e+1};
  static const G4double SL40[nL]={
    3.209841e-2,7.976722e-2,2.053207e-1,6.550332e-1,2.231758e+0,7.040054e+0,1.977368e+1,
    4.932232e+1,1.081304e+2,2.014268e+2,3.129604e+2,4.273335e+2,4.672288e+2,3.790370e+2,
    2.780145e+2,2.044508e+2,1.540079e+2,1.191075e+2,9.446637e+1,7.671257e+1,6.368984e+1,
    5.399013e+1,4.667044e+1,4.108296e+1,3.677199e+1,3.341039e+1,3.075957e+1,2.864366e+1,
    2.693231e+1,2.552882e+1,2.436152e+1,2.337743e+1,2.253741e+1,2.181251e+1,2.118116e+1,
    2.062710e+1,2.013787e+1,1.970371e+1,1.931682e+1,1.897076e+1,1.866018e+1,1.838048e+1,
    1.812768e+1,1.789832e+1,1.768932e+1,1.749797e+1,1.732189e+1,1.715894e+1,1.700725e+1,
    1.686515e+1,1.673120e+1,1.660410e+1,1.648274e+1,1.636614e+1,1.625346e+1,1.614398e+1,
    1.603707e+1,1.593219e+1,1.582891e+1,1.572684e+1,1.562568e+1,1.552516e+1,1.542508e+1,
    1.532528e+1,1.522563e+1,1.512602e+1,1.502639e+1,1.492671e+1,1.482693e+1,1.472705e+1,
    1.462708e+1,1.452704e+1,1.442695e+1,1.432685e+1,1.422678e+1,1.412679e+1,1.402694e+1,
    1.392729e+1,1.382789e+1,1.372881e+1,1.363012e+1,1.353188e+1,1.343417e+1,1.333705e+1,
    1.324062e+1,1.314493e+1,1.305008e+1,1.295615e+1,1.286323e+1,1.277140e+1,1.268077e+1,
    1.259142e+1,1.250348e+1,1.241705e+1,1.233226e+1,1.224922e+1,1.216807e+1,1.208897e+1,
    1.201207e+1,1.193753e+1,1.186553e+1,1.179627e+1,1.172995e+1,1.166680e+1,1.160705e+1};
  static const G4double SL41[nL]={
    3.651675e-2,9.110302e-2,2.379739e-1,7.729064e-1,2.661008e+0,8.435802e+0,2.376201e+1,
    5.938433e+1,1.302605e+2,2.436285e+2,3.917132e+2,5.585883e+2,5.737641e+2,4.370072e+2,
    3.113977e+2,2.250866e+2,1.673725e+2,1.281023e+2,1.007441e+2,8.125622e+1,6.710300e+1,
    5.665354e+1,4.882987e+1,4.290028e+1,3.835570e+1,3.483427e+1,3.207439e+1,2.988467e+1,
    2.812414e+1,2.668871e+1,2.550153e+1,2.450594e+1,2.366020e+1,2.293344e+1,2.230273e+1,
    2.175080e+1,2.126446e+1,2.083343e+1,2.044949e+1,2.010594e+1,1.979722e+1,1.951859e+1,
    1.926597e+1,1.903585e+1,1.882512e+1,1.863109e+1,1.845136e+1,1.828384e+1,1.812670e+1,
    1.797832e+1,1.783729e+1,1.770239e+1,1.757256e+1,1.744688e+1,1.732458e+1,1.720499e+1,
    1.708753e+1,1.697174e+1,1.685722e+1,1.674365e+1,1.663076e+1,1.651834e+1,1.640623e+1,
    1.629429e+1,1.618244e+1,1.607062e+1,1.595877e+1,1.584690e+1,1.573499e+1,1.562305e+1,
    1.551112e+1,1.539924e+1,1.528743e+1,1.517577e+1,1.506429e+1,1.495307e+1,1.484216e+1,
    1.473164e+1,1.462157e+1,1.451203e+1,1.440308e+1,1.429481e+1,1.418729e+1,1.408060e+1,
    1.397483e+1,1.387005e+1,1.376636e+1,1.366384e+1,1.356258e+1,1.346270e+1,1.336429e+1,
    1.326746e+1,1.317233e+1,1.307902e+1,1.298767e+1,1.289841e+1,1.281140e+1,1.272680e+1,
    1.264478e+1,1.256554e+1,1.248927e+1,1.241619e+1,1.234653e+1,1.228054e+1,1.221850e+1};
  static const G4double SL42[nL]={
    3.967023e-2,9.923729e-2,2.617528e-1,8.598507e-1,2.978671e+0,9.461208e+0,2.661567e+1,
    6.608688e+1,1.423929e+2,2.556972e+2,3.845421e+2,5.509335e+2,6.958427e+2,6.038492e+2,
    4.319374e+2,3.026257e+2,2.171501e+2,1.607002e+2,1.225958e+2,9.625007e+1,7.763059e+1,
    6.421712e+1,5.439270e+1,4.709187e+1,4.159468e+1,3.740346e+1,3.416762e+1,3.163650e+1,
    2.962908e+1,2.801379e+1,2.669475e+1,2.560194e+1,2.468411e+1,2.390360e+1,2.323257e+1,
    2.265018e+1,2.214060e+1,2.169159e+1,2.129351e+1,2.093858e+1,2.062040e+1,2.033365e+1,
    2.007380e+1,1.983699e+1,1.961986e+1,1.941954e+1,1.923351e+1,1.905958e+1,1.889585e+1,
    1.874067e+1,1.859262e+1,1.845046e+1,1.831315e+1,1.817977e+1,1.804956e+1,1.792188e+1,
    1.779617e+1,1.767200e+1,1.754898e+1,1.742682e+1,1.730528e+1,1.718417e+1,1.706334e+1,
    1.694269e+1,1.682214e+1,1.670165e+1,1.658119e+1,1.646076e+1,1.634038e+1,1.622006e+1,
    1.609985e+1,1.597979e+1,1.585993e+1,1.574033e+1,1.562106e+1,1.550217e+1,1.538375e+1,
    1.526585e+1,1.514856e+1,1.503195e+1,1.491609e+1,1.480108e+1,1.468698e+1,1.457387e+1,
    1.446186e+1,1.435101e+1,1.424143e+1,1.413321e+1,1.402645e+1,1.392125e+1,1.381771e+1,
    1.371597e+1,1.361614e+1,1.351835e+1,1.342275e+1,1.332949e+1,1.323872e+1,1.315063e+1,
    1.306540e+1,1.298324e+1,1.290436e+1,1.282900e+1,1.275740e+1,1.268985e+1,1.262663e+1};
  static const G4double SL43[nL]={
    4.090181e-2,1.024236e-1,2.711450e-1,8.945313e-1,3.107339e+0,9.893526e+0,2.794798e+1,
    7.002434e+1,1.540039e+2,2.909962e+2,4.857228e+2,6.871085e+2,6.405347e+2,4.648476e+2,
    3.262939e+2,2.342891e+2,1.734786e+2,1.323680e+2,1.038662e+2,8.364716e+1,6.901712e+1,
    5.825243e+1,5.021882e+1,4.414930e+1,3.951210e+1,3.593043e+1,3.313262e+1,3.092037e+1,
    2.914791e+1,2.770777e+1,2.652069e+1,2.552830e+1,2.468759e+1,2.396683e+1,2.334242e+1,
    2.279664e+1,2.231596e+1,2.188986e+1,2.150996e+1,2.116947e+1,2.086272e+1,2.058497e+1,
    2.033214e+1,2.010070e+1,1.988760e+1,1.969016e+1,1.950605e+1,1.933324e+1,1.916994e+1,
    1.901461e+1,1.886590e+1,1.872267e+1,1.858390e+1,1.844877e+1,1.831654e+1,1.818660e+1,
    1.805846e+1,1.793169e+1,1.780594e+1,1.768095e+1,1.755648e+1,1.743238e+1,1.730850e+1,
    1.718477e+1,1.706112e+1,1.693751e+1,1.681394e+1,1.669040e+1,1.656692e+1,1.644353e+1,
    1.632027e+1,1.619719e+1,1.607435e+1,1.595181e+1,1.582965e+1,1.570791e+1,1.558669e+1,
    1.546605e+1,1.534607e+1,1.522683e+1,1.510840e+1,1.499087e+1,1.487432e+1,1.475882e+1,
    1.464448e+1,1.453138e+1,1.441961e+1,1.430927e+1,1.420046e+1,1.409328e+1,1.398785e+1,
    1.388429e+1,1.378272e+1,1.368328e+1,1.358612e+1,1.349139e+1,1.339925e+1,1.330989e+1,
    1.322351e+1,1.314030e+1,1.306050e+1,1.298434e+1,1.291209e+1,1.284402e+1,1.278043e+1};
  static const G4double SL44[nL]={
    4.170472e-2,1.045035e-1,2.772928e-1,9.171279e-1,3.188391e+0,1.013541e+1,2.846282e+1,
    7.022927e+1,1.487823e+2,2.562178e+2,3.490445e+2,4.157065e+2,4.491619e+2,3.805077e+2,
    2.830064e+2,2.087710e+2,1.575446e+2,1.221104e+2,9.711794e+1,7.913014e+1,6.595149e+1,
    5.615029e+1,4.876829e+1,4.314714e+1,3.882350e+1,3.546449e+1,3.282708e+1,3.073199e+1,
    2.904621e+1,2.767098e+1,2.653307e+1,2.557832e+1,2.476670e+1,2.406860e+1,2.346198e+1,
    2.293024e+1,2.246069e+1,2.204342e+1,2.167052e+1,2.133555e+1,2.103313e+1,2.075871e+1,
    2.050838e+1,2.027875e+1,2.006688e+1,1.987017e+1,1.968637e+1,1.951349e+1,1.934981e+1,
    1.919382e+1,1.904422e+1,1.889988e+1,1.875984e+1,1.862328e+1,1.848949e+1,1.835789e+1,
    1.822798e+1,1.809938e+1,1.797173e+1,1.784478e+1,1.771833e+1,1.759220e+1,1.746629e+1,
    1.734050e+1,1.721479e+1,1.708912e+1,1.696349e+1,1.683791e+1,1.671240e+1,1.658700e+1,
    1.646175e+1,1.633671e+1,1.621193e+1,1.608749e+1,1.596344e+1,1.583986e+1,1.571683e+1,
    1.559442e+1,1.547270e+1,1.535176e+1,1.523167e+1,1.511252e+1,1.499439e+1,1.487736e+1,
    1.476153e+1,1.464698e+1,1.453380e+1,1.442210e+1,1.431197e+1,1.420353e+1,1.409689e+1,
    1.399216e+1,1.388949e+1,1.378900e+1,1.369084e+1,1.359517e+1,1.350217e+1,1.341201e+1,
    1.332489e+1,1.324103e+1,1.316065e+1,1.308400e+1,1.301134e+1,1.294296e+1,1.287917e+1};
  static const G4double SL45[nL]={
    5.274378e-2,1.333110e-1,3.641511e-1,1.243435e+0,4.396136e+0,1.406245e+1,3.947630e+1,
    9.639532e+1,1.976332e+2,3.177821e+2,3.939795e+2,4.252091e+2,4.177958e+2,3.393139e+2,
    2.539062e+2,1.908284e+2,1.468253e+2,1.158554e+2,9.364418e+1,7.742461e+1,6.539740e+1,
    5.636581e+1,4.951281e+1,4.426627e+1,4.021621e+1,3.706322e+1,3.458554e+1,3.261752e+1,
    3.103504e+1,2.974522e+1,2.867873e+1,2.778411e+1,2.702322e+1,2.636780e+1,2.579681e+1,
    2.529440e+1,2.484848e+1,2.444965e+1,2.409041e+1,2.376468e+1,2.346741e+1,2.319435e+1,
    2.294188e+1,2.270688e+1,2.248667e+1,2.227890e+1,2.208156e+1,2.189291e+1,2.171143e+1,
    2.153584e+1,2.136504e+1,2.119811e+1,2.103426e+1,2.087283e+1,2.071330e+1,2.055522e+1,
    2.039825e+1,2.024210e+1,2.008657e+1,1.993150e+1,1.977677e+1,1.962233e+1,1.946811e+1,
    1.931412e+1,1.916036e+1,1.900686e+1,1.885366e+1,1.870082e+1,1.854839e+1,1.839645e+1,
    1.824507e+1,1.809433e+1,1.794431e+1,1.779509e+1,1.764676e+1,1.749940e+1,1.735310e+1,
    1.720795e+1,1.706403e+1,1.692142e+1,1.678022e+1,1.664052e+1,1.650241e+1,1.636599e+1,
    1.623135e+1,1.609859e+1,1.596783e+1,1.583917e+1,1.571273e+1,1.558865e+1,1.546706e+1,
    1.534810e+1,1.523193e+1,1.511872e+1,1.500865e+1,1.490192e+1,1.479874e+1,1.469935e+1,
    1.460399e+1,1.451294e+1,1.442648e+1,1.434493e+1,1.426862e+1,1.419792e+1,1.413321e+1};
  static const G4double SL46[nL]={
    5.429151e-2,1.373796e-1,3.766506e-1,1.291127e+0,4.573808e+0,1.464102e+1,4.109390e+1,
    1.001858e+2,2.044221e+2,3.254972e+2,3.979352e+2,4.222069e+2,4.093963e+2,3.318361e+2,
    2.490584e+2,1.878782e+2,1.450486e+2,1.147964e+2,9.303167e+1,7.709613e+1,6.525408e+1,
    5.634628e+1,4.957827e+1,4.439170e+1,4.038515e+1,3.726461e+1,3.481177e+1,3.286321e+1,
    3.129630e+1,3.001910e+1,2.896297e+1,2.807689e+1,2.732305e+1,2.667342e+1,2.610713e+1,
    2.560848e+1,2.516547e+1,2.476879e+1,2.441101e+1,2.408612e+1,2.378911e+1,2.351580e+1,
    2.326261e+1,2.302646e+1,2.280470e+1,2.259503e+1,2.239548e+1,2.220433e+1,2.202011e+1,
    2.184156e+1,2.166762e+1,2.149737e+1,2.133007e+1,2.116508e+1,2.100189e+1,2.084008e+1,
    2.067932e+1,2.051935e+1,2.035997e+1,2.020104e+1,2.004246e+1,1.988416e+1,1.972612e+1,
    1.956833e+1,1.941080e+1,1.925358e+1,1.909671e+1,1.894024e+1,1.878425e+1,1.862880e+1,
    1.847398e+1,1.831986e+1,1.816653e+1,1.801408e+1,1.786258e+1,1.771212e+1,1.756280e+1,
    1.741470e+1,1.726790e+1,1.712250e+1,1.697859e+1,1.683625e+1,1.669557e+1,1.655667e+1,
    1.641963e+1,1.628456e+1,1.615157e+1,1.602077e+1,1.589228e+1,1.576624e+1,1.564279e+1,
    1.552206e+1,1.540424e+1,1.528948e+1,1.517797e+1,1.506993e+1,1.496557e+1,1.486512e+1,
    1.476885e+1,1.467703e+1,1.458997e+1,1.450798e+1,1.443142e+1,1.436066e+1,1.429608e+1};
  static const G4double SL47[nL]={
    5.586443e-2,1.415214e-1,3.894297e-1,1.340049e+0,4.756205e+0,1.523387e+1,4.273995e+1,
    1.039652e+2,2.107611e+2,3.306660e+2,3.922061e+2,3.966237e+2,4.015659e+2,3.735019e+2,
    2.931095e+2,2.186308e+2,1.651709e+2,1.280396e+2,1.019243e+2,8.320469e+1,6.954857e+1,
    5.943717e+1,5.185713e+1,4.611438e+1,4.172171e+1,3.832984e+1,3.568426e+1,3.359747e+1,
    3.193051e+1,3.058030e+1,2.947046e+1,2.854456e+1,2.776093e+1,2.708881e+1,2.650532e+1,
    2.599335e+1,2.553984e+1,2.513469e+1,2.476993e+1,2.443910e+1,2.413690e+1,2.385890e+1,
    2.360137e+1,2.336109e+1,2.313534e+1,2.292175e+1,2.271830e+1,2.252324e+1,2.233509e+1,
    2.215258e+1,2.197463e+1,2.180035e+1,2.162898e+1,2.145990e+1,2.129261e+1,2.112668e+1,
    2.096180e+1,2.079771e+1,2.063423e+1,2.047122e+1,2.030859e+1,2.014627e+1,1.998426e+1,
    1.982254e+1,1.966115e+1,1.950011e+1,1.933948e+1,1.917932e+1,1.901970e+1,1.886070e+1,
    1.870239e+1,1.854486e+1,1.838819e+1,1.823247e+1,1.807778e+1,1.792421e+1,1.777185e+1,
    1.762079e+1,1.747111e+1,1.732291e+1,1.717627e+1,1.703129e+1,1.688806e+1,1.674668e+1,
    1.660725e+1,1.646987e+1,1.633466e+1,1.620173e+1,1.607120e+1,1.594322e+1,1.581792e+1,
    1.569546e+1,1.557599e+1,1.545971e+1,1.534680e+1,1.523747e+1,1.513195e+1,1.503048e+1,
    1.493334e+1,1.484080e+1,1.475318e+1,1.467081e+1,1.459405e+1,1.452327e+1,1.445890e+1};
  static const G4double SL48[nL]={
    5.639434e-2,1.429184e-1,3.937524e-1,1.356644e+0,4.818271e+0,1.543710e+1,4.331539e+1,
    1.053565e+2,2.134958e+2,3.351211e+2,4.008312e+2,4.175491e+2,4.294249e+2,3.789784e+2,
    2.880226e+2,2.137731e+2,1.618452e+2,1.258669e+2,1.005048e+2,8.227129e+1,6.893285e+1,
    5.903288e+1,5.159623e+1,4.595252e+1,4.162949e+1,3.828743e+1,3.567806e+1,3.361799e+1,
    3.197104e+1,3.063601e+1,2.953782e+1,2.862094e+1,2.784434e+1,2.717773e+1,2.659856e+1,
    2.608997e+1,2.563908e+1,2.523594e+1,2.487267e+1,2.454291e+1,2.424142e+1,2.396383e+1,
    2.370644e+1,2.346608e+1,2.324005e+1,2.302601e+1,2.282197e+1,2.262619e+1,2.243721e+1,
    2.225378e+1,2.207484e+1,2.189950e+1,2.172701e+1,2.155676e+1,2.138826e+1,2.122110e+1,
    2.105496e+1,2.088961e+1,2.072484e+1,2.056055e+1,2.039663e+1,2.023303e+1,2.006974e+1,
    1.990676e+1,1.974411e+1,1.958183e+1,1.941998e+1,1.925861e+1,1.909781e+1,1.893764e+1,
    1.877819e+1,1.861953e+1,1.846176e+1,1.830496e+1,1.814922e+1,1.799463e+1,1.784127e+1,
    1.768923e+1,1.753860e+1,1.738947e+1,1.724194e+1,1.709608e+1,1.695201e+1,1.680981e+1,
    1.666958e+1,1.653144e+1,1.639550e+1,1.626187e+1,1.613067e+1,1.600205e+1,1.587614e+1,
    1.575310e+1,1.563310e+1,1.551631e+1,1.540294e+1,1.529319e+1,1.518730e+1,1.508550e+1,
    1.498808e+1,1.489531e+1,1.480751e+1,1.472502e+1,1.464821e+1,1.457745e+1,1.451316e+1};
  static const G4double SH0[nH]={
    1.718841e-5,1.912141e-5,2.128656e-5,2.372770e-5,2.651339e-5,2.976162e-5,3.369201e-5,
    3.873597e-5,4.577051e-5,5.661516e-5,7.508997e-5,1.092699e-4,1.762839e-4,3.124886e-4,
    5.948094e-4,1.184449e-3,2.411855e-3,4.923726e-3,9.871386e-3,1.894320e-2,3.373152e-2,
    5.419455e-2,7.777948e-2,1.011811e-1,1.227807e-1,1.428966e-1,1.626818e-1,1.833195e-1,
    2.057743e-1,2.307930e-1,2.589428e-1,2.906090e-1,3.259289e-1,3.646554e-1,4.059556e-1,
    4.481828e-1,4.887166e-1,5.240358e-1,5.501959e-1,5.637401e-1,5.627614e-1,5.475832e-1,
    5.206446e-1,4.856647e-1,4.465759e-1,4.067172e-1,3.684796e-1,3.333189e-1,3.019524e-1,
    2.745971e-1,2.511726e-1,2.314485e-1,2.151395e-1,2.019637e-1,1.916740e-1,1.840748e-1,
    1.790291e-1,1.764601e-1,1.763488e-1,1.787259e-1,1.836564e-1,1.912090e-1,2.014025e-1,
    2.141163e-1,2.289594e-1,2.451064e-1,2.611598e-1,2.751583e-1,2.848795e-1,2.884723e-1,
    2.851743e-1,2.756664e-1,2.618121e-1,2.459864e-1,2.304469e-1,2.170242e-1,2.071089e-1,
    2.017331e-1,2.014838e-1,2.059886e-1,2.130499e-1,2.185478e-1,2.186039e-1,2.124513e-1,
    2.023557e-1,1.911989e-1,1.808918e-1,1.722630e-1,1.654744e-1,1.603770e-1,1.567046e-1,
    1.541608e-1,1.524546e-1,1.513189e-1,1.505256e-1,1.498980e-1,1.493175e-1,1.487199e-1,
    1.480828e-1,1.474096e-1,1.467148e-1,1.460147e-1,1.453221e-1,1.446452e-1,1.439881e-1,
    1.433514e-1,1.427339e-1,1.421336e-1,1.415477e-1,1.409739e-1,1.404099e-1,1.398539e-1,
    1.393046e-1,1.387609e-1,1.382221e-1,1.376879e-1,1.371581e-1,1.366326e-1,1.361116e-1,
    1.355952e-1,1.350837e-1,1.345775e-1,1.340767e-1,1.335816e-1,1.330926e-1,1.326099e-1,
    1.321338e-1,1.316644e-1,1.312019e-1,1.307465e-1,1.302983e-1,1.298574e-1,1.294239e-1,
    1.289978e-1,1.285792e-1,1.281681e-1,1.277645e-1,1.273684e-1,1.269797e-1,1.265984e-1,
    1.262246e-1,1.258580e-1,1.254987e-1,1.251465e-1,1.248015e-1,1.244635e-1,1.241324e-1,
    1.238082e-1,1.234908e-1,1.231801e-1,1.228760e-1,1.225784e-1,1.222872e-1,1.220024e-1,
    1.217239e-1,1.214515e-1,1.211852e-1,1.209249e-1,1.206706e-1,1.204221e-1,1.201793e-1,
    1.199423e-1,1.197109e-1,1.194850e-1,1.192646e-1,1.190497e-1,1.188400e-1,1.186357e-1,
    1.184365e-1,1.182425e-1,1.180536e-1,1.178697e-1,1.176908e-1,1.175169e-1,1.173477e-1,
    1.171834e-1,1.170239e-1,1.168690e-1,1.167189e-1,1.165733e-1,1.164323e-1,1.162959e-1,
    1.161639e-1,1.160364e-1,1.159132e-1,1.157944e-1,1.156800e-1,1.155698e-1,1.154639e-1,
    1.153622e-1,1.152646e-1,1.151712e-1,1.150819e-1,1.149967e-1,1.149155e-1,1.148384e-1,
    1.147652e-1,1.146960e-1,1.146307e-1,1.145693e-1,1.145118e-1,1.144581e-1,1.144082e-1,
    1.143621e-1,1.143198e-1,1.142812e-1,1.142464e-1,1.142152e-1,1.141877e-1,1.141639e-1,
    1.141437e-1,1.141271e-1,1.141140e-1,1.141046e-1,1.140986e-1,1.140962e-1,1.140973e-1,
    1.141019e-1,1.141099e-1,1.141214e-1,1.141363e-1,1.141546e-1,1.141763e-1,1.142013e-1};
  static const G4double SH1[nH]={
    6.668702e-2,6.471599e-2,6.280838e-2,6.096276e-2,5.917858e-2,5.745696e-2,5.580253e-2,
    5.422715e-2,5.275790e-2,5.145368e-2,5.043968e-2,4.997674e-2,5.059345e-2,5.330701e-2,
    5.988460e-2,7.281273e-2,9.420106e-2,1.234283e-1,1.561378e-1,1.872520e-1,2.145181e-1,
    2.385444e-1,2.610207e-1,2.835546e-1,3.073666e-1,3.333473e-1,3.621748e-1,3.943976e-1,
    4.304707e-1,4.707493e-1,5.154414e-1,5.645200e-1,6.175936e-1,6.737415e-1,7.313391e-1,
    7.879273e-1,8.402160e-1,8.843244e-1,9.163187e-1,9.329844e-1,9.326046e-1,9.154272e-1,
    8.836029e-1,8.406235e-1,7.905189e-1,7.371246e-1,6.836113e-1,6.323070e-1,5.847347e-1,
    5.417618e-1,5.037798e-1,4.708677e-1,4.429213e-1,4.197458e-1,4.011180e-1,3.868245e-1,
    3.766827e-1,3.705478e-1,3.683091e-1,3.698743e-1,3.751405e-1,3.839482e-1,3.960151e-1,
    4.108509e-1,4.276609e-1,4.452665e-1,4.620878e-1,4.762505e-1,4.858563e-1,4.893850e-1,
    4.860977e-1,4.762589e-1,4.610644e-1,4.423195e-1,4.220327e-1,4.020934e-1,3.841080e-1,
    3.693549e-1,3.587356e-1,3.525376e-1,3.499160e-1,3.485284e-1,3.454232e-1,3.391105e-1,
    3.303977e-1,3.210857e-1,3.125479e-1,3.053688e-1,2.995858e-1,2.949875e-1,2.913061e-1,
    2.883070e-1,2.858170e-1,2.837185e-1,2.819325e-1,2.804025e-1,2.790830e-1,2.779344e-1,
    2.769203e-1,2.760078e-1,2.751677e-1,2.743752e-1,2.736094e-1,2.728539e-1,2.720958e-1,
    2.713259e-1,2.705376e-1,2.697270e-1,2.688919e-1,2.680317e-1,2.671470e-1,2.662393e-1,
    2.653104e-1,2.643629e-1,2.633994e-1,2.624226e-1,2.614353e-1,2.604402e-1,2.594400e-1,
    2.584370e-1,2.574337e-1,2.564319e-1,2.554338e-1,2.544410e-1,2.534551e-1,2.524775e-1,
    2.515094e-1,2.505518e-1,2.496056e-1,2.486718e-1,2.477508e-1,2.468434e-1,2.459499e-1,
    2.450708e-1,2.442063e-1,2.433567e-1,2.425221e-1,2.417026e-1,2.408984e-1,2.401094e-1,
    2.393357e-1,2.385771e-1,2.378336e-1,2.371052e-1,2.363916e-1,2.356928e-1,2.350087e-1,
    2.343390e-1,2.336837e-1,2.330425e-1,2.324153e-1,2.318018e-1,2.312020e-1,2.306156e-1,
    2.300425e-1,2.294824e-1,2.289352e-1,2.284007e-1,2.278786e-1,2.273690e-1,2.268715e-1,
    2.263859e-1,2.259122e-1,2.254502e-1,2.249996e-1,2.245604e-1,2.241323e-1,2.237153e-1,
    2.233091e-1,2.229137e-1,2.225289e-1,2.221545e-1,2.217905e-1,2.214366e-1,2.210928e-1,
    2.207590e-1,2.204350e-1,2.201206e-1,2.198159e-1,2.195207e-1,2.192348e-1,2.189582e-1,
    2.186908e-1,2.184324e-1,2.181831e-1,2.179426e-1,2.177109e-1,2.174878e-1,2.172735e-1,
    2.170676e-1,2.168702e-1,2.166812e-1,2.165004e-1,2.163279e-1,2.161635e-1,2.160072e-1,
    2.158589e-1,2.157185e-1,2.155860e-1,2.154613e-1,2.153443e-1,2.152351e-1,2.151334e-1,
    2.150393e-1,2.149526e-1,2.148734e-1,2.148017e-1,2.147372e-1,2.146800e-1,2.146301e-1,
    2.145873e-1,2.145516e-1,2.145230e-1,2.145015e-1,2.144869e-1,2.144793e-1,2.144786e-1,
    2.144847e-1,2.144976e-1,2.145173e-1,2.145437e-1,2.145768e-1,2.146166e-1,2.146629e-1};
  static const G4double SH2[nH]={
    1.542383e-1,1.519749e-1,1.500571e-1,1.485008e-1,1.473380e-1,1.466231e-1,1.464406e-1,
    1.469146e-1,1.482202e-1,1.505946e-1,1.543461e-1,1.598551e-1,1.675614e-1,1.779302e-1,
    1.913930e-1,2.082703e-1,2.286959e-1,2.525746e-1,2.795999e-1,3.093368e-1,3.413423e-1,
    3.752810e-1,4.110019e-1,4.485617e-1,4.882066e-1,5.303277e-1,5.754066e-1,6.239567e-1,
    6.764587e-1,7.332866e-1,7.946123e-1,8.602851e-1,9.296820e-1,1.001541e+0,1.073807e+0,
    1.143552e+0,1.207044e+0,1.260057e+0,1.298438e+0,1.318861e+0,1.319550e+0,1.300718e+0,
    1.264525e+0,1.214592e+0,1.155254e+0,1.090830e+0,1.025091e+0,9.610003e-1,9.006681e-1,
    8.454562e-1,7.961371e-1,7.530613e-1,7.163006e-1,6.857589e-1,6.612494e-1,6.425427e-1,
    6.293909e-1,6.215310e-1,6.186705e-1,6.204555e-1,6.264246e-1,6.359503e-1,6.481772e-1,
    6.619737e-1,6.759222e-1,6.883771e-1,6.976142e-1,7.020655e-1,7.005923e-1,6.927105e-1,
    6.786831e-1,6.594410e-1,6.363628e-1,6.109982e-1,5.848221e-1,5.590708e-1,5.346711e-1,
    5.122400e-1,4.921277e-1,4.744756e-1,4.592744e-1,4.464140e-1,4.357218e-1,4.269910e-1,
    4.199996e-1,4.145245e-1,4.103489e-1,4.072685e-1,4.050942e-1,4.036540e-1,4.027937e-1,
    4.023771e-1,4.022854e-1,4.024166e-1,4.026846e-1,4.030181e-1,4.033588e-1,4.036607e-1,
    4.038884e-1,4.040156e-1,4.040242e-1,4.039024e-1,4.036443e-1,4.032483e-1,4.027164e-1,
    4.020532e-1,4.012654e-1,4.003614e-1,3.993501e-1,3.982412e-1,3.970447e-1,3.957703e-1,
    3.944278e-1,3.930263e-1,3.915747e-1,3.900812e-1,3.885536e-1,3.869989e-1,3.854235e-1,
    3.838335e-1,3.822342e-1,3.806303e-1,3.790262e-1,3.774257e-1,3.758321e-1,3.742484e-1,
    3.726773e-1,3.711210e-1,3.695814e-1,3.680602e-1,3.665588e-1,3.650785e-1,3.636202e-1,
    3.621847e-1,3.607727e-1,3.593848e-1,3.580212e-1,3.566824e-1,3.553685e-1,3.540796e-1,
    3.528157e-1,3.515768e-1,3.503629e-1,3.491738e-1,3.480093e-1,3.468693e-1,3.457535e-1,
    3.446616e-1,3.435935e-1,3.425488e-1,3.415272e-1,3.405285e-1,3.395522e-1,3.385982e-1,
    3.376660e-1,3.367553e-1,3.358660e-1,3.349975e-1,3.341496e-1,3.333220e-1,3.325144e-1,
    3.317265e-1,3.309579e-1,3.302084e-1,3.294778e-1,3.287656e-1,3.280716e-1,3.273956e-1,
    3.267373e-1,3.260964e-1,3.254728e-1,3.248660e-1,3.242759e-1,3.237024e-1,3.231450e-1,
    3.226036e-1,3.220781e-1,3.215681e-1,3.210735e-1,3.205941e-1,3.201297e-1,3.196801e-1,
    3.192451e-1,3.188245e-1,3.184182e-1,3.180260e-1,3.176477e-1,3.172832e-1,3.169323e-1,
    3.165949e-1,3.162708e-1,3.159598e-1,3.156619e-1,3.153769e-1,3.151046e-1,3.148449e-1,
    3.145978e-1,3.143631e-1,3.141406e-1,3.139302e-1,3.137319e-1,3.135455e-1,3.133709e-1,
    3.132080e-1,3.130567e-1,3.129170e-1,3.127886e-1,3.126715e-1,3.125656e-1,3.124708e-1,
    3.123871e-1,3.123143e-1,3.122524e-1,3.122012e-1,3.121607e-1,3.121308e-1,3.121115e-1,
    3.121025e-1,3.121040e-1,3.121157e-1,3.121377e-1,3.121699e-1,3.122121e-1,3.122643e-1};
  static const G4double SH3[nH]={
    2.629521e-1,2.526203e-1,2.431577e-1,2.345868e-1,2.269548e-1,2.203418e-1,2.148720e-1,
    2.107267e-1,2.081600e-1,2.075143e-1,2.092319e-1,2.138553e-1,2.220073e-1,2.343402e-1,
    2.514498e-1,2.737640e-1,3.014354e-1,3.342839e-1,3.718253e-1,4.133919e-1,4.583047e-1,
    5.060360e-1,5.563139e-1,6.091512e-1,6.648144e-1,7.237579e-1,7.865441e-1,8.537601e-1,
    9.259290e-1,1.003412e+0,1.086287e+0,1.174207e+0,1.266226e+0,1.360618e+0,1.454730e+0,
    1.544923e+0,1.626692e+0,1.695047e+0,1.745152e+0,1.773146e+0,1.776925e+0,1.756615e+0,
    1.714557e+0,1.654820e+0,1.582408e+0,1.502458e+0,1.419603e+0,1.337617e+0,1.259313e+0,
    1.186603e+0,1.120659e+0,1.062094e+0,1.011127e+0,9.677216e-1,9.316827e-1,9.027247e-1,
    8.805089e-1,8.646577e-1,8.547494e-1,8.502958e-1,8.507072e-1,8.552495e-1,8.630006e-1,
    8.728198e-1,8.833474e-1,8.930526e-1,9.003428e-1,9.037330e-1,9.020462e-1,8.945963e-1,
    8.812964e-1,8.626545e-1,8.396586e-1,8.135868e-1,7.858018e-1,7.575779e-1,7.299877e-1,
    7.038504e-1,6.797291e-1,6.579596e-1,6.386918e-1,6.219342e-1,6.075947e-1,5.955137e-1,
    5.854904e-1,5.773026e-1,5.707205e-1,5.655170e-1,5.614749e-1,5.583909e-1,5.560792e-1,
    5.543723e-1,5.531226e-1,5.522016e-1,5.514995e-1,5.509243e-1,5.504007e-1,5.498680e-1,
    5.492791e-1,5.485986e-1,5.478012e-1,5.468703e-1,5.457963e-1,5.445754e-1,5.432087e-1,
    5.417006e-1,5.400585e-1,5.382913e-1,5.364098e-1,5.344250e-1,5.323486e-1,5.301922e-1,
    5.279674e-1,5.256850e-1,5.233556e-1,5.209890e-1,5.185943e-1,5.161800e-1,5.137537e-1,
    5.113224e-1,5.088925e-1,5.064694e-1,5.040582e-1,5.016634e-1,4.992887e-1,4.969376e-1,
    4.946129e-1,4.923172e-1,4.900525e-1,4.878208e-1,4.856233e-1,4.834613e-1,4.813359e-1,
    4.792477e-1,4.771973e-1,4.751850e-1,4.732112e-1,4.712759e-1,4.693792e-1,4.675209e-1,
    4.657009e-1,4.639189e-1,4.621746e-1,4.604677e-1,4.587977e-1,4.571643e-1,4.555669e-1,
    4.540052e-1,4.524784e-1,4.509863e-1,4.495282e-1,4.481036e-1,4.467119e-1,4.453527e-1,
    4.440254e-1,4.427295e-1,4.414644e-1,4.402296e-1,4.390246e-1,4.378489e-1,4.367021e-1,
    4.355835e-1,4.344928e-1,4.334295e-1,4.323931e-1,4.313831e-1,4.303992e-1,4.294410e-1,
    4.285079e-1,4.275996e-1,4.267157e-1,4.258558e-1,4.250196e-1,4.242067e-1,4.234167e-1,
    4.226493e-1,4.219041e-1,4.211809e-1,4.204793e-1,4.197990e-1,4.191397e-1,4.185012e-1,
    4.178831e-1,4.172851e-1,4.167071e-1,4.161487e-1,4.156098e-1,4.150899e-1,4.145890e-1,
    4.141068e-1,4.136430e-1,4.131974e-1,4.127699e-1,4.123602e-1,4.119681e-1,4.115935e-1,
    4.112360e-1,4.108956e-1,4.105721e-1,4.102652e-1,4.099749e-1,4.097009e-1,4.094431e-1,
    4.092013e-1,4.089753e-1,4.087651e-1,4.085704e-1,4.083912e-1,4.082272e-1,4.080784e-1,
    4.079445e-1,4.078256e-1,4.077213e-1,4.076317e-1,4.075566e-1,4.074959e-1,4.074494e-1,
    4.074170e-1,4.073987e-1,4.073943e-1,4.074037e-1,4.074268e-1,4.074635e-1,4.075137e-1};
  static const G4double SH4[nH]={
    7.405778e-2,7.529642e-2,7.695159e-2,7.911585e-2,8.191973e-2,8.554566e-2,9.024590e-2,
    9.636447e-2,1.043620e-1,1.148402e-1,1.285589e-1,1.464340e-1,1.695011e-1,1.988282e-1,
    2.353712e-1,2.797917e-1,3.322864e-1,3.925018e-1,4.595934e-1,5.324297e-1,6.098729e-1,
    6.910366e-1,7.754430e-1,8.630582e-1,9.542326e-1,1.049589e+0,1.149890e+0,1.255911e+0,
    1.368302e+0,1.487447e+0,1.613306e+0,1.745230e+0,1.881750e+0,2.020382e+0,2.157473e+0,
    2.288176e+0,2.406638e+0,2.506475e+0,2.581539e+0,2.626868e+0,2.639586e+0,2.619491e+0,
    2.569112e+0,2.493232e+0,2.398030e+0,2.290122e+0,2.175753e+0,2.060259e+0,1.947831e+0,
    1.841504e+0,1.743297e+0,1.654408e+0,1.575420e+0,1.506477e+0,1.447431e+0,1.397940e+0,
    1.357542e+0,1.325684e+0,1.301743e+0,1.285014e+0,1.274698e+0,1.269869e+0,1.269456e+0,
    1.272221e+0,1.276770e+0,1.281590e+0,1.285130e+0,1.285911e+0,1.282675e+0,1.274519e+0,
    1.260999e+0,1.242178e+0,1.218585e+0,1.191121e+0,1.160911e+0,1.129160e+0,1.097017e+0,
    1.065486e+0,1.035376e+0,1.007278e+0,9.815797e-1,9.584915e-1,9.380753e-1,9.202798e-1,
    9.049702e-1,8.919552e-1,8.810082e-1,8.718847e-1,8.643355e-1,8.581166e-1,8.529964e-1,
    8.487605e-1,8.452151e-1,8.421882e-1,8.395306e-1,8.371154e-1,8.348367e-1,8.326085e-1,
    8.303629e-1,8.280477e-1,8.256249e-1,8.230684e-1,8.203621e-1,8.174981e-1,8.144752e-1,
    8.112975e-1,8.079728e-1,8.045118e-1,8.009274e-1,7.972336e-1,7.934448e-1,7.895759e-1,
    7.856415e-1,7.816556e-1,7.776315e-1,7.735820e-1,7.695185e-1,7.654519e-1,7.613920e-1,
    7.573474e-1,7.533262e-1,7.493353e-1,7.453809e-1,7.414684e-1,7.376024e-1,7.337870e-1,
    7.300255e-1,7.263208e-1,7.226752e-1,7.190907e-1,7.155687e-1,7.121104e-1,7.087166e-1,
    7.053878e-1,7.021244e-1,6.989265e-1,6.957938e-1,6.927263e-1,6.897235e-1,6.867848e-1,
    6.839097e-1,6.810975e-1,6.783474e-1,6.756586e-1,6.730302e-1,6.704613e-1,6.679511e-1,
    6.654984e-1,6.631024e-1,6.607622e-1,6.584766e-1,6.562449e-1,6.540659e-1,6.519387e-1,
    6.498624e-1,6.478361e-1,6.458588e-1,6.439296e-1,6.420476e-1,6.402120e-1,6.384218e-1,
    6.366764e-1,6.349747e-1,6.333162e-1,6.316999e-1,6.301252e-1,6.285912e-1,6.270973e-1,
    6.256428e-1,6.242270e-1,6.228493e-1,6.215089e-1,6.202054e-1,6.189380e-1,6.177062e-1,
    6.165095e-1,6.153472e-1,6.142189e-1,6.131240e-1,6.120620e-1,6.110324e-1,6.100348e-1,
    6.090686e-1,6.081335e-1,6.072289e-1,6.063545e-1,6.055099e-1,6.046946e-1,6.039083e-1,
    6.031506e-1,6.024211e-1,6.017194e-1,6.010453e-1,6.003983e-1,5.997781e-1,5.991845e-1,
    5.986171e-1,5.980756e-1,5.975596e-1,5.970690e-1,5.966035e-1,5.961627e-1,5.957464e-1,
    5.953544e-1,5.949863e-1,5.946420e-1,5.943211e-1,5.940236e-1,5.937490e-1,5.934973e-1,
    5.932682e-1,5.930615e-1,5.928769e-1,5.927143e-1,5.925735e-1,5.924543e-1,5.923565e-1,
    5.922798e-1,5.922242e-1,5.921894e-1,5.921754e-1,5.921818e-1,5.922085e-1,5.922555e-1};
  static const G4double SH5[nH]={
    4.659776e-1,4.476902e-1,4.309775e-1,4.158946e-1,4.025449e-1,3.910970e-1,3.818062e-1,
    3.750402e-1,3.713091e-1,3.712933e-1,3.758637e-1,3.860770e-1,4.031303e-1,4.282529e-1,
    4.625322e-1,5.066943e-1,5.609019e-1,6.246586e-1,6.968852e-1,7.761692e-1,8.611016e-1,
    9.505815e-1,1.043997e+0,1.141265e+0,1.242749e+0,1.349132e+0,1.461255e+0,1.579962e+0,
    1.705950e+0,1.839605e+0,1.980828e+0,2.128829e+0,2.281907e+0,2.437240e+0,2.590731e+0,
    2.736990e+0,2.869544e+0,2.981348e+0,3.065603e+0,3.116761e+0,3.131496e+0,3.109354e+0,
    3.052863e+0,2.967066e+0,2.858644e+0,2.734889e+0,2.602808e+0,2.468500e+0,2.336846e+0,
    2.211463e+0,2.094818e+0,1.988432e+0,1.893096e+0,1.809073e+0,1.736261e+0,1.674311e+0,
    1.622712e+0,1.580840e+0,1.547984e+0,1.523343e+0,1.506019e+0,1.494995e+0,1.489116e+0,
    1.487081e+0,1.487447e+0,1.488669e+0,1.489169e+0,1.487439e+0,1.482168e+0,1.472372e+0,
    1.457487e+0,1.437429e+0,1.412571e+0,1.383670e+0,1.351749e+0,1.317954e+0,1.283431e+0,
    1.249221e+0,1.216200e+0,1.185043e+0,1.156226e+0,1.130034e+0,1.106598e+0,1.085916e+0,
    1.067888e+0,1.052343e+0,1.039061e+0,1.027798e+0,1.018293e+0,1.010290e+0,1.003540e+0,
    9.978077e-1,9.928816e-1,9.885710e-1,9.847092e-1,9.811537e-1,9.777849e-1,9.745051e-1,
    9.712366e-1,9.679192e-1,9.645088e-1,9.609745e-1,9.572969e-1,9.534658e-1,9.494785e-1,
    9.453386e-1,9.410538e-1,9.366353e-1,9.320965e-1,9.274525e-1,9.227188e-1,9.179114e-1,
    9.130459e-1,9.081377e-1,9.032012e-1,8.982499e-1,8.932966e-1,8.883526e-1,8.834286e-1,
    8.785340e-1,8.736773e-1,8.688658e-1,8.641061e-1,8.594039e-1,8.547640e-1,8.501906e-1,
    8.456872e-1,8.412565e-1,8.369010e-1,8.326224e-1,8.284221e-1,8.243011e-1,8.202600e-1,
    8.162992e-1,8.124188e-1,8.086186e-1,8.048983e-1,8.012573e-1,7.976949e-1,7.942105e-1,
    7.908030e-1,7.874715e-1,7.842150e-1,7.810324e-1,7.779224e-1,7.748840e-1,7.719159e-1,
    7.690169e-1,7.661858e-1,7.634213e-1,7.607221e-1,7.580872e-1,7.555153e-1,7.530051e-1,
    7.505555e-1,7.481653e-1,7.458334e-1,7.435587e-1,7.413400e-1,7.391763e-1,7.370666e-1,
    7.350099e-1,7.330050e-1,7.310511e-1,7.291473e-1,7.272925e-1,7.254859e-1,7.237267e-1,
    7.220140e-1,7.203469e-1,7.187247e-1,7.171466e-1,7.156118e-1,7.141196e-1,7.126694e-1,
    7.112604e-1,7.098919e-1,7.085633e-1,7.072740e-1,7.060234e-1,7.048108e-1,7.036357e-1,
    7.024976e-1,7.013959e-1,7.003300e-1,6.992995e-1,6.983039e-1,6.973426e-1,6.964153e-1,
    6.955214e-1,6.946606e-1,6.938324e-1,6.930363e-1,6.922720e-1,6.915391e-1,6.908371e-1,
    6.901658e-1,6.895248e-1,6.889136e-1,6.883320e-1,6.877797e-1,6.872562e-1,6.867613e-1,
    6.862947e-1,6.858561e-1,6.854452e-1,6.850616e-1,6.847052e-1,6.843755e-1,6.840725e-1,
    6.837958e-1,6.835452e-1,6.833203e-1,6.831211e-1,6.829472e-1,6.827983e-1,6.826744e-1,
    6.825752e-1,6.825004e-1,6.824498e-1,6.824233e-1,6.824206e-1,6.824416e-1,6.824860e-1};
  static const G4double SH6[nH]={
    5.445765e-1,5.259720e-1,5.092147e-1,4.943922e-1,4.816568e-1,4.712475e-1,4.635187e-1,
    4.589754e-1,4.583113e-1,4.624465e-1,4.725502e-1,4.900327e-1,5.164768e-1,5.534896e-1,
    6.024632e-1,6.642824e-1,7.390624e-1,8.260352e-1,9.236699e-1,1.030017e+0,1.143159e+0,
    1.261606e+0,1.384524e+0,1.511770e+0,1.643781e+0,1.781393e+0,1.925637e+0,2.077538e+0,
    2.237925e+0,2.407239e+0,2.585320e+0,2.771179e+0,2.962748e+0,3.156646e+0,3.348019e+0,
    3.530524e+0,3.696566e+0,3.837865e+0,3.946347e+0,4.015266e+0,4.040307e+0,4.020377e+0,
    3.957835e+0,3.858084e+0,3.728662e+0,3.578110e+0,3.414904e+0,3.246656e+0,3.079662e+0,
    2.918748e+0,2.767352e+0,2.627710e+0,2.501104e+0,2.388097e+0,2.288730e+0,2.202685e+0,
    2.129397e+0,2.068133e+0,2.018030e+0,1.978119e+0,1.947324e+0,1.924455e+0,1.908200e+0,
    1.897116e+0,1.889649e+0,1.884158e+0,1.878985e+0,1.872537e+0,1.863397e+0,1.850433e+0,
    1.832895e+0,1.810468e+0,1.783286e+0,1.751879e+0,1.717086e+0,1.679941e+0,1.641548e+0,
    1.602978e+0,1.565190e+0,1.528974e+0,1.494936e+0,1.463490e+0,1.434878e+0,1.409189e+0,
    1.386391e+0,1.366358e+0,1.348893e+0,1.333755e+0,1.320678e+0,1.309384e+0,1.299598e+0,
    1.291059e+0,1.283522e+0,1.276765e+0,1.270593e+0,1.264837e+0,1.259353e+0,1.254023e+0,
    1.248754e+0,1.243470e+0,1.238117e+0,1.232657e+0,1.227063e+0,1.221321e+0,1.215428e+0,
    1.209386e+0,1.203202e+0,1.196890e+0,1.190463e+0,1.183940e+0,1.177337e+0,1.170673e+0,
    1.163967e+0,1.157235e+0,1.150494e+0,1.143761e+0,1.137048e+0,1.130370e+0,1.123738e+0,
    1.117164e+0,1.110656e+0,1.104223e+0,1.097872e+0,1.091609e+0,1.085440e+0,1.079369e+0,
    1.073400e+0,1.067535e+0,1.061776e+0,1.056126e+0,1.050585e+0,1.045154e+0,1.039834e+0,
    1.034624e+0,1.029523e+0,1.024532e+0,1.019650e+0,1.014875e+0,1.010206e+0,1.005642e+0,
    1.001182e+0,9.968228e-1,9.925643e-1,9.884044e-1,9.843414e-1,9.803736e-1,9.764992e-1,
    9.727164e-1,9.690235e-1,9.654188e-1,9.619005e-1,9.584670e-1,9.551164e-1,9.518472e-1,
    9.486577e-1,9.455463e-1,9.425115e-1,9.395517e-1,9.366654e-1,9.338511e-1,9.311074e-1,
    9.284330e-1,9.258264e-1,9.232864e-1,9.208117e-1,9.184010e-1,9.160531e-1,9.137668e-1,
    9.115411e-1,9.093748e-1,9.072668e-1,9.052161e-1,9.032217e-1,9.012826e-1,8.993979e-1,
    8.975665e-1,8.957878e-1,8.940607e-1,8.923845e-1,8.907582e-1,8.891813e-1,8.876527e-1,
    8.861720e-1,8.847382e-1,8.833507e-1,8.820089e-1,8.807120e-1,8.794595e-1,8.782507e-1,
    8.770849e-1,8.759618e-1,8.748805e-1,8.738407e-1,8.728417e-1,8.718831e-1,8.709643e-1,
    8.700849e-1,8.692443e-1,8.684421e-1,8.676779e-1,8.669512e-1,8.662615e-1,8.656085e-1,
    8.649917e-1,8.644108e-1,8.638654e-1,8.633550e-1,8.628794e-1,8.624381e-1,8.620308e-1,
    8.616572e-1,8.613169e-1,8.610097e-1,8.607351e-1,8.604929e-1,8.602828e-1,8.601045e-1,
    8.599577e-1,8.598422e-1,8.597575e-1,8.597036e-1,8.596801e-1,8.596867e-1,8.597232e-1};
  static const G4double SH7[nH]={
    6.601789e-1,6.429790e-1,6.279582e-1,6.152537e-1,6.050918e-1,5.978188e-1,5.939410e-1,
    5.941728e-1,5.994893e-1,6.111752e-1,6.308524e-1,6.604603e-1,7.021523e-1,7.580774e-1,
    8.300428e-1,9.191103e-1,1.025252e+0,1.147222e+0,1.282761e+0,1.429089e+0,1.583531e+0,
    1.744040e+0,1.909484e+0,2.079667e+0,2.255167e+0,2.437080e+0,2.626743e+0,2.825481e+0,
    3.034365e+0,3.253973e+0,3.484139e+0,3.723691e+0,3.970159e+0,4.219512e+0,4.465971e+0,
    4.701980e+0,4.918458e+0,5.105403e+0,5.252875e+0,5.352255e+0,5.397546e+0,5.386392e+0,
    5.320522e+0,5.205479e+0,5.049727e+0,4.863392e+0,4.656975e+0,4.440289e+0,4.221774e+0,
    4.008169e+0,3.804492e+0,3.614203e+0,3.439460e+0,3.281400e+0,3.140395e+0,3.016263e+0,
    2.908433e+0,2.816062e+0,2.738112e+0,2.673395e+0,2.620597e+0,2.578285e+0,2.544911e+0,
    2.518816e+0,2.498248e+0,2.481393e+0,2.466431e+0,2.451611e+0,2.435343e+0,2.416298e+0,
    2.393491e+0,2.366348e+0,2.334725e+0,2.298889e+0,2.259451e+0,2.217280e+0,2.173389e+0,
    2.128833e+0,2.084617e+0,2.041630e+0,2.000598e+0,1.962072e+0,1.926419e+0,1.893839e+0,
    1.864390e+0,1.838007e+0,1.814532e+0,1.793744e+0,1.775373e+0,1.759127e+0,1.744707e+0,
    1.731819e+0,1.720184e+0,1.709544e+0,1.699669e+0,1.690356e+0,1.681435e+0,1.672762e+0,
    1.664220e+0,1.655719e+0,1.647190e+0,1.638584e+0,1.629868e+0,1.621024e+0,1.612045e+0,
    1.602933e+0,1.593696e+0,1.584347e+0,1.574905e+0,1.565388e+0,1.555818e+0,1.546216e+0,
    1.536603e+0,1.526999e+0,1.517424e+0,1.507896e+0,1.498432e+0,1.489047e+0,1.479755e+0,
    1.470568e+0,1.461496e+0,1.452549e+0,1.443734e+0,1.435059e+0,1.426529e+0,1.418148e+0,
    1.409919e+0,1.401847e+0,1.393931e+0,1.386174e+0,1.378576e+0,1.371137e+0,1.363856e+0,
    1.356733e+0,1.349767e+0,1.342956e+0,1.336298e+0,1.329791e+0,1.323434e+0,1.317223e+0,
    1.311158e+0,1.305234e+0,1.299449e+0,1.293802e+0,1.288289e+0,1.282908e+0,1.277655e+0,
    1.272529e+0,1.267527e+0,1.262647e+0,1.257885e+0,1.253239e+0,1.248707e+0,1.244286e+0,
    1.239975e+0,1.235770e+0,1.231670e+0,1.227671e+0,1.223773e+0,1.219973e+0,1.216269e+0,
    1.212660e+0,1.209142e+0,1.205714e+0,1.202375e+0,1.199123e+0,1.195955e+0,1.192871e+0,
    1.189869e+0,1.186947e+0,1.184104e+0,1.181338e+0,1.178648e+0,1.176032e+0,1.173489e+0,
    1.171019e+0,1.168619e+0,1.166288e+0,1.164026e+0,1.161831e+0,1.159702e+0,1.157638e+0,
    1.155638e+0,1.153701e+0,1.151825e+0,1.150011e+0,1.148257e+0,1.146563e+0,1.144926e+0,
    1.143348e+0,1.141826e+0,1.140360e+0,1.138949e+0,1.137593e+0,1.136290e+0,1.135041e+0,
    1.133844e+0,1.132699e+0,1.131605e+0,1.130561e+0,1.129567e+0,1.128623e+0,1.127727e+0,
    1.126879e+0,1.126079e+0,1.125326e+0,1.124620e+0,1.123959e+0,1.123344e+0,1.122775e+0,
    1.122249e+0,1.121768e+0,1.121331e+0,1.120937e+0,1.120586e+0,1.120278e+0,1.120011e+0,
    1.119786e+0,1.119603e+0,1.119460e+0,1.119358e+0,1.119296e+0,1.119274e+0,1.119291e+0};
  static const G4double SH8[nH]={
    8.326620e-1,8.187860e-1,8.074156e-1,7.987501e-1,7.931123e-1,7.909917e-1,7.931005e-1,
    8.004391e-1,8.143676e-1,8.366683e-1,8.695766e-1,9.157404e-1,9.780603e-1,1.059371e+0,
    1.161963e+0,1.287041e+0,1.434275e+0,1.601688e+0,1.785989e+0,1.983289e+0,2.189948e+0,
    2.403252e+0,2.621760e+0,2.845287e+0,3.074662e+0,3.311376e+0,3.557223e+0,3.813967e+0,
    4.083048e+0,4.365291e+0,4.660621e+0,4.967743e+0,5.283822e+0,5.604165e+0,5.921995e+0,
    6.228396e+0,6.512543e+0,6.762334e+0,6.965444e+0,7.110740e+0,7.189809e+0,7.198283e+0,
    7.136586e+0,7.009905e+0,6.827389e+0,6.600804e+0,6.343002e+0,6.066546e+0,5.782692e+0,
    5.500809e+0,5.228183e+0,4.970100e+0,4.730096e+0,4.510269e+0,4.311603e+0,4.134245e+0,
    3.977736e+0,3.841183e+0,3.723389e+0,3.622934e+0,3.538226e+0,3.467539e+0,3.409022e+0,
    3.360729e+0,3.320633e+0,3.286663e+0,3.256758e+0,3.228933e+0,3.201359e+0,3.172457e+0,
    3.140973e+0,3.106044e+0,3.067232e+0,3.024518e+0,2.978262e+0,2.929127e+0,2.877988e+0,
    2.825830e+0,2.773655e+0,2.722401e+0,2.672884e+0,2.625763e+0,2.581520e+0,2.540465e+0,
    2.502745e+0,2.468369e+0,2.437227e+0,2.409125e+0,2.383802e+0,2.360958e+0,2.340275e+0,
    2.321429e+0,2.304110e+0,2.288023e+0,2.272905e+0,2.258519e+0,2.244664e+0,2.231170e+0,
    2.217898e+0,2.204739e+0,2.191611e+0,2.178454e+0,2.165229e+0,2.151912e+0,2.138493e+0,
    2.124976e+0,2.111368e+0,2.097687e+0,2.083952e+0,2.070187e+0,2.056416e+0,2.042664e+0,
    2.028956e+0,2.015316e+0,2.001766e+0,1.988328e+0,1.975021e+0,1.961862e+0,1.948866e+0,
    1.936048e+0,1.923419e+0,1.910988e+0,1.898765e+0,1.886757e+0,1.874968e+0,1.863402e+0,
    1.852064e+0,1.840954e+0,1.830074e+0,1.819424e+0,1.809004e+0,1.798812e+0,1.788847e+0,
    1.779107e+0,1.769589e+0,1.760290e+0,1.751207e+0,1.742337e+0,1.733676e+0,1.725221e+0,
    1.716968e+0,1.708913e+0,1.701051e+0,1.693380e+0,1.685895e+0,1.678592e+0,1.671467e+0,
    1.664517e+0,1.657738e+0,1.651125e+0,1.644675e+0,1.638385e+0,1.632251e+0,1.626270e+0,
    1.620437e+0,1.614751e+0,1.609207e+0,1.603802e+0,1.598534e+0,1.593400e+0,1.588396e+0,
    1.583519e+0,1.578768e+0,1.574140e+0,1.569631e+0,1.565240e+0,1.560964e+0,1.556801e+0,
    1.552749e+0,1.548805e+0,1.544967e+0,1.541233e+0,1.537602e+0,1.534071e+0,1.530639e+0,
    1.527303e+0,1.524063e+0,1.520915e+0,1.517860e+0,1.514895e+0,1.512018e+0,1.509229e+0,
    1.506525e+0,1.503906e+0,1.501370e+0,1.498916e+0,1.496542e+0,1.494247e+0,1.492031e+0,
    1.489891e+0,1.487827e+0,1.485838e+0,1.483923e+0,1.482080e+0,1.480310e+0,1.478609e+0,
    1.476979e+0,1.475418e+0,1.473924e+0,1.472498e+0,1.471138e+0,1.469844e+0,1.468614e+0,
    1.467449e+0,1.466346e+0,1.465306e+0,1.464328e+0,1.463411e+0,1.462555e+0,1.461758e+0,
    1.461021e+0,1.460342e+0,1.459721e+0,1.459157e+0,1.458650e+0,1.458199e+0,1.457804e+0,
    1.457464e+0,1.457179e+0,1.456948e+0,1.456770e+0,1.456645e+0,1.456574e+0,1.456554e+0};
  static const G4double SH9[nH]={
    1.425410e+0,1.421381e+0,1.420064e+0,1.421846e+0,1.427350e+0,1.437522e+0,1.453729e+0,
    1.477877e+0,1.512528e+0,1.560992e+0,1.627344e+0,1.716285e+0,1.832783e+0,1.981419e+0,
    2.165508e+0,2.386174e+0,2.641752e+0,2.927832e+0,3.238097e+0,3.565696e+0,3.904657e+0,
    4.250887e+0,4.602541e+0,4.959837e+0,5.324557e+0,5.699439e+0,6.087597e+0,6.492035e+0,
    6.915215e+0,7.358674e+0,7.822640e+0,8.305618e+0,8.803962e+0,9.311461e+0,9.819000e+0,
    1.031441e+1,1.078267e+1,1.120652e+1,1.156771e+1,1.184869e+1,1.203464e+1,1.211552e+1,
    1.208752e+1,1.195370e+1,1.172359e+1,1.141180e+1,1.103616e+1,1.061560e+1,1.016840e+1,
    9.710864e+0,9.256615e+0,8.816313e+0,8.397788e+0,8.006357e+0,7.645246e+0,7.316023e+0,
    7.018988e+0,6.753498e+0,6.518232e+0,6.311380e+0,6.130788e+0,5.974058e+0,5.838619e+0,
    5.721783e+0,5.620790e+0,5.532858e+0,5.455242e+0,5.385289e+0,5.320521e+0,5.258707e+0,
    5.197938e+0,5.136693e+0,5.073880e+0,5.008852e+0,4.941393e+0,4.871670e+0,4.800166e+0,
    4.727591e+0,4.654788e+0,4.582647e+0,4.512022e+0,4.443675e+0,4.378228e+0,4.316142e+0,
    4.257712e+0,4.203068e+0,4.152199e+0,4.104970e+0,4.061152e+0,4.020445e+0,3.982508e+0,
    3.946978e+0,3.913492e+0,3.881699e+0,3.851274e+0,3.821926e+0,3.793398e+0,3.765476e+0,
    3.737984e+0,3.710781e+0,3.683763e+0,3.656854e+0,3.630006e+0,3.603191e+0,3.576400e+0,
    3.549640e+0,3.522927e+0,3.496284e+0,3.469742e+0,3.443334e+0,3.417094e+0,3.391057e+0,
    3.365257e+0,3.339727e+0,3.314496e+0,3.289591e+0,3.265039e+0,3.240859e+0,3.217072e+0,
    3.193692e+0,3.170733e+0,3.148206e+0,3.126119e+0,3.104477e+0,3.083284e+0,3.062543e+0,
    3.042254e+0,3.022416e+0,3.003026e+0,2.984082e+0,2.965578e+0,2.947509e+0,2.929870e+0,
    2.912654e+0,2.895854e+0,2.879463e+0,2.863473e+0,2.847876e+0,2.832664e+0,2.817830e+0,
    2.803363e+0,2.789258e+0,2.775505e+0,2.762096e+0,2.749024e+0,2.736279e+0,2.723855e+0,
    2.711744e+0,2.699939e+0,2.688431e+0,2.677214e+0,2.666280e+0,2.655624e+0,2.645238e+0,
    2.635115e+0,2.625250e+0,2.615637e+0,2.606269e+0,2.597141e+0,2.588247e+0,2.579582e+0,
    2.571141e+0,2.562918e+0,2.554909e+0,2.547110e+0,2.539514e+0,2.532119e+0,2.524920e+0,
    2.517912e+0,2.511092e+0,2.504456e+0,2.498001e+0,2.491721e+0,2.485615e+0,2.479679e+0,
    2.473909e+0,2.468302e+0,2.462856e+0,2.457567e+0,2.452433e+0,2.447450e+0,2.442617e+0,
    2.437930e+0,2.433387e+0,2.428985e+0,2.424723e+0,2.420598e+0,2.416608e+0,2.412750e+0,
    2.409023e+0,2.405425e+0,2.401953e+0,2.398607e+0,2.395383e+0,2.392280e+0,2.389297e+0,
    2.386432e+0,2.383683e+0,2.381049e+0,2.378528e+0,2.376118e+0,2.373819e+0,2.371629e+0,
    2.369546e+0,2.367569e+0,2.365697e+0,2.363929e+0,2.362263e+0,2.360698e+0,2.359233e+0,
    2.357868e+0,2.356600e+0,2.355429e+0,2.354353e+0,2.353373e+0,2.352486e+0,2.351691e+0,
    2.350989e+0,2.350378e+0,2.349856e+0,2.349424e+0,2.349080e+0,2.348823e+0,2.348653e+0};
  static const G4double SH10[nH]={
    3.918292e+0,3.904931e+0,3.893792e+0,3.886847e+0,3.886858e+0,3.897612e+0,3.924175e+0,
    3.973155e+0,4.052892e+0,4.173448e+0,4.346251e+0,4.583168e+0,4.894929e+0,5.289011e+0,
    5.767472e+0,6.325587e+0,6.952077e+0,7.631192e+0,8.346046e+0,9.081993e+0,9.828955e+0,
    1.058224e+1,1.134205e+1,1.211228e+1,1.289914e+1,1.370979e+1,1.455140e+1,1.543028e+1,
    1.635126e+1,1.731704e+1,1.832759e+1,1.937949e+1,2.046524e+1,2.157253e+1,2.268371e+1,
    2.377554e+1,2.481942e+1,2.578236e+1,2.662884e+1,2.732355e+1,2.783477e+1,2.813791e+1,
    2.821866e+1,2.807489e+1,2.771701e+1,2.716665e+1,2.645395e+1,2.561414e+1,2.468397e+1,
    2.369877e+1,2.269021e+1,2.168502e+1,2.070456e+1,1.976487e+1,1.887726e+1,1.804896e+1,
    1.728387e+1,1.658334e+1,1.594670e+1,1.537185e+1,1.485564e+1,1.439420e+1,1.398317e+1,
    1.361791e+1,1.329360e+1,1.300539e+1,1.274849e+1,1.251824e+1,1.231020e+1,1.212025e+1,
    1.194464e+1,1.178007e+1,1.162372e+1,1.147334e+1,1.132717e+1,1.118402e+1,1.104314e+1,
    1.090421e+1,1.076726e+1,1.063254e+1,1.050046e+1,1.037150e+1,1.024612e+1,1.012471e+1,
    1.000757e+1,9.894853e+0,9.786587e+0,9.682662e+0,9.582857e+0,9.486860e+0,9.394299e+0,
    9.304761e+0,9.217824e+0,9.133076e+0,9.050131e+0,8.968642e+0,8.888313e+0,8.808898e+0,
    8.730205e+0,8.652092e+0,8.574465e+0,8.497273e+0,8.420498e+0,8.344153e+0,8.268275e+0,
    8.192917e+0,8.118146e+0,8.044035e+0,7.970660e+0,7.898102e+0,7.826434e+0,7.755730e+0,
    7.686057e+0,7.617475e+0,7.550037e+0,7.483791e+0,7.418775e+0,7.355021e+0,7.292555e+0,
    7.231395e+0,7.171554e+0,7.113040e+0,7.055853e+0,6.999992e+0,6.945450e+0,6.892218e+0,
    6.840281e+0,6.789624e+0,6.740229e+0,6.692076e+0,6.645144e+0,6.599408e+0,6.554846e+0,
    6.511433e+0,6.469143e+0,6.427951e+0,6.387831e+0,6.348757e+0,6.310703e+0,6.273644e+0,
    6.237554e+0,6.202409e+0,6.168183e+0,6.134852e+0,6.102393e+0,6.070783e+0,6.039998e+0,
    6.010017e+0,5.980819e+0,5.952382e+0,5.924686e+0,5.897711e+0,5.871440e+0,5.845852e+0,
    5.820930e+0,5.796658e+0,5.773017e+0,5.749993e+0,5.727569e+0,5.705730e+0,5.684462e+0,
    5.663751e+0,5.643583e+0,5.623944e+0,5.604824e+0,5.586208e+0,5.568086e+0,5.550445e+0,
    5.533276e+0,5.516567e+0,5.500309e+0,5.484490e+0,5.469103e+0,5.454137e+0,5.439585e+0,
    5.425436e+0,5.411684e+0,5.398319e+0,5.385335e+0,5.372724e+0,5.360478e+0,5.348591e+0,
    5.337056e+0,5.325867e+0,5.315016e+0,5.304499e+0,5.294309e+0,5.284441e+0,5.274889e+0,
    5.265647e+0,5.256711e+0,5.248076e+0,5.239736e+0,5.231688e+0,5.223926e+0,5.216445e+0,
    5.209243e+0,5.202315e+0,5.195656e+0,5.189262e+0,5.183131e+0,5.177258e+0,5.171640e+0,
    5.166273e+0,5.161154e+0,5.156279e+0,5.151646e+0,5.147251e+0,5.143092e+0,5.139165e+0,
    5.135467e+0,5.131996e+0,5.128750e+0,5.125724e+0,5.122918e+0,5.120328e+0,5.117952e+0,
    5.115788e+0,5.113833e+0,5.112084e+0,5.110541e+0,5.109201e+0,5.108061e+0,5.107120e+0};
  static const G4double SH11[nH]={
    7.590321e+0,7.509120e+0,7.439927e+0,7.389122e+0,7.365094e+0,7.378718e+0,7.443759e+0,
    7.577042e+0,7.798111e+0,8.128047e+0,8.587142e+0,9.191438e+0,9.948696e+0,1.085507e+1,
    1.189410e+1,1.303902e+1,1.425810e+1,1.552108e+1,1.680454e+1,1.809464e+1,1.938727e+1,
    2.068629e+1,2.200102e+1,2.334386e+1,2.472823e+1,2.616717e+1,2.767219e+1,2.925247e+1,
    3.091405e+1,3.265909e+1,3.448491e+1,3.638304e+1,3.833817e+1,4.032706e+1,4.231784e+1,
    4.426974e+1,4.613359e+1,4.785354e+1,4.937004e+1,5.062412e+1,5.156263e+1,5.214377e+1,
    5.234188e+1,5.215066e+1,5.158403e+1,5.067443e+1,4.946906e+1,4.802467e+1,4.640214e+1,
    4.466144e+1,4.285783e+1,4.103934e+1,3.924553e+1,3.750731e+1,3.584751e+1,3.428176e+1,
    3.281973e+1,3.146624e+1,3.022230e+1,2.908612e+1,2.805376e+1,2.711986e+1,2.627800e+1,
    2.552116e+1,2.484197e+1,2.423289e+1,2.368645e+1,2.319532e+1,2.275246e+1,2.235122e+1,
    2.198543e+1,2.164944e+1,2.133824e+1,2.104744e+1,2.077331e+1,2.051280e+1,2.026347e+1,
    2.002347e+1,1.979144e+1,1.956643e+1,1.934784e+1,1.913527e+1,1.892847e+1,1.872725e+1,
    1.853145e+1,1.834084e+1,1.815514e+1,1.797402e+1,1.779705e+1,1.762379e+1,1.745373e+1,
    1.728639e+1,1.712128e+1,1.695797e+1,1.679605e+1,1.663522e+1,1.647522e+1,1.631586e+1,
    1.615706e+1,1.599876e+1,1.584098e+1,1.568380e+1,1.552732e+1,1.537170e+1,1.521709e+1,
    1.506367e+1,1.491163e+1,1.476115e+1,1.461241e+1,1.446558e+1,1.432082e+1,1.417827e+1,
    1.403806e+1,1.390031e+1,1.376511e+1,1.363253e+1,1.350265e+1,1.337550e+1,1.325113e+1,
    1.312955e+1,1.301078e+1,1.289481e+1,1.278164e+1,1.267124e+1,1.256360e+1,1.245866e+1,
    1.235641e+1,1.225680e+1,1.215977e+1,1.206529e+1,1.197330e+1,1.188374e+1,1.179657e+1,
    1.171172e+1,1.162914e+1,1.154878e+1,1.147057e+1,1.139446e+1,1.132039e+1,1.124831e+1,
    1.117817e+1,1.110991e+1,1.104348e+1,1.097882e+1,1.091590e+1,1.085465e+1,1.079504e+1,
    1.073701e+1,1.068053e+1,1.062555e+1,1.057202e+1,1.051991e+1,1.046918e+1,1.041979e+1,
    1.037170e+1,1.032488e+1,1.027930e+1,1.023491e+1,1.019170e+1,1.014962e+1,1.010865e+1,
    1.006876e+1,1.002992e+1,9.992112e+0,9.955302e+0,9.919468e+0,9.884586e+0,9.850632e+0,
    9.817586e+0,9.785425e+0,9.754129e+0,9.723679e+0,9.694054e+0,9.665237e+0,9.637210e+0,
    9.609955e+0,9.583456e+0,9.557697e+0,9.532663e+0,9.508339e+0,9.484709e+0,9.461761e+0,
    9.439481e+0,9.417856e+0,9.396874e+0,9.376521e+0,9.356788e+0,9.337661e+0,9.319131e+0,
    9.301186e+0,9.283817e+0,9.267013e+0,9.250765e+0,9.235064e+0,9.219900e+0,9.205264e+0,
    9.191149e+0,9.177546e+0,9.164447e+0,9.151845e+0,9.139731e+0,9.128099e+0,9.116942e+0,
    9.106252e+0,9.096024e+0,9.086249e+0,9.076924e+0,9.068040e+0,9.059593e+0,9.051576e+0,
    9.043985e+0,9.036813e+0,9.030055e+0,9.023706e+0,9.017761e+0,9.012215e+0,9.007064e+0,
    9.002302e+0,8.997926e+0,8.993930e+0,8.990311e+0,8.987065e+0,8.984186e+0,8.981672e+0};
  static const G4double SH12[nH]={
    1.274173e+1,1.261154e+1,1.253680e+1,1.253678e+1,1.263549e+1,1.286155e+1,1.324712e+1,
    1.382513e+1,1.462465e+1,1.566474e+1,1.694815e+1,1.845745e+1,2.015574e+1,2.199299e+1,
    2.391598e+1,2.587848e+1,2.784823e+1,2.980938e+1,3.176107e+1,3.371382e+1,3.568546e+1,
    3.769757e+1,3.977278e+1,4.193303e+1,4.419831e+1,4.658588e+1,4.910951e+1,5.177874e+1,
    5.459799e+1,5.756537e+1,6.067137e+1,6.389725e+1,6.721344e+1,7.057806e+1,7.393584e+1,
    7.721799e+1,8.034325e+1,8.322076e+1,8.575479e+1,8.785141e+1,8.942639e+1,9.041346e+1,
    9.077151e+1,9.048934e+1,8.958707e+1,8.811385e+1,8.614235e+1,8.376122e+1,8.106674e+1,
    7.815525e+1,7.511689e+1,7.203142e+1,6.896584e+1,6.597376e+1,6.309588e+1,6.036126e+1,
    5.778900e+1,5.538996e+1,5.316850e+1,5.112392e+1,4.925181e+1,4.754501e+1,4.599451e+1,
    4.459006e+1,4.332071e+1,4.217515e+1,4.114200e+1,4.021009e+1,3.936858e+1,3.860711e+1,
    3.791594e+1,3.728600e+1,3.670900e+1,3.617743e+1,3.568464e+1,3.522478e+1,3.479282e+1,
    3.438448e+1,3.399619e+1,3.362496e+1,3.326834e+1,3.292429e+1,3.259109e+1,3.226728e+1,
    3.195155e+1,3.164275e+1,3.133978e+1,3.104162e+1,3.074731e+1,3.045594e+1,3.016670e+1,
    2.987885e+1,2.959177e+1,2.930495e+1,2.901801e+1,2.873069e+1,2.844289e+1,2.815459e+1,
    2.786589e+1,2.757701e+1,2.728821e+1,2.699984e+1,2.671229e+1,2.642597e+1,2.614131e+1,
    2.585875e+1,2.557869e+1,2.530156e+1,2.502772e+1,2.475752e+1,2.449130e+1,2.422931e+1,
    2.397182e+1,2.371903e+1,2.347111e+1,2.322821e+1,2.299044e+1,2.275787e+1,2.253056e+1,
    2.230854e+1,2.209182e+1,2.188038e+1,2.167419e+1,2.147321e+1,2.127737e+1,2.108662e+1,
    2.090087e+1,2.072003e+1,2.054401e+1,2.037271e+1,2.020604e+1,2.004387e+1,1.988611e+1,
    1.973265e+1,1.958337e+1,1.943817e+1,1.929694e+1,1.915958e+1,1.902597e+1,1.889600e+1,
    1.876959e+1,1.864662e+1,1.852699e+1,1.841062e+1,1.829740e+1,1.818725e+1,1.808008e+1,
    1.797579e+1,1.787432e+1,1.777556e+1,1.767946e+1,1.758593e+1,1.749489e+1,1.740629e+1,
    1.732004e+1,1.723608e+1,1.715436e+1,1.707480e+1,1.699735e+1,1.692196e+1,1.684856e+1,
    1.677711e+1,1.670754e+1,1.663983e+1,1.657390e+1,1.650973e+1,1.644726e+1,1.638646e+1,
    1.632728e+1,1.626968e+1,1.621363e+1,1.615908e+1,1.610601e+1,1.605438e+1,1.600415e+1,
    1.595530e+1,1.590779e+1,1.586160e+1,1.581669e+1,1.577303e+1,1.573061e+1,1.568940e+1,
    1.564937e+1,1.561049e+1,1.557275e+1,1.553612e+1,1.550058e+1,1.546611e+1,1.543269e+1,
    1.540030e+1,1.536892e+1,1.533854e+1,1.530912e+1,1.528067e+1,1.525316e+1,1.522657e+1,
    1.520089e+1,1.517611e+1,1.515221e+1,1.512917e+1,1.510699e+1,1.508564e+1,1.506512e+1,
    1.504542e+1,1.502651e+1,1.500840e+1,1.499106e+1,1.497449e+1,1.495868e+1,1.494362e+1,
    1.492928e+1,1.491568e+1,1.490279e+1,1.489061e+1,1.487912e+1,1.486832e+1,1.485821e+1,
    1.484876e+1,1.483998e+1,1.483185e+1,1.482437e+1,1.481753e+1,1.481133e+1,1.480575e+1};
  static const G4double SH13[nH]={
    1.444282e+1,1.433200e+1,1.430197e+1,1.437852e+1,1.459264e+1,1.497946e+1,1.557548e+1,
    1.641373e+1,1.751708e+1,1.889114e+1,2.051932e+1,2.236282e+1,2.436683e+1,2.647125e+1,
    2.862233e+1,3.078116e+1,3.292713e+1,3.505686e+1,3.718040e+1,3.931668e+1,4.148940e+1,
    4.372398e+1,4.604543e+1,4.847705e+1,5.103960e+1,5.375059e+1,5.662368e+1,5.966783e+1,
    6.288632e+1,6.627543e+1,6.982292e+1,7.350621e+1,7.729059e+1,8.112759e+1,8.495381e+1,
    8.869083e+1,9.224647e+1,9.551803e+1,9.839764e+1,1.007797e+2,1.025696e+2,1.036932e+2,
    1.041041e+2,1.037897e+2,1.027723e+2,1.011066e+2,9.887372e+1,9.617299e+1,9.311251e+1,
    8.980074e+1,8.633954e+1,8.281934e+1,7.931645e+1,7.589222e+1,7.259351e+1,6.945404e+1,
    6.649625e+1,6.373321e+1,6.117052e+1,5.880801e+1,5.664119e+1,5.466241e+1,5.286187e+1,
    5.122830e+1,4.974958e+1,4.841316e+1,4.720643e+1,4.611690e+1,4.513247e+1,4.424158e+1,
    4.343327e+1,4.269734e+1,4.202440e+1,4.140591e+1,4.083420e+1,4.030251e+1,3.980491e+1,
    3.933629e+1,3.889231e+1,3.846927e+1,3.806404e+1,3.767398e+1,3.729683e+1,3.693063e+1,
    3.657365e+1,3.622434e+1,3.588127e+1,3.554314e+1,3.520875e+1,3.487703e+1,3.454699e+1,
    3.421782e+1,3.388880e+1,3.355941e+1,3.322925e+1,3.289811e+1,3.256591e+1,3.223270e+1,
    3.189868e+1,3.156414e+1,3.122945e+1,3.089506e+1,3.056146e+1,3.022918e+1,2.989873e+1,
    2.957065e+1,2.924544e+1,2.892361e+1,2.860559e+1,2.829182e+1,2.798266e+1,2.767845e+1,
    2.737947e+1,2.708598e+1,2.679819e+1,2.651624e+1,2.624029e+1,2.597040e+1,2.570666e+1,
    2.544908e+1,2.519768e+1,2.495244e+1,2.471333e+1,2.448028e+1,2.425323e+1,2.403210e+1,
    2.381680e+1,2.360721e+1,2.340323e+1,2.320475e+1,2.301165e+1,2.282379e+1,2.264106e+1,
    2.246332e+1,2.229045e+1,2.212233e+1,2.195881e+1,2.179978e+1,2.164512e+1,2.149469e+1,
    2.134838e+1,2.120607e+1,2.106765e+1,2.093300e+1,2.080201e+1,2.067458e+1,2.055060e+1,
    2.042998e+1,2.031261e+1,2.019839e+1,2.008725e+1,1.997909e+1,1.987382e+1,1.977137e+1,
    1.967165e+1,1.957458e+1,1.948010e+1,1.938812e+1,1.929859e+1,1.921144e+1,1.912659e+1,
    1.904399e+1,1.896358e+1,1.888531e+1,1.880911e+1,1.873493e+1,1.866273e+1,1.859245e+1,
    1.852404e+1,1.845746e+1,1.839266e+1,1.832961e+1,1.826826e+1,1.820857e+1,1.815050e+1,
    1.809402e+1,1.803908e+1,1.798567e+1,1.793373e+1,1.788325e+1,1.783418e+1,1.778651e+1,
    1.774019e+1,1.769521e+1,1.765154e+1,1.760914e+1,1.756801e+1,1.752810e+1,1.748940e+1,
    1.745188e+1,1.741553e+1,1.738032e+1,1.734624e+1,1.731325e+1,1.728134e+1,1.725050e+1,
    1.722071e+1,1.719194e+1,1.716418e+1,1.713742e+1,1.711164e+1,1.708682e+1,1.706295e+1,
    1.704001e+1,1.701800e+1,1.699688e+1,1.697667e+1,1.695733e+1,1.693885e+1,1.692124e+1,
    1.690446e+1,1.688852e+1,1.687340e+1,1.685909e+1,1.684558e+1,1.683285e+1,1.682090e+1,
    1.680973e+1,1.679931e+1,1.678964e+1,1.678072e+1,1.677252e+1,1.676505e+1,1.675830e+1};
  static const G4double* SL[nLA]={
    SL0,SL1,SL2,SL3,SL4,SL5,SL6,SL7,SL8,SL9,SL10,SL11,SL12,SL13,SL14,SL15,SL16,SL17,SL18,
    SL19,SL20,SL21,SL22,SL23,SL24,SL25,SL26,SL27,SL28,SL29,SL30,SL31,SL32,SL33,SL34,SL35,
    SL36,SL37,SL38,SL39,SL40,SL41,SL42,SL43,SL44,SL45,SL46,SL47,SL48};
  static const G4double* SH[nHA]={SH0,SH1,SH2,SH3,SH4,SH5,SH6,SH7,SH8,SH9,SH10,
                                  SH11,SH12,SH13};
  if(a<=.9)
  {
    G4cout<<"***G4QPhotonNuclearCS::GetFunctions: A="<<a<<"(?). No CS returned!"<<G4endl;
    return -1;
  }
  G4int r=0;                            // Low channel for GDR (filling-flag for GDR)
  for(G4int i=0; i<nLA; i++) if(std::fabs(a-LA[i])<.0005)
  {
    for(G4int k=0; k<nL; k++) y[k]=SL[i][k];
    r=1;                                // Flag of filled GDR part 
  }
  G4int h=0;
  for(G4int j=0; j<nHA; j++) if(std::fabs(a-HA[j])<.0005)
  {
    for(G4int k=0; k<nH; k++) z[k]=SH[j][k];
    h=1;                                // Flag of filled GDR part 
  }
  if(!r)                                // GDR part is not filled
  {
    G4int k=0;                          // !! To be good for different compilers !!
    for(k=1; k<nLA; k++) if(a<LA[k]) break;
    if(k<1) k=1;                        // Extrapolation from the first bin (D/He)
    if(k>=nLA) k=nLA-1;                 // Extrapolation from the last bin (U)
    G4int     k1=k-1;
    G4double  xi=LA[k1];
    G4double   b=(a-xi)/(LA[k]-xi);
    for(G4int j=0; j<nL; j++)
    {
      if(a>1.5)
      {
        G4double yi=SL[k1][j];
        y[j]=yi+(SL[k][j]-yi)*b;
#ifdef debugs
        if(y[j]<0.)G4cout<<"G4QPhotNucCS::GetF:y="<<y[j]<<",k="<<k<<",yi="<<yi<<",ya="
                         <<SL[k][j]<<",b="<<b<<",xi="<<xi<<",xa="<<LA[k]<<",a="<<a<<G4endl;
#endif
      }
      else y[j]=0.;
    }
    r=1;
  }
  if(!h)                                // High Energy part is not filled
  {
    G4int k=0;
    for(k=1; k<nHA; k++) if(a<HA[k]) break;
    if(k<1) k=1;                        // Extrapolation from the first bin (D/He)
    if(k>=nHA) k=nHA-1;                 // Extrapolation from the last bin (Pu)
    G4int     k1=k-1;
    G4double  xi=HA[k1];
    G4double   b=(a-xi)/(HA[k]-xi);
    for(G4int j=0; j<nH; j++)
    {
      G4double zi=SH[k1][j];
      z[j]=zi+(SH[k][j]-zi)*b;
    }
    h=1;
  }
  return r*h;
}
