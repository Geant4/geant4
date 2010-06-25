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
// The lust update: M.V. Kossov, CERN/ITEP(Moscow) 19-Aug-07
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Physics class: G4QIonIonCrossSection for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-03
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Feb-04
// --------------------------------------------------------------------------------
// ****************************************************************************************
// This Header is a part of the CHIPS physics package (author: M. Kosov)
// ****************************************************************************************
// Short description: CHIPS cross-sectons for Ion-Ion interactions
// ---------------------------------------------------------------
//
//#define debug
//#define pdebug
//#define debug3
//#define debugn
//#define debugs

#include "G4QIonIonCrossSection.hh"

// Initialization of the
G4double* G4QIonIonCrossSection::lastLENI=0;// Pointer to the lastArray of LowEn Inelast CS
G4double* G4QIonIonCrossSection::lastHENI=0;// Pointer to the lastArray of HighEn InelastCS
G4double* G4QIonIonCrossSection::lastLENE=0;// Pointer to the lastArray of LowEn Elastic CS
G4double* G4QIonIonCrossSection::lastHENE=0;// Pointer to the lastArray of HighEn ElasticCS
G4int     G4QIonIonCrossSection::lastPDG=0; // The last PDG code of the projectile
G4int     G4QIonIonCrossSection::lastN=0;   // The last N of calculated nucleus
G4int     G4QIonIonCrossSection::lastZ=0;   // The last Z of calculated nucleus
G4double  G4QIonIonCrossSection::lastP=0.;  // Last used in cross section Momentum
G4double  G4QIonIonCrossSection::lastTH=0.; // Last threshold momentum
G4double  G4QIonIonCrossSection::lastICS=0.;// Last value of the Inelastic Cross Section
G4double  G4QIonIonCrossSection::lastECS=0.;// Last value of the Elastic Cross Section
G4int     G4QIonIonCrossSection::lastI=0;   // The last position in the DAMDB

// Returns Pointer to the G4VQCrossSection class
G4VQCrossSection* G4QIonIonCrossSection::GetPointer()
{
  static G4QIonIonCrossSection theCrossSection; //**Static body of Cross Section**
  return &theCrossSection;
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units !(Now it is MeV): fCS=true->Inelastic, fCS=false->Elastic
G4double G4QIonIonCrossSection::GetCrossSection(G4bool fCS, G4double pMom, G4int tZ,
                                                G4int tN, G4int pPDG)
{
  static G4int j;                      // A#0f records found in DB for this projectile
  static std::vector <G4int>    colPDG;// Vector of the projectile PDG code
  static std::vector <G4int>    colN;  // Vector of N for calculated nuclei (isotops)
  static std::vector <G4int>    colZ;  // Vector of Z for calculated nuclei (isotops)
  static std::vector <G4double> colP;  // Vector of last momenta for the reaction
  static std::vector <G4double> colTH; // Vector of energy thresholds for the reaction
  static std::vector <G4double> colICS;// Vector of last inelastic cross-sections
  static std::vector <G4double> colECS;// Vector of last elastic cross-sections
  // ***---*** End of the mandatory Static Definitions of the Associative Memory ***---***
#ifdef pdebug
  G4cout<<"G4QIICS::GetCS:>>> f="<<fCS<<", Z="<<tZ<<"("<<lastZ<<"), N="<<tN<<"("<<lastN
        <<"),PDG="<<pPDG<<"("<<lastPDG<<"), p="<<pMom<<"("<<lastTH<<")"<<",Sz="
        <<colN.size()<<G4endl;
#endif
  if(!pPDG)
  {
#ifdef pdebug
    G4cout<<"G4QIonIonCS::GetCS: *** Found pPDG="<<pPDG<<" ====> CS=0"<<G4endl;
#endif
    return 0.;                         // projectile PDG=0 is a mistake (?!) @@
  }
  G4bool in=false;                     // By default the isotope must be found in the AMDB
  if(tN!=lastN || tZ!=lastZ || pPDG!=lastPDG)// The nucleus was not the last used isotope
  {
    in = false;                        // By default the isotope haven't be found in AMDB  
    lastP   = 0.;                      // New momentum history (nothing to compare with)
    lastPDG = pPDG;                    // The last PDG of the projectile
    lastN   = tN;                      // The last N of the calculated nucleus
    lastZ   = tZ;                      // The last Z of the calculated nucleus
    lastI   = colN.size();             // Size of the Associative Memory DB in the heap
    j  = 0;                            // A#0f records found in DB for this projectile
#ifdef pdebug
    G4cout<<"G4QIICS::GetCS:FindI="<<lastI<<",pPDG="<<pPDG<<",tN="<<tN<<",tZ="<<tZ<<G4endl;
#endif
    if(lastI) for(G4int i=0; i<lastI; i++) // Loop over all DB
    {                                  // The nucleus with projPDG is found in AMDB
#ifdef pdebug
      G4cout<<"G4QII::GCS:P="<<colPDG[i]<<",N="<<colN[i]<<",Z="<<colZ[i]<<",j="<<j<<G4endl;
#endif
      if(colPDG[i]==pPDG && colN[i]==tN && colZ[i]==tZ)
      {
        lastI=i;
        lastTH =colTH[i];                // Last THreshold (A-dependent)
#ifdef pdebug
        G4cout<<"G4QIICS::GetCS:*Found* P="<<pMom<<",Threshold="<<lastTH<<",j="<<j<<G4endl;
#endif
        if(pMom<=lastTH)
        {
#ifdef pdebug
          G4cout<<"G4QIICS::GetCS:Found P="<<pMom<<"<Threshold="<<lastTH<<"->XS=0"<<G4endl;
#endif
          return 0.;                     // Energy is below the Threshold value
        }
        lastP  =colP [i];                // Last Momentum  (A-dependent)
        lastICS=colICS[i];               // Last Inelastic Cross-Section (A-dependent)
        lastECS=colECS[i];               // Last Elastic Cross-Section (A-dependent)
        if(std::fabs(lastP/pMom-1.)<tolerance)
        {
#ifdef pdebug
          G4cout<<"G4QIonIonCS::GetCS:P="<<pMom<<",InXS="<<lastICS*millibarn<<",ElXS="
                <<lastECS*millibarn<<G4endl;
#endif
          CalculateCrossSection(fCS,-1,j,lastPDG,lastZ,lastN,pMom); // Update param's only
          if(fCS) return lastICS*millibarn;     // Use theLastInelasticCS
          return         lastECS*millibarn;     // Use theLastElasticCS
        }
        in = true;                       // This is the case when the isotop is found in DB
        // Momentum pMom is in IU ! @@ Units
#ifdef pdebug
        G4cout<<"G4QIICS::G:UpdatDB P="<<pMom<<",f="<<fCS<<",lI="<<lastI<<",j="<<j<<G4endl;
#endif
        lastICS=CalculateCrossSection( true,-1,j,lastPDG,lastZ,lastN,pMom);// read & update
        lastECS=CalculateCrossSection(false,-1,j,lastPDG,lastZ,lastN,pMom);// read & update
#ifdef pdebug
        G4cout<<"G4QIonIonCS::GetCS:=>New(inDB) InCS="<<lastICS<<",ElCS="<<lastECS<<G4endl;
#endif
        if((lastICS<=0. || lastECS<=0.) && pMom>lastTH) // Correct the threshold
        {
#ifdef pdebug
          G4cout<<"G4QIonIonCS::GetCS:New,T="<<pMom<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pMom;
        }
        break;                           // Go out of the LOOP
      }
#ifdef pdebug
      G4cout<<"--->G4QIonIonCrossSec::GetCrosSec: pPDG="<<pPDG<<",j="<<j<<",N="<<colN[i]
            <<",Z["<<i<<"]="<<colZ[i]<<",PDG="<<colPDG[i]<<G4endl;
#endif
      j++;                             // Increment a#0f records found in DB for this pPDG
    }
    if(!in)                            // This nucleus has not been calculated previously
    {
#ifdef pdebug
      G4cout<<"G4QIICS::GetCrosSec:CalcNew P="<<pMom<<",f="<<fCS<<",lastI="<<lastI<<G4endl;
#endif
      //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
      lastICS=CalculateCrossSection(true ,0,j,lastPDG,lastZ,lastN,pMom); //calculate&create
      lastECS=CalculateCrossSection(false,0,j,lastPDG,lastZ,lastN,pMom); //calculate&create
      if(lastICS<=0. || lastECS<=0.)
      {
        lastTH = ThresholdEnergy(tZ, tN); // Threshold Energy=Mom=0 which is now the last
#ifdef pdebug
        G4cout<<"G4QIonIonCrossSect::GetCrossSect:NewThresh="<<lastTH<<",P="<<pMom<<G4endl;
#endif
        if(pMom>lastTH)
        {
#ifdef pdebug
          G4cout<<"G4QIonIonCS::GetCS:1-st,P="<<pMom<<">Thresh="<<lastTH<<"->XS=0"<<G4endl;
#endif
          lastTH=pMom;
        }
      }
#ifdef pdebug
      G4cout<<"G4QIICS::GetCS: *New* ICS="<<lastICS<<", ECS="<<lastECS<<",N="<<lastN<<",Z="
            <<lastZ<<G4endl;
#endif
      colN.push_back(tN);
      colZ.push_back(tZ);
      colPDG.push_back(pPDG);
      colP.push_back(pMom);
      colTH.push_back(lastTH);
      colICS.push_back(lastICS);
      colECS.push_back(lastECS);
#ifdef pdebug
      G4cout<<"G4QIICS::GetCS:*1st*, P="<<pMom<<"(MeV), InCS="<<lastICS*millibarn
            <<", ElCS="<<lastECS*millibarn<<"(mb)"<<G4endl;
#endif
      if(fCS) return lastICS*millibarn;     // Use theLastInelasticCS
      return         lastECS*millibarn;     // Use theLastElasticCS
    } // End of creation of the new set of parameters
    else
    {
#ifdef pdebug
      G4cout<<"G4QIICS::GetCS: Update lastI="<<lastI<<",j="<<j<<G4endl;
#endif
      colP[lastI]=pMom;
      colPDG[lastI]=pPDG;
      colICS[lastI]=lastICS;
      colECS[lastI]=lastECS;
    }
  } // End of parameters udate
  else if(pMom<=lastTH)
  {
#ifdef pdebug
    G4cout<<"G4QIICS::GetCS: Current T="<<pMom<<" < Threshold="<<lastTH<<", CS=0"<<G4endl;
#endif
    return 0.;                         // Momentum is below the Threshold Value -> CS=0
  }
  else if(std::fabs(lastP/pMom-1.)<tolerance)
  {
#ifdef pdebug
    G4cout<<"G4QIICS::GetCS:OldCur P="<<pMom<<"="<<pMom<<", InCS="<<lastICS*millibarn
          <<", ElCS="<<lastECS*millibarn<<"(mb)"<<G4endl;
#endif
    if(fCS) return lastICS*millibarn;     // Use theLastInelasticCS
    return         lastECS*millibarn;     // Use theLastElasticCS
  }
  else
  {
#ifdef pdebug
    G4cout<<"G4QIICS::GetCS:UpdatCur P="<<pMom<<",f="<<fCS<<",I="<<lastI<<",j="<<j<<G4endl;
#endif
    lastICS=CalculateCrossSection( true,1,j,lastPDG,lastZ,lastN,pMom); // Only UpdateDB
    lastECS=CalculateCrossSection(false,1,j,lastPDG,lastZ,lastN,pMom); // Only UpdateDB
    lastP=pMom;
  }
#ifdef pdebug
  G4cout<<"G4QIICS::GetCroSec:*End*,P="<<pMom<<"(MeV), InCS="<<lastICS*millibarn<<", ElCS="
        <<lastECS*millibarn<<"(mb)"<<G4endl;
#endif
    if(fCS) return lastICS*millibarn;     // Use theLastInelasticCS
    return         lastECS*millibarn;     // Use theLastElasticCS
}

// The main member function giving the A-A cross section (Momentum in MeV, CS in mb)
G4double G4QIonIonCrossSection::CalculateCrossSection(G4bool XS,G4int F,G4int I,G4int pPDG,
                                                      G4int tZ,G4int tN, G4double TotMom)
{
  //static const G4double third=1./3.; // power for A^P->R conversion [R=1.1*A^(1/3)]
  //static const G4double conv=38.; // coeff. R2->sig=c*(pR+tR)^2, c=pi*10(mb/fm^2)*1.21
  // If change the following, please change in ::GetFunctions:
  static const G4double THmin=0.;  // @@ start from threshold (?) minimum Energy Threshold
  static const G4double dP=10.;    // step for the LEN table
  static const G4int    nL=100;    // A#of LENesonance points in E (each MeV from 2 to 106)
  static const G4double Pmin=THmin+(nL-1)*dP; // minE for the HighE part
  static const G4double Pmax=300000.;   // maxE for the HighE part
  static const G4int    nH=100;         // A#of HResonance points in lnE
  static const G4double milP=std::log(Pmin); // Low logarithm energy for the HighE part
  static const G4double malP=std::log(Pmax); // High logarithm energy (each 2.75 percent)
  static const G4double dlP=(malP-milP)/(nH-1); // Step in log energy in the HighE part
  //
  // Associative memory for acceleration
  static std::vector <G4double*> LENI;   // Vector of pointers: LowEnIneIonIonCrossSection
  static std::vector <G4double*> HENI;   // Vector of pointers: HighEnIneIonIonCrossSection
  static std::vector <G4double*> LENE;   // Vector of pointers: LowEnElaIonIonCrossSection
  static std::vector <G4double*> HENE;   // Vector of pointers: HighEnElaIonIonCrossSection
#ifdef debug
  G4cout<<"G4QIonIonCrossSection::CalcCS: Z="<<tZ<<", N="<<tN<<", P="<<TotMom<<G4endl;
#endif
  G4int dPDG=pPDG/10;                // 10SZZZAAA
  G4int zPDG=dPDG/1000;              // 10SZZZ (?)
  G4int zA=dPDG%1000;                // proj A
  G4int pZ=zPDG%1000;                // proj Z (?)
  G4int pN=zA-pZ;                    // proj N (?)
  G4double Momentum=TotMom/zA;       // Momentum per nucleon
  if (Momentum<THmin) return 0.;     // @@ This can be dangerouse for the heaviest nuc.!
  G4double sigma=0.;
  if(F&&I) sigma=0.;                 // @@ *!* Fake line *!* to use F & I !!!Temporary!!!
  G4double tA=tN+tZ;                 // Target weight
  G4double pA=zA;                    // Projectile weight
  if(F<=0)                           // This isotope was not the last used isotop
  {
    if(F<0 || !XS)                   // This isotope was found in DAMDB or Elast =>RETRIEVE
    {
      lastLENI=LENI[I];              // Pointer to Low Energy inelastic cross sections
      lastHENI=HENI[I];              // Pointer to High Energy inelastic cross sections
      lastLENE=LENE[I];              // Pointer to Low Energy inelastic cross sections
      lastHENE=HENE[I];              // Pointer to High Energy inelastic cross sections
    }
    else                             // This isotope wasn't calculated previously => CREATE
    {
      lastLENI = new G4double[nL];   // Allocate memory for the new LEN cross sections
      lastHENI = new G4double[nH];   // Allocate memory for the new HEN cross sections
      lastLENE = new G4double[nL];   // Allocate memory for the new LEN cross sections
      lastHENE = new G4double[nH];   // Allocate memory for the new HEN cross sections
      G4int er=GetFunctions(pZ,pN,tZ,tN,lastLENI,lastHENI,lastLENE,lastHENE);
      if(er<1) G4cerr<<"*W*G4QIonIonCroSec::CalcCrossSection: pA="<<tA<<",tA="<<tA<<G4endl;
#ifdef debug
      G4cout<<"G4QIonIonCrossSection::CalcCS: GetFunctions er="<<er<<",pA="<<pA<<",tA="<<tA
            <<G4endl;
#endif
      // *** The synchronization check ***
      G4int sync=LENI.size();
      if(sync!=I) G4cout<<"*W*G4IonIonCrossSec::CalcCrossSect:Sync="<<sync<<"#"<<I<<G4endl;
      LENI.push_back(lastLENI);      // added LEN Inelastic
      HENI.push_back(lastHENI);      // added HEN Inelastic
      LENE.push_back(lastLENE);      // added LEN Elastic
      HENE.push_back(lastHENE);      // added HEN Elastic
    } // End of creation of the new set of parameters
  } // End of parameters udate
  // ============================== NOW the Magic Formula =================================
  if (Momentum<lastTH) return 0.;    // It must be already checked in the interface class
  else if (Momentum<Pmin)            // LEN region (approximated in E, not in lnE)
  {
#ifdef debug
    G4cout<<"G4QIICS::CalCS:p="<<pA<<",t="<<tA<<",n="<<nL<<",T="<<THmin<<",d="<<dP<<G4endl;
#endif
    if(tA<1. || pA<1.)
    {
      G4cout<<"-Warning-G4QIICS::CalcCS: pA="<<pA<<" or tA="<<tA<<" aren't nuclei"<<G4endl;
      sigma=0.;
    }
    else
    {
      G4double dPp=dP*pA;
      if(XS) sigma=EquLinearFit(Momentum,nL,THmin,dPp,lastLENI);
      else   sigma=EquLinearFit(Momentum,nL,THmin,dPp,lastLENE);
    }
#ifdef debugn
    if(sigma<0.) G4cout<<"-Warning-G4QIICS::CalcCS:pA="<<pA<<",tA="<<tA<<",XS="<<XS<<",P="
                       <<Momentum<<", Th="<<THmin<<", dP="<<dP<<G4endl;
#endif
  }
  else if (Momentum<Pmax*pA)                     // High Energy region
  {
    G4double lP=std::log(Momentum);
#ifdef debug
    G4cout<<"G4QIonIonCS::CalcCS:before HEN nH="<<nH<<",iE="<<milP<<",dlP="<<dlP<<G4endl;
#endif
    if(tA<1. || pA<1.)
    {
      G4cout<<"-Warning-G4QIICS::CalCS:pA="<<pA<<" or tA="<<tA<<" aren't composit"<<G4endl;
      sigma=0.;
    }
    else
    {
      G4double milPp=milP+std::log(pA);
      if(XS) sigma=EquLinearFit(lP,nH,milPp,dlP,lastHENI);
      else   sigma=EquLinearFit(lP,nH,milPp,dlP,lastHENE);
    }
  }
  else                                      // UltraHighE region (not frequent)
  {
    std::pair<G4double, G4double> inelel = CalculateXS(pZ, pN, tZ, tN, Momentum);
    if(XS) sigma=inelel.first;
    else   sigma=inelel.second;
  }
#ifdef debug
  G4cout<<"G4IonIonCrossSection::CalculateCrossSection: sigma="<<sigma<<G4endl;
#endif
  if(sigma<0.) return 0.;
  return sigma;
}

// Linear fit for YN[N] tabulated (from X0 with fixed step DX) function to X point

// Calculate the functions for the log(A)
G4int G4QIonIonCrossSection::GetFunctions(G4int pZ,G4int pN,G4int tZ,G4int tN,G4double* li,
                                          G4double* hi, G4double* le, G4double* he)
{
  // If change the following, please change in ::CalculateCrossSection:
  static const G4double THmin=0.;  // @@ start from threshold (?) minimum Energy Threshold
  static const G4double dP=10.;    // step for the LEN table
  static const G4int    nL=100;    // A#of LENesonance points in E (each MeV from 2 to 106)
  static const G4double Pmin=THmin+(nL-1)*dP;   // minE for the HighE part
  static const G4double Pmax=300000.;           // maxE for the HighE part
  static const G4int    nH=100;                 // A#of HResonance points in lnE
  static const G4double milP=std::log(Pmin);    // Low logarithm energy for the HighE part
  static const G4double malP=std::log(Pmax);    // High logarithm energy
  static const G4double dlP=(malP-milP)/(nH-1); // Step in log energy in the HighE part
  static const G4double lP=std::exp(dlP);       // Multiplication factor in the HighE part
  // If the cross section approximation formula is changed - replace from file.
  if(pZ<1 || pN<0 || tZ<1 || tN<0)
  {
    G4cout<<"-W-G4QIonIonCS::GetFunct:pZ="<<pZ<<",pN="<<pN<<",tZ="<<tZ<<",tN="<<tN<<G4endl;
    return -1;
  }
  G4int pA=pN+pZ;
  G4double dPp=dP*pA;
  G4double Mom=THmin;
  for(G4int k=0; k<nL; k++)
  {
    std::pair<G4double,G4double> len = CalculateXS(pZ, pN, tZ, tN, Mom);
    li[k]=len.first;
    le[k]=len.second;
    Mom+=dPp;
  }
  G4double lMom=Pmin*pA;
  for(G4int j=0; j<nH; j++)
  {
    std::pair<G4double,G4double> len = CalculateXS(pZ, pN, pZ, pN, lMom);
    hi[j]=len.first;
    he[j]=len.second;
    lMom*=lP;
  }
#ifdef debug
  G4cout<<"G4QIonIonCS::GetFunctions: pZ="<<pZ<<", tZ="<<tZ<<" pair is calculated"<<G4endl;
#endif
  return 1;
}

// Momentum (Mom=p/A) is in MeV/c, first=InelasticXS, second=ElasticXS (mb)
std::pair<G4double,G4double> G4QIonIonCrossSection::CalculateXS(G4int pZ,G4int pN,G4int tZ,
                                                                G4int tN, G4double Mom)
{
  static G4VQCrossSection* PElCSman = G4QProtonElasticCrossSection::GetPointer();
  static G4VQCrossSection* NElCSman = G4QNeutronElasticCrossSection::GetPointer();
  static G4VQCrossSection* InelPCSman = G4QProtonNuclearCrossSection::GetPointer();
  static G4VQCrossSection* InelNCSman = G4QNeutronNuclearCrossSection::GetPointer();
  G4double pA=pZ+pN;
  G4double tA=tZ+tN;
  if(pA<.9 || tA<.9 ||pA>239. || tA>239 || Mom < 0.) return std::make_pair(0.,0.);
  G4double inCS=0.;
  G4double elCS=0.;
  if(pA<1.1 )               // nucleon-ion interaction use NA(in,el)
  {
    if     (pZ == 1 && !pN) // proton-nuclear
    {
      inCS=InelPCSman->GetCrossSection(true, Mom, tZ, tN, 2212);
      elCS=PElCSman->GetCrossSection(true, Mom, tZ, tN, 2212);
    }
    else if(pN == 1 && !pZ) // neutron-nuclear
    {
      inCS=InelNCSman->GetCrossSection(true, Mom, tZ, tN, 2112);
      elCS=NElCSman->GetCrossSection(true, Mom, tZ, tN, 2112);
    }
    else G4cerr<<"-Warn-G4QIICS::CaCS:pZ="<<pZ<<",pN="<<pN<<",tZ="<<tZ<<",tN="<<tN<<G4endl;
  }
  else
  {
    G4double T=ThresholdMomentum(pZ, pN, tZ, tN); // @@ Can be cashed as lastTH (?)
    if(Mom<=T)
    {
      elCS=0.;
      inCS=0.;
    }
    else
    {
      G4double P2=Mom*Mom;
      G4double T2=T*T;
      G4double R=1.-T2/P2;                        // @@ Very rough threshold effect
      //G4double P4=P2*P2;
      //G4double P8=P4*P4;
      //G4double T4=T2*T2;
      //G4double tot=CalculateTotal(pA, tA, Mom)*P8/(P8+T4*T4); // @@ convert to IndepUnits
      G4double tot=R*CalculateTotal(pA, tA, Mom); // @@ convert to IndepUnits
      G4double rat=CalculateElTot(pA, tA, Mom);
      elCS=tot*rat;
      inCS=tot-elCS;
    }
  }
  return std::make_pair(inCS,elCS);
}

// Total Ion-ion cross-section (mb), Momentum (Mom) here is p/A (MeV/c=IU) (No Threshold)
G4double G4QIonIonCrossSection::CalculateTotal(G4double pA, G4double tA, G4double Mom)
{
  G4double y=std::log(Mom/1000.); // Log of momentum in GeV/c
  G4double ab=pA+tA;
  G4double al=std::log(ab);
  G4double ap=std::log(pA*tA);
  G4double e=std::pow(pA,0.1)+std::pow(tA,0.1);
  G4double d=e-1.55/std::pow(al,0.2);
  G4double f=4.;
  if(pA>4. && tA>4.) f=3.3+130./ab/ab+2.25/e;
  G4double c=(385.+11.*ab/(1.+.02*ab*al)+(.5*ab+500./al/al)/al)*std::pow(d,f);
  G4double r=y-3.-4./ap/ap;
#ifdef pdebug
  G4cout<<"G4QIonIonCS::CalcTot:P="<<Mom<<", stot="<<c+d*r*r<<", d="<<d<<", r="<<r<<G4endl;
#endif
  return c+d*r*r;
}

// Ratio elastic/Total, Momentum (Mom) here is p/A (MeV/c=IU)
G4double G4QIonIonCrossSection::CalculateElTot(G4double pA, G4double tA, G4double Mom)
{
  G4double y=std::log(Mom/1000.); // Log of momentum in GeV/c
  G4double ab=pA*tA;
  G4double ap=std::log(ab);
  G4double r=y-3.92-1.73/ap/ap;
  G4double d=.00166/(1.+.002*std::pow(ab,1.33333));
  G4double al=std::log(pA+tA);
  G4double e=1.+.235*(std::fabs(ap-5.78));
  if     (std::fabs(pA-2.)<.1 && std::fabs(tA-2.)<.1) e=2.18;
  else if(std::fabs(pA-2.)<.1 && std::fabs(tA-3.)<.1) e=1.90;
  else if(std::fabs(pA-3.)<.1 && std::fabs(tA-2.)<.1) e=1.90;
  else if(std::fabs(pA-2.)<.1 && std::fabs(tA-4.)<.1) e=1.65;
  else if(std::fabs(pA-4.)<.1 && std::fabs(tA-2.)<.1) e=1.65;
  else if(std::fabs(pA-3.)<.1 && std::fabs(tA-4.)<.1) e=1.32;
  else if(std::fabs(pA-4.)<.1 && std::fabs(tA-3.)<.1) e=1.32;
  else if(std::fabs(pA-4.)<.1 && std::fabs(tA-4.)<.1) e=1.;
  G4double f=.37+.0188*al;
  G4double g=std::log(std::pow(pA,0.35)+std::pow(tA,0.35));
  G4double h=g*g;
  G4double c=f/(1.+e/h/h);
#ifdef pdebug
  G4cout<<"G4QIonIonCS::CalcElT:P="<<Mom<<",el/tot="<<c+d*r*r<<",d="<<d<<", r="<<r<<G4endl;
#endif
  return c+d*r*r;
}

// Electromagnetic momentum/A-threshold (in MeV/c) 
G4double G4QIonIonCrossSection::ThresholdMomentum(G4int pZ, G4int pN, G4int tZ, G4int tN)
{
  static const G4double third=1./3.;
  static const G4double pM = G4QPDGCode(2212).GetMass(); // Proton mass in MeV
  static const G4double tpM= pM+pM;       // Doubled proton mass (MeV)
  if(pZ<.99 || pN<0. || tZ<.99 || tN<0.) return 0.;
  G4double tA=tZ+tN;
  G4double pA=pZ+pN;
  //G4double dE=1.263*tZ/(1.+std::pow(tA,third));
  G4double dE=pZ*tZ/(std::pow(pA,third)+std::pow(tA,third))/pA; // dE/pA (per projNucleon)
  return std::sqrt(dE*(tpM+dE));
}
