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
// G4 Physics class: G4QKaonMinusNuclearCrossSection for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-03
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Feb-04
// --------------------------------------------------------------------------------
// ****************************************************************************************
// This Header is a part of the CHIPS physics package (author: M. Kosov)
// ****************************************************************************************
// Short description: CHIPS cross-sections for kaon(minus)-nuclear interactions
// -------------------------------------------------------------------------------------
//
//#define debug
//#define pdebug
//#define debug3
//#define debugn
//#define debugs

#include "G4QKaonMinusNuclearCrossSection.hh"

// Initialization of the
G4double* G4QKaonMinusNuclearCrossSection::lastLEN=0; // Pointer to lastArray of LowEn CS
G4double* G4QKaonMinusNuclearCrossSection::lastHEN=0; // Pointer to lastArray of HighEn CS
G4int     G4QKaonMinusNuclearCrossSection::lastN=0;   // The last N of calculated nucleus
G4int     G4QKaonMinusNuclearCrossSection::lastZ=0;   // The last Z of calculated nucleus
G4double  G4QKaonMinusNuclearCrossSection::lastP=0.;  // Last used in CrossSection Momentum
G4double  G4QKaonMinusNuclearCrossSection::lastTH=0.; // Last threshold momentum
G4double  G4QKaonMinusNuclearCrossSection::lastCS=0.; // Last value of the Cross Section
G4int     G4QKaonMinusNuclearCrossSection::lastI=0;   // The last position in the DAMDB
std::vector<G4double*>* G4QKaonMinusNuclearCrossSection::LEN = new std::vector<G4double*>;
std::vector<G4double*>* G4QKaonMinusNuclearCrossSection::HEN = new std::vector<G4double*>;

// Returns Pointer to the G4VQCrossSection class
G4VQCrossSection* G4QKaonMinusNuclearCrossSection::GetPointer()
{
  static G4QKaonMinusNuclearCrossSection theCrossSection;//**Static body of Cross Section**
  return &theCrossSection;
}

G4QKaonMinusNuclearCrossSection::~G4QKaonMinusNuclearCrossSection()
{
  G4int lens=LEN->size();
  for(G4int i=0; i<lens; ++i) delete[] (*LEN)[i];
  delete LEN;
  G4int hens=HEN->size();
  for(G4int i=0; i<hens; ++i) delete[] (*HEN)[i];
  delete HEN;
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4QKaonMinusNuclearCrossSection::GetCrossSection(G4bool fCS, G4double pMom,
                                                       G4int tgZ, G4int tgN, G4int PDG)
{
  static G4double tolerance=0.001;     // Tolerance (0.1%) to consider as "the same mom"
  static G4int j;                      // A#0f Z/N-records already tested in AMDB
  static std::vector <G4int>    colN;  // Vector of N for calculated nuclei (isotops)
  static std::vector <G4int>    colZ;  // Vector of Z for calculated nuclei (isotops)
  static std::vector <G4double> colP;  // Vector of last momenta for the reaction
  static std::vector <G4double> colTH; // Vector of energy thresholds for the reaction
  static std::vector <G4double> colCS; // Vector of last cross sections for the reaction
  // ***---*** End of the mandatory Static Definitions of the Associative Memory ***---***
#ifdef debug
  G4cout<<"G4QKmCS::GetCS:>>> f="<<fCS<<", p="<<pMom<<", Z="<<tgZ<<"("<<lastZ<<") ,N="<<tgN
        <<"("<<lastN<<"),PDG="<<PDG<<", thresh="<<lastTH<<",Sz="<<colN.size()<<G4endl;
#endif
  if(PDG!=-321) G4cout<<"-Warning-G4QPiMinusCS::GetCS:**Not a KMinus**,PDG="<<PDG<<G4endl;
  G4bool in=false;                     // By default the isotope must be found in the AMDB
  if(tgN!=lastN || tgZ!=lastZ)         // The nucleus was not the last used isotope
  {
    in = false;                        // By default the isotope haven't be found in AMDB  
    lastP   = 0.;                      // New momentum history (nothing to compare with)
    lastN   = tgN;                     // The last N of the calculated nucleus
    lastZ   = tgZ;                     // The last Z of the calculated nucleus
    lastI   = colN.size();             // Size of the Associative Memory DB in the heap
    j  = 0;                            // A#0f records found in DB for this projectile
#ifdef debug
    G4cout<<"G4QKmCS::GetCS: the amount of records in the AMDB lastI="<<lastI<<G4endl;
#endif
    if(lastI) for(G4int i=0; i<lastI; i++) // AMDB exists, try to find the (Z,N) isotope
    {
      if(colN[i]==tgN && colZ[i]==tgZ) // Try the record "i" in the AMDB
      {
        lastI=i;                       // Remember the index for future fast/last use
        lastTH =colTH[i];              // The last THreshold (A-dependent)
#ifdef debug
        G4cout<<"G4QKmCS::GetCS:*Found* P="<<pMom<<",Threshold="<<lastTH<<",j="<<j<<G4endl;
#endif
        if(pMom<=lastTH)
        {
#ifdef debug
          G4cout<<"G4QPCS::GetCS:Found,P="<<pMom<<" < Threshold="<<lastTH<<",CS=0"<<G4endl;
#endif
          return 0.;                   // Energy is below the Threshold value
        }
        lastP  =colP [i];              // Last Momentum  (A-dependent)
        lastCS =colCS[i];              // Last CrossSect (A-dependent)
        if(std::fabs(lastP-pMom)<tolerance*pMom)
        //if(lastP==pMom)              // VI do not use tolerance
        {
#ifdef debug
          G4cout<<"..G4QKmCS::GetCS:.DoNothing.P="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
#endif
          //CalculateCrossSection(fCS,-1,j,-321,lastZ,lastN,pMom); // Update param's only
          return lastCS*millibarn;     // Use theLastCS
        }
        in = true;                     // This is the case when the isotop is found in DB
        // Momentum pMom is in IU ! @@ Units
#ifdef debug
        G4cout<<"G4QKmCS::G:UpdatDB P="<<pMom<<",f="<<fCS<<",lI="<<lastI<<",j="<<j<<G4endl;
#endif
        lastCS=CalculateCrossSection(fCS,-1,j,-321,lastZ,lastN,pMom); // read & update
#ifdef debug
        G4cout<<"G4QKmCS::GetCrosSec: *****> New (inDB) Calculated CS="<<lastCS<<G4endl;
#endif
        if(lastCS<=0. && pMom>lastTH)  // Correct the threshold (@@ No intermediate Zeros)
        {
#ifdef debug
          G4cout<<"G4QKmCS::GetCS: New P="<<pMom<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastCS=0.;
          lastTH=pMom;
        }
        break;                         // Go out of the LOOP
      }
#ifdef debug
      G4cout<<"-->G4QKmCrossSec::GetCrosSec: pPDG=-321, j="<<j<<", N="<<colN[i]
            <<",Z["<<i<<"]="<<colZ[i]<<G4endl;
#endif
      j++;                             // Increment a#0f records found in DB
    }
#ifdef debug
    G4cout<<"-?-G4QKmCS::GetCS:RC Z="<<tgZ<<",N="<<tgN<<",in="<<in<<",j="<<j<<" ?"<<G4endl;
#endif
    if(!in)                            // This isotope has not been calculated previously
    {
#ifdef debug
      G4cout<<"^^^G4QKmCS::GetCS:CalcNew P="<<pMom<<", f="<<fCS<<", lastI="<<lastI<<G4endl;
#endif
      //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
      lastCS=CalculateCrossSection(fCS,0,j,-321,lastZ,lastN,pMom); //calculate & create
      //if(lastCS>0.)                   // It means that the AMBD was initialized
      //{

        lastTH = ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last
#ifdef debug
        G4cout<<"G4QKmCrossSection::GetCrossSect: NewThresh="<<lastTH<<",P="<<pMom<<G4endl;
#endif
        colN.push_back(tgN);
        colZ.push_back(tgZ);
        colP.push_back(pMom);
        colTH.push_back(lastTH);
        colCS.push_back(lastCS);
#ifdef debug
        G4cout<<"G4QKmCS::GetCrosSec:recCS="<<lastCS<<",lZ="<<lastN<<",lN="<<lastZ<<G4endl;
#endif
      //} // M.K. Presence of H1 with high threshold breaks the syncronization
#ifdef pdebug
      G4cout<<"G4QKmCS::GetCS:1st,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
#endif
      return lastCS*millibarn;
    } // End of creation of the new set of parameters
    else
    {
#ifdef debug
      G4cout<<"G4QKmCS::GetCS: Update lastI="<<lastI<<",j="<<j<<G4endl;
#endif
      colP[lastI]=pMom;
      colCS[lastI]=lastCS;
    }
  } // End of parameters udate
  else if(pMom<=lastTH)
  {
#ifdef debug
    G4cout<<"G4QKmCS::GetCS: Current P="<<pMom<<" < Threshold="<<lastTH<<", CS=0"<<G4endl;
#endif
    return 0.;                         // Momentum is below the Threshold Value -> CS=0
  }
  else if(std::fabs(lastP-pMom)<tolerance*pMom)
  //else if(lastP==pMom)               // VI do not use tolerance
  {
#ifdef debug
    G4cout<<"..G4QPCS::GetCS:OldNZ&P="<<lastP<<"="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
#endif
    return lastCS*millibarn;           // Use theLastCS
  }
  else                                 // It is the last used -> use the current tables
  {
#ifdef debug
    G4cout<<"-!-G4QPCS::GetCS:UseCur P="<<pMom<<",f="<<fCS<<",I="<<lastI<<",j="<<j<<G4endl;
#endif
    lastCS=CalculateCrossSection(fCS,1,j,-321,lastZ,lastN,pMom); // Only read and UpdateDB
    lastP=pMom;
  }
#ifdef debug
  G4cout<<"==>G4QKmCS::GetCroSec: P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
#endif
  return lastCS*millibarn;
}

// The main member function giving the gamma-A cross section (E in GeV, CS in mb)
G4double G4QKaonMinusNuclearCrossSection::CalculateCrossSection(G4bool, G4int F, G4int I,
                                        G4int, G4int targZ, G4int targN, G4double Momentum)
{
  static const G4double THmin=27.;     // default minimum Momentum (MeV/c) Threshold
  static const G4double THmiG=THmin*.001; // minimum Momentum (GeV/c) Threshold
  static const G4double dP=10.;        // step for the LEN (Low ENergy) table MeV/c
  static const G4double dPG=dP*.001;   // step for the LEN (Low ENergy) table GeV/c
  static const G4int    nL=105;        // A#of LEN points in E (step 10 MeV/c)
  static const G4double Pmin=THmin+(nL-1)*dP; // minP for the HighE part with safety
  static const G4double Pmax=227000.;  // maxP for the HEN (High ENergy) part 227 GeV
  static const G4int    nH=224;        // A#of HEN points in lnE
  static const G4double milP=std::log(Pmin);// Low logarithm energy for the HEN part
  static const G4double malP=std::log(Pmax);// High logarithm energy (each 2.75 percent)
  static const G4double dlP=(malP-milP)/(nH-1); // Step in log energy in the HEN part
  static const G4double milPG=std::log(.001*Pmin);// Low logarithmEnergy for HEN part GeV/c
#ifdef debug
  G4cout<<"G4QKmNCS::CalCS:N="<<targN<<",Z="<<targZ<<",P="<<Momentum<<">"<<THmin<<G4endl;
#endif
  G4double sigma=0.;
  if(F&&I) sigma=0.;                   // @@ *!* Fake line *!* to use F & I !!!Temporary!!!
  //G4double A=targN+targZ;              // A of the target
#ifdef debug
  G4cout<<"G4QKmNucCS::CalCS: A="<<A<<",F="<<F<<",I="<<I<<",nL="<<nL<<",nH="<<nH<<G4endl;
#endif
  if(F<=0)                             // This isotope was not the last used isotop
  {
    if(F<0)                            // This isotope was found in DAMDB =======> RETRIEVE
    {
      G4int sync=LEN->size();
      if(sync<=I) G4cerr<<"*!*G4QPiMinusNuclCS::CalcCrosSect:Sync="<<sync<<"<="<<I<<G4endl;
      lastLEN=(*LEN)[I];                // Pointer to prepared LowEnergy cross sections
      lastHEN=(*HEN)[I];                // Pointer to prepared High Energy cross sections
    }
    else                               // This isotope wasn't calculated before => CREATE
    {
      lastLEN = new G4double[nL];      // Allocate memory for the new LEN cross sections
      lastHEN = new G4double[nH];      // Allocate memory for the new HEN cross sections
      // --- Instead of making a separate function ---
      G4double P=THmiG;                // Table threshold in GeV/c
      for(G4int m=0; m<nL; m++)
      {
        lastLEN[m] = CrossSectionLin(targZ, targN, P);
        P+=dPG;
      }
      G4double lP=milPG;
      for(G4int n=0; n<nH; n++)
      {
        lastHEN[n] = CrossSectionLog(targZ, targN, lP);
        lP+=dlP;
      }
#ifdef debug
      G4cout<<"-*->G4QKmNucCS::CalcCS:Tab for Z="<<targZ<<",N="<<targN<<",I="<<I<<G4endl;
#endif
      // --- End of possible separate function
      // *** The synchronization check ***
      G4int sync=LEN->size();
      if(sync!=I)
      {
        G4cerr<<"***G4QPiMinusNuclCS::CalcCrossSect: Sinc="<<sync<<"#"<<I<<", Z=" <<targZ
              <<", N="<<targN<<", F="<<F<<G4endl;
        //G4Exception("G4PiMinusNuclearCS::CalculateCS:","39",FatalException,"DBoverflow");
      }
      LEN->push_back(lastLEN);         // remember the Low Energy Table
      HEN->push_back(lastHEN);         // remember the High Energy Table
    } // End of creation of the new set of parameters
  } // End of parameters udate
  // ============================== NOW the Magic Formula =================================
#ifdef debug
  G4cout<<"G4QKmNCS::CalcCS:lTH="<<lastTH<<",Pmi="<<Pmin<<",dP="<<dP<<",dlP="<<dlP<<G4endl;
#endif
  if (Momentum<lastTH) return 0.;      // It must be already checked in the interface class
  else if (Momentum<Pmin)              // High Energy region
  {
#ifdef debug
    G4cout<<"G4QKmNCS::CalcCS:bLEN nL="<<nL<<",TH="<<THmin<<",dP="<<dP<<G4endl;
#endif
    sigma=EquLinearFit(Momentum,nL,THmin,dP,lastLEN);
#ifdef debugn
    if(sigma<0.)
      G4cout<<"G4QKmNuCS::CalcCS: E="<<Momentum<<",T="<<THmin<<",dP="<<dP<<G4endl;
#endif
  }
  else if (Momentum<Pmax)              // High Energy region
  {
    G4double lP=std::log(Momentum);
#ifdef debug
    G4cout<<"G4QKmNucCS::CalcCS: before HEN nH="<<nH<<",iE="<<milP<<",dlP="<<dlP<<G4endl;
#endif
    sigma=EquLinearFit(lP,nH,milP,dlP,lastHEN);
  }
  else                                 // UHE region (calculation, not frequent)
  {
    G4double P=0.001*Momentum;         // Approximation formula is for P in GeV/c
    sigma=CrossSectionFormula(targZ, targN, P, std::log(P));
  }
#ifdef debug
  G4cout<<"G4QKaonMinusNuclearCrossSection::CalcCS: CS="<<sigma<<G4endl;
#endif
  if(sigma<0.) return 0.;
  return sigma;
}

// Calculation formula for piMinus-nuclear inelastic cross-section (mb) (P in GeV/c)
G4double G4QKaonMinusNuclearCrossSection::CrossSectionLin(G4int tZ, G4int tN, G4double P)
{
  G4double lP=std::log(P);
  return CrossSectionFormula(tZ, tN, P, lP);
}

// Calculation formula for piMinus-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4QKaonMinusNuclearCrossSection::CrossSectionLog(G4int tZ, G4int tN, G4double lP)
{
  G4double P=std::exp(lP);
  return CrossSectionFormula(tZ, tN, P, lP);
}
// Calculation formula for piMinus-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4QKaonMinusNuclearCrossSection::CrossSectionFormula(G4int tZ, G4int tN,
                                                              G4double P, G4double lP)
{
  G4double sigma=0.;
  if(tZ==1 && !tN)                        // PiMin-Proton interaction from G4QuasiElRatios
  {
    G4double ld=lP-3.5;
    G4double ld2=ld*ld;
    G4double p2=P*P;
    G4double p4=p2*p2;
    G4double sp=std::sqrt(P);
    G4double psp=P*sp;
    G4double lm=P-.39;
    G4double md=lm*lm+.000156;
    G4double lh=P-1.;
    G4double hd=lh*lh+.0156;
    G4double El=(.0557*ld2+2.23)/(1.-.7/sp+.075/p4);
    G4double To=(.3*ld2+19.5)/(1.-.21/sp+.52/p4);
    sigma=8.8/psp+(To-El)+.002/md+.15/hd;
  }
  else if(tZ==1 && tN==1)                  // kmp_tot
  {
    G4double p2=P*P;
    G4double dX=lP-3.7;
    G4double dR=P-.94;
    G4double sp=std::sqrt(P);
    sigma=(.6*dX*dX+36.)/(1.-.11/sp+.52/p2/p2)+.7/(dR*dR+.0256)+18./P/sp;
  }
  else if(tZ<97 && tN<152)                // General solution
  {
    G4double d=lP-4.2;
    G4double sp=std::sqrt(P);
    G4double p2=P*P;
    G4double a=tN+tZ;                       // A of the target
    G4double sa=std::sqrt(a);
    G4double al=std::log(a);
    G4double a2=a*a;
    G4double c=52.*std::exp(al*0.6)*(1.+97./a2)/(1.+9.8/a)/(1.+47./a2);
    G4double g=-.2-.003*a;
    G4double h=.5+.07*a;
    G4double v=P-1.;
    G4double f=.6*a*sa/(1.+.00002*a2);
    G4double u=.125+.127*al;
    sigma=(c+d*d)/(1.+g/sp+h/p2/p2)+f/(v*v+u*u)+20.*sa/P/sp;
#ifdef pdebug
    G4cout<<"G4QKMinNucCS::CSForm: A="<<a<<",P="<<P<<",CS="<<sigma<<",c="<<c<<",g="<<g
          <<",d="<<d<<",r="<<r<<",e="<<e<<",h="<<h<<G4endl;
#endif
  }
  else
  {
    G4cerr<<"-Warning-G4QKMinusNuclearCroSect::CSForm:*Bad A* Z="<<tZ<<", N="<<tN<<G4endl;
    sigma=0.;
  }
  if(sigma<0.) return 0.;
  return sigma;  
}
