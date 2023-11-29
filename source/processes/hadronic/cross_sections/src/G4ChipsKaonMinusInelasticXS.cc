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
//
//
// G4 Physics class: G4ChipsKaonMinusInelasticXS for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-03
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Feb-04
// --------------------------------------------------------------------------------
// Short description: Cross-sections extracted from the CHIPS package for 
// kaon(minus)-nuclear interactions. Author: M. Kossov
// -------------------------------------------------------------------------------------
//

#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4KaonMinus.hh"

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4ChipsKaonMinusInelasticXS);

namespace {
    const G4double THmin=27.;     // default minimum Momentum (MeV/c) Threshold
    const G4double THmiG=THmin*.001; // minimum Momentum (GeV/c) Threshold
    const G4double dP=10.;        // step for the LEN (Low ENergy) table MeV/c
    const G4double dPG=dP*.001;   // step for the LEN (Low ENergy) table GeV/c
    const G4int    nL=105;        // A#of LEN points in E (step 10 MeV/c)
    const G4double Pmin=THmin+(nL-1)*dP; // minP for the HighE part with safety
    const G4double Pmax=227000.;  // maxP for the HEN (High ENergy) part 227 GeV
    const G4int    nH=224;        // A#of HEN points in lnE
    const G4double milP=std::log(Pmin);// Low logarithm energy for the HEN part
    const G4double malP=std::log(Pmax);// High logarithm energy (each 2.75 percent)
    const G4double dlP=(malP-milP)/(nH-1); // Step in log energy in the HEN part
    const G4double milPG=std::log(.001*Pmin);// Low logarithmEnergy for HEN part GeV/c
}
// Initialization of the

G4ChipsKaonMinusInelasticXS::G4ChipsKaonMinusInelasticXS():G4VCrossSectionDataSet(Default_Name())
{
  lastLEN=0; // Pointer to lastArray of LowEn CS
  lastHEN=0; // Pointer to lastArray of HighEn CS
  lastN=0;   // The last N of calculated nucleus
  lastZ=0;   // The last Z of calculated nucleus
  lastP=0.;  // Last used in CrossSection Momentum
  lastTH=0.; // Last threshold momentum
  lastCS=0.; // Last value of the Cross Section
  lastI=0;   // The last position in the DAMDB
  LEN = new std::vector<G4double*>;
  HEN = new std::vector<G4double*>;
}

G4ChipsKaonMinusInelasticXS::~G4ChipsKaonMinusInelasticXS()
{  
  std::size_t lens=LEN->size();
  for(std::size_t i=0; i<lens; ++i) delete[] (*LEN)[i];
  delete LEN;
 
  std::size_t hens=HEN->size();
  for(std::size_t i=0; i<hens; ++i) delete[] (*HEN)[i];
  delete HEN; 
}

void
G4ChipsKaonMinusInelasticXS::CrossSectionDescription(std::ostream& outFile) const
{
    outFile << "G4ChipsKaonMinusInelasticXS provides the inelastic cross\n"
            << "section for K- nucleus scattering as a function of incident\n"
            << "momentum. The cross section is calculated using M. Kossov's\n"
            << "CHIPS parameterization of cross section data.\n";
}

G4bool G4ChipsKaonMinusInelasticXS::IsIsoApplicable(const G4DynamicParticle*, G4int, G4int,    
				 const G4Element*,
				 const G4Material*)
{
  return true;
}


// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4ChipsKaonMinusInelasticXS::GetIsoCrossSection(const G4DynamicParticle* Pt, G4int tgZ, G4int A,  
							 const G4Isotope*,
							 const G4Element*,
							 const G4Material*)
{
  G4double pMom=Pt->GetTotalMomentum();
  G4int tgN = A - tgZ;
  
  return GetChipsCrossSection(pMom, tgZ, tgN, -321);
}

G4double G4ChipsKaonMinusInelasticXS::GetChipsCrossSection(G4double pMom, G4int tgZ, G4int tgN, G4int)
{
  G4bool in=false;                     // By default the isotope must be found in the AMDB
  if(tgN!=lastN || tgZ!=lastZ)         // The nucleus was not the last used isotope
  {
    in = false;                        // By default the isotope haven't be found in AMDB  
    lastP   = 0.;                      // New momentum history (nothing to compare with)
    lastN   = tgN;                     // The last N of the calculated nucleus
    lastZ   = tgZ;                     // The last Z of the calculated nucleus
    lastI   = (G4int)colN.size();      // Size of the Associative Memory DB in the heap
    j  = 0;                            // A#0f records found in DB for this projectile
    if(lastI) for(G4int i=0; i<lastI; ++i) // AMDB exists, try to find the (Z,N) isotope
    {
      if(colN[i]==tgN && colZ[i]==tgZ) // Try the record "i" in the AMDB
      {
        lastI=i;                       // Remember the index for future fast/last use
        lastTH =colTH[i];              // The last THreshold (A-dependent)
        if(pMom<=lastTH)
        {
          return 0.;                   // Energy is below the Threshold value
        }
        lastP  =colP [i];              // Last Momentum  (A-dependent)
        lastCS =colCS[i];              // Last CrossSect (A-dependent)
        in = true;                     // This is the case when the isotop is found in DB
        // Momentum pMom is in IU ! @@ Units
        lastCS=CalculateCrossSection(-1,j,-321,lastZ,lastN,pMom); // read & update
        if(lastCS<=0. && pMom>lastTH)  // Correct the threshold (@@ No intermediate Zeros)
        {
          lastCS=0.;
          lastTH=pMom;
        }
        break;                         // Go out of the LOOP
      }
      j++;                             // Increment a#0f records found in DB
    }
    if(!in)                            // This isotope has not been calculated previously
    {
      //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
      lastCS=CalculateCrossSection(0,j,-321,lastZ,lastN,pMom); //calculate & create
      //if(lastCS>0.)                   // It means that the AMBD was initialized
      //{

      //        lastTH = ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last

      lastTH = 0; // WP - to be checked!!!
        colN.push_back(tgN);
        colZ.push_back(tgZ);
        colP.push_back(pMom);
        colTH.push_back(lastTH);
        colCS.push_back(lastCS);
      //} // M.K. Presence of H1 with high threshold breaks the syncronization
      return lastCS*millibarn;
    } // End of creation of the new set of parameters
    else
    {
      colP[lastI]=pMom;
      colCS[lastI]=lastCS;
    }
  } // End of parameters udate
  else if(pMom<=lastTH)
  {
    return 0.;                         // Momentum is below the Threshold Value -> CS=0
  }
  else                                 // It is the last used -> use the current tables
  {
    lastCS=CalculateCrossSection(1,j,-321,lastZ,lastN,pMom); // Only read and UpdateDB
    lastP=pMom;
  }
  return lastCS*millibarn;
}

// The main member function giving the gamma-A cross section (E in GeV, CS in mb)
G4double G4ChipsKaonMinusInelasticXS::CalculateCrossSection(G4int F, G4int I,
                                        G4int, G4int targZ, G4int targN, G4double Momentum)
{
  G4double sigma=0.;
  if(F&&I) sigma=0.;                   // @@ *!* Fake line *!* to use F & I !!!Temporary!!!
  //G4double A=targN+targZ;              // A of the target
  if(F<=0)                             // This isotope was not the last used isotop
  {
    if(F<0)                            // This isotope was found in DAMDB =-----=> RETRIEVE
    {
      G4int sync=(G4int)LEN->size();
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
      for(G4int k=0; k<nL; k++)
      {
        lastLEN[k] = CrossSectionLin(targZ, targN, P);
        P+=dPG;
      }
      G4double lP=milPG;
      for(G4int n=0; n<nH; n++)
      {
        lastHEN[n] = CrossSectionLog(targZ, targN, lP);
        lP+=dlP;
      }
      // --- End of possible separate function
      // *** The synchronization check ***
      G4int sync=(G4int)LEN->size();
      if(sync!=I)
      {
        G4cerr<<"***G4ChipsKaonMinusCS::CalcCrossSect: Sinc="<<sync<<"#"<<I<<", Z=" <<targZ
              <<", N="<<targN<<", F="<<F<<G4endl;
        //G4Exception("G4PiMinusNuclearCS::CalculateCS:","39",FatalException,"DBoverflow");
      }
      LEN->push_back(lastLEN);         // remember the Low Energy Table
      HEN->push_back(lastHEN);         // remember the High Energy Table
    } // End of creation of the new set of parameters
  } // End of parameters udate
  // =------------------= NOW the Magic Formula =--------------------------=
  if (Momentum<lastTH) return 0.;      // It must be already checked in the interface class
  else if (Momentum<Pmin)              // High Energy region
  {
    sigma=EquLinearFit(Momentum,nL,THmin,dP,lastLEN);
  }
  else if (Momentum<Pmax)              // High Energy region
  {
    G4double lP=std::log(Momentum);
    sigma=EquLinearFit(lP,nH,milP,dlP,lastHEN);
  }
  else                                 // UHE region (calculation, not frequent)
  {
    G4double P=0.001*Momentum;         // Approximation formula is for P in GeV/c
    sigma=CrossSectionFormula(targZ, targN, P, std::log(P));
  }
  if(sigma<0.) return 0.;
  return sigma;
}

// Calculation formula for piMinus-nuclear inelastic cross-section (mb) (P in GeV/c)
G4double G4ChipsKaonMinusInelasticXS::CrossSectionLin(G4int tZ, G4int tN, G4double P)
{
  G4double lP=std::log(P);
  return CrossSectionFormula(tZ, tN, P, lP);
}

// Calculation formula for piMinus-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4ChipsKaonMinusInelasticXS::CrossSectionLog(G4int tZ, G4int tN, G4double lP)
{
  G4double P=std::exp(lP);
  return CrossSectionFormula(tZ, tN, P, lP);
}
// Calculation formula for piMinus-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4ChipsKaonMinusInelasticXS::CrossSectionFormula(G4int tZ, G4int tN,
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
    G4double gg=-.2-.003*a;
    G4double h=.5+.07*a;
    G4double v=P-1.;
    G4double f=.6*a*sa/(1.+.00002*a2);
    G4double u=.125+.127*al;
    sigma=(c+d*d)/(1.+gg/sp+h/p2/p2)+f/(v*v+u*u)+20.*sa/P/sp;
  }
  else
  {
    G4cerr<<"-Warning-G4ChipsKMinusNuclearCroSect::CSForm:*Bad A* Z="<<tZ<<", N="<<tN<<G4endl;
    sigma=0.;
  }
  if(sigma<0.) return 0.;
  return sigma;  
}

G4double G4ChipsKaonMinusInelasticXS::EquLinearFit(G4double X, G4int N, G4double X0, G4double DX, G4double* Y)
{
  if(DX<=0. || N<2)
    {
      G4cerr<<"***G4ChipsKaonMinusInelasticXS::EquLinearFit: DX="<<DX<<", N="<<N<<G4endl;
      return Y[0];
    }
  
  G4int    N2=N-2;
  G4double d=(X-X0)/DX;
  G4int         jj=static_cast<int>(d);
  if     (jj<0)  jj=0;
  else if(jj>N2) jj=N2;
  d-=jj; // excess
  G4double yi=Y[jj];
  G4double sigma=yi+(Y[jj+1]-yi)*d;
  
  return sigma;
}
