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
// ****************************************************************************************
// Short description: Cross-sections extracted (by W.Pokorski) from the CHIPS package for 
// Hyperon-nuclear  interactions. Original author: M. Kossov
// -------------------------------------------------------------------------------------
//

#include "G4ChipsHyperonInelasticXS.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Lambda.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaMinus.hh"
#include "G4SigmaZero.hh"
#include "G4XiMinus.hh"
#include "G4XiZero.hh"
#include "G4OmegaMinus.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4ChipsHyperonInelasticXS);

G4ChipsHyperonInelasticXS::G4ChipsHyperonInelasticXS():G4VCrossSectionDataSet(Default_Name())
{
  // Initialization of the
  lastLEN=0; // Pointer to the lastArray of LowEn CS
  lastHEN=0; // Pointer to the lastArray of HighEn CS
  lastN=0;   // The last N of calculated nucleus
  lastZ=0;   // The last Z of calculated nucleus
  lastP=0.;  // Last used in cross section Momentum
  lastTH=0.; // Last threshold momentum
  lastCS=0.; // Last value of the Cross Section
  lastI=0;   // The last position in the DAMDB
  LEN = new std::vector<G4double*>;
  HEN = new std::vector<G4double*>;
}

G4ChipsHyperonInelasticXS::~G4ChipsHyperonInelasticXS()
{
    G4int lens=LEN->size();
    for(G4int i=0; i<lens; ++i) delete[] (*LEN)[i];
    delete LEN;

    G4int hens=HEN->size();
    for(G4int i=0; i<hens; ++i) delete[] (*HEN)[i];
    delete HEN;
}

void G4ChipsHyperonInelasticXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4ChipsHyperonInelasticXS provides the inelastic cross\n"
          << "section for hyperon nucleus scattering as a function of incident\n"
          << "momentum. The cross section is calculated using M. Kossov's\n"
          << "CHIPS parameterization of cross section data.\n";
}

G4bool G4ChipsHyperonInelasticXS::IsIsoApplicable(const G4DynamicParticle*, G4int, G4int,    
				 const G4Element*,
				 const G4Material*)
{
  /*
  const G4ParticleDefinition* particle = Pt->GetDefinition();
  if (particle == G4Lambda::Lambda()) 
    {
      return true;
    }
  else if(particle == G4SigmaPlus::SigmaPlus())
    {
    return true;
    }
  else if(particle == G4SigmaMinus::SigmaMinus())
    {
    return true;
    }
  else if(particle == G4SigmaZero::SigmaZero())
    {
      return true;
    }
  else if(particle == G4XiMinus::XiMinus())
    {
      return true;
    }
  else if(particle == G4XiZero::XiZero())
    {
      return true;
    }
  else if(particle == G4OmegaMinus::OmegaMinus())
    {
      return true;
    }
  */
  return true;
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4ChipsHyperonInelasticXS::GetIsoCrossSection(const G4DynamicParticle* Pt, G4int tgZ, G4int A,  
							 const G4Isotope*,
							 const G4Element*,
							 const G4Material*)
{
  G4double pMom=Pt->GetTotalMomentum();
  G4int tgN = A - tgZ;
  G4int pdg = Pt->GetDefinition()->GetPDGEncoding();
  
  return GetChipsCrossSection(pMom, tgZ, tgN, pdg);
}

G4double G4ChipsHyperonInelasticXS::GetChipsCrossSection(G4double pMom, G4int tgZ, G4int tgN, G4int PDG)
{
 
  G4bool in=false;                     // By default the isotope must be found in the AMDB
  if(tgN!=lastN || tgZ!=lastZ)         // The nucleus was not the last used isotope
  {
    in = false;                        // By default the isotope haven't be found in AMDB  
    lastP   = 0.;                      // New momentum history (nothing to compare with)
    lastN   = tgN;                     // The last N of the calculated nucleus
    lastZ   = tgZ;                     // The last Z of the calculated nucleus
    lastI   = colN.size();             // Size of the Associative Memory DB in the heap
    j  = 0;                            // A#0f records found in DB for this projectile

    if(lastI) for(G4int i=0; i<lastI; i++) // AMDB exists, try to find the (Z,N) isotope
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
        lastCS=CalculateCrossSection(-1,j,PDG,lastZ,lastN,pMom); // read & update

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
      lastCS=CalculateCrossSection(0,j,PDG,lastZ,lastN,pMom); //calculate & create
      //if(lastCS>0.)                   // It means that the AMBD was initialized
      //{

      lastTH = 0; //ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last
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
    lastCS=CalculateCrossSection(1,j,PDG,lastZ,lastN,pMom); // Only read and UpdateDB
    lastP=pMom;
  }
  return lastCS*millibarn;
}

// The main member function giving the gamma-A cross section (E in GeV, CS in mb)
G4double G4ChipsHyperonInelasticXS::CalculateCrossSection(G4int F, G4int I,
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
  static const G4double milP=G4Log(Pmin);// Low logarithm energy for the HEN part
  static const G4double malP=G4Log(Pmax);// High logarithm energy (each 2.75 percent)
  static const G4double dlP=(malP-milP)/(nH-1); // Step in log energy in the HEN part
  static const G4double milPG=G4Log(.001*Pmin);// Low logarithmEnergy for HEN part GeV/c
  G4double sigma=0.;
  if(F&&I) sigma=0.;                   // @@ *!* Fake line *!* to use F & I !!!Temporary!!!
  //G4double A=targN+targZ;              // A of the target
  if(F<=0)                             // This isotope was not the last used isotop
  {
    if(F<0)                            // This isotope was found in DAMDB =-----=> RETRIEVE
    {
      G4int sync=LEN->size();
      if(sync<=I) G4cerr<<"*!*G4QPiMinusNuclCS::CalcCrosSect:Sync="<<sync<<"<="<<I<<G4endl;
      lastLEN=(*LEN)[I];               // Pointer to prepared LowEnergy cross sections
      lastHEN=(*HEN)[I];               // Pointer to prepared High Energy cross sections
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
      G4int sync=LEN->size();
      if(sync!=I)
      {
        G4cerr<<"***G4QHyperNuclCS::CalcCrossSect: Sinc="<<sync<<"#"<<I<<", Z=" <<targZ
              <<", N="<<targN<<", F="<<F<<G4endl;
        //G4Exception("G4PiMinusNuclearCS::CalculateCS:","39",FatalException,"DBoverflow");
      }
      LEN->push_back(lastLEN);         // remember the Low Energy Table
      HEN->push_back(lastHEN);         // remember the High Energy Table
    } // End of creation of the new set of parameters
  } // End of parameters udate
  // =--------------------------= NOW the Magic Formula =------------------------------=
  if (Momentum<lastTH) return 0.;      // It must be already checked in the interface class
  else if (Momentum<Pmin)              // High Energy region
  {
    sigma=EquLinearFit(Momentum,nL,THmin,dP,lastLEN);
  }
  else if (Momentum<Pmax)              // High Energy region
  {
    G4double lP=G4Log(Momentum);
    sigma=EquLinearFit(lP,nH,milP,dlP,lastHEN);
  }
  else                                 // UHE region (calculation, not frequent)
  {
    G4double P=0.001*Momentum;         // Approximation formula is for P in GeV/c
    sigma=CrossSectionFormula(targZ, targN, P, G4Log(P));
  }
  if(sigma<0.) return 0.;
  return sigma;
}

// Calculation formula for piMinus-nuclear inelastic cross-section (mb) (P in GeV/c)
G4double G4ChipsHyperonInelasticXS::CrossSectionLin(G4int tZ, G4int tN, G4double P)
{
  G4double lP=G4Log(P);
  return CrossSectionFormula(tZ, tN, P, lP);
}

// Calculation formula for piMinus-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4ChipsHyperonInelasticXS::CrossSectionLog(G4int tZ, G4int tN, G4double lP)
{
  G4double P=G4Exp(lP);
  return CrossSectionFormula(tZ, tN, P, lP);
}
// Calculation formula for piMinus-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4ChipsHyperonInelasticXS::CrossSectionFormula(G4int tZ, G4int tN,
                                                              G4double P, G4double lP)
{
  G4double sigma=0.;

  //AR-24Apr2018 Switch to allow transuranic elements
  const G4bool isHeavyElementAllowed = true;

  if(tZ==1 && !tN)                        // Hyperon-P interaction from G4QuasiElastRatios
  {
    G4double ld=lP-3.5;
    G4double ld2=ld*ld;
    G4double p2=P*P;
    G4double p4=p2*p2;
    G4double sp=std::sqrt(P);
    G4double El=(.0557*ld2+6.72+99./p2)/(1.+2./sp+2./p4);
    G4double To=(.3*ld2+38.2+900./sp)/(1.+27./sp+3./p4);
    sigma=To-El;
  }
  else if((tZ<97 && tN<152) || isHeavyElementAllowed)                 // General solution
  {
    G4double d=lP-4.2;
    G4double p2=P*P;
    G4double p4=p2*p2;
    G4double sp=std::sqrt(P);
    G4double ssp=std::sqrt(sp);
    G4double a=tN+tZ;                      // A of the target
    G4double al=G4Log(a);
    G4double sa=std::sqrt(a);
    G4double a2=a*a;
    G4double a2s=a2*sa;
    G4double a4=a2*a2;
    G4double a8=a4*a4;
    G4double c=(170.+3600./a2s)/(1.+65./a2s);
    G4double gg=42.*(G4Exp(al*0.8)+4.E-8*a4)/(1.+28./a)/(1.+5.E-5*a2);
    G4double e=390.;                       // Defolt values for deutrons
    G4double r=0.27;
    G4double h=2.E-7;
    G4double t=0.3;
    if(tZ>1 || tN>1)
    {
      e=380.+18.*a2/(1.+a2/60.)/(1.+2.E-19*a8);
      r=0.15;
      h=1.E-8*a2/(1.+a2/17.)/(1.+3.E-20*a8);
      t=(.2+.00056*a2)/(1.+a2*.0006);
    }
    sigma=(c+d*d)/(1.+t/ssp+r/p4)+(gg+e*G4Exp(-6.*P))/(1.+h/p4/p4);
#ifdef pdebug
    G4cout<<"G4QHyperonNucCS::CSForm: A="<<a<<",P="<<P<<",CS="<<sigma<<",c="<<c<<",g="<<gg
          <<",d="<<d<<",r="<<r<<",e="<<e<<",h="<<h<<G4endl;
#endif
  }
  else
  {
    G4cerr<<"-Warning-G4QHyperonNuclearCroSect::CSForm:*Bad A* Z="<<tZ<<", N="<<tN<<G4endl;
    sigma=0.;
  }
  if(sigma<0.) return 0.;
  return sigma;  
}

G4double G4ChipsHyperonInelasticXS::EquLinearFit(G4double X, G4int N, G4double X0, G4double DX, G4double* Y)
{
  if(DX<=0. || N<2)
    {
      G4cerr<<"***G4ChipsHyperonInelasticXS::EquLinearFit: DX="<<DX<<", N="<<N<<G4endl;
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
