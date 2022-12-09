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
// G4 Physics class: G4ChipsProtonInelasticXS for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-03
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Feb-04
//
//
// ****************************************************************************************
// Short description: Cross-sections extracted (by W.Pokorski) from the CHIPS package for 
// proton-nuclear  interactions. Original author: M. Kossov
// -------------------------------------------------------------------------------------
//


#include "G4ChipsProtonInelasticXS.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"


// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4ChipsProtonInelasticXS);

G4ChipsProtonInelasticXS::G4ChipsProtonInelasticXS():G4VCrossSectionDataSet(Default_Name())
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

G4ChipsProtonInelasticXS::~G4ChipsProtonInelasticXS()
{
  std::size_t lens=LEN->size();
  for(std::size_t i=0; i<lens; ++i) delete[] (*LEN)[i];
  delete LEN;
  std::size_t hens=HEN->size();
  for(std::size_t i=0; i<hens; ++i) delete[] (*HEN)[i];
  delete HEN;
}

void
G4ChipsProtonInelasticXS::CrossSectionDescription(std::ostream& outFile) const
{
    outFile << "G4ChipsProtonInelasticXS provides the inelastic cross\n"
            << "section for proton nucleus scattering as a function of incident\n"
            << "momentum. The cross section is calculated using M. Kossov's\n"
            << "CHIPS parameterization of cross section data.\n";
}

G4bool G4ChipsProtonInelasticXS::IsIsoApplicable(const G4DynamicParticle*, G4int, G4int,    
				 const G4Element*,
				 const G4Material*)
{
  return true;
}


// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4ChipsProtonInelasticXS::GetIsoCrossSection(const G4DynamicParticle* Pt, G4int tgZ, G4int A,  
                                                      const G4Isotope*,
                                                      const G4Element*,
                                                      const G4Material*)
{
  G4double pMom=Pt->GetTotalMomentum();
  G4int tgN = A - tgZ;

  return GetChipsCrossSection(pMom, tgZ, tgN, 2212);
}

G4double G4ChipsProtonInelasticXS::GetChipsCrossSection(G4double pMom, G4int tgZ, G4int tgN, G4int)
{

  G4bool in=false;                     // By default the isotope must be found in the AMDB
  if(tgN!=lastN || tgZ!=lastZ)         // The nucleus was not the last used isotope
  {
    in      = false;                   // By default the isotope haven't been found in AMDB
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
        lastCS=CalculateCrossSection(-1,j,2212,lastZ,lastN,pMom); // read & update
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
      lastCS=CalculateCrossSection(0,j,2212,lastZ,lastN,pMom); //calculate & create
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
    lastCS=CalculateCrossSection(1,j,2212,lastZ,lastN,pMom); // Only read and UpdateDB
    lastP=pMom;
  }
  return lastCS*millibarn;
}

// The main member function giving the gamma-A cross section (E in GeV, CS in mb)
G4double G4ChipsProtonInelasticXS::CalculateCrossSection(G4int F, G4int I,
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
      G4int sync=(G4int)LEN->size();
      if(sync<=I) G4cout<<"*!*G4QProtonNuclCS::CalcCrossSect:Sync="<<sync<<"<="<<I<<G4endl;
      lastLEN=(*LEN)[I];               // Pointer to prepared LowEnergy cross sections
      lastHEN=(*HEN)[I];               // Pointer to prepared High Energy cross sections
    }
    else                               // This isotope wasn't calculated before => CREATE
    {
      lastLEN = new G4double[nL];      // Allocate memory for the new LEN cross sections
      lastHEN = new G4double[nH];      // Allocate memory for the new HEN cross sections
      // --- Instead of making a separate function ---
      G4double P=THmiG;                // Table threshold in GeV/c
      for(G4int k=0; k<nL; ++k)
      {
        lastLEN[k] = CrossSectionLin(targZ, targN, P);
        P+=dPG;
      }
      G4double lP=milPG;
      for(G4int n=0; n<nH; ++n)
      {
        lastHEN[n] = CrossSectionLog(targZ, targN, lP);
        lP+=dlP;
      }
      // --- End of possible separate function
      // *** The synchronization check ***
      G4int sync=(G4int)LEN->size();
      if(sync!=I)
      {
        G4cout<<"***G4ChipsProtonNuclCS::CalcCrossSect: Sinc="<<sync<<"#"<<I<<", Z=" <<targZ
              <<", N="<<targN<<", F="<<F<<G4endl;
        //G4Exception("G4ProtonNuclearCS::CalculateCS:","39",FatalException,"overflow DB");
      }
      LEN->push_back(lastLEN);          // remember the Low Energy Table
      HEN->push_back(lastHEN);          // remember the High Energy Table
    } // End of creation of the new set of parameters
  } // End of parameters udate
  // =------------------= NOW the Magic Formula =-----------------------=
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

// Electromagnetic momentum-threshold (in MeV/c) 
G4double G4ChipsProtonInelasticXS::ThresholdMomentum(G4int tZ, G4int tN)
{
  static const G4double third=1./3.;
  static const G4double pM = G4Proton::Proton()->Definition()->GetPDGMass(); // Projectile mass in MeV
  static const G4double tpM= pM+pM;       // Doubled projectile mass (MeV)

  G4double tA=tZ+tN;
  if(tZ<.99 || tN<0.) return 0.;
  else if(tZ==1 && tN==0) return 800.;    // A threshold on the free proton
  //G4double dE=1.263*tZ/(1.+G4Pow::GetInstance()->powA(tA,third));
  G4double dE=tZ/(1.+G4Pow::GetInstance()->powA(tA,third)); // Safety for diffused edge of the nucleus (QE)
  G4double tM=931.5*tA;
  G4double T=dE+dE*(dE/2+pM)/tM;
  return std::sqrt(T*(tpM+T));
}

// Calculation formula for proton-nuclear inelastic cross-section (mb) (P in GeV/c)
G4double G4ChipsProtonInelasticXS::CrossSectionLin(G4int tZ, G4int tN, G4double P)
{
  G4double sigma=0.;
  if(P<ThresholdMomentum(tZ,tN)*.001) return sigma;
  G4double lP=G4Log(P);
  if(tZ==1&&!tN){if(P>.35) sigma=CrossSectionFormula(tZ,tN,P,lP);}// s(pp)=0 below 350Mev/c
  else if(tZ<97 && tN<152)                // General solution
  {
    G4double pex=0.;
    G4double pos=0.;
    G4double wid=1.;
    if(tZ==13 && tN==14)                  // Excited metastable states
    {
      pex=230.;
      pos=.13;
      wid=8.e-5;
    }
    else if(tZ<7)
    {
      if(tZ==6 && tN==6)
      {
        pex=320.;
        pos=.14;
        wid=7.e-6;
      }
      else if(tZ==5 && tN==6)
      {
        pex=270.;
        pos=.17;
        wid=.002;
      }
      else if(tZ==4 && tN==5)
      {
        pex=600.;
        pos=.132;
        wid=.005;
      }
      else if(tZ==3 && tN==4)
      {
        pex=280.;
        pos=.19;
        wid=.0025;
      }
      else if(tZ==3 && tN==3)
      {
        pex=370.;
        pos=.171;
        wid=.006;
      }
      else if(tZ==2 && tN==1)
      {
        pex=30.;
        pos=.22;
        wid=.0005;
      }
    }
    sigma=CrossSectionFormula(tZ,tN,P,lP);
    if(pex>0.)
    {
      G4double dp=P-pos;
      sigma+=pex*G4Exp(-dp*dp/wid);
    }
  }
  else
  {
    G4cerr<<"-Warning-G4ChipsProtonNuclearXS::CSLin:*Bad A* Z="<<tZ<<", N="<<tN<<G4endl;
    sigma=0.;
  }
  if(sigma<0.) return 0.;
  return sigma;  
}

// Calculation formula for proton-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4ChipsProtonInelasticXS::CrossSectionLog(G4int tZ, G4int tN, G4double lP)
{
  G4double P=G4Exp(lP);
  return CrossSectionFormula(tZ, tN, P, lP);
}
// Calculation formula for proton-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4ChipsProtonInelasticXS::CrossSectionFormula(G4int tZ, G4int tN,
                                                           G4double P, G4double lP)
{
  G4double sigma=0.;
  if(tZ==1 && !tN)                        // pp interaction (from G4QuasiElasticRatios)
  {
    G4double El(0.),To(0.);               // Uzhi
    if(P<0.1)                             // Copied from G4QuasiElasticRatios Uzhi / start
    {
      G4double p2=P*P;
      El=1./(0.00012+p2*0.2);
      To=El;
    }
    else if(P>1000.)
    {
      G4double lp=G4Log(P)-3.5;
      G4double lp2=lp*lp;
      El=0.0557*lp2+6.72;
      To=0.3*lp2+38.2;
    }
    else
    {
      G4double p2=P*P;
      G4double LE=1./(0.00012+p2*0.2);
      G4double lp=G4Log(P)-3.5;
      G4double lp2=lp*lp;
      G4double rp2=1./p2;
      El=LE+(0.0557*lp2+6.72+32.6/P)/(1.+rp2/P);
      To=LE+(0.3   *lp2+38.2+52.7*rp2)/(1.+2.72*rp2*rp2);
    }                                   // Copied from G4QuasiElasticRatios Uzhi / end

/*                                                      // Uzhi 4.03.2013
    G4double p2=P*P;
    G4double lp=lP-3.5;
    G4double lp2=lp*lp;
    G4double rp2=1./p2;
    G4double El=(.0557*lp2+6.72+30./P)/(1.+.49*rp2/P);
    G4double To=(.3*lp2+38.2)/(1.+.54*rp2*rp2);
*/                                                      // Uzhi 4.03.2013

    sigma=To-El;
  }
  else if(tZ<97 && tN<152)                // General solution
  {
    //G4double lP=G4Log(P);            // Already calculated
    G4double d=lP-4.2;
    G4double p2=P*P;
    G4double p4=p2*p2;
    G4double a=tN+tZ;                       // A of the target
    G4double al=G4Log(a);
    G4double sa=std::sqrt(a);
    G4double a2=a*a;
    G4double a2s=a2*sa;
    G4double a4=a2*a2;
    G4double a8=a4*a4;
    G4double a12=a8*a4;
    G4double a16=a8*a8;
    G4double c=(170.+3600./a2s)/(1.+65./a2s);
    G4double dl=al-3.;
    G4double dl2=dl*dl;
    G4double r=.21+.62*dl2/(1.+.5*dl2);
    G4double gg=40.*G4Exp(al*0.712)/(1.+12.2/a)/(1.+34./a2);
    G4double e=318.+a4/(1.+.0015*a4/G4Exp(al*0.09))/(1.+4.e-28*a12)+
               8.e-18/(1./a16+1.3e-20)/(1.+1.e-21*a12);
    G4double ss=3.57+.009*a2/(1.+.0001*a2*a);
    G4double h=(.01/a4+2.5e-6/a)*(1.+6.e-6*a2*a)/(1.+6.e7/a12/a2);
    sigma=(c+d*d)/(1.+r/p4)+(gg+e*G4Exp(-ss*P))/(1.+h/p4/p4);
  }
  else
  {
    G4cerr<<"-Warning-G4QProtonNuclearCroSect::CSForm:*Bad A* Z="<<tZ<<", N="<<tN<<G4endl;
    sigma=0.;
  }
  if(sigma<0.) return 0.;
  return sigma;  
}

G4double G4ChipsProtonInelasticXS::EquLinearFit(G4double X, G4int N, G4double X0, G4double DX, G4double* Y)
{
  if(DX<=0. || N<2)
    {
      G4cerr<<"***G4ChipsProtonInelasticXS::EquLinearFit: DX="<<DX<<", N="<<N<<G4endl;
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
