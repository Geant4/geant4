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
// The lust update: M.V. Kossov, CERN/ITEP(Moscow) 17-May-09
//
//
// G4 Physics class: G4ChipsNeutronInelasticXS for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 17-May-2009
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 21-May-04
//
// ****************************************************************************************
// Short description: Cross-sections extracted (by W.Pokorski) from the CHIPS package for 
// neutron-nuclear  interactions. Original author: M. Kossov
// -------------------------------------------------------------------------------------
//

#include "G4ChipsNeutronInelasticXS.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Neutron.hh"
#include "G4Log.hh"
#include "G4Exp.hh"


// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4ChipsNeutronInelasticXS);

// Initialization of the

G4ChipsNeutronInelasticXS::G4ChipsNeutronInelasticXS():G4VCrossSectionDataSet(Default_Name())
{
  lastLEN=0; // Pointer to the lastArray of LowEn CS
  lastHEN=0; // Pointer to the lastArray of HighEnCS
  lastN=0;   // The last N of calculated nucleus
  lastZ=0;   // The last Z of calculated nucleus
  lastP=0.;  // Last used in cross section Momentum
  lastTH=0.; // Last threshold momentum
  lastCS=0.; // Last value of the Cross Section
  lastI=0;   // The last position in the DAMDB
  HEthresh=0.;// HE threshold for the CS calculation
  LEN = new std::vector<G4double*>;
  HEN = new std::vector<G4double*>;
}

G4ChipsNeutronInelasticXS::~G4ChipsNeutronInelasticXS()
{
  std::size_t lens=LEN->size();
  for(std::size_t i=0; i<lens; ++i) delete[] (*LEN)[i];
  delete LEN;
  std::size_t hens=HEN->size();
  for(std::size_t i=0; i<hens; ++i) delete[] (*HEN)[i];
  delete HEN;
}

void
G4ChipsNeutronInelasticXS::CrossSectionDescription(std::ostream& outFile) const
{
    outFile << "G4ChipsNeutronInelasticXS provides the inelastic cross\n"
            << "section for neutron nucleus scattering as a function of incident\n"
            << "momentum. The cross section is calculated using M. Kossov's\n"
            << "CHIPS parameterization of cross section data.\n";
}

G4bool G4ChipsNeutronInelasticXS::IsIsoApplicable(const G4DynamicParticle*, G4int, G4int,    
				 const G4Element*,
				 const G4Material*)
{
  return true;
}


G4double G4ChipsNeutronInelasticXS::GetIsoCrossSection(const G4DynamicParticle* Pt, G4int tgZ, G4int A,  
						       const G4Isotope*,
						       const G4Element*,
						       const G4Material*)
{
  G4double pMom=Pt->GetTotalMomentum();
  G4int tgN = A - tgZ;

  return GetChipsCrossSection(pMom, tgZ, tgN, 2112);
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4ChipsNeutronInelasticXS::GetChipsCrossSection(G4double pMom, G4int tgZ, G4int tgN, G4int)
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
        lastCS=CalculateCrossSection(-1,j,2112,lastZ,lastN,pMom); // read & update
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
      lastCS=CalculateCrossSection(0,j,2112,lastZ,lastN,pMom); //calculate & create
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
    lastCS=CalculateCrossSection(1,j,2112,lastZ,lastN,pMom); // Only read and UpdateDB
    lastP=pMom;
  }
  return lastCS*millibarn;
}

// The main member function giving the gamma-A cross section (E in GeV, CS in mb)
G4double G4ChipsNeutronInelasticXS::CalculateCrossSection(G4int F, G4int I,
                                        G4int, G4int targZ, G4int targN, G4double Momentum)
{
  static const G4double THmin=1.;      // default minimum Momentum (MeV/c) Threshold
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
  //
  // Associative memory for acceleration
  //static std::vector <G4double>  spA;  // shadowing coefficients (A-dependent)
  G4double sigma=0.;
  if(F&&I) sigma=0.;                   // @@ *!* Fake line *!* to use F & I !!!Temporary!!!
  //G4double A=targN+targZ;              // A of the target
  if(F<=0)                             // This isotope was not the last used isotop
  {
    if(F<0)                            // This isotope was found in DAMDB =-----=> RETRIEVE
    {
      G4int sync=(G4int)LEN->size();
      if(sync<=I) G4cerr<<"*!*G4ChipsNetronNuclCS::CalcCrossSect:Sync="<<sync<<"<="<<I<<G4endl;
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
      G4int sync=(G4int)LEN->size();
      if(sync!=I)
      {
        G4cerr<<"***G4ChipsNetronNuclearCS::CalcCrossSect: Sync="<<sync<<"#"<<I<<", Z=" <<targZ
              <<", N="<<targN<<", F="<<F<<G4endl;
        //G4Exception("G4ProtonNuclearCS::CalculateCS:","39",FatalException,"overflow DB");
      }
      LEN->push_back(lastLEN);          // remember the Low Energy Table
      HEN->push_back(lastHEN);          // remember the High Energy Table
    } // End of creation of the new set of parameters
  } // End of parameters udate
  // =------------------= NOW the Magic Formula =---------------------------=
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

// Calculation formula for proton-nuclear inelastic cross-section (mb) (P in GeV/c)
G4double G4ChipsNeutronInelasticXS::CrossSectionLin(G4int tZ, G4int tN, G4double P)
{
  //==> n (Z=0)
  static const G4int N0=1;
  static const G4double pZ0N1[4]={1., 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z0N1(1,pZ0N1);
  static const std::pair<G4int, const G4double*> Z0[N0]={Z0N1};
  //==> H (Z=1) *** no protons, which are treated separately ***
  static const G4int N1=2;
  static const G4double pZ1N1[4]={6.E-8, 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z1N1(1,pZ1N1);
  static const G4double pZ1N2[4]={9.E-8, 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z1N2(2,pZ1N2);
  static const std::pair<G4int, const G4double*> Z1[N1]={Z1N1, Z1N2};
  //==> He(Z=2)
  static const G4int N2=2;
  static const G4double pZ2N1[4]={1.E-13, 9000., 1.E-4, 2.E-4};
  static const std::pair<G4int, const G4double*> Z2N1(1,pZ2N1);
  static const G4double pZ2N2[4]={7.E-4, 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z2N2(2,pZ2N2);
  static const std::pair<G4int, const G4double*> Z2[N2]={Z2N1, Z2N2};
  //==> Li(Z=3)
  static const G4int N3=2;
  static const G4double pZ3N3[4]={1.E-9, 3200., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z3N1(3,pZ3N3);
  static const G4double pZ3N4[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z3N2(4,pZ3N4);
  static const std::pair<G4int, const G4double*> Z3[N3]={Z3N1, Z3N2};
  //==> Be(Z=4)
  static const G4int N4=1;
  static const G4double pZ4N5[4]={9.E-9, 400., .088, 4.E-4};
  static const std::pair<G4int, const G4double*> Z4N5(5,pZ4N5);
  static const std::pair<G4int, const G4double*> Z4[N4]={Z4N5};
  //==> B (Z=5)
  static const G4int N5=2;
  static const G4double pZ5N5[4]={2.E-10, 2700., .009, 4.E-4};
  static const std::pair<G4int, const G4double*> Z5N5(5,pZ5N5);
  static const G4double pZ5N6[4]={2.E-8, 110., .030, 1.E-4};
  static const std::pair<G4int, const G4double*> Z5N6(6,pZ5N6);
  static const std::pair<G4int, const G4double*> Z5[N5]={Z5N5, Z5N6};
  //==> C (Z=6)
  static const G4int N6=2;
  static const G4double pZ6N6[4]={1.5E-7, 300., .129, 5.E-4}; // *** Only Nat Mix ***
  static const std::pair<G4int, const G4double*> Z6N6(6,pZ6N6);
  static const G4double pZ6N7[4]={1.5E-7, 300., .129, 5.E-4}; // *** Only Nat Mix ***
  static const std::pair<G4int, const G4double*> Z6N7(7,pZ6N7);
  static const std::pair<G4int, const G4double*> Z6[N6]={Z6N6, Z6N7};
  //==> N (Z=7)
  static const G4int N7=2;
  static const G4double pZ7N7[4]={5.E-8, 500., .085, 2.E-4};
  static const std::pair<G4int, const G4double*> Z7N7(7,pZ7N7);
  static const G4double pZ7N8[4]={5.E-8, 140., .15, 9.E-4};
  static const std::pair<G4int, const G4double*> Z7N8(8,pZ7N8);
  static const std::pair<G4int, const G4double*> Z7[N7]={Z7N7, Z7N8};
  //==> O (Z=8)
  static const G4int N8=3;
  static const G4double pZ8N8[4]={7.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z8N8(8,pZ8N8);
  static const G4double pZ8N9[4]={2.E-8, 170., .062, 1.E-3};
  static const std::pair<G4int, const G4double*> Z8N9(9,pZ8N9);
  static const G4double pZ8N10[4]={1.E-9, 0., .051, 2.5E-4}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z8N10(10,pZ8N10);
  static const std::pair<G4int, const G4double*> Z8[N8]={Z8N8, Z8N9, Z8N10};
  //==> F (Z=9)
  static const G4int N9=1;
  static const G4double pZ9N10[4]={1.E-11, 3000., .026, 3.E-5};
  static const std::pair<G4int, const G4double*> Z9N10(10,pZ9N10);
  static const std::pair<G4int, const G4double*> Z9[N9]={Z9N10};
  //==> Ne(Z=10)
  static const G4int N10=3;
  static const G4double pZ10N10[4]={4.E-8, 0., .021, 1.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z10N10(10,pZ10N10);
  static const G4double pZ10N11[4]={4.E-8, 0., .021, 1.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z10N11(11,pZ10N11);
  static const G4double pZ10N12[4]={4.E-8, 0., .051, 2.5E-4}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z10N12(12,pZ10N12);
  static const std::pair<G4int, const G4double*> Z10[N10]={Z10N10, Z10N11, Z10N12};
  //==> Na(Z=11)
  static const G4int N11=1;
  static const G4double pZ11N12[4]={8.E-10, 500., .05, 3.E-4};
  static const std::pair<G4int, const G4double*> Z11N12(12,pZ11N12);
  static const std::pair<G4int, const G4double*> Z11[N11]={Z11N12};
  //==> Mg(Z=12)
  static const G4int N12=3;
  static const G4double pZ12N12[4]={2.E-9, 350., .065, 3.E-4};
  static const std::pair<G4int, const G4double*> Z12N12(12,pZ12N12);
  static const G4double pZ12N13[4]={2.E-9, 350., .068, 2.E-4};
  static const std::pair<G4int, const G4double*> Z12N13(13,pZ12N13);
  static const G4double pZ12N14[4]={2.E-9, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z12N14(14,pZ12N14);
  static const std::pair<G4int, const G4double*> Z12[N12]={Z12N12, Z12N13, Z12N14};
  //==> Al(Z=13)
  static const G4int N13=1;
  static const G4double pZ13N14[4]={9.E-9, 500., .075, 4.E-4};
  static const std::pair<G4int, const G4double*> Z13N14(14,pZ13N14);
  static const std::pair<G4int, const G4double*> Z13[N13]={Z13N14};
  //==> Si(Z=14)
  static const G4int N14=3;
  static const G4double pZ14N14[4]={4.E-9, 200., .076, 1.E-4};
  static const std::pair<G4int, const G4double*> Z14N14(14,pZ14N14);
  static const G4double pZ14N15[4]={6.E-9, 500., .073, 4.E-4};
  static const std::pair<G4int, const G4double*> Z14N15(15,pZ14N15);
  static const G4double pZ14N16[4]={4.E-9, 200., .076, 1.E-4};
  static const std::pair<G4int, const G4double*> Z14N16(16,pZ14N16);
  static const std::pair<G4int, const G4double*> Z14[N14]={Z14N14, Z14N15, Z14N16};
  //==> P (Z=15)
  static const G4int N15=1;
  static const G4double pZ15N16[4]={6.E-9, 550., .077, 2.E-4};
  static const std::pair<G4int, const G4double*> Z15N16(16,pZ15N16);
  static const std::pair<G4int, const G4double*> Z15[N15]={Z15N16};
  //==> S (Z=16)
  static const G4int N16=4;
  static const G4double pZ16N16[4]={1.5E-8, 500., .087, 5.E-4};
  static const std::pair<G4int, const G4double*> Z16N16(16,pZ16N16);
  static const G4double pZ16N17[4]={1.E-8, 300., .07, 4.E-3};
  static const std::pair<G4int, const G4double*> Z16N17(17,pZ16N17);
  static const G4double pZ16N18[4]={2.E-8, 300., .094, 3.E-4};
  static const std::pair<G4int, const G4double*> Z16N18(18,pZ16N18);
  static const G4double pZ16N20[4]={2.E-8, 200., .11, 3.E-4};
  static const std::pair<G4int, const G4double*> Z16N20(20,pZ16N20);
  static const std::pair<G4int, const G4double*> Z16[N16]={Z16N16, Z16N17, Z16N18, Z16N20};
  //==> Cl(Z=17)
  static const G4int N17=2;
  static const G4double pZ17N18[4]={3.E-9, 300., .072, 4.E-4};
  static const std::pair<G4int, const G4double*> Z17N18(18,pZ17N18);
  static const G4double pZ17N20[4]={5.E-9, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z17N20(20,pZ17N20);
  static const std::pair<G4int, const G4double*> Z17[N17]={Z17N18, Z17N20};
  //==> Ar(Z=18)
  static const G4int N18=3;
  static const G4double pZ18N18[4]={2.5E-9, 300., .074, 2.E-4};
  static const std::pair<G4int, const G4double*> Z18N18(18,pZ18N18);
  static const G4double pZ18N20[4]={2.E-8, 400., .084, 4.E-4};
  static const std::pair<G4int, const G4double*> Z18N20(20,pZ18N20);
  static const G4double pZ18N22[4]={1.E-9, 100., .065, 2.E-4};
  static const std::pair<G4int, const G4double*> Z18N22(22,pZ18N22);
  static const std::pair<G4int, const G4double*> Z18[N18]={Z18N18, Z18N20, Z18N22};
  //==> K (Z=19)
  static const G4int N19=3;
  static const G4double pZ19N20[4]={3.E-9, 4., .02, 2.E-4};
  static const std::pair<G4int, const G4double*> Z19N20(20,pZ19N20);
  static const G4double pZ19N21[4]={3.E-9, 500., .062, 7.E-4};
  static const std::pair<G4int, const G4double*> Z19N21(21,pZ19N21);
  static const G4double pZ19N22[4]={3.E-9, 400., .073, 3.E-4};
  static const std::pair<G4int, const G4double*> Z19N22(22,pZ19N22);
  static const std::pair<G4int, const G4double*> Z19[N19]={Z19N20, Z19N21, Z19N22};
  //==> Ca(Z=20)
  static const G4int N20=6;
  static const G4double pZ20N20[4]={3.E-9, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z20N20(20,pZ20N20);
  static const G4double pZ20N22[4]={2.E-9, 400., .072, 4.E-4};
  static const std::pair<G4int, const G4double*> Z20N22(22,pZ20N22);
  static const G4double pZ20N23[4]={.3E-9, 280., .042, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z20N23(23,pZ20N23);
  static const G4double pZ20N24[4]={1.E-9, 300., .062, 2.E-4};
  static const std::pair<G4int, const G4double*> Z20N24(24,pZ20N24);
  static const G4double pZ20N26[4]={1.5E-8, 400., .064, 2.E-4};
  static const std::pair<G4int, const G4double*> Z20N26(26,pZ20N26);
  static const G4double pZ20N28[4]={7.E-9, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z20N28(28,pZ20N28);
  static const std::pair<G4int, const G4double*> Z20[N20]={Z20N20, Z20N22, Z20N23,
                                                           Z20N24, Z20N26, Z20N28};
  //==> Sc(Z=21)
  static const G4int N21=1;
  static const G4double pZ21N24[4]={5.E-9, 1000., .068, 6.E-4};
  static const std::pair<G4int, const G4double*> Z21N24(24,pZ21N24);
  static const std::pair<G4int, const G4double*> Z21[N21]={Z21N24};
  //==> Ti(Z=22)
  static const G4int N22=5;
  static const G4double pZ22N24[4]={4.E-9, 900., .065, 6.E-4};
  static const std::pair<G4int, const G4double*> Z22N24(24,pZ22N24);
  static const G4double pZ22N25[4]={4.E-9, 1000., .065, 1.E-3};
  static const std::pair<G4int, const G4double*> Z22N25(25,pZ22N25);
  static const G4double pZ22N26[4]={4.E-9, 900., .066, 4.E-4};
  static const std::pair<G4int, const G4double*> Z22N26(26,pZ22N26);
  static const G4double pZ22N27[4]={4.E-9, 800., .021, 3.E-4};
  static const std::pair<G4int, const G4double*> Z22N27(27,pZ22N27);
  static const G4double pZ22N28[4]={4.E-9, 550., .067, 2.E-4};
  static const std::pair<G4int, const G4double*> Z22N28(28,pZ22N28);
  static const std::pair<G4int, const G4double*> Z22[N22]={Z22N24, Z22N25, Z22N26,
                                                         Z22N27, Z22N28};
  //==> V (Z=23)
  static const G4int N23=2;
  static const G4double pZ23N27[4]={4.E-9, 700., .065, 1.E-3}; // *** Only Nat mix ***
  static const std::pair<G4int, const G4double*> Z23N27(27,pZ23N27);
  static const G4double pZ23N28[4]={4.E-9, 700., .065, 1.E-3}; // *** Only Nat mix ***
  static const std::pair<G4int, const G4double*> Z23N28(28,pZ23N28);
  static const std::pair<G4int, const G4double*> Z23[N23]={Z23N27, Z23N28};
  //==> Cr(Z=24)
  static const G4int N24=4;
  static const G4double pZ24N26[4]={1.E-9, 750., .056, 2.E-4};
  static const std::pair<G4int, const G4double*> Z24N26(26,pZ24N26);
  static const G4double pZ24N28[4]={1.E-9, 350., .061, 1.E-4};
  static const std::pair<G4int, const G4double*> Z24N28(28,pZ24N28);
  static const G4double pZ24N29[4]={.4E-9, 650., .056, 1.5E-4};
  static const std::pair<G4int, const G4double*> Z24N29(29,pZ24N29);
  static const G4double pZ24N30[4]={1.E-9, 700., .054, 3.E-4};
  static const std::pair<G4int, const G4double*> Z24N30(30,pZ24N30);
  static const std::pair<G4int, const G4double*> Z24[N24]={Z24N26, Z24N28, Z24N29, Z24N30};
  //==> Mn(Z=25)
  static const G4int N25=1;
  static const G4double pZ25N30[4]={.3E-9, 650., .042, 3.5E-4};
  static const std::pair<G4int, const G4double*> Z25N30(30,pZ25N30);
  static const std::pair<G4int, const G4double*> Z25[N25]={Z25N30};
  //==> Fe(Z=26)
  static const G4int N26=4;
  static const G4double pZ26N28[4]={.9E-9, 200., .062, 1.E-4};
  static const std::pair<G4int, const G4double*> Z26N28(28,pZ26N28);
  static const G4double pZ26N30[4]={.9E-9, 1500., .055, 5.E-5};
  static const std::pair<G4int, const G4double*> Z26N30(30,pZ26N30);
  static const G4double pZ26N31[4]={.9E-9, 1100., .048, 9.E-4};
  static const std::pair<G4int, const G4double*> Z26N31(31,pZ26N31);
  static const G4double pZ26N32[4]={.9E-9, 500., .055, 2.E-4};
  static const std::pair<G4int, const G4double*> Z26N32(32,pZ26N32);
  static const std::pair<G4int, const G4double*> Z26[N26]={Z26N28, Z26N30, Z26N31, Z26N32};
  //==> Co(Z=27)
  static const G4int N27=1;
  static const G4double pZ27N32[4]={.2E-9, 21., .008, 3.E-6};
  static const std::pair<G4int, const G4double*> Z27N32(32,pZ27N32);
  static const std::pair<G4int, const G4double*> Z27[N27]={Z27N32};
  //==> Ni(Z=28)
  static const G4int N28=5;
  static const G4double pZ28N30[4]={.3E-9, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z28N30(30,pZ28N30);
  static const G4double pZ28N32[4]={.3E-9, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z28N32(32,pZ28N32);
  static const G4double pZ28N33[4]={.3E-9, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z28N33(33,pZ28N33);
  static const G4double pZ28N34[4]={.3E-9, 700., .0065, 2.E-6};
  static const std::pair<G4int, const G4double*> Z28N34(34,pZ28N34);
  static const G4double pZ28N36[4]={.3E-9, 75., .0107, 4.E-6};
  static const std::pair<G4int, const G4double*> Z28N36(36,pZ28N36);
  static const std::pair<G4int, const G4double*> Z28[N28]={Z28N30, Z28N32, Z28N33,
                                                         Z28N34, Z28N36};
  //==> Cu(Z=29)
  static const G4int N29=2;
  static const G4double pZ29N34[4]={.1E-9, 35., .005, 6.E-4};
  static const std::pair<G4int, const G4double*> Z29N34(34,pZ29N34);
  static const G4double pZ29N36[4]={.2E-9, 23., .01, 2.E-4};
  static const std::pair<G4int, const G4double*> Z29N36(36,pZ29N36);
  static const std::pair<G4int, const G4double*> Z29[N29]={Z29N34, Z29N36};
  //==> Zn(Z=30)
  static const G4int N30=5;
  static const G4double pZ30N34[4]={.2E-9, 140., .02, 8.E-6}; // *** only NAT mix ***
  static const std::pair<G4int, const G4double*> Z30N34(34,pZ30N34);
  static const G4double pZ30N36[4]={.2E-9, 140., .02, 8.E-6}; // *** only NAT mix ***
  static const std::pair<G4int, const G4double*> Z30N36(36,pZ30N36);
  static const G4double pZ30N37[4]={.2E-9, 140., .02, 8.E-6}; // *** only NAT mix ***
  static const std::pair<G4int, const G4double*> Z30N37(37,pZ30N37);
  static const G4double pZ30N38[4]={.2E-9, 140., .02, 8.E-6}; // *** only NAT mix ***
  static const std::pair<G4int, const G4double*> Z30N38(38,pZ30N38);
  static const G4double pZ30N40[4]={.2E-9, 140., .02, 8.E-6}; // *** only NAT mix ***
  static const std::pair<G4int, const G4double*> Z30N40(40,pZ30N40);
  static const std::pair<G4int, const G4double*> Z30[N30]={Z30N34, Z30N36, Z30N37,
                                                           Z30N38, Z30N40};
  //==> Ga(Z=31)
  static const G4int N31=2;
  static const G4double pZ31N38[4]={.3E-9, 450., .050, 3.E-4};
  static const std::pair<G4int, const G4double*> Z31N38(38,pZ31N38);
  static const G4double pZ31N40[4]={.3E-9, 600., .048, 2.E-4};
  static const std::pair<G4int, const G4double*> Z31N40(40,pZ31N40);
  static const std::pair<G4int, const G4double*> Z31[N31]={Z31N38, Z31N40};
  //==> Ge(Z=32)
  static const G4int N32=5;
  static const G4double pZ32N38[4]={.2E-9, 200., .05, 2.E-4};
  static const std::pair<G4int, const G4double*> Z32N38(38,pZ32N38);
  static const G4double pZ32N40[4]={.2E-9, 600., .05, 2.E-4};
  static const std::pair<G4int, const G4double*> Z32N40(40,pZ32N40);
  static const G4double pZ32N41[4]={1.5E-11, 600., .028, 3.E-4};
  static const std::pair<G4int, const G4double*> Z32N41(41,pZ32N41);
  static const G4double pZ32N42[4]={9.E-11, 400., .048, 3.E-4};
  static const std::pair<G4int, const G4double*> Z32N42(42,pZ32N42);
  static const G4double pZ32N44[4]={9.E-11, 400., .043, 3.E-4};
  static const std::pair<G4int, const G4double*> Z32N44(44,pZ32N44);
  static const std::pair<G4int, const G4double*> Z32[N32]={Z32N38, Z32N40, Z32N41,
                                                           Z32N42, Z32N44};
  //==> As(Z=33)
  static const G4int N33=1;
  static const G4double pZ33N42[4]={1.E-11, 1000., .032, 1.E-4};
  static const std::pair<G4int, const G4double*> Z33N42(42,pZ33N42);
  static const std::pair<G4int, const G4double*> Z33[N33]={Z33N42};
  //==> Se(Z=34)
  static const G4int N34=6;
  static const G4double pZ34N40[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N40(40,pZ34N40);
  static const G4double pZ34N42[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N42(42,pZ34N42);
  static const G4double pZ34N43[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N43(43,pZ34N43);
  static const G4double pZ34N44[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N44(44,pZ34N44);
  static const G4double pZ34N46[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N46(46,pZ34N46);
  static const G4double pZ34N48[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N48(48,pZ34N48);
  static const std::pair<G4int, const G4double*> Z34[N34]={Z34N40, Z34N42, Z34N43,
                                                           Z34N44, Z34N46, Z34N48};
  //==> Br(Z=35)
  static const G4int N35=2;
  static const G4double pZ35N44[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z35N44(44,pZ35N44);
  static const G4double pZ35N46[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z35N46(46,pZ35N46);
  static const std::pair<G4int, const G4double*> Z35[N35]={Z35N44, Z35N46};
  //==> Kr(Z=36)
  static const G4int N36=6;
  static const G4double pZ36N42[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N42(42,pZ36N42);
  static const G4double pZ36N44[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N44(44,pZ36N44);
  static const G4double pZ36N46[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N46(46,pZ36N46);
  static const G4double pZ36N47[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N47(47,pZ36N47);
  static const G4double pZ36N48[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N48(48,pZ36N48);
  static const G4double pZ36N50[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N50(50,pZ36N50);
  static const std::pair<G4int, const G4double*> Z36[N36]={Z36N42, Z36N44, Z36N46,
                                                           Z36N47, Z36N48, Z36N50};
  //==> Rb(Z=37)
  static const G4int N37=2;
  static const G4double pZ37N48[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z37N48(48,pZ37N48);
  static const G4double pZ37N50[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z37N50(50,pZ37N50);
  static const std::pair<G4int, const G4double*> Z37[N37]={Z37N48, Z37N50};
  //==> Sr(Z=38)
  static const G4int N38=4;
  static const G4double pZ38N46[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N46(46,pZ38N46);
  static const G4double pZ38N48[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N48(48,pZ38N48);
  static const G4double pZ38N49[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N49(49,pZ38N49);
  static const G4double pZ38N50[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N50(50,pZ38N50);
  static const std::pair<G4int, const G4double*> Z38[N38]={Z38N46, Z38N48, Z38N49, Z38N50};
  //==> Y (Z=39)
  static const G4int N39=1;
  static const G4double pZ39N50[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z39N50(50,pZ39N50);
  static const std::pair<G4int, const G4double*> Z39[N39]={Z39N50};
  //==> Zr(Z=40)
  static const G4int N40=5;
  static const G4double pZ40N50[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N50(50,pZ40N50);
  static const G4double pZ40N51[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N51(51,pZ40N51);
  static const G4double pZ40N52[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N52(52,pZ40N52);
  static const G4double pZ40N54[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N54(54,pZ40N54);
  static const G4double pZ40N56[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N56(56,pZ40N56);
  static const std::pair<G4int, const G4double*> Z40[N40]={Z40N50, Z40N51, Z40N52,
                                                           Z40N54, Z40N56};
  //==> Nb(Z=41)
  static const G4int N41=1;
  static const G4double pZ41N52[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z41N52(52,pZ41N52);
  static const std::pair<G4int, const G4double*> Z41[N41]={Z41N52};
  //==> Mo(Z=42)
  static const G4int N42=7;
  static const G4double pZ42N50[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N50(50,pZ42N50);
  static const G4double pZ42N52[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N52(52,pZ42N52);
  static const G4double pZ42N53[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N53(53,pZ42N53);
  static const G4double pZ42N54[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N54(54,pZ42N54);
  static const G4double pZ42N55[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N55(55,pZ42N55);
  static const G4double pZ42N56[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N56(56,pZ42N56);
  static const G4double pZ42N58[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N58(58,pZ42N58);
  static const std::pair<G4int, const G4double*> Z42[N42]={Z42N50, Z42N52, Z42N53, Z42N54,
                                                           Z42N55, Z42N56, Z42N58};
  //==> Mo(Z=43)
  static const G4int N43=1;
  static const G4double pZ43N0[4]={3.E-12, 500., .01, 2.5E-4}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z43N0(0,pZ43N0);
  static const std::pair<G4int, const G4double*> Z43[N43]={Z43N0};
  //==> Ru(Z=44)
  static const G4int N44=7;
  static const G4double pZ44N52[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N52(52,pZ44N52);
  static const G4double pZ44N54[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N54(54,pZ44N54);
  static const G4double pZ44N55[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N55(55,pZ44N55);
  static const G4double pZ44N56[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N56(56,pZ44N56);
  static const G4double pZ44N57[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N57(57,pZ44N57);
  static const G4double pZ44N58[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N58(58,pZ44N58);
  static const G4double pZ44N60[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N60(60,pZ44N60);
  static const std::pair<G4int, const G4double*> Z44[N44]={Z44N52, Z44N54, Z44N55, Z44N56,
                                                           Z44N57, Z44N58, Z44N60};
  //==> Rh(Z=45)
  static const G4int N45=1;
  static const G4double pZ45N58[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z45N58(58,pZ45N58);
  static const std::pair<G4int, const G4double*> Z45[N45]={Z45N58};
  //==> Pd(Z=46)
  static const G4int N46=6;
  static const G4double pZ46N56[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N56(56,pZ46N56);
  static const G4double pZ46N58[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N58(58,pZ46N58);
  static const G4double pZ46N59[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N59(59,pZ46N59);
  static const G4double pZ46N60[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N60(60,pZ46N60);
  static const G4double pZ46N62[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N62(62,pZ46N62);
  static const G4double pZ46N64[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N64(64,pZ46N64);
  static const std::pair<G4int, const G4double*> Z46[N46]={Z46N56, Z46N58, Z46N59,
                                                           Z46N60, Z46N62, Z46N64};
  //==> Ag(Z=47)
  static const G4int N47=2;
  static const G4double pZ47N60[4]={3.E-12, 500., .01, 2.7E-5};
  static const std::pair<G4int, const G4double*> Z47N60(60,pZ47N60);
  static const G4double pZ47N62[4]={3.E-12, 480., .01, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z47N62(62,pZ47N62);
  static const std::pair<G4int, const G4double*> Z47[N47]={Z47N60, Z47N62};
  //==> Cd(Z=48)
  static const G4int N48=8;
  static const G4double pZ48N58[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N58(58,pZ48N58);
  static const G4double pZ48N60[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N60(60,pZ48N60);
  static const G4double pZ48N62[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N62(62,pZ48N62);
  static const G4double pZ48N63[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N63(63,pZ48N63);
  static const G4double pZ48N64[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N64(64,pZ48N64);
  static const G4double pZ48N65[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N65(65,pZ48N65);
  static const G4double pZ48N66[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N66(66,pZ48N66);
  static const G4double pZ48N68[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N68(68,pZ48N68);
  static const std::pair<G4int, const G4double*> Z48[N48]={Z48N58, Z48N60, Z48N62, Z48N63,
                                                           Z48N64, Z48N65, Z48N66, Z48N68};
  //==> In(Z=49)
  static const G4int N49=2;
  static const G4double pZ49N64[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z49N64(64,pZ49N64);
  static const G4double pZ49N66[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z49N66(66,pZ49N66);
  static const std::pair<G4int, const G4double*> Z49[N49]={Z49N64, Z49N66};
  //==> Sn(Z=50)
  static const G4int N50=10;
  static const G4double pZ50N62[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N62(62,pZ50N62);
  static const G4double pZ50N64[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N64(64,pZ50N64);
  static const G4double pZ50N65[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N65(65,pZ50N65);
  static const G4double pZ50N66[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N66(66,pZ50N66);
  static const G4double pZ50N67[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N67(67,pZ50N67);
  static const G4double pZ50N68[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N68(68,pZ50N68);
  static const G4double pZ50N69[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N69(69,pZ50N69);
  static const G4double pZ50N70[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N70(70,pZ50N70);
  static const G4double pZ50N72[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N72(72,pZ50N72);
  static const G4double pZ50N74[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N74(74,pZ50N74);
  static const std::pair<G4int, const G4double*> Z50[N50]={Z50N62, Z50N64, Z50N65, Z50N66,
                                                           Z50N67, Z50N68, Z50N69, Z50N70,
                                                           Z50N72, Z50N74};
  //==> Sb(Z=51)
  static const G4int N51=2;
  static const G4double pZ51N70[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z51N70(70,pZ51N70);
  static const G4double pZ51N72[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z51N72(72,pZ51N72);
  static const std::pair<G4int, const G4double*> Z51[N51]={Z51N70, Z51N72};
  //==> Te(Z=52)
  static const G4int N52=8;
  static const G4double pZ52N68[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N68(68,pZ52N68);
  static const G4double pZ52N70[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N70(70,pZ52N70);
  static const G4double pZ52N71[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N71(71,pZ52N71);
  static const G4double pZ52N72[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N72(72,pZ52N72);
  static const G4double pZ52N73[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N73(73,pZ52N73);
  static const G4double pZ52N74[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N74(74,pZ52N74);
  static const G4double pZ52N76[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N76(76,pZ52N76);
  static const G4double pZ52N78[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N78(78,pZ52N78);
  static const std::pair<G4int, const G4double*> Z52[N52]={Z52N68, Z52N70, Z52N71, Z52N72,
                                                           Z52N73, Z52N74, Z52N76, Z52N78};
  //==> I (Z=53)
  static const G4int N53=1;
  static const G4double pZ53N74[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z53N74(74,pZ53N74);
  static const std::pair<G4int, const G4double*> Z53[N53]={Z53N74};
  //==> Xe(Z=54)
  static const G4int N54=9;
  static const G4double pZ54N70[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N70(70,pZ54N70);
  static const G4double pZ54N72[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N72(72,pZ54N72);
  static const G4double pZ54N74[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N74(74,pZ54N74);
  static const G4double pZ54N75[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N75(75,pZ54N75);
  static const G4double pZ54N76[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N76(76,pZ54N76);
  static const G4double pZ54N77[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N77(77,pZ54N77);
  static const G4double pZ54N78[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N78(78,pZ54N78);
  static const G4double pZ54N80[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N80(80,pZ54N80);
  static const G4double pZ54N82[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N82(82,pZ54N82);
  static const std::pair<G4int, const G4double*> Z54[N54]={Z54N70, Z54N72, Z54N74,
                                                           Z54N75, Z54N76, Z54N77,
                                                           Z54N78, Z54N80, Z54N82};
  //==> Cs(Z=55)
  static const G4int N55=1;
  static const G4double pZ55N78[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z55N78(78,pZ55N78);
  static const std::pair<G4int, const G4double*> Z55[N55]={Z55N78};
  //==> Ba(Z=56)
  static const G4int N56=7;
  static const G4double pZ56N74[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N74(74,pZ56N74);
  static const G4double pZ56N76[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N76(76,pZ56N76);
  static const G4double pZ56N78[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N78(78,pZ56N78);
  static const G4double pZ56N79[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N79(79,pZ56N79);
  static const G4double pZ56N80[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N80(80,pZ56N80);
  static const G4double pZ56N81[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N81(81,pZ56N81);
  static const G4double pZ56N82[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N82(82,pZ56N82);
  static const std::pair<G4int, const G4double*> Z56[N56]={Z56N74, Z56N76, Z56N78, Z56N79,
                                                           Z56N80, Z56N81, Z56N82};
  //==> La(Z=57)
  static const G4int N57=2;
  static const G4double pZ57N81[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z57N81(81,pZ57N81);
  static const G4double pZ57N82[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z57N82(82,pZ57N82);
  static const std::pair<G4int, const G4double*> Z57[N57]={Z57N81, Z57N82};
  //==> Ce(Z=58)
  static const G4int N58=4;
  static const G4double pZ58N78[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z58N78(78,pZ58N78);
  static const G4double pZ58N80[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z58N80(80,pZ58N80);
  static const G4double pZ58N82[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z58N82(82,pZ58N82);
  static const G4double pZ58N84[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z58N84(84,pZ58N84);
  static const std::pair<G4int, const G4double*> Z58[N58]={Z58N78, Z58N80, Z58N82, Z58N84};
  //==> Pr(Z=59)
  static const G4int N59=1;
  static const G4double pZ59N82[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z59N82(82,pZ59N82);
  static const std::pair<G4int, const G4double*> Z59[N59]={Z59N82};
  //==> Nd(Z=60)
  static const G4int N60=7;
  static const G4double pZ60N82[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N82(82,pZ60N82);
  static const G4double pZ60N83[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N83(83,pZ60N83);
  static const G4double pZ60N84[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N84(84,pZ60N84);
  static const G4double pZ60N85[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N85(85,pZ60N85);
  static const G4double pZ60N86[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N86(86,pZ60N86);
  static const G4double pZ60N88[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N88(88,pZ60N88);
  static const G4double pZ60N90[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N90(90,pZ60N90);
  static const std::pair<G4int, const G4double*> Z60[N60]={Z60N82, Z60N83, Z60N84, Z60N85,
                                                           Z60N86, Z60N88, Z60N90};
  //==> Mo(Z=61)
  static const G4int N61=1;
  static const G4double pZ61N0[4]={3.E-12, 500., .01, 2.5E-4}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z61N0(0,pZ61N0);
  static const std::pair<G4int, const G4double*> Z61[N61]={Z61N0};
  //==> Sm(Z=62)
  static const G4int N62=7;
  static const G4double pZ62N82[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N82(82,pZ62N82);
  static const G4double pZ62N85[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N85(85,pZ62N85);
  static const G4double pZ62N86[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N86(86,pZ62N86);
  static const G4double pZ62N87[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N87(87,pZ62N87);
  static const G4double pZ62N88[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N88(88,pZ62N88);
  static const G4double pZ62N90[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N90(90,pZ62N90);
  static const G4double pZ62N92[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N92(92,pZ62N92);
  static const std::pair<G4int, const G4double*> Z62[N62]={Z62N82, Z62N85, Z62N86, Z62N87,
                                                           Z62N88, Z62N90, Z62N92};
  //==> Eu(Z=63)
  static const G4int N63=2;
  static const G4double pZ63N88[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z63N88(88,pZ63N88);
  static const G4double pZ63N90[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z63N90(90,pZ63N90);
  static const std::pair<G4int, const G4double*> Z63[N63]={Z63N88, Z63N90};
  //==> Gd(Z=64)
  static const G4int N64=7;
  static const G4double pZ64N88[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N88(88,pZ64N88);
  static const G4double pZ64N90[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N90(90,pZ64N90);
  static const G4double pZ64N91[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N91(91,pZ64N91);
  static const G4double pZ64N92[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N92(92,pZ64N92);
  static const G4double pZ64N93[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N93(93,pZ64N93);
  static const G4double pZ64N94[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N94(94,pZ64N94);
  static const G4double pZ64N96[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N96(96,pZ64N96);
  static const std::pair<G4int, const G4double*> Z64[N64]={Z64N88, Z64N90, Z64N91, Z64N92,
                                                           Z64N93, Z64N94, Z64N96};
  //==> Tb(Z=65)
  static const G4int N65=1;
  static const G4double pZ65N94[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z65N94(82,pZ65N94);
  static const std::pair<G4int, const G4double*> Z65[N65]={Z65N94};
  //==> Dy(Z=66)
  static const G4int N66=7;
  static const G4double pZ66N90[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N90(90,pZ66N90);
  static const G4double pZ66N92[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N92(92,pZ66N92);
  static const G4double pZ66N94[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N94(94,pZ66N94);
  static const G4double pZ66N95[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N95(95,pZ66N95);
  static const G4double pZ66N96[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N96(96,pZ66N96);
  static const G4double pZ66N97[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N97(97,pZ66N97);
  static const G4double pZ66N98[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N98(98,pZ66N98);
  static const std::pair<G4int, const G4double*> Z66[N66]={Z66N90, Z66N92, Z66N94, Z66N95,
                                                           Z66N96, Z66N97, Z66N98};
  //==> Ho(Z=67)
  static const G4int N67=1;
  static const G4double pZ67N98[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z67N98(98,pZ67N98);
  static const std::pair<G4int, const G4double*> Z67[N67]={Z67N98};
  //==> Er(Z=68)
  static const G4int N68=6;
  static const G4double pZ68N94[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N94(94,pZ68N94);
  static const G4double pZ68N96[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N96(96,pZ68N96);
  static const G4double pZ68N98[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N98(98,pZ68N98);
  static const G4double pZ68N99[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N99(99,pZ68N99);
  static const G4double pZ68N100[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N100(100,pZ68N100);
  static const G4double pZ68N102[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N102(102,pZ68N102);
  static const std::pair<G4int, const G4double*> Z68[N68]={Z68N94, Z68N96, Z68N98,
                                                           Z68N99, Z68N100, Z68N102};
  //==> Tm(Z=69)
  static const G4int N69=1;
  static const G4double pZ69N100[4]={3.E-12, 500., .01, 2.5E-4}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z69N100(100,pZ69N100);
  static const std::pair<G4int, const G4double*> Z69[N69]={Z69N100};
  //==> Yb(Z=70)
  static const G4int N70=7;
  static const G4double pZ70N98[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N98(98,pZ70N98);
  static const G4double pZ70N100[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N100(100,pZ70N100);
  static const G4double pZ70N101[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N101(101,pZ70N101);
  static const G4double pZ70N102[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N102(102,pZ70N102);
  static const G4double pZ70N103[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N103(103,pZ70N103);
  static const G4double pZ70N104[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N104(104,pZ70N104);
  static const G4double pZ70N106[4]={3.E-12, 500., .01, 2.5E-4}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N106(106,pZ70N106);
  static const std::pair<G4int, const G4double*> Z70[N70]={Z70N98, Z70N100, Z70N101,
                                                           Z70N102, Z70N103, Z70N104,
                                                           Z70N106};
  //==> Lu(Z=71)
  static const G4int N71=2;
  static const G4double pZ71N104[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z71N104(104,pZ71N104);
  static const G4double pZ71N105[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z71N105(105,pZ71N105);
  static const std::pair<G4int, const G4double*> Z71[N71]={Z71N104, Z71N105};
  //==> Hf(Z=72)
  static const G4int N72=6;
  static const G4double pZ72N102[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N102(102,pZ72N102);
  static const G4double pZ72N104[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N104(104,pZ72N104);
  static const G4double pZ72N105[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N105(105,pZ72N105);
  static const G4double pZ72N106[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N106(106,pZ72N106);
  static const G4double pZ72N107[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N107(107,pZ72N107);
  static const G4double pZ72N108[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N108(108,pZ72N108);
  static const std::pair<G4int, const G4double*> Z72[N72]={Z72N102, Z72N104, Z72N105,
                                                           Z72N106, Z72N107, Z72N108};
  //==> Ta(Z=73)
  static const G4int N73=1;
  static const G4double pZ73N108[4]={4.E-12, 1100., .027, 1.E-3};
  static const std::pair<G4int, const G4double*> Z73N108(108,pZ73N108);
  static const std::pair<G4int, const G4double*> Z73[N73]={Z73N108};
  //==> W (Z=74)
  static const G4int N74=5;
  static const G4double pZ74N106[4]={7.E-12, 1000., .03, 2.E-4}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z74N106(106,pZ74N106);
  static const G4double pZ74N108[4]={7.E-12, 1300., .03, 1.5E-4};
  static const std::pair<G4int, const G4double*> Z74N108(108,pZ74N108);
  static const G4double pZ74N109[4]={2.E-12, 1700., .023, 2.E-4};
  static const std::pair<G4int, const G4double*> Z74N109(109,pZ74N109);
  static const G4double pZ74N110[4]={7.E-12, 1100., .03, 1.5E-4};
  static const std::pair<G4int, const G4double*> Z74N110(110,pZ74N110);
  static const G4double pZ74N112[4]={7.E-12, 1100., .03, 1.5E-4};
  static const std::pair<G4int, const G4double*> Z74N112(112,pZ74N112);
  static const std::pair<G4int, const G4double*> Z74[N74]={Z74N106, Z74N108, Z74N109,
                                                           Z74N110, Z74N112};
  //==> Re(Z=75)
  static const G4int N75=2;
  static const G4double pZ75N110[4]={5.E-12, 1000., .025, 3.E-4};
  static const std::pair<G4int, const G4double*> Z75N110(110,pZ75N110);
  static const G4double pZ75N112[4]={5.E-12, 1000., .025, 3.E-4};
  static const std::pair<G4int, const G4double*> Z75N112(112,pZ75N112);
  static const std::pair<G4int, const G4double*> Z75[N75]={Z75N110, Z75N112};
  //==> Os(Z=76)
  static const G4int N76=7;
  static const G4double pZ76N108[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N108(108,pZ76N108);
  static const G4double pZ76N110[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N110(110,pZ76N110);
  static const G4double pZ76N111[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N111(111,pZ76N111);
  static const G4double pZ76N112[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N112(112,pZ76N112);
  static const G4double pZ76N113[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N113(113,pZ76N113);
  static const G4double pZ76N114[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N114(114,pZ76N114);
  static const G4double pZ76N116[4]={3.E-12, 500., .01, 2.5E-4}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N116(116,pZ76N116);
  static const std::pair<G4int, const G4double*> Z76[N76]={Z76N108, Z76N110, Z76N111,
                                                           Z76N112, Z76N113, Z76N114,
                                                           Z76N116};
  //==> Ir(Z=77)
  static const G4int N77=2;
  static const G4double pZ77N114[4]={4.E-12, 1700., .028, 2.E-4};
  static const std::pair<G4int, const G4double*> Z77N114(114,pZ77N114);
  static const G4double pZ77N116[4]={5.E-12, 1500., .028, 2.E-4};
  static const std::pair<G4int, const G4double*> Z77N116(116,pZ77N116);
  static const std::pair<G4int, const G4double*> Z77[N77]={Z77N114, Z77N116};
  //==> Pt(Z=78)
  static const G4int N78=6;
  static const G4double pZ78N112[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N112(112,pZ78N112);
  static const G4double pZ78N114[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N114(114,pZ78N114);
  static const G4double pZ78N116[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N116(116,pZ78N116);
  static const G4double pZ78N117[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N117(117,pZ78N117);
  static const G4double pZ78N118[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N118(118,pZ78N118);
  static const G4double pZ78N120[4]={3.E-12, 500., .01, 2.5E-4}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N120(120,pZ78N120);
  static const std::pair<G4int, const G4double*> Z78[N78]={Z78N112, Z78N114, Z78N116,
                                                           Z78N117, Z78N118, Z78N120};
  //==> Au(Z=79)
  static const G4int N79=1;
  static const G4double pZ79N118[4]={.2E-9, 1600., .043, 5.E-4};
  static const std::pair<G4int, const G4double*> Z79N118(118,pZ79N118);
  static const std::pair<G4int, const G4double*> Z79[N79]={Z79N118};
  //==> Hg(Z=80)
  static const G4int N80=7;
  static const G4double pZ80N116[4]={6.E-8, 2500., .085, 2.E-3};
  static const std::pair<G4int, const G4double*> Z80N116(116,pZ80N116);
  static const G4double pZ80N118[4]={6.E-8, 2500., .083, 1.7E-3};
  static const std::pair<G4int, const G4double*> Z80N118(118,pZ80N118);
  static const G4double pZ80N119[4]={6.E-8, 2600., .073, 2.5E-3};
  static const std::pair<G4int, const G4double*> Z80N119(119,pZ80N119);
  static const G4double pZ80N120[4]={6.E-8, 2500., .084, 1.7E-3};
  static const std::pair<G4int, const G4double*> Z80N120(120,pZ80N120);
  static const G4double pZ80N121[4]={1.5E-7, 2600., .078, 4.E-3};
  static const std::pair<G4int, const G4double*> Z80N121(121,pZ80N121);
  static const G4double pZ80N122[4]={6.E-8, 2500., .083, 1.6E-3};
  static const std::pair<G4int, const G4double*> Z80N122(122,pZ80N122);
  static const G4double pZ80N124[4]={6.E-8, 2500., .083, 1.5E-3};
  static const std::pair<G4int, const G4double*> Z80N124(124,pZ80N124);
  static const std::pair<G4int, const G4double*> Z80[N80]={Z80N116, Z80N118, Z80N119,
                                                           Z80N120, Z80N121, Z80N122,
                                                           Z80N124};
  //==> Tl(Z=81)
  static const G4int N81=2;
  static const G4double pZ81N122[4]={3.E-12, 500., .01, 2.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z81N122(122,pZ81N122);
  static const G4double pZ81N124[4]={3.E-12, 500., .01, 2.5E-4}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z81N124(124,pZ81N124);
  static const std::pair<G4int, const G4double*> Z81[N81]={Z81N122, Z81N124};
  //==> Pb(Z=82)
  static const G4int N82=4;
  static const G4double pZ82N122[4]={.2E-9, 40., .002, 6.E-4};
  static const std::pair<G4int, const G4double*> Z82N122(122,pZ82N122);
  static const G4double pZ82N124[4]={6.E-9, 1700., .076, 7.E-4};
  static const std::pair<G4int, const G4double*> Z82N124(124,pZ82N124);
  static const G4double pZ82N125[4]={.2E-9, 770., .057, 4.5E-4};
  static const std::pair<G4int, const G4double*> Z82N125(125,pZ82N125);
  static const G4double pZ82N126[4]={4.E-9, 0., .051, 2.E-4};
  static const std::pair<G4int, const G4double*> Z82N126(126,pZ82N126);
  static const std::pair<G4int, const G4double*> Z82[N82]={Z82N122, Z82N124, Z82N125,
                                                           Z82N126};
  //==> Bi(Z=83)
  static const G4int N83=1;
  static const G4double pZ83N126[4]={1.5E-9, 150., .052, 5.E-5};
  static const std::pair<G4int, const G4double*> Z83N126(126,pZ83N126);
  static const std::pair<G4int, const G4double*> Z83[N83]={Z83N126};
  //==> Po(Z=84)
  static const G4int N84=1;
  static const G4double pZ84N0[4]={3.E-12, 500., .01, 2.5E-4}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z84N0(0,pZ84N0);
  static const std::pair<G4int, const G4double*> Z84[N84]={Z84N0};
  //==> At(Z=85)
  static const G4int N85=1;
  static const G4double pZ85N0[4]={3.E-12, 500., .01, 2.5E-4}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z85N0(0,pZ85N0);
  static const std::pair<G4int, const G4double*> Z85[N85]={Z85N0};
  //==> Rn(Z=86)
  static const G4int N86=1;
  static const G4double pZ86N0[4]={3.E-12, 500., .01, 2.5E-4}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z86N0(0,pZ86N0);
  static const std::pair<G4int, const G4double*> Z86[N86]={Z86N0};
  //==> Fr(Z=87)
  static const G4int N87=1;
  static const G4double pZ87N0[4]={3.E-12, 500., .01, 2.5E-4}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z87N0(0,pZ87N0);
  static const std::pair<G4int, const G4double*> Z87[N87]={Z87N0};
  //==> Ra(Z=88)
  static const G4int N88=1;
  static const G4double pZ88N138[4]={3.E-9, 2200., .057, 1.2E-3};
  static const std::pair<G4int, const G4double*> Z88N138(138,pZ88N138);
  static const std::pair<G4int, const G4double*> Z88[N88]={Z88N138};
  //==> Ac(Z=89)
  static const G4int N89=1;
  static const G4double pZ89N0[4]={3.E-12, 500., .01, 2.5E-4}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z89N0(0,pZ89N0);
  static const std::pair<G4int, const G4double*> Z89[N89]={Z89N0};
  //==> Th(Z=90)
  static const G4int N90=1;
  static const G4double pZ90N142[4]={1.E-11, 1200., .028, 3.E-4};
  static const std::pair<G4int, const G4double*> Z90N142(142,pZ90N142);
  static const std::pair<G4int, const G4double*> Z90[N90]={Z90N142};
  //==> Pa(Z=91)
  static const G4int N91=1;
  static const G4double pZ91N0[4]={3.E-12, 500., .01, 2.5E-4}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z91N0(0,pZ91N0);
  static const std::pair<G4int, const G4double*> Z91[N91]={Z91N0};
  //==> U (Z=92)
  static const G4int N92=2;
  static const G4double pZ92N143[4]={2.E-11, 2700., .026, 6.E-4};
  static const std::pair<G4int, const G4double*> Z92N143(143,pZ92N143);
  static const G4double pZ92N146[4]={1.E-11, 1700., .029, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z92N146(146,pZ92N146);
  static const std::pair<G4int, const G4double*> Z92[N92]={Z92N143, Z92N146};
  //==> Np(Z=93)
  static const G4int N93=1;
  static const G4double pZ93N144[4]={4.E-8, 3700., .066, 3.5E-3};
  static const std::pair<G4int, const G4double*> Z93N144(144,pZ93N144);
  static const std::pair<G4int, const G4double*> Z93[N93]={Z93N144};
  //==> Pu(Z=94)
  static const G4int N94=3;
  static const G4double pZ94N145[4]={8.E-11, 2900., .029, 1.3E-3}; // *** Artificial ***
  static const std::pair<G4int, const G4double*> Z94N145(145,pZ94N145);
  static const G4double pZ94N148[4]={9.E-12, 1400., .025, 3.E-4}; // *** Artificial ***
  static const std::pair<G4int, const G4double*> Z94N148(148,pZ94N148);
  static const G4double pZ94N150[4]={4.E-12, 1500., .023, 1.2E-4};
  static const std::pair<G4int, const G4double*> Z94N150(150,pZ94N150);
  static const std::pair<G4int, const G4double*> Z94[N94]={Z94N145, Z94N148, Z94N150};
  //==> Am(Z=95)
  static const G4int N95=1;
  static const G4double pZ95N0[4]={3.E-12, 500., .01, 2.5E-4}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z95N0(0,pZ95N0);
  static const std::pair<G4int, const G4double*> Z95[N95]={Z95N0};
  //==> Cm(Z=96)
  static const G4int N96=1;
  static const G4double pZ96N151[4]={1.5E-8, 3700., .055, 2.E-3};
  static const std::pair<G4int, const G4double*> Z96N151(151,pZ96N151);
  static const std::pair<G4int, const G4double*> Z96[N96]={Z96N151};
 
  static const G4int NZ=97; // #of Elements covered by CHIPS
  static const std::pair<G4int, const G4double*>* Pars[NZ]={Z0,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,
    Z10,Z11,Z12,Z13,Z14,Z15,Z16,Z17,Z18,Z19,Z20,Z21,Z22,Z23,Z24,Z25,Z26,Z27,Z28,Z29,Z30,
    Z31,Z32,Z33,Z34,Z35,Z36,Z37,Z38,Z39,Z40,Z41,Z42,Z43,Z44,Z45,Z46,Z47,Z48,Z49,Z50,Z51,
    Z52,Z53,Z54,Z55,Z56,Z57,Z58,Z59,Z60,Z61,Z62,Z63,Z64,Z65,Z66,Z67,Z68,Z69,Z70,Z71,Z72,
    Z73,Z74,Z75,Z76,Z77,Z78,Z79,Z80,Z81,Z82,Z83,Z84,Z85,Z86,Z87,Z88,Z89,Z90,Z91,Z92,Z93,
    Z94,Z95,Z96};
  static const G4int NIso[NZ]={N0,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16,
    N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,N36,N37,
    N38,N39,N40,N41,N42,N43,N44,N45,N46,N47,N48,N49,N50,N51,N52,N53,N54,N55,N56,N57,N58,
    N59,N60,N61,N62,N63,N64,N65,N66,N67,N68,N69,N70,N71,N72,N73,N74,N75,N76,N77,N78,N79,
    N80,N81,N82,N83,N84,N85,N86,N87,N88,N89,N90,N91,N92,N93,N94,N95,N96};
  //G4int curN=Pars[1][0].first;
  //G4double par=Pars[1][0].second[1];
  //G4cout<<"-Warning-G4ChipsNeutronInelasticXS::CSLin: N="<<curN<<", P="<<par<<G4endl;
  G4double sigma=0.;
  G4double lP=G4Log(P);
  if( (tZ==1 && !tN) || (!tZ && tN==1)){if(P>.35) sigma=CrossSectionFormula(tZ,tN,P,lP);}
  else if(tZ<97 && tN<152)                // General solution (*** Z/A limits ***)
  {
    HEthresh=1.E-4; // Default guess
    G4double pex=0.;
    G4double pos=0.;
    G4double wid=1.;
    G4int nn=NIso[tZ];
    G4bool nfound=true;
    if(nn) for (G4int in=0; in<nn; in++)
    {
      std::pair<G4int, const G4double*> curIs=Pars[tZ][in];
      if(curIs.first == tN)
      {
        const G4double* curT=curIs.second;
        HEthresh= curT[0];
        pex     = curT[1];
        pos     = curT[2];
        wid     = curT[3];
        nfound  = false;
        break;
      }
    }
    if(nfound) G4cout<<"-Warning-G4ChipsNeutronInelasticXS::CSLin: Z="<<tZ<<", N="
                     <<tN<<" isotope is not implemented in CHIPS"<<G4endl;
    sigma=CrossSectionFormula(tZ,tN,P,lP);
    if(pex>0.)
    {
      G4double dp=P-pos;
      sigma+=pex*G4Exp(-dp*dp/wid);
    }
  }
  else
  {
    G4cerr<<"-Warning-G4ChipsNeutronNuclearCroSect::CSLin:*Bad A* Z="<<tZ<<", N="<<tN<<G4endl;
    sigma=0.;
  }
  if(sigma<0.) return 0.;
  return sigma;  
}

// Calculation formula for proton-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4ChipsNeutronInelasticXS::CrossSectionLog(G4int tZ, G4int tN, G4double lP)
{
  G4double P=G4Exp(lP);
  return CrossSectionFormula(tZ, tN, P, lP);
}
// Calculation formula for proton-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4ChipsNeutronInelasticXS::CrossSectionFormula(G4int tZ, G4int tN,
                                                           G4double P, G4double lP)
{
  G4double sigma=0.;
  if(tZ==1 && !tN)                        // np interaction from G4QuasiElasticRatios
  {

    G4double El(0.), To(0.);              // Uzhi
    if(P<0.1)                             // Copied from G4QuasiElasticRatios Uzhi / start
    {
      G4double p2=P*P;
      El=1./(0.00012+p2*(0.051+0.1*p2));
      To=El;
    }
    else if(P>1000.)
    {
      G4double lp=G4Log(P)-3.5;
      G4double lp2=lp*lp;
      El=0.0557*lp2+6.72;
      To=0.3   *lp2+38.2;
    }
    else
    {
      G4double p2=P*P;
      G4double LE=1./(0.00012+p2*(0.051+0.1*p2));
      G4double lp=G4Log(P)-3.5;
      G4double lp2=lp*lp;
      G4double rp2=1./p2;
      El=LE+(0.0557*lp2+6.72+30./P)/(1.+0.49*rp2/P);
      To=LE+(0.3   *lp2+38.2)/(1.+0.54*rp2*rp2);
    }                                   // Copied from G4QuasiElasticRatios Uzhi / end

/*                                                          // Uzhi 4.03.2013
    G4double p2=P*P;
    G4double lp=lP-3.5;
    G4double lp2=lp*lp;
    G4double rp2=1./p2;
    G4double El=(.0557*lp2+6.72+32.6/P)/(1.+rp2/P);
    G4double To=(.3*lp2+38.2+52.7*rp2)/(1.+2.72*rp2*rp2);
*/                                                          // Uzhi 4.03.2013
    sigma=To-El;
  }
  else if(tZ<97 && tN<152)                // General solution
  {
    //G4double lP=G4Log(P);            // Already calculated
    G4double d=lP-4.2;        //
    G4double p2=P*P;          //
    G4double p4=p2*p2;        //
    G4double a=tN+tZ;                     // A of the target
    G4double al=G4Log(a);  //
    G4double sa=std::sqrt(a); //
    G4double a2=a*a;          //
    G4double sa2=sa*a2;       //
    G4double a3=a2*a;         //
    G4double a4=a2*a2;        //
    //G4double a5=a4*a;
    G4double a6=a4*a2;        //
    G4double a7=a6*a;         //
    G4double a8=a4*a4;        //
    //G4double a12=a8*a4;
    //G4double a16=a8*a8;
    G4double c=(170.+3600./sa2)/(1.+65./sa2);
    G4double dl=al-3.;
    G4double dl2=dl*dl;
    G4double r=.21+.62*dl2/(1.+.5*dl2);
    G4double gg=42.*(G4Exp(al*0.8)+4.E-8*a4)/(1.+28./a)/(1.+5.E-5*a2);
    G4double e=5.*((a6+.021*a8)/(1.+.0013*a7)+.001*a3)/(1.+.0007*a2);
    G4double ss=5./(1.+144./a8);
				G4double h=HEthresh; // Individual

    //G4double h=(.01/a4+2.5e-6/a)*(1.+7.e-8*a4)/(1.+6.e7/a12/a2);
    //sigma=(c+d*d)/(1.+r/p4)+(gg+e*G4Exp(-ss*P))/(1.+h/p4/p4);
    sigma=(c+d*d)/(1+r/p4)+(gg+e*G4Exp(-ss*P))/(1+h/p4/p4);
  }
  else
  {
    G4cerr<<"-Warning-G4ChipsNeutronNuclearCroSect::CSForm:*Bad A* Z="<<tZ<<", N="<<tN<<G4endl;
    sigma=0.;
  }
  if(sigma<0.) return 0.;
  return sigma;  
}

G4double G4ChipsNeutronInelasticXS::EquLinearFit(G4double X, G4int N, G4double X0, G4double DX, G4double* Y)
{
  if(DX<=0. || N<2)
    {
      G4cerr<<"***G4ChipsNeutronInelasticXS::EquLinearFit: DX="<<DX<<", N="<<N<<G4endl;
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
