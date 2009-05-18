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
// G4 Physics class: G4QNeutronNuclearCrossSection for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-03
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Feb-04
// --------------------------------------------------------------------------------
// ****************************************************************************************
// ***** This HEADER is a property of the CHIPS hadronic package in Geant4 (M. Kosov) *****
// *********** DO NOT MAKE ANY CHANGE without approval of Mikhail.Kossov@cern.ch **********
// ****************************************************************************************
// Short description: CHIPS cross-sections for proton-nuclear interactions
// -----------------------------------------------------------------------
//
//#define debug
//#define pdebug
//#define debug3
//#define debugn
//#define debugs

#include "G4QNeutronNuclearCrossSection.hh"

// Initialization of the
G4double* G4QNeutronNuclearCrossSection::lastLEN=0; // Pointer to the lastArray of LowEn CS
G4double* G4QNeutronNuclearCrossSection::lastHEN=0; // Pointer to the lastArray of HighEnCS
G4int     G4QNeutronNuclearCrossSection::lastN=0;   // The last N of calculated nucleus
G4int     G4QNeutronNuclearCrossSection::lastZ=0;   // The last Z of calculated nucleus
G4double  G4QNeutronNuclearCrossSection::lastP=0.;  // Last used in cross section Momentum
G4double  G4QNeutronNuclearCrossSection::lastTH=0.; // Last threshold momentum
G4double  G4QNeutronNuclearCrossSection::lastCS=0.; // Last value of the Cross Section
G4int     G4QNeutronNuclearCrossSection::lastI=0;   // The last position in the DAMDB

// Returns Pointer to the G4VQCrossSection class
G4VQCrossSection* G4QNeutronNuclearCrossSection::GetPointer()
{
  static G4QNeutronNuclearCrossSection theCrossSection; //**Static body of Cross Section**
  return &theCrossSection;
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4QNeutronNuclearCrossSection::GetCrossSection(G4bool fCS, G4double pMom,
                                                        G4int tgZ, G4int tgN, G4int)
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
  G4cout<<"G4QPrCS::GetCS:>>> f="<<fCS<<", p="<<pMom<<", Z="<<tgZ<<"("<<lastZ<<") ,N="<<tgN
        <<"("<<lastN<<"),PDG=2112, thresh="<<lastTH<<",Sz="<<colN.size()<<G4endl;
#endif
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
    G4cout<<"G4QPrCS::GetCS: the amount of records in the AMDB lastI="<<lastI<<G4endl;
#endif
    if(lastI) for(G4int i=0; i<lastI; i++) // AMDB exists, try to find the (Z,N) isotope
    {
      if(colN[i]==tgN && colZ[i]==tgZ) // Try the record "i" in the AMDB
      {
        lastI=i;                       // Remember the index for future fast/last use
        lastTH =colTH[i];              // The last THreshold (A-dependent)
#ifdef debug
        G4cout<<"G4QPrCS::GetCS:*Found* P="<<pMom<<",Threshold="<<lastTH<<",j="<<j<<G4endl;
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
          G4cout<<"..G4QPrCS::GetCS:.DoNothing.P="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
#endif
          //CalculateCrossSection(fCS,-1,j,2112,lastZ,lastN,pMom); // Update param's only
          return lastCS*millibarn;     // Use theLastCS
        }
        in = true;                     // This is the case when the isotop is found in DB
        // Momentum pMom is in IU ! @@ Units
#ifdef debug
        G4cout<<"G4QPrCS::G:UpdatDB P="<<pMom<<",f="<<fCS<<",lI="<<lastI<<",j="<<j<<G4endl;
#endif
        lastCS=CalculateCrossSection(fCS,-1,j,2112,lastZ,lastN,pMom); // read & update
#ifdef debug
        G4cout<<"G4QPrCS::GetCrosSec: *****> New (inDB) Calculated CS="<<lastCS<<G4endl;
#endif
        if(lastCS<=0. && pMom>lastTH)  // Correct the threshold (@@ No intermediate Zeros)
        {
#ifdef debug
          G4cout<<"G4QPrCS::GetCS: New P="<<pMom<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastCS=0.;
          lastTH=pMom;
        }
        break;                         // Go out of the LOOP
      }
#ifdef debug
      G4cout<<"-->G4QPrCrossSec::GetCrosSec: pPDG=2112, j="<<j<<", N="<<colN[i]
            <<",Z["<<i<<"]="<<colZ[i]<<G4endl;
#endif
      j++;                             // Increment a#0f records found in DB
    }
#ifdef debug
    G4cout<<"-?-G4QPrCS::GetCS:RC Z="<<tgZ<<",N="<<tgN<<",in="<<in<<",j="<<j<<" ?"<<G4endl;
#endif
    if(!in)                            // This isotope has not been calculated previously
    {
#ifdef debug
      G4cout<<"^^^G4QPrCS::GetCS:CalcNew P="<<pMom<<", f="<<fCS<<", lastI="<<lastI<<G4endl;
#endif
      //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
      lastCS=CalculateCrossSection(fCS,0,j,2112,lastZ,lastN,pMom); //calculate & create
      if(lastCS>0.)                   // It means that the AMBD was initialized
      {

        lastTH = ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last
#ifdef debug
        G4cout<<"G4QPrCrossSection::GetCrossSect: NewThresh="<<lastTH<<",P="<<pMom<<G4endl;
#endif
        colN.push_back(tgN);
        colZ.push_back(tgZ);
        colP.push_back(pMom);
        colTH.push_back(lastTH);
        colCS.push_back(lastCS);
#ifdef debug
        G4cout<<"G4QPrCS::GetCrosSec:recCS="<<lastCS<<",lZ="<<lastN<<",lN="<<lastZ<<G4endl;
#endif
      }
#ifdef pdebug
      G4cout<<"G4QPrCS::GetCS:1st,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
#endif
      return lastCS*millibarn;
    } // End of creation of the new set of parameters
    else
    {
#ifdef debug
      G4cout<<"G4QPrCS::GetCS: Update lastI="<<lastI<<",j="<<j<<G4endl;
#endif
      colP[lastI]=pMom;
      colCS[lastI]=lastCS;
    }
  } // End of parameters udate
  else if(pMom<=lastTH)
  {
#ifdef debug
    G4cout<<"G4QPrCS::GetCS: Current P="<<pMom<<" < Threshold="<<lastTH<<", CS=0"<<G4endl;
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
    lastCS=CalculateCrossSection(fCS,1,j,2112,lastZ,lastN,pMom); // Only read and UpdateDB
    lastP=pMom;
  }
#ifdef debug
  G4cout<<"==>G4QPrCS::GetCroSec: P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
#endif
  return lastCS*millibarn;
}

// The main member function giving the gamma-A cross section (E in GeV, CS in mb)
G4double G4QNeutronNuclearCrossSection::CalculateCrossSection(G4bool, G4int F, G4int I,
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
  //
  // Associative memory for acceleration
  //static std::vector <G4double>  spA;  // shadowing coefficients (A-dependent)
  static std::vector <G4double*> LEN;  // Vector of pointers to LowEnProtonCrossSection
  static std::vector <G4double*> HEN;  // Vector of pointers to HighEnProtonCrossSection
#ifdef debug
  G4cout<<"G4QProtNCS::CalCS:N="<<targN<<",Z="<<targZ<<",P="<<Momentum<<">"<<THmin<<G4endl;
#endif
  if (Momentum<THmin) return 0.;       // @@ This can be dangerouse for the heaviest nuc.?!
  G4double sigma=0.;
  if(F&&I) sigma=0.;                   // @@ *!* Fake line *!* to use F & I !!!Temporary!!!
  G4double A=targN+targZ;              // A of the target
#ifdef debug
  G4cout<<"G4QProtNucCS::CalCS: A="<<A<<",F="<<F<<",I="<<I<<",nL="<<nL<<",nH="<<nH<<G4endl;
#endif
  if(F<=0)                             // This isotope was not the last used isotop
  {
    if(F<0)                            // This isotope was found in DAMDB =======> RETRIEVE
    {
      G4int sync=LEN.size();
      if(sync<=I) G4cerr<<"*!*G4QProtonNuclCS::CalcCrossSect:Sync="<<sync<<"<="<<I<<G4endl;
      lastLEN=LEN[I];                  // Pointer to prepared LowEnergy cross sections
      lastHEN=HEN[I];                  // Pointer to prepared High Energy cross sections
    }
    else if(Momentum<ThresholdMomentum(targZ,targN)) return 0.; // BelowThreshold -> NotIni
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
      G4cout<<"-*->G4QPr0tNucCS::CalcCS:Tab for Z="<<targZ<<",N="<<targN<<",I="<<I<<G4endl;
#endif
      // --- End of possible separate function
      // *** The synchronization check ***
      G4int sync=LEN.size();
      if(sync!=I)
      {
        G4cerr<<"***G4QProtonNuclCS::CalcCrossSect: Sinc="<<sync<<"#"<<I<<", Z=" <<targZ
              <<", N="<<targN<<", F="<<F<<G4endl;
        //G4Exception("G4ProtonNuclearCS::CalculateCS:","39",FatalException,"overflow DB");
      }
      LEN.push_back(lastLEN);          // remember the Low Energy Table
      HEN.push_back(lastHEN);          // remember the High Energy Table
    } // End of creation of the new set of parameters
  } // End of parameters udate
  // ============================== NOW the Magic Formula =================================
#ifdef debug
  G4cout<<"G4QPrNCS::CalcCS:lTH="<<lastTH<<",Pmi="<<Pmin<<",dP="<<dP<<",dlP="<<dlP<<G4endl;
#endif
  if (Momentum<lastTH) return 0.;      // It must be already checked in the interface class
  else if (Momentum<Pmin)              // High Energy region
  {
#ifdef debug
    G4cout<<"G4QPrNCS::CalcCS:bLEN A="<<A<<", nL="<<nL<<",TH="<<THmin<<",dP="<<dP<<G4endl;
#endif
    if(A<=1.) sigma=0.;
    else      sigma=EquLinearFit(Momentum,nL,THmin,dP,lastLEN);
#ifdef debugn
    if(sigma<0.)
      G4cout<<"G4QPrNuCS::CalcCS:A="<<A<<",E="<<Momentum<<",T="<<THmin<<",dP="<<dP<<G4endl;
#endif
  }
  else if (Momentum<Pmax)              // High Energy region
  {
    G4double lP=std::log(Momentum);
#ifdef debug
    G4cout<<"G4QProtNucCS::CalcCS: before HEN nH="<<nH<<",iE="<<milP<<",dlP="<<dlP<<G4endl;
#endif
    sigma=EquLinearFit(lP,nH,milP,dlP,lastHEN);
  }
  else                                 // UHE region (calculation, not frequent)
  {
    G4double P=0.001*Momentum;         // Approximation formula is for P in GeV/c
    sigma=CrossSectionFormula(targZ, targN, P, std::log(P));
  }
#ifdef debug
  G4cout<<"G4QNeutronNuclearCrossSection::CalcCS: CS="<<sigma<<G4endl;
#endif
  if(sigma<0.) return 0.;
  return sigma;
}

// Electromagnetic momentum-threshold (in MeV/c) 
G4double G4QNeutronNuclearCrossSection::ThresholdMomentum(G4int tZ, G4int tN)
{
  static const G4double third=1./3.;
  static const G4double pM = G4QPDGCode(2112).GetMass(); // Proton mass in MeV
  static const G4double tpM= pM+pM;       // Doubled proton mass (MeV)
  G4double tA=tZ+tN;
  if(tZ<.99 || tN<0.) return 0.;
  else if(tZ==1 && tN==0) return 800.;    // A threshold on the free proton
  //G4double dE=1.263*tZ/(1.+std::pow(tA,third));
  G4double dE=tZ/(1.+std::pow(tA,third)); // Safety for diffused edge of the nucleus (QE)
  return std::sqrt(dE*(tpM+dE));
}

// Calculation formula for proton-nuclear inelastic cross-section (mb) (P in GeV/c)
G4double G4QNeutronNuclearCrossSection::CrossSectionLin(G4int tZ, G4int tN, G4double P)
{
  //==> n (Z=0)
  static const G4double pZ0N1[4]={1., 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z0N1=std::make_pair(1,pZ0N1);
  static const std::pair<G4int, const G4double*> Z0[1]={Z0N1};
  //==> H (Z=1)
  static const G4double pZ1N1[4]={6.E-8, 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z1N1=std::make_pair(1,pZ1N1);
  static const G4double pZ1N2[4]={9.E-8, 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z1N2=std::make_pair(2,pZ1N2);
  static const std::pair<G4int, const G4double*> Z1[2]={Z1N1, Z1N2};
  //==> He(Z=2)
  static const G4double pZ2N1[4]={1.E-13, 9000., 1.E-4, 2.E-4};
  static const std::pair<G4int, const G4double*> Z2N1=std::make_pair(1,pZ2N1);
  static const G4double pZ2N2[4]={7.E-4, 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z2N2=std::make_pair(2,pZ2N2);
  static const std::pair<G4int, const G4double*> Z2[2]={Z2N1, Z2N2};
  //==> Li(Z=3)
  static const G4double pZ3N3[4]={1.E-9, 3200., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z3N1=std::make_pair(3,pZ3N3);
  static const G4double pZ3N4[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z3N2=std::make_pair(4,pZ3N4);
  static const std::pair<G4int, const G4double*> Z3[2]={Z3N1, Z3N2};
  //==> Be(Z=4)
  static const G4double pZ4N5[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z4N5=std::make_pair(5,pZ4N5);
  static const std::pair<G4int, const G4double*> Z4[1]={Z4N5};
  //==> B (Z=5)
  static const G4double pZ5N5[4]={1.E-9, 3200., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z5N5=std::make_pair(5,pZ5N5);
  static const G4double pZ5N6[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z5N6=std::make_pair(6,pZ5N6);
  static const std::pair<G4int, const G4double*> Z5[2]={Z5N5, Z5N6};
  //==> C (Z=6)
  static const G4double pZ6N6[4]={1.E-9, 3200., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z6N6=std::make_pair(6,pZ6N6);
  static const G4double pZ6N7[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z6N7=std::make_pair(7,pZ6N7);
  static const std::pair<G4int, const G4double*> Z6[2]={Z6N6, Z6N7};
  //==> N (Z=7)
  static const G4double pZ7N7[4]={1.E-9, 3200., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z7N7=std::make_pair(7,pZ7N7);
  static const G4double pZ7N8[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z7N8=std::make_pair(8,pZ7N8);
  static const std::pair<G4int, const G4double*> Z7[2]={Z7N7, Z7N8};
  //==> O (Z=8)
  static const G4double pZ8N8[4]={1.E-9, 3200., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z8N8=std::make_pair(8,pZ8N8);
  static const G4double pZ8N9[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z8N9=std::make_pair(9,pZ8N9);
  static const G4double pZ8N10[4]={3.E-9, 200., .051, 2.5E-4}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z8N10=std::make_pair(10,pZ8N10);
  static const std::pair<G4int, const G4double*> Z8[3]={Z8N8, Z8N9, Z8N10};
  //==> F (Z=9)
  static const G4double pZ9N10[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z9N10=std::make_pair(10,pZ9N10);
  static const std::pair<G4int, const G4double*> Z9[1]={Z9N10};
  //==> Ne(Z=10)
  static const G4double pZ10N10[4]={4.E-8, 0., .021, 1.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z10N10=std::make_pair(10,pZ10N10);
  static const G4double pZ10N11[4]={4.E-8, 0., .021, 1.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z10N11=std::make_pair(11,pZ10N11);
  static const G4double pZ10N12[4]={4.E-8, 0., .051, 2.5E-4}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z10N12=std::make_pair(12,pZ10N12);
  static const std::pair<G4int, const G4double*> Z10[3]={Z10N10, Z10N11, Z10N12};
  //==> Na(Z=11)
  static const G4double pZ11N12[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z11N12=std::make_pair(12,pZ11N12);
  static const std::pair<G4int, const G4double*> Z11[1]={Z11N12};
  //==> Mg(Z=12)
  static const G4double pZ12N12[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z12N12=std::make_pair(12,pZ12N12);
  static const G4double pZ12N13[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z12N13=std::make_pair(13,pZ12N13);
  static const G4double pZ12N14[4]={4.E-8, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z12N14=std::make_pair(14,pZ12N14);
  static const std::pair<G4int, const G4double*> Z12[3]={Z12N12, Z12N13, Z12N14};
  //==> Al(Z=13)
  static const G4double pZ13N14[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z13N14=std::make_pair(14,pZ13N14);
  static const std::pair<G4int, const G4double*> Z13[1]={Z13N14};
  //==> Si(Z=14)
  static const G4double pZ14N14[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z14N14=std::make_pair(14,pZ14N14);
  static const G4double pZ14N15[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z14N15=std::make_pair(15,pZ14N15);
  static const G4double pZ14N16[4]={4.E-8, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z14N16=std::make_pair(16,pZ14N16);
  static const std::pair<G4int, const G4double*> Z14[3]={Z14N14, Z14N15, Z14N16};
  //==> P (Z=15)
  static const G4double pZ15N16[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z15N16=std::make_pair(16,pZ15N16);
  static const std::pair<G4int, const G4double*> Z15[1]={Z15N16};
  //==> S (Z=16)
  static const G4double pZ16N16[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z16N16=std::make_pair(16,pZ16N16);
  static const G4double pZ16N17[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z16N17=std::make_pair(17,pZ16N17);
  static const G4double pZ16N18[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z16N18=std::make_pair(18,pZ16N18);
  static const G4double pZ16N20[4]={4.E-8, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z16N20=std::make_pair(20,pZ16N20);
  static const std::pair<G4int, const G4double*> Z16[4]={Z16N16, Z16N17, Z16N18, Z16N20};
  //==> Cl(Z=17)
  static const G4double pZ17N18[4]={1.E-9, 3200., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z17N18=std::make_pair(18,pZ17N18);
  static const G4double pZ17N20[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z17N20=std::make_pair(20,pZ17N20);
  static const std::pair<G4int, const G4double*> Z17[2]={Z17N18, Z17N20};
  //==> Ar(Z=18)
  static const G4double pZ18N18[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z18N18=std::make_pair(18,pZ18N18);
  static const G4double pZ18N19[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z18N19=std::make_pair(19,pZ18N19);
  static const G4double pZ18N20[4]={4.E-8, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z18N20=std::make_pair(20,pZ18N20);
  static const std::pair<G4int, const G4double*> Z18[3]={Z18N18, Z18N19, Z18N20};
  //==> K (Z=19)
  static const G4double pZ19N20[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z19N20=std::make_pair(20,pZ19N20);
  static const G4double pZ19N21[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z19N21=std::make_pair(21,pZ19N21);
  static const G4double pZ19N22[4]={4.E-8, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z19N22=std::make_pair(22,pZ19N22);
  static const std::pair<G4int, const G4double*> Z19[3]={Z19N20, Z19N21, Z19N22};
  //==> Ca(Z=20)
  static const G4double pZ20N20[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z20N20=std::make_pair(20,pZ20N20);
  static const G4double pZ20N22[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z20N22=std::make_pair(22,pZ20N22);
  static const G4double pZ20N23[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z20N23=std::make_pair(23,pZ20N23);
  static const G4double pZ20N24[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z20N24=std::make_pair(24,pZ20N24);
  static const G4double pZ20N26[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z20N26=std::make_pair(26,pZ20N26);
  static const G4double pZ20N28[4]={4.E-8, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z20N28=std::make_pair(28,pZ20N28);
  static const std::pair<G4int, const G4double*> Z20[6]={Z20N20, Z20N22, Z20N23,
                                                         Z20N24, Z20N26, Z20N28};
  //==> Sc(Z=21)
  static const G4double pZ21N24[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z21N24=std::make_pair(24,pZ21N24);
  static const std::pair<G4int, const G4double*> Z21[1]={Z21N24};
  //==> Ti(Z=22)
  static const G4double pZ22N24[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z22N24=std::make_pair(24,pZ22N24);
  static const G4double pZ22N25[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z22N25=std::make_pair(25,pZ22N25);
  static const G4double pZ22N26[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z22N26=std::make_pair(26,pZ22N26);
  static const G4double pZ22N27[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z22N27=std::make_pair(27,pZ22N27);
  static const G4double pZ22N28[4]={4.E-8, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z22N28=std::make_pair(28,pZ22N28);
  static const std::pair<G4int, const G4double*> Z22[5]={Z22N24, Z22N25, Z22N26,
                                                         Z22N27, Z22N28};
  //==> V (Z=23)
  static const G4double pZ23N27[4]={1.E-9, 3200., .021, 1.5E-5}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z23N27=std::make_pair(27,pZ23N27);
  static const G4double pZ23N28[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z23N28=std::make_pair(28,pZ23N28);
  static const std::pair<G4int, const G4double*> Z23[2]={Z23N27, Z23N28};
  //==> Cr(Z=24)
  static const G4double pZ24N26[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z24N26=std::make_pair(26,pZ24N26);
  static const G4double pZ24N28[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z24N28=std::make_pair(28,pZ24N28);
  static const G4double pZ24N29[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z24N29=std::make_pair(29,pZ24N29);
  static const G4double pZ24N30[4]={4.E-8, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z24N30=std::make_pair(30,pZ24N30);
  static const std::pair<G4int, const G4double*> Z24[4]={Z24N26, Z24N28, Z24N29, Z24N30};
  //==> Mn(Z=25)
  static const G4double pZ25N30[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z25N30=std::make_pair(30,pZ25N30);
  static const std::pair<G4int, const G4double*> Z25[1]={Z25N30};
  //==> Fe(Z=26)
  static const G4double pZ26N28[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z26N28=std::make_pair(28,pZ26N28);
  static const G4double pZ26N30[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z26N30=std::make_pair(30,pZ26N30);
  static const G4double pZ26N31[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z26N31=std::make_pair(31,pZ26N31);
  static const G4double pZ26N32[4]={4.E-8, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z26N32=std::make_pair(32,pZ26N32);
  static const std::pair<G4int, const G4double*> Z26[4]={Z26N28, Z26N30, Z26N31, Z26N32};
  //==> Co(Z=27)
  static const G4double pZ27N32[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z27N32=std::make_pair(32,pZ27N32);
  static const std::pair<G4int, const G4double*> Z27[1]={Z27N32};
  //==> Ni(Z=28)
  static const G4double pZ28N30[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z28N30=std::make_pair(30,pZ28N30);
  static const G4double pZ28N32[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z28N32=std::make_pair(32,pZ28N32);
  static const G4double pZ28N33[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z28N33=std::make_pair(33,pZ28N33);
  static const G4double pZ28N34[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z28N34=std::make_pair(34,pZ28N34);
  static const G4double pZ28N36[4]={4.E-8, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z28N36=std::make_pair(36,pZ28N36);
  static const std::pair<G4int, const G4double*> Z28[5]={Z28N30, Z28N32, Z28N33,
                                                         Z28N34, Z28N36};
  //==> Cu(Z=29)
  static const G4double pZ29N34[4]={1.E-9, 3200., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z29N34=std::make_pair(34,pZ29N34);
  static const G4double pZ29N36[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z29N36=std::make_pair(36,pZ29N36);
  static const std::pair<G4int, const G4double*> Z29[2]={Z29N34, Z29N34};
  //==> Zn(Z=30)
  static const G4double pZ30N34[4]={2.E-10, 140., .02, 8.E-6}; // *** only NAT ***
  static const std::pair<G4int, const G4double*> Z30N34=std::make_pair(34,pZ30N34);
  static const G4double pZ30N36[4]={2.E-10, 140., .02, 8.E-6}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z30N36=std::make_pair(36,pZ30N36);
  static const G4double pZ30N37[4]={2.E-10, 140., .02, 8.E-6}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z30N37=std::make_pair(37,pZ30N37);
  static const G4double pZ30N38[4]={2.E-10, 140., .02, 8.E-6}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z30N38=std::make_pair(38,pZ30N38);
  static const G4double pZ30N40[4]={2.E-10, 140., .02, 8.E-6}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z30N40=std::make_pair(40,pZ30N40);
  static const std::pair<G4int, const G4double*> Z30[5]={Z30N34, Z30N36, Z30N37,
                                                         Z30N38, Z30N40};
  //==> Ga(Z=31)
  static const G4double pZ31N38[4]={1.E-9, 3200., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z31N38=std::make_pair(38,pZ31N38);
  static const G4double pZ31N40[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z31N40=std::make_pair(40,pZ31N40);
  static const std::pair<G4int, const G4double*> Z31[2]={Z31N38, Z31N40};
  //==> Ge(Z=32)
  static const G4double pZ32N38[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z32N38=std::make_pair(38,pZ32N38);
  static const G4double pZ32N40[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z32N40=std::make_pair(40,pZ32N40);
  static const G4double pZ32N41[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z32N41=std::make_pair(41,pZ32N41);
  static const G4double pZ32N42[4]={4.E-8, 0., .021, 1.5E-5};
  static const std::pair<G4int, const G4double*> Z32N42=std::make_pair(42,pZ32N42);
  static const G4double pZ32N44[4]={4.E-8, 0., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z32N44=std::make_pair(44,pZ32N44);
  static const std::pair<G4int, const G4double*> Z32[5]={Z32N38, Z32N40, Z32N41,
                                                         Z32N42, Z32N44};
  //==> As(Z=33)
  static const G4double pZ33N42[4]={3.E-9, 200., .051, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z33N42=std::make_pair(42,pZ33N42);
  static const std::pair<G4int, const G4double*> Z33[1]={Z33N42};
  //==> Se(Z=34)
  static const G4double pZ34N40[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N40=std::make_pair(40,pZ34N40);
  static const G4double pZ34N42[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N42=std::make_pair(42,pZ34N42);
  static const G4double pZ34N43[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N43=std::make_pair(43,pZ34N43);
  static const G4double pZ34N44[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N44=std::make_pair(44,pZ34N44);
  static const G4double pZ34N46[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N46=std::make_pair(46,pZ34N46);
  static const G4double pZ34N48[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N48=std::make_pair(48,pZ34N48);
  static const std::pair<G4int, const G4double*> Z34[6]={Z34N40, Z34N42, Z34N43,
                                                         Z34N44, Z34N46, Z34N48};
  //==> Br(Z=35)
  static const G4double pZ35N44[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z35N44=std::make_pair(44,pZ35N44);
  static const G4double pZ35N46[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z35N46=std::make_pair(46,pZ35N46);
  static const std::pair<G4int, const G4double*> Z35[2]={Z35N44, Z35N46};
  //==> Kr(Z=36)
  static const G4double pZ36N42[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N42=std::make_pair(42,pZ36N42);
  static const G4double pZ36N44[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N44=std::make_pair(44,pZ36N44);
  static const G4double pZ36N46[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N46=std::make_pair(46,pZ36N46);
  static const G4double pZ36N47[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N47=std::make_pair(47,pZ36N47);
  static const G4double pZ36N48[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N48=std::make_pair(48,pZ36N48);
  static const G4double pZ36N50[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N50=std::make_pair(50,pZ36N50);
  static const std::pair<G4int, const G4double*> Z36[6]={Z36N42, Z36N44, Z36N46,
                                                         Z36N47, Z36N48, Z36N50};
  //==> Rb(Z=37)
  static const G4double pZ37N48[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z37N48=std::make_pair(48,pZ37N48);
  static const G4double pZ37N50[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z37N50=std::make_pair(50,pZ37N50);
  static const std::pair<G4int, const G4double*> Z37[2]={Z37N48, Z37N50};
  //==> Sr(Z=38)
  static const G4double pZ38N46[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N46=std::make_pair(46,pZ38N46);
  static const G4double pZ38N48[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N48=std::make_pair(48,pZ38N48);
  static const G4double pZ38N49[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N49=std::make_pair(49,pZ38N49);
  static const G4double pZ38N50[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N50=std::make_pair(50,pZ38N50);
  static const std::pair<G4int, const G4double*> Z38[4]={Z38N46, Z38N48, Z38N49, Z38N50};
  //==> Y (Z=39)
  static const G4double pZ39N50[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z39N50=std::make_pair(50,pZ39N50);
  static const std::pair<G4int, const G4double*> Z39[1]={Z39N50};
  //==> Zr(Z=40)
  static const G4double pZ40N50[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N50=std::make_pair(50,pZ40N50);
  static const G4double pZ40N51[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N51=std::make_pair(51,pZ40N51);
  static const G4double pZ40N52[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N52=std::make_pair(52,pZ40N52);
  static const G4double pZ40N54[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N54=std::make_pair(54,pZ40N54);
  static const G4double pZ40N56[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N56=std::make_pair(56,pZ40N56);
  static const std::pair<G4int, const G4double*> Z40[5]={Z40N50, Z40N51, Z40N52,
                                                         Z40N54, Z40N56};
  //==> Nb(Z=41)
  static const G4double pZ41N52[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z41N52=std::make_pair(52,pZ41N52);
  static const std::pair<G4int, const G4double*> Z41[1]={Z41N52};
  //==> Mo(Z=42)
  static const G4double pZ42N50[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N50=std::make_pair(50,pZ42N50);
  static const G4double pZ42N52[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N52=std::make_pair(52,pZ42N52);
  static const G4double pZ42N53[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N53=std::make_pair(53,pZ42N53);
  static const G4double pZ42N54[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N54=std::make_pair(54,pZ42N54);
  static const G4double pZ42N55[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N55=std::make_pair(55,pZ42N55);
  static const G4double pZ42N56[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N56=std::make_pair(56,pZ42N56);
  static const G4double pZ42N58[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N58=std::make_pair(58,pZ42N58);
  static const std::pair<G4int, const G4double*> Z42[7]={Z42N50, Z42N52, Z42N53, Z42N54,
                                                         Z42N55, Z42N56, Z42N58};
  //==> Mo(Z=43)
  static const G4double pZ43N0[4]={3.E-12, 500., .01, 2.5E-4}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z43N0=std::make_pair(0,pZ43N0);
  static const std::pair<G4int, const G4double*> Z43[1]={Z43N0};
  //==> Ru(Z=44)
  static const G4double pZ44N52[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N52=std::make_pair(52,pZ44N52);
  static const G4double pZ44N54[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N54=std::make_pair(54,pZ44N54);
  static const G4double pZ44N55[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N55=std::make_pair(55,pZ44N55);
  static const G4double pZ44N56[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N56=std::make_pair(56,pZ44N56);
  static const G4double pZ44N57[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N57=std::make_pair(57,pZ44N57);
  static const G4double pZ44N58[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N58=std::make_pair(58,pZ44N58);
  static const G4double pZ44N60[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N60=std::make_pair(60,pZ44N60);
  static const std::pair<G4int, const G4double*> Z44[7]={Z44N52, Z44N54, Z44N55, Z44N56,
                                                         Z44N57, Z44N58, Z44N60};
  //==> Rh(Z=45)
  static const G4double pZ45N58[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z45N58=std::make_pair(58,pZ45N58);
  static const std::pair<G4int, const G4double*> Z45[1]={Z45N58};
  //==> Pd(Z=46)
  static const G4double pZ46N56[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N56=std::make_pair(56,pZ46N56);
  static const G4double pZ46N58[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N58=std::make_pair(58,pZ46N58);
  static const G4double pZ46N59[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N59=std::make_pair(59,pZ46N59);
  static const G4double pZ46N60[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N60=std::make_pair(60,pZ46N60);
  static const G4double pZ46N62[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N62=std::make_pair(62,pZ46N62);
  static const G4double pZ46N64[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N64=std::make_pair(64,pZ46N64);
  static const std::pair<G4int, const G4double*> Z46[6]={Z46N56, Z46N58, Z46N59,
                                                         Z46N60, Z46N62, Z46N64};
  //==> Ag(Z=47)
  static const G4double pZ47N60[4]={3.E-12, 500., .01, 2.5E-5};
  static const std::pair<G4int, const G4double*> Z47N60=std::make_pair(48,pZ47N60);
  static const G4double pZ47N62[4]={3.E-12, 500., .01, 2.5E-4};
  static const std::pair<G4int, const G4double*> Z47N62=std::make_pair(50,pZ47N62);
  static const std::pair<G4int, const G4double*> Z47[2]={Z47N60, Z47N62};
  //==> Cd(Z=48)
  static const G4double pZ48N58[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N58=std::make_pair(58,pZ48N58);
  static const G4double pZ48N60[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N60=std::make_pair(60,pZ48N60);
  static const G4double pZ48N62[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N62=std::make_pair(62,pZ48N62);
  static const G4double pZ48N63[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N63=std::make_pair(63,pZ48N63);
  static const G4double pZ48N64[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N64=std::make_pair(64,pZ48N64);
  static const G4double pZ48N65[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N65=std::make_pair(65,pZ48N65);
  static const G4double pZ48N66[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N66=std::make_pair(66,pZ48N66);
  static const G4double pZ48N68[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N68=std::make_pair(68,pZ48N68);
  static const std::pair<G4int, const G4double*> Z48[8]={Z48N58, Z48N60, Z48N62, Z48N63,
                                                         Z48N64, Z48N65, Z48N66, Z48N68};
  //==> In(Z=49)
  static const G4double pZ49N64[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z49N64=std::make_pair(64,pZ49N64);
  static const G4double pZ49N66[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z49N66=std::make_pair(66,pZ49N66);
  static const std::pair<G4int, const G4double*> Z49[2]={Z49N64, Z49N66};
  //==> Sn(Z=50)
  static const G4double pZ50N62[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N62=std::make_pair(62,pZ50N62);
  static const G4double pZ50N64[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N64=std::make_pair(64,pZ50N64);
  static const G4double pZ50N65[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N65=std::make_pair(65,pZ50N65);
  static const G4double pZ50N66[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N66=std::make_pair(66,pZ50N66);
  static const G4double pZ50N67[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N67=std::make_pair(67,pZ50N67);
  static const G4double pZ50N68[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N68=std::make_pair(68,pZ50N68);
  static const G4double pZ50N69[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N69=std::make_pair(69,pZ50N69);
  static const G4double pZ50N70[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N70=std::make_pair(70,pZ50N70);
  static const G4double pZ50N72[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N72=std::make_pair(72,pZ50N72);
  static const G4double pZ50N74[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N74=std::make_pair(74,pZ50N74);
  static const std::pair<G4int, const G4double*> Z50[10]={Z50N62, Z50N64, Z50N65, Z50N66,
                                                          Z50N67, Z50N68, Z50N69, Z50N70,
                                                          Z50N72, Z50N74};
  //==> Sb(Z=51)
  static const G4double pZ51N70[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z51N70=std::make_pair(70,pZ51N70);
  static const G4double pZ51N72[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z51N72=std::make_pair(72,pZ51N72);
  static const std::pair<G4int, const G4double*> Z51[2]={Z51N70, Z51N72};
  //==> Te(Z=52)
  static const G4double pZ52N68[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N68=std::make_pair(68,pZ52N68);
  static const G4double pZ52N70[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N70=std::make_pair(70,pZ52N70);
  static const G4double pZ52N71[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N71=std::make_pair(71,pZ52N71);
  static const G4double pZ52N72[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N72=std::make_pair(72,pZ52N72);
  static const G4double pZ52N73[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N73=std::make_pair(73,pZ52N73);
  static const G4double pZ52N74[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N74=std::make_pair(74,pZ52N74);
  static const G4double pZ52N76[4]={3.E-12, 500., .01, 2.5E-5}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N76=std::make_pair(76,pZ52N76);
  static const G4double pZ52N78[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N78=std::make_pair(78,pZ52N78);
  static const std::pair<G4int, const G4double*> Z52[8]={Z52N68, Z52N70, Z52N71, Z52N72,
                                                         Z52N73, Z52N74, Z52N76, Z52N78};
  //==> I(Z=53)
  static const G4double pZ53N74[4]={3.E-12, 500., .01, 2.5E-4}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z53N74=std::make_pair(74,pZ53N74);
  static const std::pair<G4int, const G4double*> Z53[1]={Z53N74};
 
  
  static const std::pair<G4int, const G4double*>* Pars[54]={Z0,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,
    Z10,Z11,Z12,Z13,Z14,Z15,Z16,Z17,Z18,Z19,Z20,Z21,Z22,Z23,Z24,Z25,Z26,Z27,Z28,Z29,Z30,
    Z31,Z32,Z33,Z34,Z35,Z36,Z37,Z38,Z39,Z40,Z41,Z42,Z43,Z44,Z45,Z46,Z47,Z48,Z49,Z50,Z51,
    Z52,Z53};
  //Z52,Z53,Z54,Z55,Z56,Z57,Z58,Z59,Z60,Z61,Z62,Z63,Z64,Z65,Z66,Z67,Z68,Z69,Z70,Z71,Z72,
  //Z73,Z74,Z75,Z76,Z77,Z78,Z79,Z80,Z81,Z82,Z83,Z84,Z85,Z86,Z87,Z88,Z89,Z90,Z91,Z92,Z93,
  //Z94,Z95,Z96};
  G4int curN=Pars[1][0].first;
  G4double par=Pars[1][0].second[1];
  G4cout<<"-Warning-G4QNeutronNuclearCroSect::CSLin: N="<<curN<<", P="<<par<<G4endl;
  G4double sigma=0.;
  if(P<ThresholdMomentum(tZ,tN)*.001) return sigma;
  G4double lP=std::log(P);
  if(tZ==1&&!tN){if(P>.35) sigma=CrossSectionFormula(tZ,tN,P,lP);}// s(np)=0 below 350Mev/c
  else if(tZ<97 && tN<248)                // General solution
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
      sigma+=pex*std::exp(-dp*dp/wid);
    }
  }
  else
  {
    G4cerr<<"-Warning-G4QNeutronNuclearCroSect::CSLin:*Bad A* Z="<<tZ<<", N="<<tN<<G4endl;
    sigma=0.;
  }
  if(sigma<0.) return 0.;
  return sigma;  
}

// Calculation formula for proton-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4QNeutronNuclearCrossSection::CrossSectionLog(G4int tZ, G4int tN, G4double lP)
{
  G4double P=std::exp(lP);
  return CrossSectionFormula(tZ, tN, P, lP);
}
// Calculation formula for proton-nuclear inelastic cross-section (mb) log(P in GeV/c)
G4double G4QNeutronNuclearCrossSection::CrossSectionFormula(G4int tZ, G4int tN,
                                                           G4double P, G4double lP)
{
  G4double sigma=0.;
  if(tZ==1 && !tN)                        // pp interaction
  {
    G4double sp=std::sqrt(P);
    G4double ds=lP-4.2;
    G4double dp=P-.35;
    G4double d3=dp*dp*dp;
    sigma=(33.+.2*ds*ds)/(1.+.4/sp)/(1.+.5/d3/d3);
  }
  else if(tZ<93 && tN<146)                // General solution
  {
    //G4double lP=std::log(P);            // Already calculated
    G4double d=lP-4.2;
    G4double p2=P*P;
    G4double p4=p2*p2;
    G4double a=tN+tZ;                       // A of the target
    G4double al=std::log(a);
    G4double a2=a*a;
    G4double a4=a2*a2;
    G4double a8=a4*a4;
    G4double a12=a8*a4;
    G4double a16=a8*a8;
    G4double c=170./(1.+5./std::exp(al*1.6));
    G4double dl=al-2.8;
    G4double r=.3+.2*dl*dl;
    G4double g=40.*std::exp(al*0.712)/(1.+12.2/a)/(1.+34./a2);
    G4double e=318.+a4/(1.+.0015*a4/std::exp(al*0.09))/(1.+4.e-28*a12)+
               8.e-18/(1./a16+1.3e-20)/(1.+1.e-21*a12);
    G4double s=3.57+.009*a2/(1.+.0001*a2*a);
    G4double h=(.01/a4+2.5e-6/a)*(1.+7.e-8*a4)/(1.+6.e7/a12/a2);
    sigma=(c+d*d)/(1.+r/p4)+(g+e*std::exp(-s*P))/(1.+h/p4/p4);
#ifdef pdebug
  G4cout<<"G4QProtNucCS::CSForm: A="<<a<<",P="<<P<<",CS="<<sigma<<",c="<<c<<",g="<<g<<",d="
        <<d<<",r="<<r<<",e="<<e<<",h="<<h<<G4endl;
#endif
  }
  else
  {
    G4cerr<<"-Warning-G4QProtonNuclearCroSect::CSForm:*Bad A* Z="<<tZ<<", N="<<tN<<G4endl;
    sigma=0.;
  }
  if(sigma<0.) return 0.;
  return sigma;  
}
