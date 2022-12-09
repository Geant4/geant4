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
//
//
// G4 Physics class: G4ChipsNeutronElasticXS for nA elastic cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 12-Jan-10 (from G4QElCrSect)
//
// -------------------------------------------------------------------------------
// Short description: Interaction cross-sections for the elastic process. 
// Class extracted from CHIPS and integrated in Geant4 by W.Pokorski
// -------------------------------------------------------------------------------


#include "G4ChipsNeutronElasticXS.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Neutron.hh"
#include "G4Nucleus.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4IonTable.hh"

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4ChipsNeutronElasticXS);

namespace {
    G4double mNeut;
    G4double mProt;
    G4double mNeut2;

}

G4ChipsNeutronElasticXS::G4ChipsNeutronElasticXS():G4VCrossSectionDataSet(Default_Name()), nPoints(128), nLast(nPoints-1) 
{
  lPMin=-8.;  // Min tabulated log Momentum      (D)
  lPMax= 8.;  // Max tabulated log Momentum      (D)
  dlnP=(lPMax-lPMin)/nLast;// Log step in table  (D)
  onlyCS=true;// Flag to calc only CS (not Si/Bi)(L)
  lastSIG=0.; // Last calculated cross section   (L)
  lastLP=-10.;// Last log(momOfIncidentHadron)   (L)
  lastTM=0.;  // Last t_maximum                  (L)
  theSS=0.;   // The Last sq.slope of 1st difMax (L)
  theS1=0.;   // The Last mantissa of 1st difMax (L)
  theB1=0.;   // The Last slope of 1st difr. Max (L)
  theS2=0.;   // The Last mantissa of 2nd difMax (L)
  theB2=0.;   // The Last slope of 2nd difr. Max (L)
  theS3=0.;   // The Last mantissa of 3d difrMax (L)
  theB3=0.;   // The Last slope of 3d difructMax (L)
  theS4=0.;   // The Last mantissa of 4th difMax (L)
  theB4=0.;   // The Last slope of 4th difr. Max (L)
  lastTZ=0;   // Last atomic number of the target
  lastTN=0;   // Last # of neutrons in the target
  lastPIN=0.; // Last initialized max momentum
  lastCST=0;  // Elastic cross-section table
  lastPAR=0;  // Parameters of FunctionalCalculation
  lastSST=0;  // E-dep of sq.slope of the 1st difMax
  lastS1T=0;  // E-dep of mantissa of the 1st difMax
  lastB1T=0;  // E-dep of theSlope of the 1st difMax
  lastS2T=0;  // E-dep of mantissa of the 2nd difMax
  lastB2T=0;  // E-dep of theSlope of the 2nd difMax
  lastS3T=0;  // E-dep of mantissa of the 3d difrMax
  lastB3T=0;  // E-dep of the slope of the 3d difMax
  lastS4T=0;  // E-dep of mantissa of the 4th difMax
  lastB4T=0;  // E-dep of theSlope of the 4th difMax
  lastN=0;    // The last N of calculated nucleus
  lastZ=0;    // The last Z of calculated nucleus
  lastP=0.;   // Last used in cross section Momentum
  lastTH=0.;  // Last threshold momentum
  lastCS=0.;  // Last value of the Cross Section
  lastI=0;    // The last position in the DAMDB

    mNeut= G4Neutron::Neutron()->GetPDGMass()*.001;// MeV to GeV
    mProt= G4Proton::Proton()->GetPDGMass()*.001;// MeV to GeV
    mNeut2= mNeut*mNeut;
}

G4ChipsNeutronElasticXS::~G4ChipsNeutronElasticXS()
{
  std::vector<G4double*>::iterator pos;
  for (pos=CST.begin(); pos<CST.end(); pos++)
  { delete [] *pos; }
  CST.clear();
  for (pos=PAR.begin(); pos<PAR.end(); pos++)
  { delete [] *pos; }
  PAR.clear();
  for (pos=SST.begin(); pos<SST.end(); pos++)
  { delete [] *pos; }
  SST.clear();
  for (pos=S1T.begin(); pos<S1T.end(); pos++)
  { delete [] *pos; }
  S1T.clear();
  for (pos=B1T.begin(); pos<B1T.end(); pos++)
  { delete [] *pos; }
  B1T.clear();
  for (pos=S2T.begin(); pos<S2T.end(); pos++)
  { delete [] *pos; }
  S2T.clear();
  for (pos=B2T.begin(); pos<B2T.end(); pos++)
  { delete [] *pos; }
  B2T.clear();
  for (pos=S3T.begin(); pos<S3T.end(); pos++)
  { delete [] *pos; }
  S3T.clear();
  for (pos=B3T.begin(); pos<B3T.end(); pos++)
  { delete [] *pos; }
  B3T.clear();
  for (pos=S4T.begin(); pos<S4T.end(); pos++)
  { delete [] *pos; }
  S4T.clear();
  for (pos=B4T.begin(); pos<B4T.end(); pos++)
  { delete [] *pos; }
  B4T.clear(); 
}

void
G4ChipsNeutronElasticXS::CrossSectionDescription(std::ostream& outFile) const
{
    outFile << "G4ChipsNeutronElasticXS provides the elastic cross\n"
            << "section for neutron nucleus scattering as a function of incident\n"
            << "momentum. The cross section is calculated using M. Kossov's\n"
            << "CHIPS parameterization of cross section data.\n";
}

G4bool G4ChipsNeutronElasticXS::IsIsoApplicable(const G4DynamicParticle*, G4int, G4int,    
				 const G4Element*,
				 const G4Material*)
{
  return true;
}

G4double G4ChipsNeutronElasticXS::GetIsoCrossSection(const G4DynamicParticle* Pt, G4int tgZ, G4int A,  
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
G4double G4ChipsNeutronElasticXS::GetChipsCrossSection(G4double pMom, G4int tgZ, G4int tgN, G4int)
{

  G4double pEn=pMom;
  onlyCS=false;

  G4bool in=false;                   // By default the isotope must be found in the AMDB
  lastP   = 0.;                      // New momentum history (nothing to compare with)
  lastN   = tgN;                     // The last N of the calculated nucleus
  lastZ   = tgZ;                     // The last Z of the calculated nucleus
  lastI   = (G4int)colN.size();      // Size of the Associative Memory DB in the heap
  if(lastI) for(G4int i=0; i<lastI; ++i) // Loop over proj/tgZ/tgN lines of DB
  {                                  // The nucleus with projPDG is found in AMDB
    if(colN[i]==tgN && colZ[i]==tgZ) // Isotope is foind in AMDB
    {
      lastI=i;
      lastTH =colTH[i];              // Last THreshold (A-dependent)
      if(pEn<=lastTH)
      {
        return 0.;                   // Energy is below the Threshold value
      }
      lastP  =colP [i];                // Last Momentum  (A-dependent)
      lastCS =colCS[i];                // Last CrossSect (A-dependent)
      //  if(std::fabs(lastP/pMom-1.)<tolerance) //VI (do not use tolerance)
      if(lastP == pMom)              // Do not recalculate
      {
        CalculateCrossSection(false,-1,i,2112,lastZ,lastN,pMom); // Update param's only
        return lastCS*millibarn;     // Use theLastCS
      }
      in = true;                       // This is the case when the isotop is found in DB
      // Momentum pMom is in IU ! @@ Units
      lastCS=CalculateCrossSection(false,-1,i,2112,lastZ,lastN,pMom); // read & update
      if(lastCS<=0. && pEn>lastTH)    // Correct the threshold
      {
        lastTH=pEn;
      }
      break;                           // Go out of the LOOP with found lastI
    }
  }
  if(!in)                            // This nucleus has not been calculated previously
  {
    //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
    lastCS=CalculateCrossSection(false,0,lastI,2112,lastZ,lastN,pMom);//calculate&create
    if(lastCS<=0.)
    {
      lastTH = 0; // ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last
      if(pEn>lastTH)
      {
        lastTH=pEn;
      }
    }
    colN.push_back(tgN);
    colZ.push_back(tgZ);
    colP.push_back(pMom);
    colTH.push_back(lastTH);
    colCS.push_back(lastCS);
    return lastCS*millibarn;
  } // End of creation of the new set of parameters
  else
  {
    colP[lastI]=pMom;
    colCS[lastI]=lastCS;
  }
  return lastCS*millibarn;
}

// Calculation of total elastic cross section (p in IU, CS in mb) @@ Units (?)
// F=0 - create AMDB, F=-1 - read&update AMDB, F=1 - update AMDB (sinchro with higher AMDB)
G4double G4ChipsNeutronElasticXS::CalculateCrossSection(G4bool CS, G4int F,G4int I,
                                             G4int PDG, G4int tgZ, G4int tgN, G4double pIU)
{

  G4double pMom=pIU/GeV;                // All calculations are in GeV
  onlyCS=CS;                            // Flag to calculate only CS (not Si/Bi)
  lastLP=std::log(pMom);                // Make a logarithm of the momentum for calculation
  if(F)                                 // This isotope was found in AMDB =>RETRIEVE/UPDATE
  {
    if(F<0)                             // the AMDB must be loded
    {
      lastPIN = PIN[I];                 // Max log(P) initialised for this table set
      lastPAR = PAR[I];                 // Pointer to the parameter set
      lastCST = CST[I];                 // Pointer to the total sross-section table
      lastSST = SST[I];                 // Pointer to the first squared slope
      lastS1T = S1T[I];                 // Pointer to the first mantissa
      lastB1T = B1T[I];                 // Pointer to the first slope
      lastS2T = S2T[I];                 // Pointer to the second mantissa
      lastB2T = B2T[I];                 // Pointer to the second slope
      lastS3T = S3T[I];                 // Pointer to the third mantissa
      lastB3T = B3T[I];                 // Pointer to the rhird slope
      lastS4T = S4T[I];                 // Pointer to the 4-th mantissa
      lastB4T = B4T[I];                 // Pointer to the 4-th slope
    }
    if(lastLP>lastPIN && lastLP<lPMax)
    {
      lastPIN=GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);// Can update upper logP-Limit in tabs
      PIN[I]=lastPIN;                   // Remember the new P-Limit of the tables
    }
  }
  else                                  // This isotope wasn't initialized => CREATE
  {
    lastPAR = new G4double[nPoints];    // Allocate memory for parameters of CS function
    lastPAR[nLast]=0;                   // Initialization for VALGRIND
    lastCST = new G4double[nPoints];    // Allocate memory for Tabulated CS function    
    lastSST = new G4double[nPoints];    // Allocate memory for Tabulated first sqaredSlope 
    lastS1T = new G4double[nPoints];    // Allocate memory for Tabulated first mantissa 
    lastB1T = new G4double[nPoints];    // Allocate memory for Tabulated first slope    
    lastS2T = new G4double[nPoints];    // Allocate memory for Tabulated second mantissa
    lastB2T = new G4double[nPoints];    // Allocate memory for Tabulated second slope   
    lastS3T = new G4double[nPoints];    // Allocate memory for Tabulated third mantissa 
    lastB3T = new G4double[nPoints];    // Allocate memory for Tabulated third slope    
    lastS4T = new G4double[nPoints];    // Allocate memory for Tabulated 4-th mantissa 
    lastB4T = new G4double[nPoints];    // Allocate memory for Tabulated 4-th slope    
    lastPIN = GetPTables(lastLP,lPMin,PDG,tgZ,tgN); // Returns the new P-limit for tables
    PIN.push_back(lastPIN);             // Fill parameters of CS function to AMDB
    PAR.push_back(lastPAR);             // Fill parameters of CS function to AMDB
    CST.push_back(lastCST);             // Fill Tabulated CS function to AMDB    
    SST.push_back(lastSST);             // Fill Tabulated first sq.slope to AMDB 
    S1T.push_back(lastS1T);             // Fill Tabulated first mantissa to AMDB 
    B1T.push_back(lastB1T);             // Fill Tabulated first slope to AMDB    
    S2T.push_back(lastS2T);             // Fill Tabulated second mantissa to AMDB 
    B2T.push_back(lastB2T);             // Fill Tabulated second slope to AMDB    
    S3T.push_back(lastS3T);             // Fill Tabulated third mantissa to AMDB 
    B3T.push_back(lastB3T);             // Fill Tabulated third slope to AMDB    
    S4T.push_back(lastS4T);             // Fill Tabulated 4-th mantissa to AMDB 
    B4T.push_back(lastB4T);             // Fill Tabulated 4-th slope to AMDB    
  } // End of creation/update of the new set of parameters and tables
  // =-------= NOW Update (if necessary) and Calculate the Cross Section =---------=
  if(lastLP>lastPIN && lastLP<lPMax)
  {
    lastPIN = GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);
  }
  if(!onlyCS) lastTM=GetQ2max(PDG, tgZ, tgN, pMom); // Calculate (-t)_max=Q2_max (GeV2)
  if(lastLP>lPMin && lastLP<=lastPIN)   // Linear fit is made using precalculated tables
  {
    if(lastLP==lastPIN)
    {
      G4double shift=(lastLP-lPMin)/dlnP+.000001; // Log distance from lPMin
      G4int    blast=static_cast<int>(shift); // this is a bin number of the lower edge (0)
      if(blast<0 || blast>=nLast) G4cout<<"G4QNeutElCS::CCS:b="<<blast<<","<<nLast<<G4endl;
      lastSIG = lastCST[blast];
      if(!onlyCS)                       // Skip the differential cross-section parameters
      {
        theSS  = lastSST[blast];
        theS1  = lastS1T[blast];
        theB1  = lastB1T[blast];
        theS2  = lastS2T[blast];
        theB2  = lastB2T[blast];
        theS3  = lastS3T[blast];
        theB3  = lastB3T[blast];
        theS4  = lastS4T[blast];
        theB4  = lastB4T[blast];
      }
    }
    else
    {
      G4double shift=(lastLP-lPMin)/dlnP;        // a shift from the beginning of the table
      G4int    blast=static_cast<int>(shift);    // the lower bin number
      if(blast<0)   blast=0;
      if(blast>=nLast) blast=nLast-1;            // low edge of the last bin
      shift-=blast;                              // step inside the unit bin
      G4int lastL=blast+1;                       // the upper bin number
      G4double SIGL=lastCST[blast];              // the basic value of the cross-section
      lastSIG= SIGL+shift*(lastCST[lastL]-SIGL); // calculated total elastic cross-section
      if(!onlyCS)                       // Skip the differential cross-section parameters
      {
        G4double SSTL=lastSST[blast];           // the low bin of the first squared slope
        theSS=SSTL+shift*(lastSST[lastL]-SSTL); // the basic value of the first sq.slope
        G4double S1TL=lastS1T[blast];           // the low bin of the first mantissa
        theS1=S1TL+shift*(lastS1T[lastL]-S1TL); // the basic value of the first mantissa
        G4double B1TL=lastB1T[blast];           // the low bin of the first slope
        theB1=B1TL+shift*(lastB1T[lastL]-B1TL); // the basic value of the first slope
        G4double S2TL=lastS2T[blast];           // the low bin of the second mantissa
        theS2=S2TL+shift*(lastS2T[lastL]-S2TL); // the basic value of the second mantissa
        G4double B2TL=lastB2T[blast];           // the low bin of the second slope
        theB2=B2TL+shift*(lastB2T[lastL]-B2TL); // the basic value of the second slope
        G4double S3TL=lastS3T[blast];           // the low bin of the third mantissa
        theS3=S3TL+shift*(lastS3T[lastL]-S3TL); // the basic value of the third mantissa
        G4double B3TL=lastB3T[blast];           // the low bin of the third slope
        theB3=B3TL+shift*(lastB3T[lastL]-B3TL); // the basic value of the third slope
        G4double S4TL=lastS4T[blast];           // the low bin of the 4-th mantissa
        theS4=S4TL+shift*(lastS4T[lastL]-S4TL); // the basic value of the 4-th mantissa
        G4double B4TL=lastB4T[blast];           // the low bin of the 4-th slope
        theB4=B4TL+shift*(lastB4T[lastL]-B4TL); // the basic value of the 4-th slope
      }
    }
  }
  else lastSIG=GetTabValues(lastLP, PDG, tgZ, tgN); // Direct calculation beyond the table
  if(lastSIG<0.) lastSIG = 0.;                   // @@ a Warning print can be added
  return lastSIG;
}

// It has parameter sets for all tZ/tN/PDG, using them the tables can be created/updated
G4double G4ChipsNeutronElasticXS::GetPTables(G4double LP, G4double ILP, G4int PDG,
                                                   G4int tgZ, G4int tgN)
{
  // @@ At present all nA==pA ---------> Each neucleus can have not more than 51 parameters
  static const G4double pwd=2727;
  const G4int n_npel=24;                // #of parameters for np-elastic (<nPoints=128)
  const G4int n_ppel=32;                // #of parameters for pp-elastic (<nPoints=128)
  //                      -0- -1-  -2- -3- -4-  -5- -6- -7- -8- -9--10--11--12--13- -14-
  G4double np_el[n_npel]={12.,.05,.0001,5.,.35,6.75,.14,19.,.6,6.75,.14,13.,.14,.6,.00013,
                          75.,.001,7.2,4.32,.012,2.5,0.0,12.,.34};
  //                      -15--16--17- -18- -19--20--21--22--23-
  //                       -0-   -1-  -2- -3- -4- -5-  -6-  -7-  -8--9--10--11--12--13-
  G4double pp_el[n_ppel]={2.865,18.9,.6461,3.,9.,.425,.4276,.0022,5.,74.,3.,3.4,.2,.17,
                          .001,8.,.055,3.64,5.e-5,4000.,1500.,.46,1.2e6,3.5e6,5.e-5,1.e10,
                          8.5e8,1.e10,1.1,3.4e6,6.8e6,0.};
  //                      -14--15- -16- -17- -18-  -19- -20- -21- -22-  -23-   -24-  -25-
  //                       -26- -27-  -28- -29- -30- -31-
  //==> n (Z=0)
  static const G4int N0=1; // *** Not used (fake)***
  static const G4double pZ0N1[7]={0., 0., 0., 0., 0., 0., 0.}; // Not used (fake)
  static const std::pair<G4int, const G4double*> Z0N1(1,pZ0N1);
  static const std::pair<G4int, const G4double*> Z0[N0]={Z0N1};
  //==> H (Z=1) *** protons are treated separately ***
  static const G4int N1=3;
  static const G4double pZ1N0[7]={0., 0., 0., 0., 0., 0., 0.}; // Not used (fake)
  static const std::pair<G4int, const G4double*> Z1N0(0,pZ1N0);
  static const G4double pZ1N1[7]={6.E-5, 4., .055, 1.1E-8, .008, 1.2E-8, .019};
  static const std::pair<G4int, const G4double*> Z1N1(1,pZ1N1);
  static const G4double pZ1N2[7]={6.E-5, 2.2, .051, 1.E-8, .04, 9.E-8, .0075};
  static const std::pair<G4int, const G4double*> Z1N2(2,pZ1N2);
  static const std::pair<G4int, const G4double*> Z1[N1]={Z1N0, Z1N1, Z1N2};
  //==> He(Z=2)
  static const G4int N2=2;
  static const G4double pZ2N1[7]={6.E-5, 3., .06, 4.E-9, .03, 7.E-8, .015};
  static const std::pair<G4int, const G4double*> Z2N1(1,pZ2N1);
  static const G4double pZ2N2[7]={3.E-4, .23, 1., 1.5E-9, 2.E-02, 1.E-8, .003};
  static const std::pair<G4int, const G4double*> Z2N2(2,pZ2N2);
  static const std::pair<G4int, const G4double*> Z2[N2]={Z2N1, Z2N2};
  //==> Li(Z=3)
  static const G4int N3=2;
  static const G4double pZ3N3[7]={3.1E-7, 1.7, 1.3E-4, 1.E-8, .02, 1.1E-7, .0023};
  static const std::pair<G4int, const G4double*> Z3N1(3,pZ3N3);
  static const G4double pZ3N4[7]={1.3E-6, 1.8, 7.6E-4, 9.E-9, .03, 1.E-7, .0029};
  static const std::pair<G4int, const G4double*> Z3N2(4,pZ3N4);
  static const std::pair<G4int, const G4double*> Z3[N3]={Z3N1, Z3N2};
  //==> Be(Z=4)
  static const G4int N4=2;
  static const G4double pZ4N3[7]={2.E-4, 1.4, 2.7, 0., .02, 5.E-8, 0.};
  static const std::pair<G4int, const G4double*> Z4N3(3,pZ4N3);
  static const G4double pZ4N5[7]={1.E-6, 5.7, .0011, 3.E-9, .007, 2.E-8, .016};
  static const std::pair<G4int, const G4double*> Z4N5(5,pZ4N5);
  static const std::pair<G4int, const G4double*> Z4[N4]={Z4N3,Z4N5};
  //==> B (Z=5)
  static const G4int N5=2;
  static const G4double pZ5N5[7]={8.E-7, 5., 3.4E-4, 7.E-9, 1.E-02, 1.E-07, .0053};
  static const std::pair<G4int, const G4double*> Z5N5(5,pZ5N5);
  static const G4double pZ5N6[7]={4.8E-6, 6.6, .0035, 4.E-9, .003, 1.E-8, .012};
  static const std::pair<G4int, const G4double*> Z5N6(6,pZ5N6);
  static const std::pair<G4int, const G4double*> Z5[N5]={Z5N5, Z5N6};
  //==> C (Z=6) *** Only nat (C13=C12=C_nat) ***
  static const G4int N6=2;
  static const G4double pZ6N6[7]={4.9E-6, 6.6, .0035, 4.E-9, .002, 6.E-9, .011};
  static const std::pair<G4int, const G4double*> Z6N6(6,pZ6N6);
  static const G4double pZ6N7[7]={4.9E-6, 6.6, .0035, 4.E-9, .002, 6.E-9, .011};
  static const std::pair<G4int, const G4double*> Z6N7(7,pZ6N7);
  static const std::pair<G4int, const G4double*> Z6[N6]={Z6N6, Z6N7};
  //==> N (Z=7)
  static const G4int N7=2;
  static const G4double pZ7N7[7]={4.9E-6, 1.6, .03, .4E-9, .02, 6.E-8, .021};
  static const std::pair<G4int, const G4double*> Z7N7(7,pZ7N7);
  static const G4double pZ7N8[7]={2.5E-6, 5., .0021, 2.5E-9, .015, 5.E-8, .009};
  static const std::pair<G4int, const G4double*> Z7N8(8,pZ7N8);
  static const std::pair<G4int, const G4double*> Z7[N7]={Z7N7, Z7N8};
  //==> O (Z=8) (O18=O17, No data)
  static const G4int N8=3;
  static const G4double pZ8N8[7]={2.5E-6, 5.3, .0018, 3.E-9, .01, 1.5E-8, .0075};
  static const std::pair<G4int, const G4double*> Z8N8(8,pZ8N8);
  static const G4double pZ8N9[7]={1.4E-6, 2.1, .0025, 1.3E-9, .02, 1.3E-7, .0072};
  static const std::pair<G4int, const G4double*> Z8N9(9,pZ8N9);
  static const G4double pZ8N10[7]={1.4E-6, 2.1, .0025, 1.3E-9, .02, 1.3E-7, .0072};
  static const std::pair<G4int, const G4double*> Z8N10(10,pZ8N10);
  static const std::pair<G4int, const G4double*> Z8[N8]={Z8N8, Z8N9, Z8N10};
  //==> F (Z=9)
  static const G4int N9=1;
  static const G4double pZ9N10[7]={1.4E-6, 7.5, 6.7E-4, 4.E-9, 2.E-5, 7.E-12, .0066};
  static const std::pair<G4int, const G4double*> Z9N10(10,pZ9N10);
  static const std::pair<G4int, const G4double*> Z9[N9]={Z9N10};
  //==> Ne(Z=10) *** No data *** (Ne20=Na22, Ne21=F19, Ne22=Na22)
  static const G4int N10=3;
  static const G4double pZ10N10[7]={1.4E-5, 7., .01, 0., 3.E-11, 7.E-24, .12};
  static const std::pair<G4int, const G4double*> Z10N10(10,pZ10N10);
  static const G4double pZ10N11[7]={1.4E-6, 7.5, 6.7E-4, 4.E-9, 2.E-5, 7.E-12, .0066};
  static const std::pair<G4int, const G4double*> Z10N11(11,pZ10N11);
  static const G4double pZ10N12[7]={1.4E-5, 7., .01, 0., 3.E-11, 7.E-24, .12};
  static const std::pair<G4int, const G4double*> Z10N12(12,pZ10N12);
  static const std::pair<G4int, const G4double*> Z10[N10]={Z10N10, Z10N11, Z10N12};
  //==> Na(Z=11)
  static const G4int N11=2;
  static const G4double pZ11N11[7]={1.4E-5, 7., .01, 0., 3.E-11, 7.E-24, .12};
  static const std::pair<G4int, const G4double*> Z11N11(11,pZ11N11);
  static const G4double pZ11N12[7]={1.4E-6, 7.6, 6.E-4, 5.E-9, 7.E-9, 3.E-18, .0056};
  static const std::pair<G4int, const G4double*> Z11N12(12,pZ11N12);
  static const std::pair<G4int, const G4double*> Z11[N11]={Z11N11, Z11N12};
  //==> Mg(Z=12)
  static const G4int N12=3;
  static const G4double pZ12N12[7]={8.E-7, 3., .001, 1.8E-9, .0015, .2E-9, .006};
  static const std::pair<G4int, const G4double*> Z12N12(12,pZ12N12);
  static const G4double pZ12N13[7]={8.E-7, 7., 3.E-4, 6.E-9, .006, 4.E-8, .0042};
  static const std::pair<G4int, const G4double*> Z12N13(13,pZ12N13);
  static const G4double pZ12N14[7]={1.2E-6, 6.8, 5.E-4, 5.E-9, .007, 2.E-8, .0044};
  static const std::pair<G4int, const G4double*> Z12N14(14,pZ12N14);
  static const std::pair<G4int, const G4double*> Z12[N12]={Z12N12, Z12N13, Z12N14};
  //==> Al(Z=13)
  static const G4int N13=1;
  static const G4double pZ13N14[7]={3.E-7, 5., 8.4E-5, 7.E-9, .008, 2.E-8, .0022};
  static const std::pair<G4int, const G4double*> Z13N14(14,pZ13N14);
  static const std::pair<G4int, const G4double*> Z13[N13]={Z13N14};
  //==> Si(Z=14)
  static const G4int N14=3;
  static const G4double pZ14N14[7]={1.2E-6, 6., 4.E-4, 6.E-9, .012, 8.E-8, .0029};
  static const std::pair<G4int, const G4double*> Z14N14(14,pZ14N14);
  static const G4double pZ14N15[7]={2.4E-6, 4., .0016, 3.E-9, .018, 6.E-8, .0037};
  static const std::pair<G4int, const G4double*> Z14N15(15,pZ14N15);
  static const G4double pZ14N16[7]={6.E-7, 4., 3.7E-4, 3.E-9, .018, 6.E-8, .0036};
  static const std::pair<G4int, const G4double*> Z14N16(16,pZ14N16);
  static const std::pair<G4int, const G4double*> Z14[N14]={Z14N14, Z14N15, Z14N16};
  //==> P (Z=15)
  static const G4int N15=1;
  static const G4double pZ15N16[7]={6.E-7, 3., 8.2E-4, 1.4E-9, .03, 8.E-8, .0059};
  static const std::pair<G4int, const G4double*> Z15N16(16,pZ15N16);
  static const std::pair<G4int, const G4double*> Z15[N15]={Z15N16};
  //==> S (Z=16)
  static const G4int N16=4;
  static const G4double pZ16N16[7]={6.E-7, 3., 1.9E-4, 5.E-9, .03, 6.E-8, .0013};
  static const std::pair<G4int, const G4double*> Z16N16(16,pZ16N16);
  static const G4double pZ16N17[7]={2.4E-6, 3., .0023, 2.E-9, .03, 6.5E-8, .004};
  static const std::pair<G4int, const G4double*> Z16N17(17,pZ16N17);
  static const G4double pZ16N18[7]={2.4E-6, 1.6, .0031, 1.4E-9, .03, 4.E-08, .0028};
  static const std::pair<G4int, const G4double*> Z16N18(18,pZ16N18);
  static const G4double pZ16N20[7]={2.4E-6, 3.1, .0017, 2.5E-9, .03, 5.E-08, .0029};
  static const std::pair<G4int, const G4double*> Z16N20(20,pZ16N20);
  static const std::pair<G4int, const G4double*> Z16[N16]={Z16N16, Z16N17, Z16N18, Z16N20};
  //==> Cl(Z=17)
  static const G4int N17=2;
  static const G4double pZ17N18[7]={1.2E-7, .04, .062, 3.E-12, 3.E-02, 3.E-08, .027};
  static const std::pair<G4int, const G4double*> Z17N18(18,pZ17N18);
  static const G4double pZ17N20[7]={1.2E-7, 2., 6.8E-5, 2.7E-9, .03, 4.E-8, .0015};
  static const std::pair<G4int, const G4double*> Z17N20(20,pZ17N20);
  static const std::pair<G4int, const G4double*> Z17[N17]={Z17N18, Z17N20};
  //==> Ar(Z=18)
  static const G4int N18=3;
  static const G4double pZ18N18[7]={1.2E-7, .52, .017, 1.1E-11, .03, 3.E-8, .095};
  static const std::pair<G4int, const G4double*> Z18N18(18,pZ18N18);
  static const G4double pZ18N20[7]={1.2E-07, .09, .012, 1.8E-11, .03, 3.E-8, .011};
  static const std::pair<G4int, const G4double*> Z18N20(20,pZ18N20);
  static const G4double pZ18N22[7]={1.2E-7, .65, 1.2E-4, 1.5E-9, .03, 5.E-8, 8.E-4};
  static const std::pair<G4int, const G4double*> Z18N22(22,pZ18N22);
  static const std::pair<G4int, const G4double*> Z18[N18]={Z18N18, Z18N20, Z18N22};
  //==> K (Z=19)
  static const G4int N19=3;
  static const G4double pZ19N20[7]={1.2E-7, 1.3, 1.9E-4, .9E-9, .04, 5.5E-8, .0026};
  static const std::pair<G4int, const G4double*> Z19N20(20,pZ19N20);
  static const G4double pZ19N21[7]={1.6E-7, 1.2, 3.7E-4, .8E-9, .04, 6.5E-8, .0034};
  static const std::pair<G4int, const G4double*> Z19N21(21,pZ19N21);
  static const G4double pZ19N22[7]={6.E-8, 1.3, 1.2E-4, .9E-9, .04, 6.E-8, .0031};
  static const std::pair<G4int, const G4double*> Z19N22(22,pZ19N22);
  static const std::pair<G4int, const G4double*> Z19[N19]={Z19N20, Z19N21, Z19N22};
  //==> Ca(Z=20)
  static const G4int N20=6;
  static const G4double pZ20N20[7]={2.4E-7, 3.4, 2.1E-4, 1.5E-9, .035, 6.E-8, .0037};
  static const std::pair<G4int, const G4double*> Z20N20(20,pZ20N20);
  static const G4double pZ20N22[7]={6.E-8, 2.7, 2.7E-5, 3.E-9, .035, 6.E-8, .0014};
  static const std::pair<G4int, const G4double*> Z20N22(22,pZ20N22);
  static const G4double pZ20N23[7]={1.5E-8, 1.8, 3.4E-5, .6E-9, .04, 6.E-8, .0049};
  static const std::pair<G4int, const G4double*> Z20N23(23,pZ20N23);
  static const G4double pZ20N24[7]={3.E-6, 5., .002, 2.E-9, .03, 7.E-8, .0038};
  static const std::pair<G4int, const G4double*> Z20N24(24,pZ20N24);
  static const G4double pZ20N26[7]={1.7E-5, 18., .0027, 1.E-8, 2.E-7, 7.E-17, .0047};
  static const std::pair<G4int, const G4double*> Z20N26(26,pZ20N26);
  static const G4double pZ20N28[7]={7.6E-6, .4, .07, .13E-9, .05, 4.E-8, .0042};
  static const std::pair<G4int, const G4double*> Z20N28(28,pZ20N28);
  static const std::pair<G4int, const G4double*> Z20[N20]={Z20N20, Z20N22, Z20N23,
                                                           Z20N24, Z20N26, Z20N28};
  //==> Sc(Z=21)
  static const G4int N21=1;
  static const G4double pZ21N24[7]={3.6E-9, 1.5, 5.2E-5, .1E-9, .05, 1.E-7, .025};
  static const std::pair<G4int, const G4double*> Z21N24(24,pZ21N24);
  static const std::pair<G4int, const G4double*> Z21[N21]={Z21N24};
  //==> Ti(Z=22)
  static const G4int N22=5;
  static const G4double pZ22N24[7]={2.8E-8, 1.8, 5.6E-5, .6E-9, .05, 8.E-8, .0042};
  static const std::pair<G4int, const G4double*> Z22N24(24,pZ22N24);
  static const G4double pZ22N25[7]={3.1E-9, 1.6, 6.E-6, .8E-9, .04, 8.E-8, .0036};
  static const std::pair<G4int, const G4double*> Z22N25(25,pZ22N25);
  static const G4double pZ22N26[7]={3.E-9, 4., 3.2E-6, 1.4E-9, .05, 2.E-7, .0048};
  static const std::pair<G4int, const G4double*> Z22N26(26,pZ22N26);
  static const G4double pZ22N27[7]={1.E-8, 2., 3.4E-6, 4.5E-9, .05, 8.E-8, 7.7E-4};
  static const std::pair<G4int, const G4double*> Z22N27(27,pZ22N27);
  static const G4double pZ22N28[7]={4.E-7, 4., 3.7E-4, 1.E-09, .05, 1.E-7, .0041};
  static const std::pair<G4int, const G4double*> Z22N28(28,pZ22N28);
  static const std::pair<G4int, const G4double*> Z22[N22]={Z22N24, Z22N25, Z22N26,
                                                         Z22N27, Z22N28};
  //==> V (Z=23) *** Only nat *** (v50=v51=v_nat)
  static const G4int N23=2;
  static const G4double pZ23N27[7]={.3E-9, 2., 8.E-7, .55E-9, .07, 1.7E-7, .0055};
  static const std::pair<G4int, const G4double*> Z23N27(27,pZ23N27);
  static const G4double pZ23N28[7]={.3E-9, 2., 8.E-7, .55E-9, .07, 1.7E-7, .0055};
  static const std::pair<G4int, const G4double*> Z23N28(28,pZ23N28);
  static const std::pair<G4int, const G4double*> Z23[N23]={Z23N27, Z23N28};
  //==> Cr(Z=24)
  static const G4int N24=4;
  static const G4double pZ24N26[7]={1.2E-9, 2.8, 1.E-6, 1.7E-9, .07, 1.7E-7, .0026};
  static const std::pair<G4int, const G4double*> Z24N26(26,pZ24N26);
  static const G4double pZ24N28[7]={4.4E-6, 11., .0012, 5.E-9, .04, 3.E-7, .0032};
  static const std::pair<G4int, const G4double*> Z24N28(28,pZ24N28);
  static const G4double pZ24N29[7]={1.8E-9, 2.4, 6.3E-6, .5E-9, .07, 2.E-7, .0085};
  static const std::pair<G4int, const G4double*> Z24N29(29,pZ24N29);
  static const G4double pZ24N30[7]={4.8E-8, 2.8, 4.4E-5, 1.4E-9, .07, 2.E-7, .0027};
  static const std::pair<G4int, const G4double*> Z24N30(30,pZ24N30);
  static const std::pair<G4int, const G4double*> Z24[N24]={Z24N26, Z24N28, Z24N29, Z24N30};
  //==> Mn(Z=25)
  static const G4int N25=1;
  static const G4double pZ25N30[7]={6.5E-11, 1.4, 1.E-7, .8E-9, .07, 1.7E-7, .0022};
  static const std::pair<G4int, const G4double*> Z25N30(30,pZ25N30);
  static const std::pair<G4int, const G4double*> Z25[N25]={Z25N30};
  //==> Fe(Z=26)
  static const G4int N26=4;
  static const G4double pZ26N28[7]={3.9E-8, 5., 1.7E-5, 3.E-9, .07, 3.E-7, .0023};
  static const std::pair<G4int, const G4double*> Z26N28(28,pZ26N28);
  static const G4double pZ26N30[7]={5.E-9, .4, 1.5E-4, 4.E-11, .1, 3.E-7, .012};
  static const std::pair<G4int, const G4double*> Z26N30(30,pZ26N30);
  static const G4double pZ26N31[7]={.5E-9, .5, 2.6E-6, .3E-9, .11, 5.E-7, .0027};
  static const std::pair<G4int, const G4double*> Z26N31(31,pZ26N31);
  static const G4double pZ26N32[7]={1.E-7, 3.1, 1.E-4, 1.3E-9, .11, 5.E-7, .0031};
  static const std::pair<G4int, const G4double*> Z26N32(32,pZ26N32);
  static const std::pair<G4int, const G4double*> Z26[N26]={Z26N28, Z26N30, Z26N31, Z26N32};
  //==> Co(Z=27)
  static const G4int N27=2;
  static const G4double pZ27N31[7]={4.E-7, 3., .004, 0., .11, 4.5E-7, .07};
  static const std::pair<G4int, const G4double*> Z27N31(31,pZ27N31);
  static const G4double pZ27N32[7]={4.E-7, 5., 5.E-4, 1.2E-9, .13, 6.E-7, .006};
  static const std::pair<G4int, const G4double*> Z27N32(32,pZ27N32);
  static const std::pair<G4int, const G4double*> Z27[N27]={Z27N32, Z27N31};
  //==> Ni(Z=28)
  static const G4int N28=6;
  static const G4double pZ28N30[7]={1.E-7, 2.5, .001, .14E-9, .13, 6.E-7, .025};
  static const std::pair<G4int, const G4double*> Z28N30(30,pZ28N30);
  static const G4double pZ28N31[7]={1.E-7, 19., 1.2E-5, 1.E-8, 4.E-12, 3.E-22, .0024};
  static const std::pair<G4int, const G4double*> Z28N31(31,pZ28N31);
  static const G4double pZ28N32[7]={1.E-8, 2.5, 3.9E-6, 3.5E-9, .13, 6.E-7, .001};
  static const std::pair<G4int, const G4double*> Z28N32(32,pZ28N32);
  static const G4double pZ28N33[7]={5.E-9, 2.6, 1.5E-5, .42E-9, .13, 7.E-7, .008};
  static const std::pair<G4int, const G4double*> Z28N33(33,pZ28N33);
  static const G4double pZ28N34[7]={.24E-9, 2., 1.2E-6, .25E-9, .13, 6.E-7, .0094};
  static const std::pair<G4int, const G4double*> Z28N34(34,pZ28N34);
  static const G4double pZ28N36[7]={1.E-8, 3., 5.5E-8, 2.8E-7, .12, 6.E-7, 1.6E-5};
  static const std::pair<G4int, const G4double*> Z28N36(36,pZ28N36);
  static const std::pair<G4int, const G4double*> Z28[N28]={Z28N30, Z28N31, Z28N32, Z28N33,
                                                           Z28N34, Z28N36};
  //==> Cu(Z=29)
  static const G4int N29=2;
  static const G4double pZ29N34[7]={1.1E-7, 3.5, 1.6E-4, .9E-9, .13, 7.E-7, .005};
  static const std::pair<G4int, const G4double*> Z29N34(34,pZ29N34);
  static const G4double pZ29N36[7]={1.1E-7, 3.5, 4.3E-4, .3E-9, .13, 7.E-7, .013};
  static const std::pair<G4int, const G4double*> Z29N36(36,pZ29N36);
  static const std::pair<G4int, const G4double*> Z29[N29]={Z29N34, Z29N36};
  //==> Zn(Z=30) *** Only nat *** (zn64=zn66=zn67=zn68=zn70=zn_nat)
  static const G4int N30=5;
  static const G4double pZ30N34[7]={1.1E-7, 4., 1.2E-4, 1.2E-9, .17, 7.E-7, .004};
  static const std::pair<G4int, const G4double*> Z30N34(34,pZ30N34);
  static const G4double pZ30N36[7]={1.1E-7, 4., 1.2E-4, 1.2E-9, .17, 7.E-7, .004};
  static const std::pair<G4int, const G4double*> Z30N36(36,pZ30N36);
  static const G4double pZ30N37[7]={1.1E-7, 4., 1.2E-4, 1.2E-9, .17, 7.E-7, .004};
  static const std::pair<G4int, const G4double*> Z30N37(37,pZ30N37);
  static const G4double pZ30N38[7]={1.1E-7, 4., 1.2E-4, 1.2E-9, .17, 7.E-7, .004};
  static const std::pair<G4int, const G4double*> Z30N38(38,pZ30N38);
  static const G4double pZ30N40[7]={1.1E-7, 4., 1.2E-4, 1.2E-9, .17, 7.E-7, .004};
  static const std::pair<G4int, const G4double*> Z30N40(40,pZ30N40);
  static const std::pair<G4int, const G4double*> Z30[N30]={Z30N34, Z30N36, Z30N37,
                                                           Z30N38, Z30N40};
  //==> Ga(Z=31)
  static const G4int N31=2;
  static const G4double pZ31N38[7]={5.E-8, 3.7, 1.1E-4, .55E-9, .17, 8.4E-7, .0076};
  static const std::pair<G4int, const G4double*> Z31N38(38,pZ31N38);
  static const G4double pZ31N40[7]={1.E-8, 3.1, 1.7E-5, .7E-9, .17, 9.E-7, .0048};
  static const std::pair<G4int, const G4double*> Z31N40(40,pZ31N40);
  static const std::pair<G4int, const G4double*> Z31[N31]={Z31N38, Z31N40};
  //==> Ge(Z=32)
  static const G4int N32=5;
  static const G4double pZ32N38[7]={5.E-5, 4., .17, .35E-9, .17, 9.E-7, .013};
  static const std::pair<G4int, const G4double*> Z32N38(38,pZ32N38);
  static const G4double pZ32N40[7]={5.E-7, 4.4, .001, .6E-9, .17, 9.E-7, .008};
  static const std::pair<G4int, const G4double*> Z32N40(40,pZ32N40);
  static const G4double pZ32N41[7]={5.E-9, 3., 8.E-6, .7E-9, .17, 1.E-6, .0043};
  static const std::pair<G4int, const G4double*> Z32N41(41,pZ32N41);
  static const G4double pZ32N42[7]={1.E-7, 4.2, 1.7E-4, .7E-9, .17, 1.E-6, .0065};
  static const std::pair<G4int, const G4double*> Z32N42(42,pZ32N42);
  static const G4double pZ32N44[7]={1.E-6, 4.6, .0018, .6E-9, .17, 1.E-6, .0073};
  static const std::pair<G4int, const G4double*> Z32N44(44,pZ32N44);
  static const std::pair<G4int, const G4double*> Z32[N32]={Z32N38, Z32N40, Z32N41,
                                                           Z32N42, Z32N44};
  //==> As(Z=33)
  static const G4int N33=2;
  static const G4double pZ33N41[7]={1.E-8, 3.4, 1.5E-5, .72E-9, .17, 1.E-6, .0045};
  static const std::pair<G4int, const G4double*> Z33N41(41,pZ33N41);
  static const G4double pZ33N42[7]={1.E-8, 4.1, 1.3E-5, .75E-9, .2, 1.2E-6, .0048};
  static const std::pair<G4int, const G4double*> Z33N42(42,pZ33N42);
  static const std::pair<G4int, const G4double*> Z33[N33]={Z33N41, Z33N42};
  //==> Se(Z=34)
  static const G4int N34=7;
  static const G4double pZ34N40[7]={6.E-8, 7.2, 6.E-5, 1.E-9, .32, 2.E-6, .0063};
  static const std::pair<G4int, const G4double*> Z34N40(40,pZ34N40);
  static const G4double pZ34N42[7]={4.E-5, 7.4, .1, .43E-9, .34, 2.1E-6, .016};
  static const std::pair<G4int, const G4double*> Z34N42(42,pZ34N42);
  static const G4double pZ34N43[7]={1.E-7, 6.2, 1.4E-4, .9E-9, .34, 2.1E-6, .0075};
  static const std::pair<G4int, const G4double*> Z34N43(43,pZ34N43);
  static const G4double pZ34N44[7]={1.E-7, 6.6, 1.3E-4, .9E-9, .34, 2.1E-6, .0075};
  static const std::pair<G4int, const G4double*> Z34N44(44,pZ34N44);
  static const G4double pZ34N45[7]={5.E-8, 6.6, 4.8E-5, 1.2E-9, .4, 2.6E-6, .0055};
  static const std::pair<G4int, const G4double*> Z34N45(45,pZ34N45);
  static const G4double pZ34N46[7]={2.E-7, 7.7, 1.3E-4, 1.7E-9, .34, 2.1E-6, .0043};
  static const std::pair<G4int, const G4double*> Z34N46(46,pZ34N46);
  static const G4double pZ34N48[7]={2.E-7, 8.3, 1.2E-4, 1.7E-9, .34, 2.1E-6, .0043};
  static const std::pair<G4int, const G4double*> Z34N48(48,pZ34N48);
  static const std::pair<G4int, const G4double*> Z34[N34]={Z34N40, Z34N42, Z34N43, Z34N44,
                                                           Z34N45, Z34N46, Z34N48};
  //==> Br(Z=35)
  static const G4int N35=2;
  static const G4double pZ35N44[7]={5.E-8, 6., 2.8E-5, 2.E-9, .34, 2.1E-6, .0028};
  static const std::pair<G4int, const G4double*> Z35N44(44,pZ35N44);
  static const G4double pZ35N46[7]={4.E-8, 6.2, 3.7E-5, 1.1E-9, .34, 2.1E-6, .0049};
  static const std::pair<G4int, const G4double*> Z35N46(46,pZ35N46);
  static const std::pair<G4int, const G4double*> Z35[N35]={Z35N44, Z35N46};
  //==> Kr(Z=36)
  static const G4int N36=7;
  static const G4double pZ36N42[7]={1.6E-7, 6.8, 2.E-4, .8E-9, .35, 2.1E-6, .0076};
  static const std::pair<G4int, const G4double*> Z36N42(42,pZ36N42);
  static const G4double pZ36N44[7]={1.6E-7, 7.3, 1.6E-4, 1.E-9, .35, 2.1E-6, .0062};
  static const std::pair<G4int, const G4double*> Z36N44(44,pZ36N44);
  static const G4double pZ36N46[7]={1.6E-7, 7.3, 3.3E-4, .7E-9, .35, 2.1E-6, .013};
  static const std::pair<G4int, const G4double*> Z36N46(46,pZ36N46);
  static const G4double pZ36N47[7]={1.6E-6, 7.3, .003, .6E-9, .35, 2.1E-6, .011};
  static const std::pair<G4int, const G4double*> Z36N47(47,pZ36N47);
  static const G4double pZ36N48[7]={1.6E-7, 7.8, 7.6E-5, 2.E-9, .35, 2.1E-6, .0031};
  static const std::pair<G4int, const G4double*> Z36N48(48,pZ36N48);
  static const G4double pZ36N49[7]={6.E-7, 8., 4.8E-4, 1.4E-9, .27, 2.1E-6, .0053};
  static const std::pair<G4int, const G4double*> Z36N49(49,pZ36N49);
  static const G4double pZ36N50[7]={4.E-7, 8.1, 2.7E-4, 1.6E-9, .35, 2.1E-6, .0045};
  static const std::pair<G4int, const G4double*> Z36N50(50,pZ36N50);
  static const std::pair<G4int, const G4double*> Z36[N36]={Z36N42, Z36N44, Z36N46,
                                                           Z36N47, Z36N48, Z36N49, Z36N50};
  //==> Rb(Z=37)
  static const G4int N37=3;
  static const G4double pZ37N48[7]={1.6E-7, 7.2, 1.4E-4, 1.2E-9, .35, 2.1E-6, .0052};
  static const std::pair<G4int, const G4double*> Z37N48(48,pZ37N48);
  static const G4double pZ37N49[7]={8.E-8, 7.1, 4.7E-5, 1.6E-9, .27, 2.1E-6, .0034};
  static const std::pair<G4int, const G4double*> Z37N49(49,pZ37N49);
  static const G4double pZ37N50[7]={1.E-7, 8., 5.5E-5, 1.9E-9, .27, 1.5E-6, .0036};
  static const std::pair<G4int, const G4double*> Z37N50(50,pZ37N50);
  static const std::pair<G4int, const G4double*> Z37[N37]={Z37N48, Z37N49, Z37N50};
  //==> Sr(Z=38)
  static const G4int N38=6;
  static const G4double pZ38N46[7]={8.E-8, 7.3, 6.E-5, 1.3E-9, .27, 2.E-6, .0045};
  static const std::pair<G4int, const G4double*> Z38N46(46,pZ38N46);
  static const G4double pZ38N48[7]={8.E-8, 9.7, 2.3E-5, 3.5E-9, .4, 2.7E-6, .0023};
  static const std::pair<G4int, const G4double*> Z38N48(48,pZ38N48);
  static const G4double pZ38N49[7]={2.6E-7, 9.5, 1.9E-4, 1.5E-9, .4, 2.7E-6, .0057};
  static const std::pair<G4int, const G4double*> Z38N49(49,pZ38N49);
  static const G4double pZ38N50[7]={2.6E-7, 9.5, 2.E-4, 1.4E-9, .37, 3.2E-6, .0059};
  static const std::pair<G4int, const G4double*> Z38N50(50,pZ38N50);
  static const G4double pZ38N51[7]={1.3E-7, 9.9, 7.5E-5, 1.7E-9, .37, 3.2E-6, .0046};
  static const std::pair<G4int, const G4double*> Z38N51(51,pZ38N51);
  static const G4double pZ38N52[7]={2.6E-7, 9.6, 1.6E-4, 1.8E-9, .37, 2.7E06, .0047};
  static const std::pair<G4int, const G4double*> Z38N52(52,pZ38N52);
  static const std::pair<G4int, const G4double*> Z38[N38]={Z38N46, Z38N48, Z38N49, Z38N50,
                                                           Z38N51, Z38N52};
  //==> Y (Z=39)
  static const G4int N39=3;
  static const G4double pZ39N50[7]={2.6E-7, 9.9, 2.E-4, 1.1E-9, .37, 3.2E-6, .0062};
  static const std::pair<G4int, const G4double*> Z39N50(50,pZ39N50);
  static const G4double pZ39N51[7]={2.7E-5, 20., .013, 2.2E-9, .37, 1.9E-6, .0078};
  static const std::pair<G4int, const G4double*> Z39N51(51,pZ39N51);
  static const G4double pZ39N52[7]={2.E-7, 9.6, 1.2E-4, 1.7E-9, .37, 3.E-6, .0046};
  static const std::pair<G4int, const G4double*> Z39N52(52,pZ39N52);
  static const std::pair<G4int, const G4double*> Z39[N39]={Z39N50, Z39N51, Z39N52};
  //==> Zr(Z=40)
  static const G4int N40=7;
  static const G4double pZ40N50[7]={1.E-7, 9., 6.2E-5, 1.7E-9, .3, 3.E-6, .0044};
  static const std::pair<G4int, const G4double*> Z40N50(50,pZ40N50);
  static const G4double pZ40N51[7]={5.E-7, 9.8, 5.E-4, 1.E-9, .25, 2.1E-6, .0079};
  static const std::pair<G4int, const G4double*> Z40N51(51,pZ40N51);
  static const G4double pZ40N52[7]={3.E-7, 9.6, 2.2E-4, 1.2E-9, .25, 2.2E-6, .0056};
  static const std::pair<G4int, const G4double*> Z40N52(52,pZ40N52);
  static const G4double pZ40N53[7]={2.E-7, 9.6, 1.2E-4, 1.6E-9, .38, 2.9E-6, .0046};
  static const std::pair<G4int, const G4double*> Z40N53(53,pZ40N53);
  static const G4double pZ40N54[7]={2.E-7, 9.6, 1.8E-4, 1.1E-9, .25, 2.1E-6, .0067};
  static const std::pair<G4int, const G4double*> Z40N54(54,pZ40N54);
  static const G4double pZ40N55[7]={1.8E-7, 9.4, 1.1E-4, 1.6E-9, .33, 2.3E-6, .0045};
  static const std::pair<G4int, const G4double*> Z40N55(55,pZ40N55);
  static const G4double pZ40N56[7]={1.8E-7, 9.4, 1.1E-4, 1.6E-9, .2, 1.5E-6, .0045};
  static const std::pair<G4int, const G4double*> Z40N56(56,pZ40N56);
  static const std::pair<G4int, const G4double*> Z40[N40]={Z40N50, Z40N51, Z40N52, Z40N53,
                                                           Z40N54, Z40N55, Z40N56};
  //==> Nb(Z=41)
  static const G4int N41=3;
  static const G4double pZ41N52[7]={2.6E-7, 8.3, 2.E-4, 1.2E-9, .4, 4.E-6, .0051};
  static const std::pair<G4int, const G4double*> Z41N52(52,pZ41N52);
  //static const G4double pZ41N53[7]={2.E-7, 8.3, 1.6E-4, 1.4E-9, .35, 2.5E-6, .0051};
  //static const std::pair<G4int, const G4double*> Z41N53(53,pZ41N53);
  //static const G4double pZ41N54[7]={1.5E-7, 8.6, 1.E-4, 1.5E-9, .35, 2.5E-6, .0045};
  //static const std::pair<G4int, const G4double*> Z41N54(54,pZ41N54);
  static const std::pair<G4int, const G4double*> Z41[N41]={Z41N52, Z41N52, Z41N52};
  //==> Mo(Z=42)
  static const G4int N42=8;
  static const G4double pZ42N50[7]={2.E-7, 10., 1.1E-4, 1.8E-9, .3, 2.7E-6, .0044};
  static const std::pair<G4int, const G4double*> Z42N50(50,pZ42N50);
  static const G4double pZ42N52[7]={2.1E-7, 10., 1.2E-4, 1.7E-9, .3, 2.8E-6, .0046};
  static const std::pair<G4int, const G4double*> Z42N52(52,pZ42N52);
  static const G4double pZ42N53[7]={3.E-7, 10., 1.9E-4, 1.5E-9, .29, 3.E-6, .005};
  static const std::pair<G4int, const G4double*> Z42N53(53,pZ42N53);
  static const G4double pZ42N54[7]={1.5E-7, 10., 7.1E-5, 2.1E-9, .29, 2.9E-6, .0037};
  static const std::pair<G4int, const G4double*> Z42N54(54,pZ42N54);
  static const G4double pZ42N55[7]={1.9E-7, 9., 1.4E-4, 1.3E-9, .29, 2.8E-6, .0052};
  static const std::pair<G4int, const G4double*> Z42N55(55,pZ42N55);
  static const G4double pZ42N56[7]={1.9E-7, 9.9, 1.1E-4, 1.8E-9, .29, 2.4E-6, .0044};
  static const std::pair<G4int, const G4double*> Z42N56(56,pZ42N56);
  static const G4double pZ42N57[7]={1.4E-7, 8., 1.E-4, 1.4E-9, .34, 2.5E-6, .0044};
  static const std::pair<G4int, const G4double*> Z42N57(57,pZ42N57);
  static const G4double pZ42N58[7]={1.8E-7, 9.5, 1.E-4, 1.7E-9, .27, 2.2E-6, .0041};
  static const std::pair<G4int, const G4double*> Z42N58(58,pZ42N58);
  static const std::pair<G4int, const G4double*> Z42[N42]={Z42N50, Z42N52, Z42N53, Z42N54,
                                                           Z42N55, Z42N56, Z42N57, Z42N58};
  //==> Tc(Z=43)
  static const G4int N43=1;
  static const G4double pZ43N56[7]={1.E-7, 8., 7.2E-5, 1.4E-9, .24, 2.5E-6, .0044};
  static const std::pair<G4int, const G4double*> Z43N56(56,pZ43N56);
  static const std::pair<G4int, const G4double*> Z43[N43]={Z43N56};
  //==> Ru(Z=44)
  static const G4int N44=9;
  static const G4double pZ44N52[7]={1.9E-7, 10., 1.E-4, 2.1E-9, .4, 4.E-6, .004};
  static const std::pair<G4int, const G4double*> Z44N52(52,pZ44N52);
  static const G4double pZ44N54[7]={1.5E-7, 10., 7.7E-5, 2.1E-9, .29, 2.9E-6, .004};
  static const std::pair<G4int, const G4double*> Z44N54(54,pZ44N54);
  static const G4double pZ44N55[7]={1.8E-7, 10., 6.6E-5, 2.6E-9, .47, 4.6E-6, .0028};
  static const std::pair<G4int, const G4double*> Z44N55(55,pZ44N55);
  static const G4double pZ44N56[7]={1.8E-6, 10., .0017, 1.1E-9, .47, 6.E-6, .0073};
  static const std::pair<G4int, const G4double*> Z44N56(56,pZ44N56);
  static const G4double pZ44N57[7]={1.8E-7, 7.8, 1.3E-4, 1.4E-9, .42, 5.E-6, .0043};
  static const std::pair<G4int, const G4double*> Z44N57(57,pZ44N57);
  static const G4double pZ44N58[7]={1.7E-6, 9.8, .0015, 1.2E-9, .32, 4.E-6, .0065};
  static const std::pair<G4int, const G4double*> Z44N58(58,pZ44N58);
  static const G4double pZ44N59[7]={3.3E-7, 8.7, 1.9E-4, 1.6E-9, .32, 3.8E-6, .0038};
  static const std::pair<G4int, const G4double*> Z44N59(59,pZ44N59);
  static const G4double pZ44N60[7]={3.E-7, 8.7, 1.8E-4, 1.6E-9, .3, 3.2E-6, .004};
  static const std::pair<G4int, const G4double*> Z44N60(60,pZ44N60);
  static const G4double pZ44N61[7]={3.E-7, 8.8, 1.4E-4, 1.9E-9, .3, 3.2E-6, .003};
  static const std::pair<G4int, const G4double*> Z44N61(61,pZ44N61);
  static const std::pair<G4int, const G4double*> Z44[N44]={Z44N52, Z44N54, Z44N55, Z44N56,
                                                           Z44N57, Z44N58, Z44N59, Z44N60,
                                                           Z44N61};
  //==> Rh(Z=45)
  static const G4int N45=2;
  static const G4double pZ45N58[7]={8.E-8, 8.7, 4.E-5, 1.8E-9, .29, 2.9E-6, .0033};
  static const std::pair<G4int, const G4double*> Z45N58(58,pZ45N58);
  static const G4double pZ45N60[7]={8.E-8, 8.7, .09, 1.3E-12, .29, 2.9E-6, 7.};
  static const std::pair<G4int, const G4double*> Z45N60(60,pZ45N60);
  static const std::pair<G4int, const G4double*> Z45[N45]={Z45N58, Z45N60};
  //==> Pd(Z=46)
  static const G4int N46=7;
  static const G4double pZ46N56[7]={2.E-7, 9.9, 1.2E-4, 1.5E-9, .35, 3.3E-6, .0045};
  static const std::pair<G4int, const G4double*> Z46N56(56,pZ46N56);
  static const G4double pZ46N58[7]={2.E-7, 9.9, 9.5E-5, 1.8E-9, .4, 4.E-6, .0036};
  static const std::pair<G4int, const G4double*> Z46N58(58,pZ46N58);
  static const G4double pZ46N59[7]={5.6E-7, 9., 4.6E-4, 1.2E-9, .5, 4.8E-6, .0056};
  static const std::pair<G4int, const G4double*> Z46N59(59,pZ46N59);
  static const G4double pZ46N60[7]={2.4E-7, 9.2, 1.2E-4, 1.8E-9, .47, 4.6E-6, .0035};
  static const std::pair<G4int, const G4double*> Z46N60(60,pZ46N60);
  static const G4double pZ46N61[7]={1.2E-7, 9.2, 4.4E-5, 2.8E-9, .5, 4.3E-6, .0025};
  static const std::pair<G4int, const G4double*> Z46N61(61,pZ46N61);
  static const G4double pZ46N62[7]={1.2E-7, 9.2, 3.2E-5, 3.4E-9, .48, 4.5E-6, .0018};
  static const std::pair<G4int, const G4double*> Z46N62(62,pZ46N62);
  static const G4double pZ46N64[7]={4.E-7, 9.1, 2.5E-4, 1.5E-9, .48, 4.7E-6, .0042};
  static const std::pair<G4int, const G4double*> Z46N64(64,pZ46N64);
  static const std::pair<G4int, const G4double*> Z46[N46]={Z46N56, Z46N58, Z46N59,
                                                           Z46N60, Z46N61, Z46N62, Z46N64};
  //==> Ag(Z=47)
  static const G4int N47=4;
  static const G4double pZ47N60[7]={1.4E-6, 9.7, .0011, 1.4E-9, .55, 5.E-6, .0056};
  static const std::pair<G4int, const G4double*> Z47N60(60,pZ47N60);
  static const G4double pZ47N62[7]={3.E-8, 8.7, 8.5E-6, 3.5E-9, .6, 5.6E-6, .0018};
  static const std::pair<G4int, const G4double*> Z47N62(62,pZ47N62);
  static const G4double pZ47N63[7]={3.E-6, 9.5, .002, 1.5E-9, .58, 5.E-6, .0047};
  static const std::pair<G4int, const G4double*> Z47N63(63,pZ47N63);
  static const G4double pZ47N64[7]={1.5E-7, 9., 9.E-5, 1.7E-9, .58, 5.6E-6, .0039};
  static const std::pair<G4int, const G4double*> Z47N64(64,pZ47N64);
  static const std::pair<G4int, const G4double*> Z47[N47]={Z47N60, Z47N62, Z47N63, Z47N64};
  //==> Cd(Z=48)
  static const G4int N48=9;
  static const G4double pZ48N58[7]={2.9E-7, 10., 1.3E-4, 1.9E-9, .4, 3.8E-6, .0034};
  static const std::pair<G4int, const G4double*> Z48N58(58,pZ48N58);
  static const G4double pZ48N60[7]={2.3E-7, 10., 8.2E-5, 2.5E-9, .5, 4.7E-6, .0026};
  static const std::pair<G4int, const G4double*> Z48N60(60,pZ48N60);
  static const G4double pZ48N62[7]={2.3E-7, 10., 9.9E-5, 2.5E-9, .5, 4.7E-6, .0031};
  static const std::pair<G4int, const G4double*> Z48N62(62,pZ48N62);
  static const G4double pZ48N63[7]={8.4E-7, 11., 4.3E-4, 1.8E-9, .5, 4.5E-6, .0042};
  static const std::pair<G4int, const G4double*> Z48N63(63,pZ48N63);
  static const G4double pZ48N64[7]={4.E-7, 11., 1.8E-4, 1.8E-9, .5, 4.6E-6, .0036};
  static const std::pair<G4int, const G4double*> Z48N64(64,pZ48N64);
  static const G4double pZ48N65[7]={1.6E-6, 12., .001, 0., .5, 4.6E-6, .013};
  static const std::pair<G4int, const G4double*> Z48N65(65,pZ48N65);
  static const G4double pZ48N66[7]={3.E-7, 11., 1.2E-4, 1.9E-9, .5, 4.6E-6, .0031};
  static const std::pair<G4int, const G4double*> Z48N66(66,pZ48N66);
  static const G4double pZ48N67[7]={3.8E-7, 11., 1.7E-4, 2.E-9, .6, 6.6E-6, .0035};
  static const std::pair<G4int, const G4double*> Z48N67(67,pZ48N67);
  static const G4double pZ48N68[7]={6.E-7, 11., 3.3E-4, 1.9E-9, .5, 4.6E-6, .0043};
  static const std::pair<G4int, const G4double*> Z48N68(68,pZ48N68);
  static const std::pair<G4int, const G4double*> Z48[N48]={Z48N58, Z48N60, Z48N62, Z48N63,
                                                           Z48N64, Z48N65, Z48N66, Z48N67,
                                                           Z48N68};
  //==> In(Z=49)
  static const G4int N49=2;
  static const G4double pZ49N64[7]={2.7E-7, 12., 8.1E-5, 2.7E-9, .5, 5.E-6, .0026};
  static const std::pair<G4int, const G4double*> Z49N64(64,pZ49N64);
  static const G4double pZ49N66[7]={2.7E-7, 12., 5.5E-5, 4.E-9, .5, 5.E-6, .0018};
  static const std::pair<G4int, const G4double*> Z49N66(66,pZ49N66);
  static const std::pair<G4int, const G4double*> Z49[N49]={Z49N64, Z49N66};
  //==> Sn(Z=50)
  static const G4int N50=14;
  static const G4double pZ50N62[7]={4.E-7, 11., 1.6E-4, 2.2E-9, .5, 4.5E-6, .0032};
  static const std::pair<G4int, const G4double*> Z50N62(62,pZ50N62);
  static const G4double pZ50N63[7]={4.1E-7, 11., 1.6E-4, 2.4E-9, .54, 6.E-6, .0031};
  static const std::pair<G4int, const G4double*> Z50N63(63,pZ50N63);
  static const G4double pZ50N64[7]={5.E-7, 12., 1.9E-4, 2.2E-9, .5, 4.4E-6, .0032};
  static const std::pair<G4int, const G4double*> Z50N64(64,pZ50N64);
  static const G4double pZ50N65[7]={1.E-5, 12., .0077, 1.4E-9, .5, 5.E-6, .0066};
  static const std::pair<G4int, const G4double*> Z50N65(65,pZ50N65);
  static const G4double pZ50N66[7]={5.E-7, 12., 1.8E-4, 2.4E-9, .5, 5.E-6, .0031};
  static const std::pair<G4int, const G4double*> Z50N66(66,pZ50N66);
  static const G4double pZ50N67[7]={1.E-6, 12., 4.4E-4, 1.8E-9, .5, 5.E-6, .0037};
  static const std::pair<G4int, const G4double*> Z50N67(67,pZ50N67);
  static const G4double pZ50N68[7]={5.E-7, 12., 2.E-4, 2.4E-9, .5, 5.E-6, .0033};
  static const std::pair<G4int, const G4double*> Z50N68(68,pZ50N68);
  static const G4double pZ50N69[7]={6.E-7, 12., 2.5E-4, 2.E-9, .5, 5.E-6, .0035};
  static const std::pair<G4int, const G4double*> Z50N69(69,pZ50N69);
  static const G4double pZ50N70[7]={1.E-6, 12., 4.7E-4, 2.E-9, .5, 5.E-6, .0039};
  static const std::pair<G4int, const G4double*> Z50N70(70,pZ50N70);
  static const G4double pZ50N72[7]={1.E-6, 12., 3.7E-4, 2.2E-9, .5, 5.E-6, .0031};
  static const std::pair<G4int, const G4double*> Z50N72(72,pZ50N72);
  static const G4double pZ50N73[7]={5.E-7, 12., 1.7E-4, 2.8E-9, .5, 5.E-6, .0028};
  static const std::pair<G4int, const G4double*> Z50N73(73,pZ50N73);
  static const G4double pZ50N74[7]={5.E-7, 12., 2.E-4, 2.E-9, .5, 5.E-6, .0033};
  static const std::pair<G4int, const G4double*> Z50N74(74,pZ50N74);
  static const G4double pZ50N75[7]={5.E-7, 12., 1.9E-4, 2.8E-9, .5, 5.E-6, .003};
  static const std::pair<G4int, const G4double*> Z50N75(75,pZ50N75);
  static const G4double pZ50N76[7]={5.E-7, 12., 1.7E-4, 2.8E-9, .5, 5.E-6, .0028};
  static const std::pair<G4int, const G4double*> Z50N76(76,pZ50N76);
  static const std::pair<G4int, const G4double*> Z50[N50]={Z50N62, Z50N63, Z50N64, Z50N65,
                                                           Z50N66, Z50N67, Z50N68, Z50N69,
                                                           Z50N70, Z50N72, Z50N73, Z50N74,
                                                           Z50N75, Z50N76};
  //==> Sb(Z=51)
  static const G4int N51=5;
  static const G4double pZ51N70[7]={6.E-7, 12., 2.E-4, 2.8E-9, .5, 5.E-6, .0028};
  static const std::pair<G4int, const G4double*> Z51N70(70,pZ51N70);
  static const G4double pZ51N72[7]={6.E-7, 12., 1.9E-4, 3.E-9, .5, 5.E-6, .0025};
  static const std::pair<G4int, const G4double*> Z51N72(72,pZ51N72);
  static const G4double pZ51N73[7]={1.1E-6, 12., 3.5E-4, 2.9E-9, .5, 5.E-6, .0026};
  static const std::pair<G4int, const G4double*> Z51N73(73,pZ51N73);
  static const G4double pZ51N74[7]={5.5E-7, 12., 1.9E-4, 2.9E-9, .5, 5.E-6, .0027};
  static const std::pair<G4int, const G4double*> Z51N74(74,pZ51N74);
  static const G4double pZ51N75[7]={6.E-7, 12., 2.E-4, 2.9E-9, .5, 5.E-6, .0027};
  static const std::pair<G4int, const G4double*> Z51N75(75,pZ51N75);
  static const std::pair<G4int, const G4double*> Z51[N51]={Z51N70, Z51N72, Z51N73, Z51N74,
                                                           Z51N75};
  //==> Te(Z=52)
  static const G4int N52=11;
  static const G4double pZ52N68[7]={2.7E-7, 12., 8.4E-5, 3.2E-9, 1., 8.E-6, .0026};
  static const std::pair<G4int, const G4double*> Z52N68(68,pZ52N68);
  static const G4double pZ52N70[7]={2.7E-7, 12., 3.8E-5, 6.E-9, 1., 8.E-6, .0012};
  static const std::pair<G4int, const G4double*> Z52N70(70,pZ52N70);
  static const G4double pZ52N71[7]={2.7E-8, 12., 1.8E-6, 2.E-8, 1., 8.E-6, 4.8E-4};
  static const std::pair<G4int, const G4double*> Z52N71(71,pZ52N71);
  static const G4double pZ52N72[7]={2.6E-6, 14., .0014, 2.E-9, 1., 9.E-6, .005};
  static const std::pair<G4int, const G4double*> Z52N72(72,pZ52N72);
  static const G4double pZ52N73[7]={1.E-6, 14., 2.4E-4, 3.9E-9, 1., 9.E-6, .0022};
  static const std::pair<G4int, const G4double*> Z52N73(73,pZ52N73);
  static const G4double pZ52N74[7]={8.E-7, 14., 2.3E-4, 3.6E-9, 1.4, 1.3E-5, .0028};
  static const std::pair<G4int, const G4double*> Z52N74(74,pZ52N74);
  static const G4double pZ52N75[7]={8.E-7, 14., 2.1E-4, 3.6E-9, 1.4, 1.3E-5, .0025};
  static const std::pair<G4int, const G4double*> Z52N75(75,pZ52N75);
  static const G4double pZ52N76[7]={8.E-7, 14., 2.5E-4, 3.E-9, 1.4, 1.3E-5, .003};
  static const std::pair<G4int, const G4double*> Z52N76(76,pZ52N76);
  static const G4double pZ52N77[7]={5.E-7, 15., 1.2E-4, 4.3E-9, 1.4, 1.4E-5, .0023};
  static const std::pair<G4int, const G4double*> Z52N77(77,pZ52N77);
  static const G4double pZ52N78[7]={8.E-7, 14., 2.7E-4, 2.7E-9, 1.4, 1.3E-5, .0031};
  static const std::pair<G4int, const G4double*> Z52N78(78,pZ52N78);
  static const G4double pZ52N80[7]={4.7E-7, 14., 1.8E-4, 2.2E-9, .83, 1.E-5, .0036};
  static const std::pair<G4int, const G4double*> Z52N80(80,pZ52N80);
  static const std::pair<G4int, const G4double*> Z52[N52]={Z52N68, Z52N70, Z52N71, Z52N72,
                                                           Z52N73, Z52N74, Z52N75, Z52N76,
                                                           Z52N77, Z52N78, Z52N80};
  //==> I (Z=53)
  static const G4int N53=5;
  static const G4double pZ53N74[7]={9.4E-7, 14., 2.5E-4, 3.E-9, .7, 7.3E-6, .0025};
  static const std::pair<G4int, const G4double*> Z53N74(74,pZ53N74);
  static const G4double pZ53N76[7]={2.1E-5, 14., .015, 1.1E-9, 1.1, 1.E-5, .007};
  static const std::pair<G4int, const G4double*> Z53N76(76,pZ53N76);
  static const G4double pZ53N77[7]={1.1E-6, 14., 2.4E-4, 3.3E-9, .9, 1.E-5, .0021};
  static const std::pair<G4int, const G4double*> Z53N77(77,pZ53N77);
  static const G4double pZ53N78[7]={5.5E-7, 14., 1.5E-4, 3.7E-9, 1.2, 1.2E-5, .0024};
  static const std::pair<G4int, const G4double*> Z53N78(78,pZ53N78);
  static const G4double pZ53N82[7]={3.2E-6, 14., .0017, 1.8E-9, .8, 8.E-6, .0024};
  static const std::pair<G4int, const G4double*> Z53N82(82,pZ53N82);
  static const std::pair<G4int, const G4double*> Z53[N53]={Z53N74, Z53N76, Z53N77, Z53N78,
                                                           Z53N82};
  //==> Xe(Z=54)
  static const G4int N54=12;
  static const G4double pZ54N69[7]={3.E-6, 14., 8.E-4, 3.7E-9, .9, 1.1E-5, .15};
  static const std::pair<G4int, const G4double*> Z54N69(69,pZ54N69);
  static const G4double pZ54N70[7]={1.5E-7, 14., 1.4E-6, 9.E-8, .7, 8.E-6, 9.5E-5};
  static const std::pair<G4int, const G4double*> Z54N70(70,pZ54N70);
  static const G4double pZ54N72[7]={1.5E-6, 14., 5.6E-4, 3.E-9, 1.2, 1.1E-5, .0036};
  static const std::pair<G4int, const G4double*> Z54N72(72,pZ54N72);
  static const G4double pZ54N74[7]={1.8E-6, 14., 8.8E-4, 2.E-9, 1.3, 1.2E-5, .0047};
  static const std::pair<G4int, const G4double*> Z54N74(74,pZ54N74);
  static const G4double pZ54N75[7]={1.E-6, 14., 2.6E-4, 3.7E-9, 1.5, 1.4E-5, .0024};
  static const std::pair<G4int, const G4double*> Z54N75(75,pZ54N75);
  static const G4double pZ54N76[7]={1.8E-6, 14., 8.E-4, 2.E-9, 1.2, 1.4E-5, .0042};
  static const std::pair<G4int, const G4double*> Z54N76(76,pZ54N76);
  static const G4double pZ54N77[7]={2.3E-7, 14., 1.9E-5, 9.E-9, 1.2, 1.4E-5, 7.7E-4};
  static const std::pair<G4int, const G4double*> Z54N77(77,pZ54N77);
  static const G4double pZ54N78[7]={6.E-7, 14., 1.6E-4, 3.E-9, 1.2, 1.4E-5, .0025};
  static const std::pair<G4int, const G4double*> Z54N78(78,pZ54N78);
  static const G4double pZ54N79[7]={6.E-7, 14., 1.6E-4, 3.3E-9, 1.6, 1.5E-5, .0024};
  static const std::pair<G4int, const G4double*> Z54N79(79,pZ54N79);
  static const G4double pZ54N80[7]={6.6E-7, 14., 2.1E-4, 2.5E-9, 1.2, 1.4E-5, .003};
  static const std::pair<G4int, const G4double*> Z54N80(80,pZ54N80);
  static const G4double pZ54N81[7]={.03, 40., 2.1, 2.5E-9, 1.E-16, 6.E-36, 140.};
  static const std::pair<G4int, const G4double*> Z54N81(81,pZ54N81);
  static const G4double pZ54N82[7]={3.1E-6, 14., .0019, 1.6E-9, 1., 1.3E-5, .0054};
  static const std::pair<G4int, const G4double*> Z54N82(82,pZ54N82);
  static const std::pair<G4int, const G4double*> Z54[N54]={Z54N69, Z54N70, Z54N72, Z54N74,
                                                           Z54N75, Z54N76, Z54N77, Z54N78,
                                                           Z54N79, Z54N80, Z54N81, Z54N82};
  //==> Cs(Z=55)
  static const G4int N55=5;
  static const G4double pZ55N78[7]={1.4E-6, 14., 4.E-4, 3.E-9, 1.2, 1.4E-5, .0026};
  static const std::pair<G4int, const G4double*> Z55N78(78,pZ55N78);
  static const G4double pZ55N79[7]={.028, 14., 44., .5E-9, 1.2, 1.3E-5, .015};
  static const std::pair<G4int, const G4double*> Z55N79(79,pZ55N79);
  static const G4double pZ55N80[7]={2.E-6, 14., 9.5E-4, 2.E-9, 1.2, 1.4E-5, .0042};
  static const std::pair<G4int, const G4double*> Z55N80(80,pZ55N80);
  static const G4double pZ55N81[7]={2.5E-7, 14., 6.5E-5, 3.8E-9, 1.2, 1.4E-5, .0023};
  static const std::pair<G4int, const G4double*> Z55N81(81,pZ55N81);
  static const G4double pZ55N82[7]={2.5E-7, 14., 6.5E-5, 3.8E-9, 1.2, 1.4E-5, .0023};
  static const std::pair<G4int, const G4double*> Z55N82(82,pZ55N82);
  static const std::pair<G4int, const G4double*> Z55[N55]={Z55N78, Z55N79, Z55N80, Z55N81,
                                                           Z55N82};
  //==> Ba(Z=56)
  static const G4int N56=9;
  static const G4double pZ56N74[7]={4.E-7, 14., 2.8E-5, 1.2E-8, 1.2, 1.5E-5, 6.6E-4};
  static const std::pair<G4int, const G4double*> Z56N74(74,pZ56N74);
  static const G4double pZ56N76[7]={4.E-6, 14., .0022, 1.4E-9, 1.3, 1.6E-5, .0053};
  static const std::pair<G4int, const G4double*> Z56N76(76,pZ56N76);
  static const G4double pZ56N77[7]={2.E-7, 14., 3.7E-5, 5.E-9, 1.1, 1.5E-5, .0016};
  static const std::pair<G4int, const G4double*> Z56N77(77,pZ56N77);
  static const G4double pZ56N78[7]={1.6E-6, 14., 6.E-4, 2.E-9, 1.3, 1.6E-5, .0033};
  static const std::pair<G4int, const G4double*> Z56N78(78,pZ56N78);
  static const G4double pZ56N79[7]={5.E-7, 17., 8.E-5, 4.5E-9, 1.3, 1.6E-5, .0018};
  static const std::pair<G4int, const G4double*> Z56N79(79,pZ56N79);
  static const G4double pZ56N80[7]={2.E-6, 20., 3.E-4, 6.E-9, 1.3, 1.8E-5, .0019};
  static const std::pair<G4int, const G4double*> Z56N80(80,pZ56N80);
  static const G4double pZ56N81[7]={5.8E-6, 20., .0018, 3.E-9, 1.3, 1.7E-5, .0041};
  static const std::pair<G4int, const G4double*> Z56N81(81,pZ56N81);
  static const G4double pZ56N82[7]={2.7E-6, 20., 5.5E-4, 4.E-9, 1.4, 2.E-5, .0027};
  static const std::pair<G4int, const G4double*> Z56N82(82,pZ56N82);
  static const G4double pZ56N84[7]={1.1E-6, 21., 1.E-4, 9.E-9, 2., 2.7E-5, .0012};
  static const std::pair<G4int, const G4double*> Z56N84(84,pZ56N84);
  static const std::pair<G4int, const G4double*> Z56[N56]={Z56N74, Z56N76, Z56N77, Z56N78,
                                                           Z56N79, Z56N80, Z56N81, Z56N82,
                                                           Z56N84};
  //==> La(Z=57)
  static const G4int N57=3;
  static const G4double pZ57N81[7]={2.7E-6, 20., .0017, 1.1E-9, 1.4, 2.E-5, .0083};
  static const std::pair<G4int, const G4double*> Z57N81(81,pZ57N81);
  static const G4double pZ57N82[7]={5.4E-6, 20., .0027, 1.5E-9, 1., 1.5E-5, .0065};
  static const std::pair<G4int, const G4double*> Z57N82(82,pZ57N82);
  static const G4double pZ57N83[7]={2.7E-6, 20., 2.6E-4, 6.E-9, 1.4, 2.E-5, .0012};
  static const std::pair<G4int, const G4double*> Z57N83(83,pZ57N83);
  static const std::pair<G4int, const G4double*> Z57[N57]={Z57N81, Z57N82, Z57N83};
  //==> Ce(Z=58)
  static const G4int N58=8;
  static const G4double pZ58N78[7]={1.8E-6, 20., 3.7E-4, 4.E-9, 1.4, 2.E-5, .0028};
  static const std::pair<G4int, const G4double*> Z58N78(78,pZ58N78);
  static const G4double pZ58N80[7]={1.8E-6, 18., 2.6E-4, 6.E-9, 1.3, 2.1E-5, .0017};
  static const std::pair<G4int, const G4double*> Z58N80(80,pZ58N80);
  static const G4double pZ58N81[7]={.0018, 18., 2.9, .6E-9, 1.3, 2.E-5, .02};
  static const std::pair<G4int, const G4double*> Z58N81(81,pZ58N81);
  static const G4double pZ58N82[7]={1.8E-6, 18., 3.7E-4, 5.E-9, 1.3, 2.E-5, .0024};
  static const std::pair<G4int, const G4double*> Z58N82(82,pZ58N82);
  static const G4double pZ58N83[7]={7.2E-6, 20., .0025, 2.3E-9, 1.1, 1.7E-5, .0045};
  static const std::pair<G4int, const G4double*> Z58N83(83,pZ58N83);
  static const G4double pZ58N84[7]={1.2E-6, 18., 2.E-4, 6.E-9, 1.5, 1.9E-5, .0018};
  static const std::pair<G4int, const G4double*> Z58N84(84,pZ58N84);
  static const G4double pZ58N85[7]={1.2E-6, 16., 6.E-4, 1.7E-9, 1.4, 2.E-5, .0053};
  static const std::pair<G4int, const G4double*> Z58N85(85,pZ58N85);
  static const G4double pZ58N86[7]={6.E-7, 18., 1.E-4, 6.E-9, 1.5, 1.9E-5, .0018};
  static const std::pair<G4int, const G4double*> Z58N86(86,pZ58N86);
  static const std::pair<G4int, const G4double*> Z58[N58]={Z58N78, Z58N80, Z58N81, Z58N82,
                                                           Z58N83, Z58N84, Z58N85, Z58N86};
  //==> Pr(Z=59)
  static const G4int N59=3;
  static const G4double pZ59N82[7]={9.5E-7, 16., 1.6E-4, 4.E-9, 1.4, 2.E-5, .0017};
  static const std::pair<G4int, const G4double*> Z59N82(82,pZ59N82);
  static const G4double pZ59N83[7]={9.5E-7, 16., 1.9E-4, 4.E-9, 1.4, 1.9E-5, .0021};
  static const std::pair<G4int, const G4double*> Z59N83(83,pZ59N83);
  static const G4double pZ59N84[7]={9.5E-6, 16., .019, .4E-9, 2., 2.4E-5, .021};
  static const std::pair<G4int, const G4double*> Z59N84(84,pZ59N84);
  static const std::pair<G4int, const G4double*> Z59[N59]={Z59N82, Z59N83, Z59N84};
  //==> Nd(Z=60)
  static const G4int N60=8;
  static const G4double pZ60N82[7]={9.6E-6, 21., .0036, 2.3E-9, 1.4, 2.E-5, .0052};
  static const std::pair<G4int, const G4double*> Z60N82(82,pZ60N82);
  static const G4double pZ60N83[7]={9.6E-4, 20., 4., .25E-9, 1.4, 2.E-5, .052};
  static const std::pair<G4int, const G4double*> Z60N83(83,pZ60N83);
  static const G4double pZ60N84[7]={4.8E-7, 21., 3.3E-5, 1.E-08, 1.3, 2.2E-5, 9.E-4};
  static const std::pair<G4int, const G4double*> Z60N84(84,pZ60N84);
  static const G4double pZ60N85[7]={.0048, 20., 4.5, .9E-9, 1.3, 2.E-5, .012};
  static const std::pair<G4int, const G4double*> Z60N85(85,pZ60N85);
  static const G4double pZ60N86[7]={1.2E-6, 16., 7.7E-4, 1.5E-9, 1.3, 1.8E-5, .0066};
  static const std::pair<G4int, const G4double*> Z60N86(86,pZ60N86);
  static const G4double pZ60N87[7]={.0012, 15., 8.4, .1E-9, 1.3, 1.5E-5, .071};
  static const std::pair<G4int, const G4double*> Z60N87(87,pZ60N87);
  static const G4double pZ60N88[7]={1.5E-7, 16., 4.1E-5, 2.5E-9, 1.3, 1.6E-5, .0027};
  static const std::pair<G4int, const G4double*> Z60N88(88,pZ60N88);
  static const G4double pZ60N90[7]={1.5E-7, 16., 4.3E-5, 2.5E-9, 1.3, 1.6E-5, .0029};
  static const std::pair<G4int, const G4double*> Z60N90(90,pZ60N90);
  static const std::pair<G4int, const G4double*> Z60[N60]={Z60N82, Z60N83, Z60N84, Z60N85,
                                                           Z60N86, Z60N87, Z60N88, Z60N90};
  //==> Pm(Z=61)
  static const G4int N61=3;
  static const G4double pZ61N86[7]={6.E-7, 16., 8.E-4, .6E-9, 3., 2.8E-5, .014};
  static const std::pair<G4int, const G4double*> Z61N86(86,pZ61N86);
  static const G4double pZ61N87[7]={6.2E-8, 16., 1.2E-5, 4.E-9, 2.2, 2.5E-5, .0019};
  static const std::pair<G4int, const G4double*> Z61N87(87,pZ61N87);
  static const G4double pZ61N88[7]={3.2E-8, 16., 6.4E-6, 4.E-9, 2.2, 2.5E-5, .002};
  static const std::pair<G4int, const G4double*> Z61N88(88,pZ61N88);
  static const std::pair<G4int, const G4double*> Z61[N61]={Z61N86, Z61N87, Z61N88};
  //==> Sm(Z=62)
  static const G4int N62=9;
  static const G4double pZ62N82[7]={1.2E-7, 16., 2.1E-5, 5.E-9, 1.4, 2.E-5, .0017};
  static const std::pair<G4int, const G4double*> Z62N82(82,pZ62N82);
  static const G4double pZ62N85[7]={1.2E-7, 16., 5.3E-5, 1.5E-9, 1.3, 1.7E-5, .0045};
  static const std::pair<G4int, const G4double*> Z62N85(85,pZ62N85);
  static const G4double pZ62N86[7]={6.E-8, 16., 1.7E-5, 3.E-9, 1.3, 1.7E-5, .0028};
  static const std::pair<G4int, const G4double*> Z62N86(86,pZ62N86);
  static const G4double pZ62N87[7]={6.E-8, 15., 5.2E-4, .11E-9, 1.3, 1.5E-5, .074};
  static const std::pair<G4int, const G4double*> Z62N87(87,pZ62N87);
  static const G4double pZ62N88[7]={6.E-7, 16., 8.6E-4, .7E-9, 1.3, 1.7E-5, .015};
  static const std::pair<G4int, const G4double*> Z62N88(88,pZ62N88);
  static const G4double pZ62N89[7]={6.E-7, 16., 5.E-4, 0., 1.3, 1.6E-5, .053};
  static const std::pair<G4int, const G4double*> Z62N89(89,pZ62N89);
  static const G4double pZ62N90[7]={6.E-8, 15., 1.3E-5, 4.5E-9, 1.3, 1.7E-5, .0019};
  static const std::pair<G4int, const G4double*> Z62N90(90,pZ62N90);
  static const G4double pZ62N91[7]={6.E-8, 15., 1.5E-5, 2.E-9, 1.3, 1.6E-5, .0024};
  static const std::pair<G4int, const G4double*> Z62N91(91,pZ62N91);
  static const G4double pZ62N92[7]={1.2E-7, 15., 8.6E-5, 1.2E-9, 1.3, 1.6E-5, .007};
  static const std::pair<G4int, const G4double*> Z62N92(92,pZ62N92);
  static const std::pair<G4int, const G4double*> Z62[N62]={Z62N82, Z62N85, Z62N86, Z62N87,
                                                           Z62N88, Z62N89, Z62N90, Z62N91,
                                                           Z62N92};
  //==> Eu(Z=63)
  static const G4int N63=7;
  static const G4double pZ63N88[7]={6.E-8, 15., 2.8E-5, 2.E-9, 1.3, 1.5E-5, .0046};
  static const std::pair<G4int, const G4double*> Z63N88(88,pZ63N88);
  static const G4double pZ63N89[7]={6.E-7, 15., .0011, .5E-9, 2.4, 2.4E-5, .017};
  static const std::pair<G4int, const G4double*> Z63N89(89,pZ63N89);
  static const G4double pZ63N90[7]={3.E-7, 15., 1.8E-4, 1.1E-9, 1., 1.2E-5, .0054};
  static const std::pair<G4int, const G4double*> Z63N90(90,pZ63N90);
  static const G4double pZ63N91[7]={4.1E-7, 15., 1.4E-4, 1.9E-9, 1., 1.4E-5, .0032};
  static const std::pair<G4int, const G4double*> Z63N91(91,pZ63N91);
  static const G4double pZ63N92[7]={5.E-8, 15., 2.4E-5, 2.8E-9, 1., 1.3E-5, .0037};
  static const std::pair<G4int, const G4double*> Z63N92(92,pZ63N92);
  static const G4double pZ63N93[7]={4.1E-8, 17., 1.6E-5, 2.E-9, 3.3, 3.4E-5, .004};
  static const std::pair<G4int, const G4double*> Z63N93(93,pZ63N93);
  static const G4double pZ63N94[7]={4.2E-8, 17., 1.6E-5, 2.E-9, 1.2, 1.6E-5, .004};
  static const std::pair<G4int, const G4double*> Z63N94(94,pZ63N94);
  static const std::pair<G4int, const G4double*> Z63[N63]={Z63N88, Z63N89, Z63N90, Z63N91,
                                                           Z63N92, Z63N93, Z63N94};
  //==> Gd(Z=64)
  static const G4int N64=8;
  static const G4double pZ64N88[7]={4.2E-8, 14., 2.E-4, 0., 1.2, 1.3E-5, .19};
  static const std::pair<G4int, const G4double*> Z64N88(88,pZ64N88);
  static const G4double pZ64N89[7]={2.E-6, 14., .0016, 1.4E-9, 1.6, 1.6E-5, .0057};
  static const std::pair<G4int, const G4double*> Z64N89(89,pZ64N89);
  static const G4double pZ64N90[7]={1.7E-7, 12., 8.4E-5, 2.E-9, 1.8, 2.2E-5, .0035};
  static const std::pair<G4int, const G4double*> Z64N90(90,pZ64N90);
  static const G4double pZ64N91[7]={1.7E-7, 13., 5.E-4, .3E-9, 1.7, 1.9E-5, .026};
  static const std::pair<G4int, const G4double*> Z64N91(91,pZ64N91);
  static const G4double pZ64N92[7]={1.7E-7, 13., 7.E-5, 2.5E-9, 1.8, 2.2E-5, .003};
  static const std::pair<G4int, const G4double*> Z64N92(92,pZ64N92);
  static const G4double pZ64N93[7]={1.7E-6, 12., .002, 0., 1.7, 1.8E-5, .47};
  static const std::pair<G4int, const G4double*> Z64N93(93,pZ64N93);
  static const G4double pZ64N94[7]={3.4E-7, 13., 1.5E-4, 2.E-9, 1.8, 2.3E-5, .0034};
  static const std::pair<G4int, const G4double*> Z64N94(94,pZ64N94);
  static const G4double pZ64N96[7]={2.6E-6, 13., .0019, 1.2E-9, 1., 1.2E-5, .0056};
  static const std::pair<G4int, const G4double*> Z64N96(96,pZ64N96);
  static const std::pair<G4int, const G4double*> Z64[N64]={Z64N88, Z64N89, Z64N90, Z64N91,
                                                           Z64N92, Z64N93, Z64N94, Z64N96};
  //==> Tb(Z=65)
  static const G4int N65=2;
  static const G4double pZ65N94[7]={9.E-7, 16., 3.9E-4, 1.7E-9, 2., 2.2E-5, .0042};
  static const std::pair<G4int, const G4double*> Z65N94(94,pZ65N94);
  //static const G4double pZ65N95[7]={4.5E-7, 16., 1.1E-4, 3.E-9, 1.7, 2.2E-5, .0024};
  //static const std::pair<G4int, const G4double*> Z65N95(95,pZ65N95);
  static const std::pair<G4int, const G4double*> Z65[N65]={Z65N94, Z65N94};
  //==> Dy(Z=66)
  static const G4int N66=7;
  static const G4double pZ66N90[7]={1.2E-7, 13., 4.E-5, 3.E-9, 1.2, 1.4E-5, .0025};
  static const std::pair<G4int, const G4double*> Z66N90(90,pZ66N90);
  static const G4double pZ66N92[7]={1.2E-7, 13., 6.7E-5, 2.E-9, 1., 1.1E-5, .004};
  static const std::pair<G4int, const G4double*> Z66N92(92,pZ66N92);
  static const G4double pZ66N94[7]={1.2E-7, 13., 5.3E-5, 1.6E-9, 1., 1.1E-5, .0034};
  static const std::pair<G4int, const G4double*> Z66N94(94,pZ66N94);
  static const G4double pZ66N95[7]={1.2E-6, 13., .0017, .7E-9, 1.3, 1.3E-5, .011};
  static const std::pair<G4int, const G4double*> Z66N95(95,pZ66N95);
  static const G4double pZ66N96[7]={1.2E-7, 13., 8.E-6, 1.5E-7, 1., 1.1E-5, 1.E-4};
  static const std::pair<G4int, const G4double*> Z66N96(96,pZ66N96);
  static const G4double pZ66N97[7]={1.5E-7, 13., 4.E-5, 4.E-9, 1.3, 1.3E-5, .002};
  static const std::pair<G4int, const G4double*> Z66N97(97,pZ66N97);
  static const G4double pZ66N98[7]={3.E-7, 13., .001, 4.E-9, 1.3, 1.3E-5, .23};
  static const std::pair<G4int, const G4double*> Z66N98(98,pZ66N98);
  static const std::pair<G4int, const G4double*> Z66[N66]={Z66N90, Z66N92, Z66N94, Z66N95,
                                                           Z66N96, Z66N97, Z66N98};
  //==> Ho(Z=67)
  static const G4int N67=2;
  static const G4double pZ67N98[7]={3.E-7, 13., 2.2E-4, 1.5E-9, 1., 1.E-5, .0054};
  static const std::pair<G4int, const G4double*> Z67N98(98,pZ67N98);
  static const G4double pZ67N99[7]={7.5E-8, 13., 2.6E-5, 4.5E-9, 1.5, 1.5E-5, .0021};
  static const std::pair<G4int, const G4double*> Z67N99(99,pZ67N99);
  static const std::pair<G4int, const G4double*> Z67[N67]={Z67N98, Z67N99};
  //==> Er(Z=68)
  static const G4int N68=6;
  static const G4double pZ68N94[7]={1.2E-7, 13., 7.8E-5, 1.6E-9, .9, 9.E-6, .005};
  static const std::pair<G4int, const G4double*> Z68N94(94,pZ68N94);
  static const G4double pZ68N96[7]={1.2E-7, 13., 8.5E-5, 1.2E-9, .9, 8.E-6, .0055};
  static const std::pair<G4int, const G4double*> Z68N96(96,pZ68N96);
  static const G4double pZ68N98[7]={1.E-6, 13., .0011, .8E-9, .9, 8.E-6, .0087};
  static const std::pair<G4int, const G4double*> Z68N98(98,pZ68N98);
  static const G4double pZ68N99[7]={1.2E-7, 13., 2.5E-5, 4.5E-9, .9, 9.E-6, .0015};
  static const std::pair<G4int, const G4double*> Z68N99(99,pZ68N99);
  static const G4double pZ68N100[7]={2.E-6, 13., .0015, 1.1E-9, .9, 8.E-6, .0058};
  static const std::pair<G4int, const G4double*> Z68N100(100,pZ68N100);
  static const G4double pZ68N102[7]={2.E-6, 13., .0018, 1.E-9, .9, 8.E-6, .007};
  static const std::pair<G4int, const G4double*> Z68N102(102,pZ68N102);
  static const std::pair<G4int, const G4double*> Z68[N68]={Z68N94, Z68N96, Z68N98,
                                                           Z68N99, Z68N100, Z68N102};
  //==> Tm(Z=69) *** No data *** (Tm169=Er167)
  static const G4int N69=1;
  static const G4double pZ69N100[7]={1.2E-7, 13., 2.5E-5, 4.5E-9, .9, 9.E-6, .0015};
  static const std::pair<G4int, const G4double*> Z69N100(100,pZ69N100);
  static const std::pair<G4int, const G4double*> Z69[N69]={Z69N100};
  //==> Yb(Z=70) *** No data *** (Yb168=Er166, Yb170=Er168, Yb171=Er167, Yb172=Er170,
  //                              Yb173=Hf177, Yb174=Hf176, Yb176=Hf178)
  static const G4int N70=7;
  static const G4double pZ70N98[7]={1.E-6, 13., .0011, .8E-9, .9, 8.E-6, .0087};
  static const std::pair<G4int, const G4double*> Z70N98(98,pZ70N98);
  static const G4double pZ70N100[7]={2.E-6, 13., .0015, 1.1E-9, .9, 8.E-6, .0058};
  static const std::pair<G4int, const G4double*> Z70N100(100,pZ70N100);
  static const G4double pZ70N101[7]={1.2E-7, 13., 2.5E-5, 4.5E-9, .9, 9.E-6, .0015};
  static const std::pair<G4int, const G4double*> Z70N101(101,pZ70N101);
  static const G4double pZ70N102[7]={2.E-6, 13., .0018, 1.E-9, .9, 8.E-6, .007};
  static const std::pair<G4int, const G4double*> Z70N102(102,pZ70N102);
  static const G4double pZ70N103[7]={5.E-7, 18., .01, 1.7E-6, 1.2, 1.4E-5, 1.4E-5};
  static const std::pair<G4int, const G4double*> Z70N103(103,pZ70N103);
  static const G4double pZ70N104[7]={5.E-7, 18., 1.9E-4, 2.5E-9, 1.2, 1.4E-5, .004};
  static const std::pair<G4int, const G4double*> Z70N104(104,pZ70N104);
  static const G4double pZ70N106[7]={5.E-7, 18., 1.3E-4, 2.E-9, 1.2, 1.4E-5, .0027};
  static const std::pair<G4int, const G4double*> Z70N106(106,pZ70N106);
  static const std::pair<G4int, const G4double*> Z70[N70]={Z70N98, Z70N100, Z70N101,
                                                           Z70N102, Z70N103, Z70N104,
                                                           Z70N106};
  //==> Lu(Z=71)
  static const G4int N71=2;
  static const G4double pZ71N104[7]={5.E-7, 18., 1.8E-4, 2.E-9, .9, 9.E-6, .0036};
  static const std::pair<G4int, const G4double*> Z71N104(104,pZ71N104);
  static const G4double pZ71N105[7]={2.5E-7, 18., 9.E-5, 1.E-8, .9, 9.E-6, .0016};
  static const std::pair<G4int, const G4double*> Z71N105(105,pZ71N105);
  static const std::pair<G4int, const G4double*> Z71[N71]={Z71N104, Z71N105};
  //==> Hf(Z=72)
  static const G4int N72=6;
  static const G4double pZ72N102[7]={1.E-6, 18., 8.8E-4, .8E-9, 1., 1.1E-5, .0092};
  static const std::pair<G4int, const G4double*> Z72N102(102,pZ72N102);
  static const G4double pZ72N104[7]={5.E-7, 18., 1.9E-4, 2.5E-9, 1.2, 1.4E-5, .004};
  static const std::pair<G4int, const G4double*> Z72N104(104,pZ72N104);
  static const G4double pZ72N105[7]={5.E-7, 18., .01, 1.7E-6, 1.2, 1.4E-5, 1.4E-5};
  static const std::pair<G4int, const G4double*> Z72N105(105,pZ72N105);
  static const G4double pZ72N106[7]={5.E-7, 18., 1.3E-4, 2.E-9, 1.2, 1.4E-5, .0027};
  static const std::pair<G4int, const G4double*> Z72N106(106,pZ72N106);
  static const G4double pZ72N107[7]={2.5E-7, 18., 1.E-4, 2.E-9, 1.2, 1.5E-5, .0041};
  static const std::pair<G4int, const G4double*> Z72N107(107,pZ72N107);
  static const G4double pZ72N108[7]={1.E-6, 18., .0012, .6E-9, 1.2, 1.5E-5, .012};
  static const std::pair<G4int, const G4double*> Z72N108(108,pZ72N108);
  static const std::pair<G4int, const G4double*> Z72[N72]={Z72N102, Z72N104, Z72N105,
                                                           Z72N106, Z72N107, Z72N108};
  //==> Ta(Z=73)
  static const G4int N73=2;
  static const G4double pZ73N108[7]={5.E-7, 18., 1.7E-4, 2.E-9, 1.2, 1.4E-5, .0035};
  static const std::pair<G4int, const G4double*> Z73N108(108,pZ73N108);
  //static const G4double pZ73N109[7]={1.E-6, 14., .002, .3E-9, 1.3, 1.5E-5, .016};
  //static const std::pair<G4int, const G4double*> Z73N109(109,pZ73N109);
  static const std::pair<G4int, const G4double*> Z73[N73]={Z73N108, Z73N108};
  //==> W (Z=74) *** W180 only bad TENDL-2008 *** (W180=Hf178)
  static const G4int N74=5;
  static const G4double pZ74N106[7]={5.E-7, 18., 1.3E-4, 2.E-9, 1.2, 1.4E-5, .0027};
  static const std::pair<G4int, const G4double*> Z74N106(106,pZ74N106);
  static const G4double pZ74N108[7]={4.E-6, 14., .0034, .9E-9, 1.3, 1.4E-5, .0067};
  static const std::pair<G4int, const G4double*> Z74N108(108,pZ74N108);
  static const G4double pZ74N109[7]={1.2E-7, 14., 3.E-5, 3.E-9, 1.3, 1.4E-5, .0019};
  static const std::pair<G4int, const G4double*> Z74N109(109,pZ74N109);
  static const G4double pZ74N110[7]={1.2E-7, 14., 3.6E-5, 2.E-9, 1.3, 1.3E-5, .0024};
  static const std::pair<G4int, const G4double*> Z74N110(110,pZ74N110);
  static const G4double pZ74N112[7]={1.2E-7, 14., 3.E-5, 1.3E-7, 1.3, 1.3E-5, 1.4E-4};
  static const std::pair<G4int, const G4double*> Z74N112(112,pZ74N112);
  static const std::pair<G4int, const G4double*> Z74[N74]={Z74N106, Z74N108, Z74N109,
                                                           Z74N110, Z74N112};
  //==> Re(Z=75)
  static const G4int N75=2;
  static const G4double pZ75N110[7]={1.2E-7, 14., 8.E-5, 1.2E-9, 1.3, 1.5E-5, .005};
  static const std::pair<G4int, const G4double*> Z75N110(110,pZ75N110);
  static const G4double pZ75N112[7]={1.2E-7, 14., 8.8E-5, 1.1E-9, 1.3, 1.5E-5, .0055};
  static const std::pair<G4int, const G4double*> Z75N112(112,pZ75N112);
  static const std::pair<G4int, const G4double*> Z75[N75]={Z75N110, Z75N112};
  //==> Os(Z=76) *** No data *** (Os184=W182, Os186=W184, Os187=Re187, Os188=W186,
  //                              Os189=Re187, Os190=Os192=W186)
  static const G4int N76=7;
  static const G4double pZ76N108[7]={4.E-6, 14., .0034, .9E-9, 1.3, 1.4E-5, .0067};
  static const std::pair<G4int, const G4double*> Z76N108(108,pZ76N108);
  static const G4double pZ76N110[7]={1.2E-7, 14., 3.6E-5, 2.E-9, 1.3, 1.3E-5, .0024};
  static const std::pair<G4int, const G4double*> Z76N110(110,pZ76N110);
  static const G4double pZ76N111[7]={1.2E-7, 14., 8.8E-5, 1.1E-9, 1.3, 1.5E-5, .0055};
  static const std::pair<G4int, const G4double*> Z76N111(111,pZ76N111);
  static const G4double pZ76N112[7]={1.2E-7, 14., 3.E-5, 1.3E-7, 1.3, 1.3E-5, 1.4E-4};
  static const std::pair<G4int, const G4double*> Z76N112(112,pZ76N112);
  static const G4double pZ76N113[7]={1.2E-7, 14., 8.8E-5, 1.1E-9, 1.3, 1.5E-5, .0055};
  static const std::pair<G4int, const G4double*> Z76N113(113,pZ76N113);
  static const G4double pZ76N114[7]={1.2E-7, 14., 3.E-5, 1.3E-7, 1.3, 1.3E-5, 1.4E-4};
  static const std::pair<G4int, const G4double*> Z76N114(114,pZ76N114);
  static const G4double pZ76N116[7]={1.2E-7, 14., 3.E-5, 1.3E-7, 1.3, 1.3E-5, 1.4E-4};
  static const std::pair<G4int, const G4double*> Z76N116(116,pZ76N116);
  static const std::pair<G4int, const G4double*> Z76[N76]={Z76N108, Z76N110, Z76N111,
                                                           Z76N112, Z76N113, Z76N114,
                                                           Z76N116};
  //==> Ir(Z=77)
  static const G4int N77=2;
  static const G4double pZ77N114[7]={4.8E-7, 14., 5.2E-4, .7E-9, 1.5, 1.7E-5, .0082};
  static const std::pair<G4int, const G4double*> Z77N114(114,pZ77N114);
  static const G4double pZ77N116[7]={4.8E-7, 14., 4.5E-4, .8E-9, 1.8, 2.3E-5, .0073};
  static const std::pair<G4int, const G4double*> Z77N116(116,pZ77N116);
  static const std::pair<G4int, const G4double*> Z77[N77]={Z77N114, Z77N116};
  //==> Pt(Z=78) *** No data *** (Pt190=Pt192=Pt194=Hg196, Pt195=Hg199, Pt196=Hg198,
  //                              Pt198=Hg200)
  static const G4int N78=6;
  static const G4double pZ78N112[7]={6.E-8, 19., 3.3E-4, .1E-9, 1.6, 1.8E-5, .06};
  static const std::pair<G4int, const G4double*> Z78N112(112,pZ78N112);
  static const G4double pZ78N114[7]={6.E-8, 19., 3.3E-4, .1E-9, 1.6, 1.8E-5, .06};
  static const std::pair<G4int, const G4double*> Z78N114(114,pZ78N114);
  static const G4double pZ78N116[7]={6.E-8, 19., 3.3E-4, .1E-9, 1.6, 1.8E-5, .06};
  static const std::pair<G4int, const G4double*> Z78N116(116,pZ78N116);
  static const G4double pZ78N117[7]={9.6E-7, 20., .001, .2E-9, 1.6, 2.E-5, .037};
  static const std::pair<G4int, const G4double*> Z78N117(117,pZ78N117);
  static const G4double pZ78N118[7]={2.4E-7, 20., 1.6E-4, 1.3E-9, 1.6, 1.8E-5, .007};
  static const std::pair<G4int, const G4double*> Z78N118(118,pZ78N118);
  static const G4double pZ78N120[7]={2.E-6, 19., .0015, .9E-9, 1.6, 1.8E-5, .0078};
  static const std::pair<G4int, const G4double*> Z78N120(120,pZ78N120);
  static const std::pair<G4int, const G4double*> Z78[N78]={Z78N112, Z78N114, Z78N116,
                                                           Z78N117, Z78N118, Z78N120};
  //==> Au(Z=79)
  static const G4int N79=1;
  static const G4double pZ79N118[7]={2.4E-7, 19., 1.E-4, 1.4E-9, 1.3, 1.7E-5, .0042};
  static const std::pair<G4int, const G4double*> Z79N118(118,pZ79N118);
  static const std::pair<G4int, const G4double*> Z79[N79]={Z79N118};
  //==> Hg(Z=80)
  static const G4int N80=7;
  static const G4double pZ80N116[7]={6.E-8, 19., 3.3E-4, .1E-9, 1.6, 1.8E-5, .06};
  static const std::pair<G4int, const G4double*> Z80N116(116,pZ80N116);
  static const G4double pZ80N118[7]={2.4E-7, 20., 1.6E-4, 1.3E-9, 1.6, 1.8E-5, .007};
  static const std::pair<G4int, const G4double*> Z80N118(118,pZ80N118);
  static const G4double pZ80N119[7]={9.6E-7, 20., .001, .2E-9, 1.6, 2.E-5, .037};
  static const std::pair<G4int, const G4double*> Z80N119(119,pZ80N119);
  static const G4double pZ80N120[7]={2.E-6, 19., .0015, .9E-9, 1.6, 1.8E-5, .0078};
  static const std::pair<G4int, const G4double*> Z80N120(120,pZ80N120);
  static const G4double pZ80N121[7]={1.E-6, 20., 7.E-4, 1.E-9, 1.6, 1.8E-5, .0076};
  static const std::pair<G4int, const G4double*> Z80N121(121,pZ80N121);
  static const G4double pZ80N122[7]={2.E-6, 18., .0016, .8E-9, 1.6, 1.8E-5, .0078};
  static const std::pair<G4int, const G4double*> Z80N122(122,pZ80N122);
  static const G4double pZ80N124[7]={2.0E-6, 18., .0032, .4E-9, 1.6, 1.8E-5, .016};
  static const std::pair<G4int, const G4double*> Z80N124(124,pZ80N124);
  static const std::pair<G4int, const G4double*> Z80[N80]={Z80N116, Z80N118, Z80N119,
                                                           Z80N120, Z80N121, Z80N122,
                                                           Z80N124};
  //==> Tl(Z=81) *** No data *** (Tl203=Au196, Tl198=Bi209)
  static const G4int N81=2;
  static const G4double pZ81N122[7]={2.4E-7, 19., 1.E-4, 1.4E-9, 1.3, 1.7E-5, .0042};
  static const std::pair<G4int, const G4double*> Z81N122(122,pZ81N122);
  static const G4double pZ81N124[7]={};
  static const std::pair<G4int, const G4double*> Z81N124(124,pZ81N124);
  static const std::pair<G4int, const G4double*> Z81[N81]={Z81N122, Z81N124};
  //==> Pb(Z=82)
  static const G4int N82=4;
  static const G4double pZ82N122[7]={4.E-6, 20., .0022, 1.E-9, 1.6, 1.8E-5, .0058};
  static const std::pair<G4int, const G4double*> Z82N122(122,pZ82N122);
  static const G4double pZ82N124[7]={4.E-6, 20., .0022, 1.E-9, 1.6, 1.8E-5, .0058};
  static const std::pair<G4int, const G4double*> Z82N124(124,pZ82N124);
  static const G4double pZ82N125[7]={2.E-6, 20., .0011, 1.2E-9, 1.6, 1.8E-5, .0056};
  static const std::pair<G4int, const G4double*> Z82N125(125,pZ82N125);
  static const G4double pZ82N126[7]={4.E-6, 20., .0023, 1.2E-9, 1.6, 1.8E-5, .0058};
  static const std::pair<G4int, const G4double*> Z82N126(126,pZ82N126);
  static const std::pair<G4int, const G4double*> Z82[N82]={Z82N122, Z82N124, Z82N125,
                                                           Z82N126};
  //==> Bi(Z=83)
  static const G4int N83=1;
  static const G4double pZ83N126[7]={8.E-7, 23., 3.3E-4, 1.8E-9, 1.6, 1.8E-5, .005};
  static const std::pair<G4int, const G4double*> Z83N126(126,pZ83N126);
  static const std::pair<G4int, const G4double*> Z83[N83]={Z83N126};
  //==> Po(Z=84) *** No data *** (Po209=Pb207)
  static const G4int N84=1;
  static const G4double pZ84N125[7]={2.E-6, 20., .0011, 1.2E-9, 1.6, 1.8E-5, .0056};
  static const std::pair<G4int, const G4double*> Z84N125(125,pZ84N125);
  static const std::pair<G4int, const G4double*> Z84[N84]={Z84N125};
  //==> At(Z=85) *** No data *** (At210=Pb207)
  static const G4int N85=1;
  static const G4double pZ85N125[7]={2.E-6, 20., .0011, 1.2E-9, 1.6, 1.8E-5, .0056};
  static const std::pair<G4int, const G4double*> Z85N125(125,pZ85N125);
  static const std::pair<G4int, const G4double*> Z85[N85]={Z85N125};
  //==> Rn(Z=86) *** No data *** (Rn222=Ra224)
  static const G4int N86=1;
  static const G4double pZ86N136[7]={1.E-7, 23., 5.5E-5, 1.2E-9, 1.6, 1.8E-5, .0062};
  static const std::pair<G4int, const G4double*> Z86N136(136,pZ86N136);
  static const std::pair<G4int, const G4double*> Z86[N86]={Z86N136};
  //==> Fr(Z=87) *** No data *** (Fr223=Ac225)
  static const G4int N87=1;
  static const G4double pZ87N136[7]={2.E-7, 23., 1.1E-4, 1.2E-9, 1.6, 1.8E-5, .0062};
  static const std::pair<G4int, const G4double*> Z87N136(136,pZ87N136);
  static const std::pair<G4int, const G4double*> Z87[N87]={Z87N136};
  //==> Ra(Z=88)
  static const G4int N88=4;
  static const G4double pZ88N135[7]={1.E-7, 23., 5.5E-5, 1.2E-9, 1.6, 1.8E-5, .0062};
  static const std::pair<G4int, const G4double*> Z88N135(135,pZ88N135);
  static const G4double pZ88N136[7]={1.E-7, 23., 5.5E-5, 1.2E-9, 1.6, 1.8E-5, .0062};
  static const std::pair<G4int, const G4double*> Z88N136(136,pZ88N136);
  static const G4double pZ88N137[7]={1.E-7, 23., 5.5E-5, 1.2E-9, 1.6, 1.8E-5, .0062};
  static const std::pair<G4int, const G4double*> Z88N137(137,pZ88N137);
  static const G4double pZ88N138[7]={4.E-7, 23., 1.7E-4, 1.5E-9, 1.6, 1.8E-5, .005};
  static const std::pair<G4int, const G4double*> Z88N138(138,pZ88N138);
  static const std::pair<G4int, const G4double*> Z88[N88]={Z88N135, Z88N136, Z88N137,
                                                           Z88N138};
  //==> Ac(Z=89)
  static const G4int N89=3;
  static const G4double pZ89N136[7]={2.E-7, 23., 1.1E-4, 1.2E-9, 1.6, 1.8E-5, .0062};
  static const std::pair<G4int, const G4double*> Z89N136(136,pZ89N136);
  static const G4double pZ89N137[7]={4.E-7, 23., 2.2E-4, 1.2E-9, 1.6, 1.8E-5, .0062};
  static const std::pair<G4int, const G4double*> Z89N137(137,pZ89N137);
  static const G4double pZ89N138[7]={1.E-7, 23., 5.5E-5, 1.2E-9, 1.6, 1.8E-5, .0062};
  static const std::pair<G4int, const G4double*> Z89N138(138,pZ89N138);
  static const std::pair<G4int, const G4double*> Z89[N89]={Z89N136, Z89N137, Z89N138};
  //==> Th(Z=90)
  static const G4int N90=7;
  static const G4double pZ90N137[7]={4.E-7, 23., 2.2E-4, 1.2E-9, 1.6, 1.8E-5, .0062};
  static const std::pair<G4int, const G4double*> Z90N137(137,pZ90N137);
  static const G4double pZ90N138[7]={1.E-6, 23., .0016, .4E-9, 3., 3.E-5, .019};
  static const std::pair<G4int, const G4double*> Z90N138(138,pZ90N138);
  static const G4double pZ90N139[7]={2.5E-7, 23., 1.1E-4, 1.4E-9, 2.4, 2.7E-5, .0049};
  static const std::pair<G4int, const G4double*> Z90N139(139,pZ90N139);
  static const G4double pZ90N140[7]={1.2E-7, 23., 3.E-5, 2.E-9, 3., 3.E-5, .003};
  static const std::pair<G4int, const G4double*> Z90N140(140,pZ90N140);
  static const G4double pZ90N142[7]={4.E-6, 23., .0023, 1.1E-9, 1.8, 2.3E-5, .0064};
  static const std::pair<G4int, const G4double*> Z90N142(142,pZ90N142);
  static const G4double pZ90N143[7]={9.4E-7, 23., 5.4E-4, 1.1E-9, 3., 3.E-5, .0066};
  static const std::pair<G4int, const G4double*> Z90N143(143,pZ90N143);
  static const G4double pZ90N144[7]={2.5E-7, 23., 1.4E-4, 1.1E-9, 3., 3.E-5, .0066};
  static const std::pair<G4int, const G4double*> Z90N144(144,pZ90N144);
  static const std::pair<G4int, const G4double*> Z90[N90]={Z90N137, Z90N138, Z90N139,
                                                           Z90N140, Z90N142, Z90N143,
                                                           Z90N144};
  //==> Pa(Z=91)
  static const G4int N91=3;
  static const G4double pZ91N140[7]={1.E-5, 23., .0052, 1.6E-9, 1.8, 2.3E-5, .0057};
  static const std::pair<G4int, const G4double*> Z91N140(140,pZ91N140);
  static const G4double pZ91N141[7]={8.E-6, 23., .006, 0., 3.5, 3.5E-5, .021};
  static const std::pair<G4int, const G4double*> Z91N141(141,pZ91N141);
  static const G4double pZ91N142[7]={8.E-6, 23., .0042, 1.E-9, 2., 2.5E-5, .006};
  static const std::pair<G4int, const G4double*> Z91N142(142,pZ91N142);
  static const std::pair<G4int, const G4double*> Z91[N91]={Z91N140, Z91N141, Z91N142};
  //==> U (Z=92)
  static const G4int N92=10;
  static const G4double pZ92N140[7]={1.4E-6, 20., 8.E-4, 1.5E-9, 2.5, 2.8E-5, .0055};
  static const std::pair<G4int, const G4double*> Z92N140(140,pZ92N140);
  static const G4double pZ92N141[7]={5.6E-6, 20., .0033, 1.E-9, 2.5, 2.8E-5, .006};
  static const std::pair<G4int, const G4double*> Z92N141(141,pZ92N141);
  static const G4double pZ92N142[7]={5.6E-6, 20., .0034, 0., 2.5, 2.8E-5, .0072};
  static const std::pair<G4int, const G4double*> Z92N142(142,pZ92N142);
  static const G4double pZ92N143[7]={5.6E-6, 20., .0032, 0., 2., 2.3E-5, .006};
  static const std::pair<G4int, const G4double*> Z92N143(143,pZ92N143);
  static const G4double pZ92N144[7]={3.6E-7, 20., 1.6E-4, 1.3E-9, 2.2, 2.7E-5, .0043};
  static const std::pair<G4int, const G4double*> Z92N144(144,pZ92N144);
  static const G4double pZ92N145[7]={3.6E-6, 20., .003, 0., 2.2, 2.7E-5, .045};
  static const std::pair<G4int, const G4double*> Z92N145(145,pZ92N145);
  static const G4double pZ92N146[7]={3.6E-7, 20., 1.6E-4, 1.3E-9, 2.2, 2.7E-5, .0043};
  static const std::pair<G4int, const G4double*> Z92N146(146,pZ92N146);
  static const G4double pZ92N147[7]={3.6E-6, 20., .0014, 1.3E-9, 2.2, 2.7E-5, 12.};
  static const std::pair<G4int, const G4double*> Z92N147(147,pZ92N147);
  static const G4double pZ92N148[7]={3.4E-7, 20., 1.3E-4, 1.3E-9, 2.2, 2.8E-5, .0036};
  static const std::pair<G4int, const G4double*> Z92N148(148,pZ92N148);
  //static const G4double pZ92N149[7]={3.3E-7, 20., 1.5E-4, 1.2E-9, 3., 3.4E-5, .0044};
  //static const std::pair<G4int, const G4double*> Z92N149(149,pZ92N149);
  static const std::pair<G4int, const G4double*> Z92[N92]={Z92N140, Z92N141, Z92N142,
                                                           Z92N143, Z92N144, Z92N145,
                                                           Z92N146, Z92N147, Z92N148,
                                                           Z92N146};
  //==> Np(Z=93)
  static const G4int N93=5;
  static const G4double pZ93N142[7]={3.4E-6, 20., .002, 1.3E-9, 3., 3.3E-5, .0056};
  static const std::pair<G4int, const G4double*> Z93N142(142,pZ93N142);
  static const G4double pZ93N143[7]={3.4E-6, 20., .002, 1.6E-9, 3.5, 3.6E-5, .005};
  static const std::pair<G4int, const G4double*> Z93N143(143,pZ93N143);
  static const G4double pZ93N144[7]={6.8E-6, 18., .0052, .8E-9, 2.4, 3.E-5, .0072};
  static const std::pair<G4int, const G4double*> Z93N144(144,pZ93N144);
  static const G4double pZ93N145[7]={3.4E-6, 20., .002, 1.E-9, 3.5, 3.6E-5, .006};
  static const std::pair<G4int, const G4double*> Z93N145(145,pZ93N145);
  static const G4double pZ93N146[7]={3.4E-6, 20., .002, 1.5E-9, 3.5, 3.6E-5, .0053};
  static const std::pair<G4int, const G4double*> Z93N146(146,pZ93N146);
  static const std::pair<G4int, const G4double*> Z93[N93]={Z93N142, Z93N143, Z93N144,
                                                           Z93N145, Z93N146};
  //==> Pu(Z=94)
  static const G4int N94=10;
  static const G4double pZ94N142[7]={6.8E-7, 16., 4.5E-4, 1.7E-9, 2.6, 3.E-5, .0047};
  static const std::pair<G4int, const G4double*> Z94N142(142,pZ94N142);
  static const G4double pZ94N143[7]={6.8E-6, 18., .0044, .9E-9, 3.3, 3.5E-5, .0058};
  static const std::pair<G4int, const G4double*> Z94N143(143,pZ94N143);
  static const G4double pZ94N144[7]={6.8E-7, 16., 6.E-4, 0., 2.7, 2.6E-5, .0082};
  static const std::pair<G4int, const G4double*> Z94N144(144,pZ94N144);
  static const G4double pZ94N145[7]={2.6E-6, 16., .0017, 1.8E-9, 1.8, 2.E-5, .004};
  static const std::pair<G4int, const G4double*> Z94N145(145,pZ94N145);
  static const G4double pZ94N146[7]={2.5E-7, 20., 9.E-5, 3.6E-8, 3.4, 3.8E-5, 5.4E-4};
  static const std::pair<G4int, const G4double*> Z94N146(146,pZ94N146);
  static const G4double pZ94N147[7]={1.4E-5, 16., .01, .8E-9, 2.7, 2.6E-5, .0055};
  static const std::pair<G4int, const G4double*> Z94N147(147,pZ94N147);
  static const G4double pZ94N148[7]={3.4E-7, 20., 1.3E-4, 1.2E-9, 3.2, 3.E-5, .0036};
  static const std::pair<G4int, const G4double*> Z94N148(148,pZ94N148);
  static const G4double pZ94N149[7]={5.2E-6, 20., .0035, .4E-9, 2.3, 3.E-5, .0095};
  static const std::pair<G4int, const G4double*> Z94N149(149,pZ94N149);
  static const G4double pZ94N150[7]={3.3E-7, 20., 1.6E-4, 1.2E-9, 3., 3.E-5, .0046};
  static const std::pair<G4int, const G4double*> Z94N150(150,pZ94N150);
  static const G4double pZ94N152[7]={2.5E-6, 16., .0018, 1.2E-9, 3., 3.1E-5, .0052};
  static const std::pair<G4int, const G4double*> Z94N152(152,pZ94N152);
  static const std::pair<G4int, const G4double*> Z94[N94]={Z94N142, Z94N143, Z94N144,
                                                           Z94N145, Z94N146, Z94N147,
                                                           Z94N148, Z94N149, Z94N150,
                                                           Z94N152};
  //==> Am(Z=95)
  static const G4int N95=4;
  static const G4double pZ95N156[7]={2.5E-6, 18., .0016, .9E-9, 2., 2.3E-5, .0058};
  static const std::pair<G4int, const G4double*> Z95N156(156,pZ95N156);
  static const G4double pZ95N157[7]={5.E-6, 18., .003, 2.7E-9, 2., 2.3E-5, .0039};
  static const std::pair<G4int, const G4double*> Z95N157(157,pZ95N157);
  static const G4double pZ95N158[7]={5.E-6, 19., .0033, 2.6E-9, 2., 2.3E-5, .0044};
  static const std::pair<G4int, const G4double*> Z95N158(158,pZ95N158);
  static const G4double pZ95N159[7]={5.E-5, 20., .029, 1.1E-9, 2., 2.3E-5, .0057};
  static const std::pair<G4int, const G4double*> Z95N159(159,pZ95N159);
  static const std::pair<G4int, const G4double*> Z95[N95]={Z95N156, Z95N157, Z95N158,
                                                           Z95N159};
  //==> Cm(Z=96)
  static const G4int N96=10;
  static const G4double pZ96N145[7]={5.E-5, 22., .027, 1.1E-9, 2.2, 2.2E-5, .006};
  static const std::pair<G4int, const G4double*> Z96N145(145,pZ96N145);
  static const G4double pZ96N146[7]={5.E-5, 24., .027, 2.E-9, 2.2, 2.2E-5, .0055};
  static const std::pair<G4int, const G4double*> Z96N146(146,pZ96N146);
  static const G4double pZ96N147[7]={5.E-5, 22., .025, 2.5E-9, 2.2, 2.4E-5, .0044};
  static const std::pair<G4int, const G4double*> Z96N147(147,pZ96N147);
  static const G4double pZ96N148[7]={5.E-5, 23., .028, 1.9E-9, 2.2, 3.E-5, .0055};
  static const std::pair<G4int, const G4double*> Z96N148(148,pZ96N148);
  static const G4double pZ96N149[7]={5.E-5, 23., .025, 1.6E-9, 3., 3.5E-5, .0054};
  static const std::pair<G4int, const G4double*> Z96N149(149,pZ96N149);
  static const G4double pZ96N150[7]={5.E-5, 24., .026, 2.E-9, 3., 3.6E-5, .0045};
  static const std::pair<G4int, const G4double*> Z96N150(150,pZ96N150);
  static const G4double pZ96N151[7]={5.E-5, 24., .022, 2.4E-9, 3., 3.6E-5, .0039};
  static const std::pair<G4int, const G4double*> Z96N151(151,pZ96N151);
  static const G4double pZ96N152[7]={6.5E-7, 25., 2.E-4, 3.4E-9, 3., 3.6E-5, .003};
  static const std::pair<G4int, const G4double*> Z96N152(152,pZ96N152);
  static const G4double pZ96N153[7]={1.6E-6, 21., 7.E-4, 1.4E-9, 3., 3.6E-5, .0045};
  static const std::pair<G4int, const G4double*> Z96N153(153,pZ96N153);
  static const G4double pZ96N154[7]={1.3E-5, 16., .016, 0., 3., 3.6E-5, .017};
  static const std::pair<G4int, const G4double*> Z96N154(154,pZ96N154);
  static const std::pair<G4int, const G4double*> Z96[N96]={Z96N145, Z96N146, Z96N147,
                                                           Z96N148, Z96N149, Z96N150,
                                                           Z96N151, Z96N152, Z96N153,
                                                           Z96N154};
  //==> Bk(Z=97)
  static const G4int N97=2;
  static const G4double pZ97N152[7]={6.5E-7, 22., 3.5E-4, 2.7E-9, 3., 4.E-5, .004};
  static const std::pair<G4int, const G4double*> Z97N152(152,pZ97N152);
  static const G4double pZ97N153[7]={6.5E-6, 22., .0036, 1.E-9, 2.7, 4.E-5, .006};
  static const std::pair<G4int, const G4double*> Z97N153(153,pZ97N153);
  static const std::pair<G4int, const G4double*> Z97[N97]={Z97N152, Z97N153};
  //==> Cf(Z=98)
  static const G4int N98=6;
  static const G4double pZ98N151[7]={6.5E-6, 22., .0035, .9E-9, 3., 4.E-5, .0068};
  static const std::pair<G4int, const G4double*> Z98N151(151,pZ98N151);
  static const G4double pZ98N152[7]={1.3E-6, 22., 7.E-4, 2.E-9, 2.7, 4.E-5, .0045};
  static const std::pair<G4int, const G4double*> Z98N152(152,pZ98N152);
  static const G4double pZ98N153[7]={2.6E-6, 22., .0014, 2.1E-9, 2.7, 4.E-5, .0044};
  static const std::pair<G4int, const G4double*> Z98N153(153,pZ98N153);
  static const G4double pZ98N154[7]={2.6E-6, 22., .0014, 1.3E-9, 2.7, 4.E-5, .0054};
  static const std::pair<G4int, const G4double*> Z98N154(154,pZ98N154);
  static const G4double pZ98N155[7]={2.6E-5, 22., .03, 0., 2.7, 4.E-5, .03};
  static const std::pair<G4int, const G4double*> Z98N155(155,pZ98N155);
  static const G4double pZ98N156[7]={5.2E-7, 22., 2.6E-4, 1.3E-9, 2.7, 4.E-5, .005};
  static const std::pair<G4int, const G4double*> Z98N156(156,pZ98N156);
  static const std::pair<G4int, const G4double*> Z98[N98]={Z98N151, Z98N152, Z98N153,
                                                           Z98N154, Z98N155, Z98N156};
 
  static const G4int NZ=99; // #of Elements covered by CHIPS elastic
  static const std::pair<G4int, const G4double*>* Pars[NZ]={Z0,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,
    Z10,Z11,Z12,Z13,Z14,Z15,Z16,Z17,Z18,Z19,Z20,Z21,Z22,Z23,Z24,Z25,Z26,Z27,Z28,Z29,Z30,
    Z31,Z32,Z33,Z34,Z35,Z36,Z37,Z38,Z39,Z40,Z41,Z42,Z43,Z44,Z45,Z46,Z47,Z48,Z49,Z50,Z51,
    Z52,Z53,Z54,Z55,Z56,Z57,Z58,Z59,Z60,Z61,Z62,Z63,Z64,Z65,Z66,Z67,Z68,Z69,Z70,Z71,Z72,
    Z73,Z74,Z75,Z76,Z77,Z78,Z79,Z80,Z81,Z82,Z83,Z84,Z85,Z86,Z87,Z88,Z89,Z90,Z91,Z92,Z93,
    Z94,Z95,Z96,Z97,Z98};
  static const G4int NIso[NZ]={N0,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16,
    N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,N36,N37,
    N38,N39,N40,N41,N42,N43,N44,N45,N46,N47,N48,N49,N50,N51,N52,N53,N54,N55,N56,N57,N58,
    N59,N60,N61,N62,N63,N64,N65,N66,N67,N68,N69,N70,N71,N72,N73,N74,N75,N76,N77,N78,N79,
    N80,N81,N82,N83,N84,N85,N86,N87,N88,N89,N90,N91,N92,N93,N94,N95,N96,N97,N98};
  if(PDG==2112)
  {
    // --- Total np elastic cross section cs & s1/b1 (t), s2/b2 (u) --- NotTuned for highE
    //p2=p*p;p3=p2*p;sp=sqrt(p);p2s=p2*sp;lp=log(p);dl1=lp-(5.=par(3));p4=p2*p2; p=|3-mom|
    //CS=12./(p2s+.05*p+.0001/sqrt(sp))+.35/p+(6.75+.14*dl1*dl1+19./p)/(1.+.6/p4);
    //  par(0)   par(1) par(2)        par(4) par(5) par(6)     par(7)     par(8)
    //s1=(6.75+.14*dl2*dl2+13./p)/(1.+.14/p4)+.6/(p4+.00013), s2=(75.+.001/p4/p)/p3
    //  par(9) par(10)   par(11)     par(12) par(13) par(14)  par(15) par(16)
    //b1=(7.2+4.32/(p4*p4+.012*p3))/(1.+2.5/p4), ss=0., b2=12./(p*sp+.34)
    //par(17) par(18)    par(19)       par(20)  par(21)   par(22)   par(23)
    //
    if(lastPAR[nLast]!=pwd) // A unique flag to avoid the repeatable definition
    {
      if ( tgZ == 1 && tgN == 0 )
      {
        for (G4int ip=0; ip<n_npel; ip++) lastPAR[ip]=np_el[ip]; // np
      }
      else if ( tgZ == 0 && tgN == 1 )
      {
        for (G4int ip=0; ip<n_ppel; ip++) lastPAR[ip]=pp_el[ip]; // nn
      }
      else
      {
        G4double a=tgZ+tgN;
        G4double ala=std::log(a); // for powers of a
        G4double sa=std::sqrt(a);
        G4double ssa=std::sqrt(sa);
        G4double asa=a*sa;
        G4double a2=a*a;
        G4double a3=a2*a;
        G4double a4=a3*a;
        G4double a5=a4*a;
        G4double a6=a4*a2;
        G4double a7=a6*a;
        G4double a8=a7*a;
        G4double a9=a8*a;
        G4double a10=a5*a5;
        G4double a12=a6*a6;
        G4double a14=a7*a7;
        G4double a16=a8*a8;
        G4double a17=a16*a;
        G4double a32=a16*a16;
        // Reaction cross-section parameters (na_el.f)
        lastPAR[ 0]=5./(1.+22./asa);                                          // p1
        lastPAR[ 1]=4.8*std::exp(ala*1.14)/(1.+3.6/a3);                       // p2
        lastPAR[ 2]=1./(1.+.004*a4)+2.E-6*a3/(1.+1.3E-6*a3);                  // p3
        lastPAR[ 3]=.07*asa/(1.+.009*a2);                                     // o4
        lastPAR[ 5]=1.7*a;                                                    // p5
        lastPAR[ 6]=5.5E-6*std::exp(ala*1.3);                                 // p6
        lastPAR[13]=0.;                                                       // reserved
        lastPAR[14]=0.;                                                       // reserved
        G4int nn=NIso[tgZ];
        G4bool nfound=true;
        if(nn) for (G4int in=0; in<nn; in++)
        {
          std::pair<G4int, const G4double*> curIs=Pars[tgZ][in];
          G4int cn=curIs.first;
          if(cn == tgN)
          {
            const G4double* curT=curIs.second;
            lastPAR[ 4]=curT[0];                                           // p4
            lastPAR[ 7]=curT[1];                                           // p7
            lastPAR[ 8]=curT[2];                                           // p8
            lastPAR[ 9]=curT[3];                                           // p9
            lastPAR[10]=curT[4];                                           // p10
            lastPAR[11]=curT[5];                                           // p11
            lastPAR[12]=curT[6];                                           // p12
            nfound = false;
            break;
          }
        }
        if(nfound) 
        {
          //AR-24Apr2018 Switch to allow transuranic elements (in this case to avoid a harmless warning)
          const G4bool isHeavyElementAllowed = true; 
          if ( ! isHeavyElementAllowed ) {
            G4cout<<"-Warning-G4ChipsNeutronElasticXS::CalcCS: Z="<<tgZ<<", N="<<tgN
                  <<" isotope is not implemented in CHIPS"<<G4endl; 
          }
          // Put default values:
          lastPAR[ 4]=5.2E-7;                                              // p4
          lastPAR[ 7]=22.;                                                 // p7
          lastPAR[ 8]=.00026;                                              // p8
          lastPAR[ 9]=1.3E-9;                                              // p9
          lastPAR[10]=2.7;                                                 // p10
          lastPAR[11]=4.E-5;                                               // p11
          lastPAR[12]=.005;                                                // p12          
        }
        // @@ the differential cross-section is parameterized separately for A>6 & A<7
        if(a<6.5)
        {
          G4double a28=a16*a12;
          // The main pre-exponent      (pel_sg)
          lastPAR[15]=4000*a;                                // p1
          lastPAR[16]=1.2e7*a8+380*a17;                      // p2
          lastPAR[17]=.7/(1.+4.e-12*a16);                    // p3
          lastPAR[18]=2.5/a8/(a4+1.e-16*a32);                // p4
          lastPAR[19]=.28*a;                                 // p5
          lastPAR[20]=1.2*a2+2.3;                            // p6
          lastPAR[21]=3.8/a;                                 // p7
          // The main slope             (pel_sl)
          lastPAR[22]=.01/(1.+.0024*a5);                     // p1
          lastPAR[23]=.2*a;                                  // p2
          lastPAR[24]=9.e-7/(1.+.035*a5);                    // p3
          lastPAR[25]=(42.+2.7e-11*a16)/(1.+.14*a);          // p4
          // The main quadratic         (pel_sh)
          lastPAR[26]=2.25*a3;                               // p1
          lastPAR[27]=18.;                                   // p2
          lastPAR[28]=.0024*a8/(1.+2.6e-4*a7);              // p3
          lastPAR[29]=3.5e-36*a32*a8/(1.+5.e-15*a32/a);      // p4
          lastPAR[30]=1.e5/(a8+2.5e12/a16);                  // p1
          lastPAR[31]=8.e7/(a12+1.e-27*a28*a28);             // p2 
          lastPAR[32]=.0006*a3;                              // p3
          // The 1st max slope          (pel_qs)
          lastPAR[33]=10.+4.e-8*a12*a;                       // p1
          lastPAR[34]=.114;                                  // p2
          lastPAR[35]=.003;                                  // p3
          lastPAR[36]=2.e-23;                                // p4
          // The effective pre-exponent (pel_ss)
          lastPAR[37]=1./(1.+.0001*a8);                      // p1
          lastPAR[38]=1.5e-4/(1.+5.e-6*a12);                 // p2
          lastPAR[39]=.03;                                   // p3
          // The effective slope        (pel_sb)
          lastPAR[40]=a/2;                                   // p1
          lastPAR[41]=2.e-7*a4;                              // p2
          lastPAR[42]=4.;                                    // p3
          lastPAR[43]=64./a3;                                // p4
          // The gloria pre-exponent    (pel_us)
          lastPAR[44]=1.e8*std::exp(.32*asa);                // p1
          lastPAR[45]=20.*std::exp(.45*asa);                 // p2
          lastPAR[46]=7.e3+2.4e6/a5;                         // p3
          lastPAR[47]=2.5e5*std::exp(.085*a3);               // p4
          lastPAR[48]=2.5*a;                                 // p5
          // The gloria slope           (pel_ub)
          lastPAR[49]=920.+.03*a8*a3;                        // p1
          lastPAR[50]=93.+.0023*a12;                         // p2
        }
        else
        {
          G4double p1a10=2.2e-28*a10;
          G4double r4a16=6.e14/a16;
          G4double s4a16=r4a16*r4a16;
          // a24
          // a36
          // The main pre-exponent      (peh_sg)
          lastPAR[15]=4.5*std::pow(a,1.15);                  // p1
          lastPAR[16]=.06*std::pow(a,.6);                    // p2
          lastPAR[17]=.6*a/(1.+2.e15/a16);                   // p3
          lastPAR[18]=.17/(a+9.e5/a3+1.5e33/a32);            // p4
          lastPAR[19]=(.001+7.e-11*a5)/(1.+4.4e-11*a5);      // p5
          lastPAR[20]=(p1a10*p1a10+2.e-29)/(1.+2.e-22*a12);  // p6
          // The main slope             (peh_sl)
          lastPAR[21]=400./a12+2.e-22*a9;                    // p1
          lastPAR[22]=1.e-32*a12/(1.+5.e22/a14);             // p2
          lastPAR[23]=1000./a2+9.5*sa*ssa;                   // p3
          lastPAR[24]=4.e-6*a*asa+1.e11/a16;                 // p4
          lastPAR[25]=(120./a+.002*a2)/(1.+2.e14/a16);       // p5
          lastPAR[26]=9.+100./a;                             // p6
          // The main quadratic         (peh_sh)
          lastPAR[27]=.002*a3+3.e7/a6;                       // p1
          lastPAR[28]=7.e-15*a4*asa;                         // p2
          lastPAR[29]=9000./a4;                              // p3
          // The 1st max pre-exponent   (peh_qq)
          lastPAR[30]=.0011*asa/(1.+3.e34/a32/a4);           // p1
          lastPAR[31]=1.e-5*a2+2.e14/a16;                    // p2
          lastPAR[32]=1.2e-11*a2/(1.+1.5e19/a12);            // p3
          lastPAR[33]=.016*asa/(1.+5.e16/a16);               // p4
          // The 1st max slope          (peh_qs)
          lastPAR[34]=.002*a4/(1.+7.e7/std::pow(a-6.83,14)); // p1
          lastPAR[35]=2.e6/a6+7.2/std::pow(a,.11);           // p2
          lastPAR[36]=11.*a3/(1.+7.e23/a16/a8);              // p3
          lastPAR[37]=100./asa;                              // p4
          // The 2nd max pre-exponent   (peh_ss)
          lastPAR[38]=(.1+4.4e-5*a2)/(1.+5.e5/a4);           // p1
          lastPAR[39]=3.5e-4*a2/(1.+1.e8/a8);                // p2
          lastPAR[40]=1.3+3.e5/a4;                           // p3
          lastPAR[41]=500./(a2+50.)+3;                       // p4
          lastPAR[42]=1.e-9/a+s4a16*s4a16;                   // p5
          // The 2nd max slope          (peh_sb)
          lastPAR[43]=.4*asa+3.e-9*a6;                       // p1
          lastPAR[44]=.0005*a5;                              // p2
          lastPAR[45]=.002*a5;                               // p3
          lastPAR[46]=10.;                                   // p4
          // The effective pre-exponent (peh_us)
          lastPAR[47]=.05+.005*a;                            // p1
          lastPAR[48]=7.e-8/sa;                              // p2
          lastPAR[49]=.8*sa;                                 // p3
          lastPAR[50]=.02*sa;                                // p4
          lastPAR[51]=1.e8/a3;                               // p5
          lastPAR[52]=3.e32/(a32+1.e32);                     // p6
          // The effective slope        (peh_ub)
          lastPAR[53]=24.;                                   // p1
          lastPAR[54]=20./sa;                                // p2
          lastPAR[55]=7.e3*a/(sa+1.);                        // p3
          lastPAR[56]=900.*sa/(1.+500./a3);                  // p4
        }
        // Parameter for lowEnergyNeutrons
        lastPAR[57]=1.e15+2.e27/a4/(1.+2.e-18*a16);
      }
      lastPAR[nLast]=pwd;
      // and initialize the zero element of the table
      G4double lp=lPMin;                                      // ln(momentum)
      G4bool memCS=onlyCS;                                    // ??
      onlyCS=false;
      lastCST[0]=GetTabValues(lp, PDG, tgZ, tgN); // Calculate AMDB tables
      onlyCS=memCS;
      lastSST[0]=theSS;
      lastS1T[0]=theS1;
      lastB1T[0]=theB1;
      lastS2T[0]=theS2;
      lastB2T[0]=theB2;
      lastS3T[0]=theS3;
      lastB3T[0]=theB3;
      lastS4T[0]=theS4;
      lastB4T[0]=theB4;
    }
    if(LP>ILP)
    {
      G4int ini = static_cast<int>((ILP-lPMin+.000001)/dlnP)+1; // already inited till this
      if(ini<0) ini=0;
      if(ini<nPoints)
      {
        G4int fin = static_cast<int>((LP-lPMin)/dlnP)+1; // final bin of initialization
        if(fin>=nPoints) fin=nLast;               // Limit of the tabular initialization
        if(fin>=ini)
        {
          G4double lp=0.;
          for(G4int ip=ini; ip<=fin; ip++)        // Calculate tabular CS,S1,B1,S2,B2,S3,B3
          {
            lp=lPMin+ip*dlnP;                     // ln(momentum)
            G4bool memCS=onlyCS;
            onlyCS=false;
            lastCST[ip]=GetTabValues(lp, PDG, tgZ, tgN); // Calculate AMDB tables (ret CS)
            onlyCS=memCS;
            lastSST[ip]=theSS;
            lastS1T[ip]=theS1;
            lastB1T[ip]=theB1;
            lastS2T[ip]=theS2;
            lastB2T[ip]=theB2;
            lastS3T[ip]=theS3;
            lastB3T[ip]=theB3;
            lastS4T[ip]=theS4;
            lastB4T[ip]=theB4;
          }
          return lp;
        }
        else G4cout<<"*Warning*G4ChipsNeutronElasticXS::GetPTables: PDG="<<PDG
                   <<", Z="<<tgZ<<", N="<<tgN<<", i="<<ini<<" > fin="<<fin<<", LP="<<LP
                   <<" > ILP="<<ILP<<" nothing is done!"<<G4endl;
      }
      else G4cout<<"*Warning*G4ChipsNeutronElasticXS::GetPTables: PDG="<<PDG<<", Z="
                 <<tgZ<<", N="<<tgN<<", i="<<ini<<">= max="<<nPoints<<", LP="<<LP
                 <<" > ILP="<<ILP<<", lPMax="<<lPMax<<" nothing is done!"<<G4endl;
    }
  }
  else
  {
    // G4cout<<"*Error*G4ChipsNeutronElasticXS::GetPTables: PDG="<<PDG<<", Z="<<tgZ
    //       <<", N="<<tgN<<", while it is defined only for PDG=2112(n)"<<G4endl;
    // throw G4QException("G4ChipsNeutronElasticXS::GetPTables:only nA're implemented");
    G4ExceptionDescription ed;
    ed << "PDG = " << PDG << ", Z = " << tgZ <<", N = " << tgN
       << ", while it is defined only for PDG=2112 (n)" << G4endl;
    G4Exception("G4ChipsNeutronElasticXS::GetPTables()", "HAD_CHPS_0000",
                FatalException, ed);
  }
  return ILP;
}

// Returns Q2=-t in independent units (MeV^2) (all internal calculations are in GeV)
G4double G4ChipsNeutronElasticXS::GetExchangeT(G4int tgZ, G4int tgN, G4int PDG)
{
  static const G4double GeVSQ=gigaelectronvolt*gigaelectronvolt;
  static const G4double third=1./3.;
  static const G4double fifth=1./5.;
  static const G4double sevth=1./7.;
  if(PDG!=2112) G4cout<<"*Warning*G4ChipsNeutronElasticXS::GetExT:PDG="<<PDG<<G4endl;
  if(onlyCS) G4cout<<"*Warning*G4ChipsNeutronElasticXS::GetExchangeT:onCS=1"<<G4endl;
  if(lastLP<-4.3) return lastTM*GeVSQ*G4UniformRand();// S-wave for p<14 MeV/c (kinE<.1MeV)
  G4double q2=0.;
  if(tgZ==1 && tgN==0)                // ===> n+p=n+p
  {
    G4double E1=lastTM*theB1;
    G4double R1=(1.-std::exp(-E1));
    G4double E2=lastTM*theB2;
    G4double R2=(1.-std::exp(-E2));
    G4double I1=R1*theS1;
    G4double I2=R2*theS2/theB2;
    //G4double I3=R3*theS3/theB3;
    G4double I12=I1+I2;
    //G4double rand=(I12+I3)*G4UniformRand();
    G4double rand=I12*G4UniformRand();
    if     (rand<I1 )
    {
      G4double ran=R1*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB1;       // t-chan
    }
    else
    {
      G4double ran=R2*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=lastTM+std::log(1.-ran)/theB2; // u-chan (ChEx)
    }
  }
  else
  {
    G4double a=tgZ+tgN;
    G4double E1=lastTM*(theB1+lastTM*theSS);
    G4double R1=(1.-std::exp(-E1));
    G4double tss=theSS+theSS; // for future solution of quadratic equation (imediate check)
    G4double tm2=lastTM*lastTM;
    G4double E2=lastTM*tm2*theB2;                   // power 3 for lowA, 5 for HighA (1st)
    if(a>6.5)E2*=tm2;                               // for heavy nuclei
    G4double R2=(1.-std::exp(-E2));
    G4double E3=lastTM*theB3;
    if(a>6.5)E3*=tm2*tm2*tm2;                       // power 1 for lowA, 7 (2nd) for HighA
    G4double R3=(1.-std::exp(-E3));
    G4double E4=lastTM*theB4;
    G4double R4=(1.-std::exp(-E4));
    G4double I1=R1*theS1;
    G4double I2=R2*theS2;
    G4double I3=R3*theS3;
    G4double I4=R4*theS4;
    G4double I12=I1+I2;
    G4double I13=I12+I3;
    G4double rand=(I13+I4)*G4UniformRand();
    if(rand<I1)
    {
      G4double ran=R1*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB1;
      if(std::fabs(tss)>1.e-7) q2=(std::sqrt(theB1*(theB1+(tss+tss)*q2))-theB1)/tss;
    }
    else if(rand<I12)
    {
      G4double ran=R2*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB2;
      if(q2<0.) q2=0.;
      if(a<6.5) q2=std::pow(q2,third);
      else      q2=std::pow(q2,fifth);
    }
    else if(rand<I13)
    {
      G4double ran=R3*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB3;
      if(q2<0.) q2=0.;
      if(a>6.5) q2=std::pow(q2,sevth);
    }
    else
    {
      G4double ran=R4*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB4;
      if(a<6.5) q2=lastTM-q2;                    // u reduced for lightA (starts from 0)
    }
  }
  if(q2<0.) q2=0.;
  if(!(q2>=-1.||q2<=1.)) G4cout<<"*NAN*G4QNeutronElCroSect::GetExchangeT: -t="<<q2<<G4endl;
  if(q2>lastTM)
  {
    q2=lastTM;
  }
  return q2*GeVSQ;
}

// Returns B in independent units (MeV^-2) (all internal calculations are in GeV) see ExT
G4double G4ChipsNeutronElasticXS::GetSlope(G4int tgZ, G4int tgN, G4int PDG)
{
  static const G4double GeVSQ=gigaelectronvolt*gigaelectronvolt;

  if(onlyCS) G4cout<<"Warning*G4ChipsNeutronElasticXS::GetSlope:onlyCS=true"<<G4endl;
  if(lastLP<-4.3) return 0.;          // S-wave for p<14 MeV/c (kinE<.1MeV)
  if(PDG!=2112)
  {
    // G4cout<<"*Error*G4ChipsNeutronElasticXS::GetSlope: PDG="<<PDG<<", Z="<<tgZ
    //       <<", N="<<tgN<<", while it is defined only for PDG=2112"<<G4endl;
    // throw G4QException("G4ChipsNeutronElasticXS::GetSlope: only nA are implemented");
    G4ExceptionDescription ed;
    ed << "PDG = " << PDG << ", Z = " << tgZ << ", N = " << tgN
       <<", while it is defined only for PDG=2112 (n) " << G4endl;
    G4Exception("G4ChipsNeutronElasticXS::GetSlope()", "HAD_CHPS_0000",
                FatalException, ed);
  }
  if(theB1<0.) theB1=0.;
  if(!(theB1>=-1.||theB1<=1.))G4cout<<"*NAN*G4QNeutElasticCrosS::Getslope:"<<theB1<<G4endl;
  return theB1/GeVSQ;
}

// Returns half max(Q2=-t) in independent units (MeV^2)
G4double G4ChipsNeutronElasticXS::GetHMaxT()
{
  static const G4double HGeVSQ=gigaelectronvolt*gigaelectronvolt/2.;
  return lastTM*HGeVSQ;
}

// lastLP is used, so calculating tables, one need to remember and then recover lastLP
G4double G4ChipsNeutronElasticXS::GetTabValues(G4double lp, G4int PDG, G4int tgZ,
                                                     G4int tgN)
{
  if(PDG!=2112) G4cout<<"*Warning*G4ChipsNeutronElasticXS::GetTaV:PDG="<<PDG<<G4endl;

  //AR-24Apr2018 Switch to allow transuranic elements
  const G4bool isHeavyElementAllowed = true;
  if(tgZ<0 || ( !isHeavyElementAllowed && tgZ>92))
  {
    G4cout<<"*Warning*G4QNElasticCrS::GetTabValue: (1-92) No isotopes for Z="<<tgZ<<G4endl;
    return 0.;
  }
  G4int iZ=tgZ-1; // Z index
  if(iZ<0)
  {
    iZ=0;         // conversion of the neutron target to the proton target
    tgZ=1;
    tgN=0;
  }
  G4double p=std::exp(lp);              // momentum
  G4double sp=std::sqrt(p);             // sqrt(p)
  G4double p2=p*p;            
  G4double p3=p2*p;
  G4double p4=p3*p;
  if ( tgZ == 1 && tgN == 0 )  // np
  {
    G4double ssp=std::sqrt(sp);           // sqrt(sqrt(p))=p^.25
    G4double p2s=p2*sp;
    G4double dl1=lp-lastPAR[3];
    theSS=lastPAR[27];
    theS1=(lastPAR[9]+lastPAR[10]*dl1*dl1+lastPAR[11]/p)/(1.+lastPAR[12]/p4)
          +lastPAR[13]/(p4+lastPAR[14]);
    theB1=(lastPAR[17]+lastPAR[18]/(p4*p4+lastPAR[19]*p3))/(1.+lastPAR[20]/p4);
    theS2=(lastPAR[15]+lastPAR[16]/p4/p)/p3;
    theB2=lastPAR[22]/(p*sp+lastPAR[23]); 
    theS3=0.;
    theB3=0.; 
    theS4=0.;
    theB4=0.; 
    // Returns the total elastic pp cross-section (to avoid spoiling lastSIG)
    return lastPAR[0]/(p2s+lastPAR[1]*p+lastPAR[2]/ssp)+lastPAR[4]/p
           +(lastPAR[5]+lastPAR[6]*dl1*dl1+lastPAR[7]/p)/(1.+lastPAR[8]/p4);

  }
  else
  {
    G4double p5=p4*p;
    G4double p6=p5*p;
    G4double p8=p6*p2;
    G4double p10=p8*p2;
    G4double p12=p10*p2;
    G4double p16=p8*p8;
    G4double dl=lp-5.;
    G4double a=tgZ+tgN;
    if(a<6.5)
    {
    G4double pah=std::pow(p,a/2);
    G4double pa=pah*pah;
    G4double pa2=pa*pa;

      theS1=lastPAR[15]/(1.+lastPAR[16]*p4*pa)+lastPAR[17]/(p4+lastPAR[18]*p4/pa2)+
            (lastPAR[19]*dl*dl+lastPAR[20])/(1.+lastPAR[21]/p2);
      theB1=(lastPAR[22]+lastPAR[23]*p2)/(p4+lastPAR[24]/pah)+lastPAR[25];
      theSS=lastPAR[26]/(1.+lastPAR[27]/p2)+lastPAR[28]/(p6/pa+lastPAR[29]/p16);
      theS2=lastPAR[30]/(pa/p2+lastPAR[31]/p4)+lastPAR[32];
      theB2=lastPAR[33]*std::pow(p,lastPAR[34])+lastPAR[35]/(p8+lastPAR[36]/p16);
      theS3=lastPAR[37]/(pa*p+lastPAR[38]/pa)+lastPAR[39];
      theB3=lastPAR[40]/(p3+lastPAR[41]/p6)+lastPAR[42]/(1.+lastPAR[43]/p2);
      theS4=p2*(pah*lastPAR[44]*std::exp(-pah*lastPAR[45])+
                lastPAR[46]/(1.+lastPAR[47]*std::pow(p,lastPAR[48])));
      theB4=lastPAR[49]*pa/p2/(1.+pa*lastPAR[50]);
    }
    else
    {
      theS1=lastPAR[15]/(1.+lastPAR[16]/p4)+lastPAR[17]/(p4+lastPAR[18]/p2)+
            lastPAR[19]/(p5+lastPAR[20]/p16);
      theB1=(lastPAR[21]/p8+lastPAR[25])/(p+lastPAR[22]/std::pow(p,lastPAR[26]))+
            lastPAR[23]/(1.+lastPAR[24]/p4);
      theSS=lastPAR[27]/(p4/std::pow(p,lastPAR[29])+lastPAR[28]/p4);
      theS2=lastPAR[30]/p4/(std::pow(p,lastPAR[31])+lastPAR[32]/p12)+lastPAR[33];
      theB2=lastPAR[34]/std::pow(p,lastPAR[35])+lastPAR[36]/std::pow(p,lastPAR[37]);
      theS3=lastPAR[38]/std::pow(p,lastPAR[41])/(1.+lastPAR[42]/p12)+
            lastPAR[39]/(1.+lastPAR[40]/p6);
      theB3=lastPAR[43]/p8+lastPAR[44]/p2+lastPAR[45]/(1.+lastPAR[46]/p8);
      theS4=(lastPAR[47]/p4+lastPAR[52]/p)/(1.+lastPAR[48]/p10)+
            (lastPAR[49]+lastPAR[50]*dl*dl)/(1.+lastPAR[51]/p12);
      theB4=lastPAR[53]/(1.+lastPAR[54]/p)+lastPAR[55]*p4/(1.+lastPAR[56]*p5);
    }
    // Returns the total elastic (n/p)A cross-section (to avoid spoiling lastSIG)
    //         p1(p6)          p2(p7)          p3(p4)       o4(p8)       (p9)p5
    return (lastPAR[0]*dl*dl+lastPAR[1])/(1.+lastPAR[2]/p+lastPAR[3]/p4)+lastPAR[5]/
	     (p3+lastPAR[6]/p3)+lastPAR[7]/(p2+lastPAR[4]/(p2+lastPAR[8])+lastPAR[9]/p)+
           lastPAR[10]/(p5+lastPAR[11]/p2)+lastPAR[12]/p;
    //         p10             p11             p12
  }
  return 0.;
} // End of GetTableValues

// Returns max -t=Q2 (GeV^2) for the momentum pP(GeV) and the target nucleus (tgN,tgZ)
G4double G4ChipsNeutronElasticXS::GetQ2max(G4int PDG, G4int tgZ, G4int tgN,
                                                 G4double pP)
{

  G4double pP2=pP*pP;                                 // squared momentum of the projectile
  if(tgZ==0 && tgN==1)
  {
    G4double tMid=std::sqrt(pP2+mNeut2)*mNeut-mNeut2;  // CMS 90deg value of -t=Q2 (GeV^2)
    return tMid+tMid;
  }
  else if(tgZ || tgN)                   // ---> nA
  {
    G4double mt=mProt;                                 // Target mass in GeV
    if(tgN||tgZ>1) mt=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(tgZ,tgZ+tgN,0)->GetPDGMass()*.001; // Target mass in GeV
    G4double dmt=mt+mt;
    G4double mds=dmt*std::sqrt(pP2+mNeut2)+mNeut2+mt*mt; // Mondelstam mds (in GeV^2)
    return dmt*dmt*pP2/mds;
  }
  else
  {
    // G4cout<<"*Error*G4ChipsNeutronElasticXS::GetQ2max:PDG="<<PDG<<", Z="<<tgZ<<", N="
    //       <<tgN<<", while it is defined only for n projectiles & Z_target>0"<<G4endl;
    // throw G4QException("G4ChipsNeutronElasticXS::GetQ2max: only nA implemented");
    G4ExceptionDescription ed;
    ed << "PDG = " << PDG << ", Z = " << tgZ << ", N =" << tgN
       <<", while it is defined only for n projectiles & Z_target>0" << G4endl;
    G4Exception("G4ChipsNeutronElasticXS::GetQ2max()", "HAD_CHPS_0000",
                FatalException, ed);
    return 0;
  }
}
