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
// $Id: G4QAntiBaryonElasticCrossSection.cc,v 1.2 2010-02-16 07:53:05 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Physics class: G4QAntiBaryonElasticCrossSection for pA elastic cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 5-Feb-2010
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 5-Feb-2010
// 
//================================================================================
// Short description: Interaction cross-sections for the G4QElastic process
// -------------------------------------------------------------------------------

//#define debug
//#define isodebug
//#define pdebug
//#define ppdebug
//#define tdebug
//#define sdebug

#include "G4QAntiBaryonElasticCrossSection.hh"

// Initialization of the static parameters
const G4int G4QAntiBaryonElasticCrossSection::nPoints=128;//#ofPt in AMDB table(>anyPar)(D)
const G4int G4QAntiBaryonElasticCrossSection::nLast=nPoints-1;// theLastElement inTable (D)
G4double  G4QAntiBaryonElasticCrossSection::lPMin=-8.;  //Min tabulatedLogarithmMomentum(D)
G4double  G4QAntiBaryonElasticCrossSection::lPMax= 8.;  //Max tabulatedLogarithmMomentum(D)
G4double  G4QAntiBaryonElasticCrossSection::dlnP=(lPMax-lPMin)/nLast;// LogStep inTable (D)
G4bool    G4QAntiBaryonElasticCrossSection::onlyCS=true;//Flag toCalculOnlyCS(not Si/Bi)(L)
G4double  G4QAntiBaryonElasticCrossSection::lastSIG=0.; //Last calculated cross section (L)
G4double  G4QAntiBaryonElasticCrossSection::lastLP=-10.;//LastLog(mom_of IncidentHadron)(L)
G4double  G4QAntiBaryonElasticCrossSection::lastTM=0.;  //Last t_maximum                (L)
G4double  G4QAntiBaryonElasticCrossSection::theSS=0.;   //TheLastSqSlope of 1st difr.Max(L)
G4double  G4QAntiBaryonElasticCrossSection::theS1=0.;   //TheLastMantissa of 1st difrMax(L)
G4double  G4QAntiBaryonElasticCrossSection::theB1=0.;   //TheLastSlope of 1st difructMax(L)
G4double  G4QAntiBaryonElasticCrossSection::theS2=0.;   //TheLastMantissa of 2nd difrMax(L)
G4double  G4QAntiBaryonElasticCrossSection::theB2=0.;   //TheLastSlope of 2nd difructMax(L)
G4double  G4QAntiBaryonElasticCrossSection::theS3=0.;   //TheLastMantissa of 3d difr.Max(L)
G4double  G4QAntiBaryonElasticCrossSection::theB3=0.;   //TheLastSlope of 3d difruct.Max(L)
G4double  G4QAntiBaryonElasticCrossSection::theS4=0.;   //TheLastMantissa of 4th difrMax(L)
G4double  G4QAntiBaryonElasticCrossSection::theB4=0.;   //TheLastSlope of 4th difructMax(L)
G4int     G4QAntiBaryonElasticCrossSection::lastTZ=0;   // Last atomic number of the target
G4int     G4QAntiBaryonElasticCrossSection::lastTN=0;   // Last # of neutrons in the target
G4double  G4QAntiBaryonElasticCrossSection::lastPIN=0.; // Last initialized max momentum
G4double* G4QAntiBaryonElasticCrossSection::lastCST=0;  // Elastic cross-section table
G4double* G4QAntiBaryonElasticCrossSection::lastPAR=0;  // ParametersForFunctionCalculation
G4double* G4QAntiBaryonElasticCrossSection::lastSST=0;  // E-dep ofSqardSlope of 1st difMax
G4double* G4QAntiBaryonElasticCrossSection::lastS1T=0;  // E-dep of mantissa of 1st dif.Max
G4double* G4QAntiBaryonElasticCrossSection::lastB1T=0;  // E-dep of the slope of 1st difMax
G4double* G4QAntiBaryonElasticCrossSection::lastS2T=0;  // E-dep of mantissa of 2nd difrMax
G4double* G4QAntiBaryonElasticCrossSection::lastB2T=0;  // E-dep of the slope of 2nd difMax
G4double* G4QAntiBaryonElasticCrossSection::lastS3T=0;  // E-dep of mantissa of 3d difr.Max
G4double* G4QAntiBaryonElasticCrossSection::lastB3T=0;  // E-dep of the slope of 3d difrMax
G4double* G4QAntiBaryonElasticCrossSection::lastS4T=0;  // E-dep of mantissa of 4th difrMax
G4double* G4QAntiBaryonElasticCrossSection::lastB4T=0;  // E-dep of the slope of 4th difMax
G4int     G4QAntiBaryonElasticCrossSection::lastN=0;    // The last N of calculated nucleus
G4int     G4QAntiBaryonElasticCrossSection::lastZ=0;    // The last Z of calculated nucleus
G4double  G4QAntiBaryonElasticCrossSection::lastP=0.;   // LastUsed inCrossSection Momentum
G4double  G4QAntiBaryonElasticCrossSection::lastTH=0.;  // Last threshold momentum
G4double  G4QAntiBaryonElasticCrossSection::lastCS=0.;  // Last value of the Cross Section
G4int     G4QAntiBaryonElasticCrossSection::lastI=0;    // The last position in the DAMDB

std::vector<G4double*> G4QAntiBaryonElasticCrossSection::PAR;// Vector ofParsForFunctCalcul
std::vector<G4double*> G4QAntiBaryonElasticCrossSection::CST;// Vector ofCrossSection table
std::vector<G4double*> G4QAntiBaryonElasticCrossSection::SST;// Vector ofThe1st SquardSlope
std::vector<G4double*> G4QAntiBaryonElasticCrossSection::S1T;// Vector of the 1st mantissa
std::vector<G4double*> G4QAntiBaryonElasticCrossSection::B1T;// Vector of the1st slope
std::vector<G4double*> G4QAntiBaryonElasticCrossSection::S2T;// Vector of the2nd mantissa
std::vector<G4double*> G4QAntiBaryonElasticCrossSection::B2T;// Vector of the2nd slope
std::vector<G4double*> G4QAntiBaryonElasticCrossSection::S3T;// Vector of the3d mantissa
std::vector<G4double*> G4QAntiBaryonElasticCrossSection::B3T;// Vector of the3d slope
std::vector<G4double*> G4QAntiBaryonElasticCrossSection::S4T;// Vector ofThe4th mantissa(g)
std::vector<G4double*> G4QAntiBaryonElasticCrossSection::B4T;// Vector ofThe4th slope(glor)

G4QAntiBaryonElasticCrossSection::G4QAntiBaryonElasticCrossSection()
{
}

G4QAntiBaryonElasticCrossSection::~G4QAntiBaryonElasticCrossSection()
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

// Returns Pointer to the G4VQCrossSection class
G4VQCrossSection* G4QAntiBaryonElasticCrossSection::GetPointer()
{
  static G4QAntiBaryonElasticCrossSection theCrossSection;//StaticBody ofTheQElCrossSection
  return &theCrossSection;
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4QAntiBaryonElasticCrossSection::GetCrossSection(G4bool fCS, G4double pMom,
                                                         G4int tgZ, G4int tgN, G4int pPDG)
{
  static std::vector <G4int>    colN;  // Vector of N for calculated nuclei (isotops)
  static std::vector <G4int>    colZ;  // Vector of Z for calculated nuclei (isotops)
  static std::vector <G4double> colP;  // Vector of last momenta for the reaction
  static std::vector <G4double> colTH; // Vector of energy thresholds for the reaction
  static std::vector <G4double> colCS; // Vector of last cross sections for the reaction
  // ***---*** End of the mandatory Static Definitions of the Associative Memory ***---***
  G4double pEn=pMom;
  onlyCS=fCS;
#ifdef pdebug
  G4cout<<"G4QaBaElCS::GetCS:>>> f="<<fCS<<", p="<<pMom<<", Z="<<tgZ<<"("<<lastZ<<") ,N="
        <<tgN<<"("<<lastN<<"), T="<<pEn<<"("<<lastTH<<")"<<",Sz="<<colN.size()<<G4endl;
  //CalculateCrossSection(fCS,-27,j,pPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  if(pPDG<-3334 || pPDG>-1111)
  {
    G4cout<<"*Warning*G4QAntiBaryElCS::GetCS:**> Found pPDG="<<pPDG<<" ====> CS=0"<<G4endl;
    //CalculateCrossSection(fCS,-27,j,pPDG,lastZ,lastN,pMom); // DUMMY TEST
    return 0.;                         // projectile PDG=0 is a mistake (?!) @@
  }
  G4bool in=false;                   // By default the isotope must be found in the AMDB
  lastP   = 0.;                      // New momentum history (nothing to compare with)
  lastN   = tgN;                     // The last N of the calculated nucleus
  lastZ   = tgZ;                     // The last Z of the calculated nucleus
  lastI   = colN.size();             // Size of the Associative Memory DB in the heap
  if(lastI) for(G4int i=0; i<lastI; i++) // Loop over proj/tgZ/tgN lines of DB
  {                                  // The nucleus with projPDG is found in AMDB
    if(colN[i]==tgN && colZ[i]==tgZ) // Isotope is foind in AMDB
    {
      lastI=i;
      lastTH =colTH[i];              // Last THreshold (A-dependent)
#ifdef pdebug
      G4cout<<"G4QaBElCS::GetCS:*Found* P="<<pMom<<",Threshold="<<lastTH<<",i="<<i<<G4endl;
      //CalculateCrossSection(fCS,-27,i,pPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      if(pEn<=lastTH)
      {
#ifdef pdebug
        G4cout<<"G4QaBElCS::GetCS:Found T="<<pEn<<" < Threshold="<<lastTH<<",CS=0"<<G4endl;
        //CalculateCrossSection(fCS,-27,i,pPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        return 0.;                   // Energy is below the Threshold value
      }
      lastP  =colP [i];              // Last Momentum  (A-dependent)
      lastCS =colCS[i];              // Last CrossSect (A-dependent)
      //  if(std::fabs(lastP/pMom-1.)<tolerance) //VI (do not use tolerance)
      if(lastP == pMom)              // Do not recalculate
      {
#ifdef pdebug
        G4cout<<"G4QAntiBaryonElasticCS::GetCS:P="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
#endif
        CalculateCrossSection(fCS,-1,i,pPDG,lastZ,lastN,pMom); // Update param's only
        return lastCS*millibarn;     // Use theLastCS
      }
      in = true;                       // This is the case when the isotop is found in DB
      // Momentum pMom is in IU ! @@ Units
#ifdef pdebug
      G4cout<<"G4QaBElCS::G:UpdateDB P="<<pMom<<",f="<<fCS<<",I="<<lastI<<",i="<<i<<G4endl;
#endif
      lastCS=CalculateCrossSection(fCS,-1,i,pPDG,lastZ,lastN,pMom); // read & update
#ifdef pdebug
      G4cout<<"G4QAntiBElCS::GetCrosSec: *****> New (inDB) Calculated CS="<<lastCS<<G4endl;
      //CalculateCrossSection(fCS,-27,i,pPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      if(lastCS<=0. && pEn>lastTH)    // Correct the threshold
      {
#ifdef pdebug
        G4cout<<"G4QAntiBarElCS::GetCS: NewT="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
        lastTH=pEn;
      }
      break;                           // Go out of the LOOP with found lastI
    }
#ifdef pdebug
    G4cout<<"---G4QAntiBarElCrossSection::GetCrosSec:pPDG="<<pPDG<<",i="<<i<<",N="<<colN[i]
          <<",Z["<<i<<"]="<<colZ[i]<<G4endl;
    //CalculateCrossSection(fCS,-27,i,pPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  } // End of attampt to find the nucleus in DB
  if(!in)                            // This nucleus has not been calculated previously
  {
#ifdef pdebug
    G4cout<<"G4QaBElCS::GetCrosSec:CalcNew P="<<pMom<<",f="<<fCS<<",lastI="<<lastI<<G4endl;
#endif
    //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
    lastCS=CalculateCrossSection(fCS,0,lastI,pPDG,lastZ,lastN,pMom);//calculate&create
    if(lastCS<=0.)
    {
      lastTH = ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last
#ifdef pdebug
      G4cout<<"G4QaBaElCrossSection::GetCrossSect: NewThresh="<<lastTH<<",T="<<pEn<<G4endl;
#endif
      if(pEn>lastTH)
      {
#ifdef pdebug
        G4cout<<"G4QaBarElCS::GetCS: First T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
        lastTH=pEn;
      }
    }
#ifdef pdebug
    G4cout<<"G4QaBElCS::GetCrosSec: New CS="<<lastCS<<",lZ="<<lastN<<",lN="<<lastZ<<G4endl;
    //CalculateCrossSection(fCS,-27,lastI,pPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    colN.push_back(tgN);
    colZ.push_back(tgZ);
    colP.push_back(pMom);
    colTH.push_back(lastTH);
    colCS.push_back(lastCS);
#ifdef pdebug
    G4cout<<"G4QaBElCS::GetCS:1st,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
    //CalculateCrossSection(fCS,-27,lastI,pPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return lastCS*millibarn;
  } // End of creation of the new set of parameters
  else
  {
#ifdef pdebug
    G4cout<<"G4QAntiBaryonElasticCrossSection::GetCS: Update lastI="<<lastI<<G4endl;
#endif
    colP[lastI]=pMom;
    colCS[lastI]=lastCS;
  }
#ifdef pdebug
  G4cout<<"G4QaBElCS::GetCSec:End,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
  //CalculateCrossSection(fCS,-27,lastI,pPDG,lastZ,lastN,pMom); // DUMMY TEST
  G4cout<<"G4QAntiBaryonElasticCrossSection::GetCrosSec:***End***,onlyCS="<<onlyCS<<G4endl;
#endif
  return lastCS*millibarn;
}

// Calculation of total elastic cross section (p in IU, CS in mb) @@ Units (?)
// F=0 - create AMDB, F=-1 - read&update AMDB, F=1 - update AMDB (sinchro with higher AMDB)
G4double G4QAntiBaryonElasticCrossSection::CalculateCrossSection(G4bool CS,G4int F,G4int I,
                                             G4int PDG, G4int tgZ, G4int tgN, G4double pIU)
{
  // *** Begin of Associative Memory DB for acceleration of the cross section calculations
  static std::vector <G4double>  PIN;   // Vector of max initialized log(P) in the table
  // *** End of Static Definitions (Associative Memory Data Base) ***
  G4double pMom=pIU/GeV;                // All calculations are in GeV
  onlyCS=CS;                            // Flag to calculate only CS (not Si/Bi)
#ifdef pdebug
  G4cout<<"G4QAntiBaryonElasticCS::CalcCS:->onlyCS="<<onlyCS<<",F="<<F<<",p="<<pIU<<G4endl;
#endif
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
#ifdef pdebug
      G4cout<<"G4QAntiBarElCS::CalcCS: DB's updated for I="<<I<<",*,PIN4="<<PIN[4]<<G4endl;
#endif
    }
#ifdef pdebug
    G4cout<<"G4QAntiBaryonElasticCS::CalcCS:*read*, LP="<<lastLP<<",PIN="<<lastPIN<<G4endl;
#endif
    if(lastLP>lastPIN && lastLP<lPMax)
    {
      lastPIN=GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);// Can update upper logP-Limit in tabs
#ifdef pdebug
      G4cout<<"G4QaBElCS::CalcCS:updated(I),LP="<<lastLP<<"<IN["<<I<<"]="<<lastPIN<<G4endl;
#endif
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
#ifdef pdebug
    G4cout<<"G4QAntiBaryonElasticCS::CalcCS:*ini*,lastLP="<<lastLP<<",min="<<lPMin<<G4endl;
#endif
    lastPIN = GetPTables(lastLP,lPMin,PDG,tgZ,tgN); // Returns the new P-limit for tables
#ifdef pdebug
    G4cout<<"G4QAnBaElCS::CCS:i,Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<",LP"<<lastPIN<<G4endl;
#endif
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
  // ============= NOW Update (if necessary) and Calculate the Cross Section ===========
#ifdef pdebug
  G4cout<<"G4QaBElCS::CalcCS:?update?,LP="<<lastLP<<",IN="<<lastPIN<<",ML="<<lPMax<<G4endl;
#endif
  if(lastLP>lastPIN && lastLP<lPMax)
  {
    lastPIN = GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);
#ifdef pdebug
    G4cout<<"G4QAntiBarElCS::CalcCS: *updated(O)*, LP="<<lastLP<<" < IN="<<lastPIN<<G4endl;
#endif
  }
#ifdef pdebug
  G4cout<<"G4QaBaElCS::CalcCS: lastLP="<<lastLP<<",lPM="<<lPMin<<",lPIN="<<lastPIN<<G4endl;
#endif
  if(!onlyCS) lastTM=GetQ2max(PDG, tgZ, tgN, pMom); // Calculate (-t)_max=Q2_max (GeV2)
#ifdef pdebug
  G4cout<<"G4QAnBarElCrosSec::CalcCS:oCS="<<onlyCS<<",-t="<<lastTM<<", p="<<lastLP<<G4endl;
#endif
  if(lastLP>lPMin && lastLP<=lastPIN)   // Linear fit is made using precalculated tables
  {
    if(lastLP==lastPIN)
    {
      G4double shift=(lastLP-lPMin)/dlnP+.000001; // Log distance from lPMin
      G4int    blast=static_cast<int>(shift); // this is a bin number of the lower edge (0)
      if(blast<0 || blast>=nLast) G4cout<<"G4QaBarElCS::CCS:b="<<blast<<","<<nLast<<G4endl;
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
#ifdef pdebug
      G4cout<<"G4QAntiBaryonElasticCS::CalculateCS:(E) S1="<<theS1<<", B1="<<theB1<<G4endl;
#endif
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
#ifdef pdebug
      G4cout<<"G4QAntiBaryonElasticCrossSection::CalcCrossSection:Sig="<<lastSIG<<",P="
            <<pMom<<",Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<",onlyCS="<<onlyCS<<G4endl;
#endif
      if(!onlyCS)                       // Skip the differential cross-section parameters
      {
        G4double SSTL=lastSST[blast];           // the low bin of the first squared slope
        theSS=SSTL+shift*(lastSST[lastL]-SSTL); // the basic value of the first sq.slope
        G4double S1TL=lastS1T[blast];           // the low bin of the first mantissa
        theS1=S1TL+shift*(lastS1T[lastL]-S1TL); // the basic value of the first mantissa
        G4double B1TL=lastB1T[blast];           // the low bin of the first slope
#ifdef pdebug
        G4cout<<"G4QAntiBaryonElasticCS::CalcCrossSection:b="<<blast<<",ls="<<lastL<<",SL="
              <<S1TL<<",SU="<<lastS1T[lastL]<<",BL="<<B1TL<<",BU="<<lastB1T[lastL]<<G4endl;
#endif
        theB1=B1TL+shift*(lastB1T[lastL]-B1TL); // the basic value of the first slope
        G4double S2TL=lastS2T[blast];           // the low bin of the second mantissa
        theS2=S2TL+shift*(lastS2T[lastL]-S2TL); // the basic value of the second mantissa
        G4double B2TL=lastB2T[blast];           // the low bin of the second slope
        theB2=B2TL+shift*(lastB2T[lastL]-B2TL); // the basic value of the second slope
        G4double S3TL=lastS3T[blast];           // the low bin of the third mantissa
        theS3=S3TL+shift*(lastS3T[lastL]-S3TL); // the basic value of the third mantissa
#ifdef pdebug
        G4cout<<"G4QAntiBaryonElasticCrossSection::CCS:s3l="<<S3TL<<",sh3="<<shift<<",s3h="
              <<lastS3T[lastL]<<",b="<<blast<<",l="<<lastL<<G4endl;
#endif
        G4double B3TL=lastB3T[blast];           // the low bin of the third slope
        theB3=B3TL+shift*(lastB3T[lastL]-B3TL); // the basic value of the third slope
        G4double S4TL=lastS4T[blast];           // the low bin of the 4-th mantissa
        theS4=S4TL+shift*(lastS4T[lastL]-S4TL); // the basic value of the 4-th mantissa
#ifdef pdebug
        G4cout<<"G4QAntiBaryonElasticCrossSection::CCS:s4l="<<S4TL<<",sh4="<<shift<<",s4h="
              <<lastS4T[lastL]<<",b="<<blast<<",l="<<lastL<<G4endl;
#endif
        G4double B4TL=lastB4T[blast];           // the low bin of the 4-th slope
        theB4=B4TL+shift*(lastB4T[lastL]-B4TL); // the basic value of the 4-th slope
      }
#ifdef pdebug
      G4cout<<"G4QAntiBaryonElasticCS::CalculateCS:(I) S1="<<theS1<<", B1="<<theB1<<G4endl;
#endif
    }
  }
  else lastSIG=GetTabValues(lastLP, PDG, tgZ, tgN); // Direct calculation beyond the table
  if(lastSIG<0.) lastSIG = 0.;                   // @@ a Warning print can be added
#ifdef pdebug
  G4cout<<"G4QAntiBaryonElasticCrossSection::CalculateCS: END, onlyCS="<<onlyCS<<G4endl;
#endif
  return lastSIG;
}

// It has parameter sets for all tZ/tN/PDG, using them the tables can be created/updated
G4double G4QAntiBaryonElasticCrossSection::GetPTables(G4double LP, G4double ILP, G4int PDG,
                                                  G4int tgZ, G4int tgN)
{
  // @@ At present all nA==pA ---------> Each neucleus can have not more than 51 parameters
  static const G4double pwd=2727;
  const G4int n_appel=30;                // #of parameters for app-elastic (<nPoints=128)
  //                          -0- -1- -2- -3-  -4- -5- -6- -7- -8--9--10--11--12--13--14-
  G4double app_el[n_appel]={1.25,3.5,80.,1.,.0557,6.72,5.,74.,3.,3.4,.2,.17,.001,8.,.055,
                            3.64,5.e-5,4000.,1500.,.46,1.2e6,3.5e6,5.e-5,1.e10,8.5e8,
                            1.e10,1.1,3.4e6,6.8e6,0.};
  //                        -15-  -16-  -17-  -18- -19- -20-  -21-  -22-  -23- -24-
  //                        -25-  -26- -27- -28- -29- 
  if(PDG>-3334 && PDG<-1111)
  {
    // -- Total pp elastic cross section cs & s1/b1 (main), s2/b2 (tail1), s3/b3 (tail2) --
    //p2=p*p;p3=p2*p;sp=sqrt(p);p2s=p2*sp;lp=log(p);dl1=lp-(3.=par(3));p4=p2*p2; p=|3-mom|
    //CS=2.865/p2s/(1+.0022/p2s)+(18.9+.6461*dl1*dl1+9./p)/(1.+.425*lp)/(1.+.4276/p4);
    //   par(0)       par(7)     par(1) par(2)      par(4)      par(5)         par(6)
    //dl2=lp-5., s1=(74.+3.*dl2*dl2)/(1+3.4/p4/p)+(.2/p2+17.*p)/(p4+.001*sp),
    //     par(8) par(9) par(10)        par(11)   par(12)par(13)    par(14)
    // b1=8.*p**.055/(1.+3.64/p3); s2=5.e-5+4000./(p4+1500.*p); b2=.46+1.2e6/(p4+3.5e6/sp);
    // par(15) par(16)  par(17)     par(18) par(19)  par(20)   par(21) par(22)  par(23)
    // s3=5.e-5+1.e10/(p4*p4+8.5e8*p2+1.e10); b3=1.1+3.4e6/(p4+6.8e6); ss=0.
    //  par(24) par(25)     par(26)  par(27) par(28) par(29)  par(30)   par(31)
    //
    if(lastPAR[nLast]!=pwd) // A unique flag to avoid the repeatable definition
    {
      if ( tgZ == 1 && tgN == 0 )
      {
        for (G4int ip=0; ip<n_appel; ip++) lastPAR[ip]=app_el[ip]; // PiMinus+P
      }
      else
      {
        G4double a=tgZ+tgN;
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
        //G4double a20=a16*a4;
        G4double a32=a16*a16;
        // Reaction cross-section parameters (pel=peh_fit.f)
        lastPAR[0]=.23*asa/(1.+a*.15);                                       // p1
        lastPAR[1]=2.8*asa/(1.+a*(.015+.05/ssa));                            // p2
        lastPAR[2]=15.*a/(1.+.005*a2);                                       // p3
        lastPAR[3]=.013*a2/(1.+a3*(.006+a*.00001));                          // p4
        lastPAR[4]=5.;                                                       // p5
        lastPAR[5]=0.;                                                       // p6 not used
        lastPAR[6]=0.;                                                       // p7 not used
        lastPAR[7]=0.;                                                       // p8 not used
        lastPAR[8]=0.;                                                       // p9 not used
        // @@ the differential cross-section is parameterized separately for A>6 & A<7
        if(a<6.5)
        {
          G4double a28=a16*a12;
          // The main pre-exponent      (pel_sg)
          lastPAR[ 9]=4000*a;                                // p1
          lastPAR[10]=1.2e7*a8+380*a17;                      // p2
          lastPAR[11]=.7/(1.+4.e-12*a16);                    // p3
          lastPAR[12]=2.5/a8/(a4+1.e-16*a32);                // p4
          lastPAR[13]=.28*a;                                 // p5
          lastPAR[14]=1.2*a2+2.3;                            // p6
          lastPAR[15]=3.8/a;                                 // p7
          // The main slope             (pel_sl)
          lastPAR[16]=.01/(1.+.0024*a5);                     // p1
          lastPAR[17]=.2*a;                                  // p2
          lastPAR[18]=9.e-7/(1.+.035*a5);                    // p3
          lastPAR[19]=(42.+2.7e-11*a16)/(1.+.14*a);          // p4
          // The main quadratic         (pel_sh)
          lastPAR[20]=2.25*a3;                               // p1
          lastPAR[21]=18.;                                   // p2
          lastPAR[22]=2.4e-3*a8/(1.+2.6e-4*a7);              // p3
          lastPAR[23]=3.5e-36*a32*a8/(1.+5.e-15*a32/a);      // p4
          // The 1st max pre-exponent   (pel_qq)
          lastPAR[24]=1.e5/(a8+2.5e12/a16);                  // p1
          lastPAR[25]=8.e7/(a12+1.e-27*a28*a28);             // p2 
          lastPAR[26]=.0006*a3;                              // p3
          // The 1st max slope          (pel_qs)
          lastPAR[27]=10.+4.e-8*a12*a;                       // p1
          lastPAR[28]=.114;                                  // p2
          lastPAR[29]=.003;                                  // p3
          lastPAR[30]=2.e-23;                                // p4
          // The effective pre-exponent (pel_ss)
          lastPAR[31]=1./(1.+.0001*a8);                      // p1
          lastPAR[32]=1.5e-4/(1.+5.e-6*a12);                 // p2
          lastPAR[33]=.03;                                   // p3
          // The effective slope        (pel_sb)
          lastPAR[34]=a/2;                                   // p1
          lastPAR[35]=2.e-7*a4;                              // p2
          lastPAR[36]=4.;                                    // p3
          lastPAR[37]=64./a3;                                // p4
          // The gloria pre-exponent    (pel_us)
          lastPAR[38]=1.e8*std::exp(.32*asa);                // p1
          lastPAR[39]=20.*std::exp(.45*asa);                 // p2
          lastPAR[40]=7.e3+2.4e6/a5;                         // p3
          lastPAR[41]=2.5e5*std::exp(.085*a3);               // p4
          lastPAR[42]=2.5*a;                                 // p5
          // The gloria slope           (pel_ub)
          lastPAR[43]=920.+.03*a8*a3;                        // p1
          lastPAR[44]=93.+.0023*a12;                         // p2
#ifdef pdebug
         G4cout<<"G4QaBElCS::CalcCS:la "<<lastPAR[38]<<", "<<lastPAR[39]<<", "<<lastPAR[40]
               <<", "<<lastPAR[42]<<", "<<lastPAR[43]<<", "<<lastPAR[44]<<G4endl;
#endif
        }
        else
        {
          G4double p1a10=2.2e-28*a10;
          G4double r4a16=6.e14/a16;
          G4double s4a16=r4a16*r4a16;
          // a24
          // a36
          // The main pre-exponent      (peh_sg)
          lastPAR[ 9]=4.5*std::pow(a,1.15);                  // p1
          lastPAR[10]=.06*std::pow(a,.6);                    // p2
          lastPAR[11]=.6*a/(1.+2.e15/a16);                   // p3
          lastPAR[12]=.17/(a+9.e5/a3+1.5e33/a32);            // p4
          lastPAR[13]=(.001+7.e-11*a5)/(1.+4.4e-11*a5);      // p5
          lastPAR[14]=(p1a10*p1a10+2.e-29)/(1.+2.e-22*a12);  // p6
          // The main slope             (peh_sl)
          lastPAR[15]=400./a12+2.e-22*a9;                    // p1
          lastPAR[16]=1.e-32*a12/(1.+5.e22/a14);             // p2
          lastPAR[17]=1000./a2+9.5*sa*ssa;                   // p3
          lastPAR[18]=4.e-6*a*asa+1.e11/a16;                 // p4
          lastPAR[19]=(120./a+.002*a2)/(1.+2.e14/a16);       // p5
          lastPAR[20]=9.+100./a;                             // p6
          // The main quadratic         (peh_sh)
          lastPAR[21]=.002*a3+3.e7/a6;                       // p1
          lastPAR[22]=7.e-15*a4*asa;                         // p2
          lastPAR[23]=9000./a4;                              // p3
          // The 1st max pre-exponent   (peh_qq)
          lastPAR[24]=.0011*asa/(1.+3.e34/a32/a4);           // p1
          lastPAR[25]=1.e-5*a2+2.e14/a16;                    // p2
          lastPAR[26]=1.2e-11*a2/(1.+1.5e19/a12);            // p3
          lastPAR[27]=.016*asa/(1.+5.e16/a16);               // p4
          // The 1st max slope          (peh_qs)
          lastPAR[28]=.002*a4/(1.+7.e7/std::pow(a-6.83,14)); // p1
          lastPAR[29]=2.e6/a6+7.2/std::pow(a,.11);           // p2
          lastPAR[30]=11.*a3/(1.+7.e23/a16/a8);              // p3
          lastPAR[31]=100./asa;                              // p4
          // The 2nd max pre-exponent   (peh_ss)
          lastPAR[32]=(.1+4.4e-5*a2)/(1.+5.e5/a4);           // p1
          lastPAR[33]=3.5e-4*a2/(1.+1.e8/a8);                // p2
          lastPAR[34]=1.3+3.e5/a4;                           // p3
          lastPAR[35]=500./(a2+50.)+3;                       // p4
          lastPAR[36]=1.e-9/a+s4a16*s4a16;                   // p5
          // The 2nd max slope          (peh_sb)
          lastPAR[37]=.4*asa+3.e-9*a6;                       // p1
          lastPAR[38]=.0005*a5;                              // p2
          lastPAR[39]=.002*a5;                               // p3
          lastPAR[40]=10.;                                   // p4
          // The effective pre-exponent (peh_us)
          lastPAR[41]=.05+.005*a;                            // p1
          lastPAR[42]=7.e-8/sa;                              // p2
          lastPAR[43]=.8*sa;                                 // p3
          lastPAR[44]=.02*sa;                                // p4
          lastPAR[45]=1.e8/a3;                               // p5
          lastPAR[46]=3.e32/(a32+1.e32);                     // p6
          // The effective slope        (peh_ub)
          lastPAR[47]=24.;                                   // p1
          lastPAR[48]=20./sa;                                // p2
          lastPAR[49]=7.e3*a/(sa+1.);                        // p3
          lastPAR[50]=900.*sa/(1.+500./a3);                  // p4
#ifdef pdebug
         G4cout<<"G4QaBElCS::CalcCS:ha "<<lastPAR[41]<<", "<<lastPAR[42]<<", "<<lastPAR[43]
               <<", "<<lastPAR[44]<<", "<<lastPAR[45]<<", "<<lastPAR[46]<<G4endl;
#endif
        }
        // Parameter for lowEnergyNeutrons
        lastPAR[51]=1.e15+2.e27/a4/(1.+2.e-18*a16);
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
#ifdef pdebug
      G4cout<<"G4QAntiBaryonElasticCrossSection::GetPTables:ip=0(init), lp="<<lp<<",S1="
            <<theS1<<",B1="<<theB1<<",S2="<<theS2<<",B2="<<theB3<<",S3="<<theS3
            <<",B3="<<theB3<<",S4="<<theS4<<",B4="<<theB4<<G4endl;
#endif
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
#ifdef pdebug
            G4cout<<"G4QAntiBaryonElasticCrossSection::GetPTables:ip="<<ip<<",lp="<<lp
                  <<",S1="<<theS1<<",B1="<<theB1<<",S2="<<theS2<<",B2="<<theB2<<",S3="
                  <<theS3<<",B3="<<theB3<<",S4="<<theS4<<",B4="<<theB4<<G4endl;
#endif
          }
          return lp;
        }
        else G4cout<<"*Warning*G4QAntiBaryonElasticCrossSection::GetPTables: PDG="<<PDG
                   <<", Z="<<tgZ<<", N="<<tgN<<", i="<<ini<<" > fin="<<fin<<", LP="<<LP
                   <<" > ILP="<<ILP<<" nothing is done!"<<G4endl;
      }
      else G4cout<<"*Warning*G4QAntiBaryonElasticCrossSection::GetPTables: PDG="<<PDG
                 <<", Z="<<tgZ<<", N="<<tgN<<", i="<<ini<<">= max="<<nPoints<<", LP="<<LP
                 <<" > ILP="<<ILP<<", lPMax="<<lPMax<<" nothing is done!"<<G4endl;
    }
#ifdef pdebug
    else G4cout<<"*Warning*G4QAntiBaryonElasticCrossSection::GetPTab:PDG="<<PDG<<",Z="<<tgZ
               <<", N="<<tgN<<", LP="<<LP<<" <= ILP="<<ILP<<" nothing is done!"<<G4endl;
#endif
  }
  else
  {
    G4cout<<"*Error*G4QAntiBaryonElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ
          <<", N="<<tgN<<", while it is defined only for Anti Baryons"<<G4endl;
    throw G4QException("G4QAntiBaryonElasticCrossSection::GetPTables:onlyaBA implemented");
  }
  return ILP;
}

// Returns Q2=-t in independent units (MeV^2) (all internal calculations are in GeV)
G4double G4QAntiBaryonElasticCrossSection::GetExchangeT(G4int tgZ, G4int tgN, G4int PDG)
{
  static const G4double GeVSQ=gigaelectronvolt*gigaelectronvolt;
  static const G4double third=1./3.;
  static const G4double fifth=1./5.;
  static const G4double sevth=1./7.;
#ifdef tdebug
  G4cout<<"G4QAnBaElCS::GetExcT: F="<<onlyCS<<",Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<G4endl;
#endif
  if(PDG<-3334 || PDG>-1111)G4cout<<"*Warning*G4QAntiBaryonElCS::GetExT:PDG="<<PDG<<G4endl;
  if(onlyCS)G4cout<<"WarningG4QAntiBaryonElasticCrossSection::GetExchanT:onlyCS=1"<<G4endl;
  if(lastLP<-4.3) return lastTM*GeVSQ*G4UniformRand();// S-wave for p<14 MeV/c (kinE<.1MeV)
  G4double q2=0.;
  if(tgZ==1 && tgN==0)                // ===> p+p=p+p
  {
#ifdef tdebug
    G4cout<<"G4QAnBarElCS::GetExchangeT: TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1<<",S2="
          <<theS2<<",B2="<<theB2<<",S3="<<theS3<<",B3="<<theB3<<",GeV2="<<GeVSQ<<G4endl;
#endif
    G4double E1=lastTM*theB1;
    G4double R1=(1.-std::exp(-E1));
#ifdef tdebug
    G4double ts1=-std::log(1.-R1)/theB1;
    G4double ds1=std::fabs(ts1-lastTM)/lastTM;
    if(ds1>.0001)
      G4cout<<"*Warning*G4QAntiBaryonElasticCrossSection::GetExT:1p "<<ts1<<"#"<<lastTM
            <<",d="<<ds1<<",R1="<<R1<<",E1="<<E1<<G4endl;
#endif
    G4double E2=lastTM*theB2;
    G4double R2=(1.-std::exp(-E2*E2*E2));
#ifdef tdebug
    G4double ts2=std::pow(-std::log(1.-R2),.333333333)/theB2;
    G4double ds2=std::fabs(ts2-lastTM)/lastTM;
    if(ds2>.0001)
      G4cout<<"*Warning*G4QAntiBaryonElasticCrossSection::GetExT:2p "<<ts2<<"#"<<lastTM
            <<",d="<<ds2<<",R2="<<R2<<",E2="<<E2<<G4endl;
#endif
    G4double E3=lastTM*theB3;
    G4double R3=(1.-std::exp(-E3));
#ifdef tdebug
    G4double ts3=-std::log(1.-R3)/theB3;
    G4double ds3=std::fabs(ts3-lastTM)/lastTM;
    if(ds3>.0001)
      G4cout<<"*Warning*G4QAntiBaryonElasticCrossSection::GetExT:3p "<<ts3<<"#"<<lastTM
            <<",d="<<ds3<<",R3="<<R1<<",E3="<<E3<<G4endl;
#endif
    G4double I1=R1*theS1/theB1;
    G4double I2=R2*theS2;
    G4double I3=R3*theS3;
    G4double I12=I1+I2;
    G4double rand=(I12+I3)*G4UniformRand();
    if     (rand<I1 )
    {
      G4double ran=R1*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB1;
    }
    else if(rand<I12)
    {
      G4double ran=R2*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran);
      if(q2<0.) q2=0.;
      q2=std::pow(q2,third)/theB2;
    }
    else
    {
      G4double ran=R3*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB3;
    }
  }
  else
  {
    G4double a=tgZ+tgN;
#ifdef tdebug
    G4cout<<"G4QAntiBaryonElasticCrossSection::GetExT:a="<<a<<",t="<<lastTM<<",S1="<<theS1
          <<",B1="<<theB1<<",SS="<<theSS<<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS3
          <<",B3="<<theB3<<",S4="<<theS4<<",B4="<<theB4<<G4endl;
#endif
    G4double E1=lastTM*(theB1+lastTM*theSS);
    G4double R1=(1.-std::exp(-E1));
    G4double tss=theSS+theSS; // for future solution of quadratic equation (imediate check)
#ifdef tdebug
    G4double ts1=-std::log(1.-R1)/theB1;
    if(std::fabs(tss)>1.e-7) ts1=(std::sqrt(theB1*(theB1+(tss+tss)*ts1))-theB1)/tss;
    G4double ds1=(ts1-lastTM)/lastTM;
    if(ds1>.0001)
      G4cout<<"*Warning*G4QAntiBaryonElasticCrossSection::GetExT:1a "<<ts1<<"#"<<lastTM
            <<",d="<<ds1<<",R1="<<R1<<",E1="<<E1<<G4endl;
#endif
    G4double tm2=lastTM*lastTM;
    G4double E2=lastTM*tm2*theB2;                   // power 3 for lowA, 5 for HighA (1st)
    if(a>6.5)E2*=tm2;                               // for heavy nuclei
    G4double R2=(1.-std::exp(-E2));
#ifdef tdebug
    G4double ts2=-std::log(1.-R2)/theB2;
    if(a<6.5)ts2=std::pow(ts2,third);
    else     ts2=std::pow(ts2,fifth);
    G4double ds2=std::fabs(ts2-lastTM)/lastTM;
    if(ds2>.0001)
      G4cout<<"*Warning*G4QAntiBaryonElasticCrossSection::GetExT:2a "<<ts2<<"#"<<lastTM
            <<",d="<<ds2<<",R2="<<R2<<",E2="<<E2<<G4endl;
#endif
    G4double E3=lastTM*theB3;
    if(a>6.5)E3*=tm2*tm2*tm2;                       // power 1 for lowA, 7 (2nd) for HighA
    G4double R3=(1.-std::exp(-E3));
#ifdef tdebug
    G4double ts3=-std::log(1.-R3)/theB3;
    if(a>6.5)ts3=std::pow(ts3,sevth);
    G4double ds3=std::fabs(ts3-lastTM)/lastTM;
    if(ds3>.0001)
      G4cout<<"*Warning*G4QAntiBaryonElasticCrossSection::GetExT:3a "<<ts3<<"#"<<lastTM
            <<",d="<<ds3<<",R3="<<R3<<",E3="<<E3<<G4endl;
#endif
    G4double E4=lastTM*theB4;
    G4double R4=(1.-std::exp(-E4));
#ifdef tdebug
    G4double ts4=-std::log(1.-R4)/theB4;
    G4double ds4=std::fabs(ts4-lastTM)/lastTM;
    if(ds4>.0001)
      G4cout<<"*Warning*G4QAntiBaryonElasticCrossSection::GetExT:4a "<<ts4<<"#"<<lastTM
            <<",d="<<ds4<<",R4="<<R4<<",E4="<<E4<<G4endl;
#endif
    G4double I1=R1*theS1;
    G4double I2=R2*theS2;
    G4double I3=R3*theS3;
    G4double I4=R4*theS4;
    G4double I12=I1+I2;
    G4double I13=I12+I3;
    G4double rand=(I13+I4)*G4UniformRand();
#ifdef tdebug
    G4cout<<"G4QaBElCS::GExT:1="<<I1<<",2="<<I2<<",3="<<I3<<",4="<<I4<<",r="<<rand<<G4endl;
#endif
    if(rand<I1)
    {
      G4double ran=R1*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB1;
      if(std::fabs(tss)>1.e-7) q2=(std::sqrt(theB1*(theB1+(tss+tss)*q2))-theB1)/tss;
#ifdef tdebug
      G4cout<<"G4QaBElCS::GExT:Q2="<<q2<<",ss="<<tss/2<<",b1="<<theB1<<",t1="<<ts1<<G4endl;
#endif
    }
    else if(rand<I12)
    {
      G4double ran=R2*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB2;
      if(q2<0.) q2=0.;
      if(a<6.5) q2=std::pow(q2,third);
      else      q2=std::pow(q2,fifth);
#ifdef tdebug
      G4cout<<"G4QaBaElCS::GetExT:Q2="<<q2<<",r2="<<R2<<",b2="<<theB2<<",t2="<<ts2<<G4endl;
#endif
    }
    else if(rand<I13)
    {
      G4double ran=R3*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB3;
      if(q2<0.) q2=0.;
      if(a>6.5) q2=std::pow(q2,sevth);
#ifdef tdebug
      G4cout<<"G4QaBaElCS::GetExT:Q2="<<q2<<",r3="<<R2<<",b3="<<theB2<<",t3="<<ts2<<G4endl;
#endif
    }
    else
    {
      G4double ran=R4*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB4;
      if(a<6.5) q2=lastTM-q2;                    // u reduced for lightA (starts from 0)
#ifdef tdebug
      G4cout<<"G4QaBElCS::GExT:Q2="<<q2<<",m="<<lastTM<<",b4="<<theB3<<",t4="<<ts3<<G4endl;
#endif
    }
  }
  if(q2<0.) q2=0.;
  if(!(q2>=-1.||q2<=1.))G4cout<<"*NAN*G4QaBElasticCrossSect::GetExchangeT:-t="<<q2<<G4endl;
  if(q2>lastTM)
  {
#ifdef tdebug
    G4cout<<"*Warning*G4QAntiBaryonElasticCrossSection::GExT:-t="<<q2<<">"<<lastTM<<G4endl;
#endif
    q2=lastTM;
  }
  return q2*GeVSQ;
}

// Returns B in independent units (MeV^-2) (all internal calculations are in GeV) see ExT
G4double G4QAntiBaryonElasticCrossSection::GetSlope(G4int tgZ, G4int tgN, G4int PDG)
{
  static const G4double GeVSQ=gigaelectronvolt*gigaelectronvolt;
#ifdef tdebug
  G4cout<<"G4QAnBarElCS::GetSlope:"<<onlyCS<<", Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<G4endl;
#endif
  if(onlyCS)G4cout<<"WarningG4QAntiBaryonElasticCrossSection::GetSlope:onlCS=true"<<G4endl;
  if(lastLP<-4.3) return 0.;          // S-wave for p<14 MeV/c (kinE<.1MeV)
  if(PDG<-3334 || PDG>-1111)
  {
    G4cout<<"*Error*G4QAntiBaryonElasticCrossSection::GetSlope: PDG="<<PDG<<", Z="<<tgZ
          <<", N="<<tgN<<", while it is defined only for Anti Baryons"<<G4endl;
    throw G4QException("G4QAntiBaryonElasticCrossSection::GetSlope: AnBa are implemented");
  }
  if(theB1<0.) theB1=0.;
  if(!(theB1>=-1.||theB1<=1.))G4cout<<"*NAN*G4QaBaElasticCrossS::Getslope:"<<theB1<<G4endl;
  return theB1/GeVSQ;
}

// Returns half max(Q2=-t) in independent units (MeV^2)
G4double G4QAntiBaryonElasticCrossSection::GetHMaxT()
{
  static const G4double HGeVSQ=gigaelectronvolt*gigaelectronvolt/2.;
  return lastTM*HGeVSQ;
}

// lastLP is used, so calculating tables, one need to remember and then recover lastLP
G4double G4QAntiBaryonElasticCrossSection::GetTabValues(G4double lp, G4int PDG, G4int tgZ,
                                                    G4int tgN)
{
  if(PDG<-3334 || PDG>-1111) G4cout<<"*Warning*G4QAntiBaryElCS::GetTabV:PDG="<<PDG<<G4endl;
  if(tgZ<0 || tgZ>92)
  {
    G4cout<<"*Warning*G4QAntiBaryonElCS::GetTabValue:(1-92) NoIsotopesFor Z="<<tgZ<<G4endl;
    return 0.;
  }
  G4int iZ=tgZ-1; // Z index
  if(iZ<0)
  {
    iZ=0;         // conversion of the neutron target to the proton target
    tgZ=1;
    tgN=0;
  }
  //if(nN[iZ][0] < 0)
  //{
#ifdef isodebug
  //  G4cout<<"*Warning*G4QAntiBaryonElastCS::GetTabValue: NoIsotopes for Z="<<tgZ<<G4endl;
#endif
  //  return 0.;
  //}
#ifdef pdebug
  G4cout<<"G4QaBElasticCS::GetTabVal:l="<<lp<<",Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<G4endl;
#endif
  G4double p=std::exp(lp);              // momentum
  G4double sp=std::sqrt(p);             // sqrt(p)
  G4double p2=p*p;            
  G4double p3=p2*p;
  G4double p4=p3*p;
  if ( tgZ == 1 && tgN == 0 ) // PiMin+P
  {
    G4double dl2=lp-lastPAR[6]; // ld ?
    theSS=lastPAR[29];
    theS1=(lastPAR[7]+lastPAR[8]*dl2*dl2)/(1.+lastPAR[9]/p4/p)+
          (lastPAR[10]/p2+lastPAR[11]*p)/(p4+lastPAR[12]*sp);
    theB1=lastPAR[13]*std::pow(p,lastPAR[14])/(1.+lastPAR[15]/p3);
    theS2=lastPAR[16]+lastPAR[17]/(p4+lastPAR[18]*p);
    theB2=lastPAR[19]+lastPAR[20]/(p4+lastPAR[21]/sp); 
    theS3=lastPAR[22]+lastPAR[23]/(p4*p4+lastPAR[24]*p2+lastPAR[25]);
    theB3=lastPAR[26]+lastPAR[27]/(p4+lastPAR[28]); 
    theS4=0.;
    theB4=0.; 
#ifdef tdebug
    G4cout<<"G4QAntiBarElasticCS::GetTableValues:TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1
          <<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS1<<",B3="<<theB1<<G4endl;
#endif
    // Returns the total elastic pim-p cross-section (to avoid spoiling lastSIG)
    G4double ye=std::exp(lp*lastPAR[0]);
    G4double dp=lp-lastPAR[1];
    return lastPAR[2]/(ye+lastPAR[3])+lastPAR[4]*dp*dp+lastPAR[5];
  }
  else
  {
    G4double p5=p4*p;
    G4double p6=p5*p;
    G4double p8=p6*p2;
    G4double p10=p8*p2;
    G4double p12=p10*p2;
    G4double p16=p8*p8;
    //G4double p24=p16*p8;
    G4double dl=lp-5.;
    G4double a=tgZ+tgN;
    G4double pah=std::pow(p,a/2);
    G4double pa=pah*pah;
    G4double pa2=pa*pa;
    if(a<6.5)
    {
      theS1=lastPAR[9]/(1.+lastPAR[10]*p4*pa)+lastPAR[11]/(p4+lastPAR[12]*p4/pa2)+
            (lastPAR[13]*dl*dl+lastPAR[14])/(1.+lastPAR[15]/p2);
      theB1=(lastPAR[16]+lastPAR[17]*p2)/(p4+lastPAR[18]/pah)+lastPAR[19];
      theSS=lastPAR[20]/(1.+lastPAR[21]/p2)+lastPAR[22]/(p6/pa+lastPAR[23]/p16);
      theS2=lastPAR[24]/(pa/p2+lastPAR[25]/p4)+lastPAR[26];
      theB2=lastPAR[27]*std::pow(p,lastPAR[28])+lastPAR[29]/(p8+lastPAR[30]/p16);
      theS3=lastPAR[31]/(pa*p+lastPAR[32]/pa)+lastPAR[33];
      theB3=lastPAR[34]/(p3+lastPAR[35]/p6)+lastPAR[36]/(1.+lastPAR[37]/p2);
      theS4=p2*(pah*lastPAR[38]*std::exp(-pah*lastPAR[39])+
                lastPAR[40]/(1.+lastPAR[41]*std::pow(p,lastPAR[42])));
      theB4=lastPAR[43]*pa/p2/(1.+pa*lastPAR[44]);
#ifdef tdebug
      G4cout<<"G4QaBaElCS::GetTabV: lA, p="<<p<<",S1="<<theS1<<",B1="<<theB1<<",SS="<<theSS
            <<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS3<<",B3="<<theB3<<",S4="<<theS4
            <<",B4="<<theB4<<G4endl;
#endif
    }
    else
    {
      theS1=lastPAR[9]/(1.+lastPAR[10]/p4)+lastPAR[11]/(p4+lastPAR[12]/p2)+
            lastPAR[13]/(p5+lastPAR[14]/p16);
      theB1=(lastPAR[15]/p8+lastPAR[19])/(p+lastPAR[16]/std::pow(p,lastPAR[20]))+
            lastPAR[17]/(1.+lastPAR[18]/p4);
      theSS=lastPAR[21]/(p4/std::pow(p,lastPAR[23])+lastPAR[22]/p4);
      theS2=lastPAR[24]/p4/(std::pow(p,lastPAR[25])+lastPAR[26]/p12)+lastPAR[27];
      theB2=lastPAR[28]/std::pow(p,lastPAR[29])+lastPAR[30]/std::pow(p,lastPAR[31]);
      theS3=lastPAR[32]/std::pow(p,lastPAR[35])/(1.+lastPAR[36]/p12)+
            lastPAR[33]/(1.+lastPAR[34]/p6);
      theB3=lastPAR[37]/p8+lastPAR[38]/p2+lastPAR[39]/(1.+lastPAR[40]/p8);
      theS4=(lastPAR[41]/p4+lastPAR[46]/p)/(1.+lastPAR[42]/p10)+
            (lastPAR[43]+lastPAR[44]*dl*dl)/(1.+lastPAR[45]/p12);
      theB4=lastPAR[47]/(1.+lastPAR[48]/p)+lastPAR[49]*p4/(1.+lastPAR[50]*p5);
#ifdef tdebug
      G4cout<<"G4QaBaElCS::GetTabV: hA, p="<<p<<",S1="<<theS1<<",B1="<<theB1<<",SS="<<theSS
            <<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS3<<",B3="<<theB3<<",S4="<<theS4
            <<",B4="<<theB4<<G4endl;
#endif
    }
    // Returns the total elastic (n/p)A cross-section (to avoid spoiling lastSIG)
#ifdef tdebug
    G4cout<<"G4QAntiBarElCS::GetTabV: PDG="<<PDG<<",P="<<p<<",N="<<tgN<<",Z="<<tgZ<<G4endl;
#endif
    G4double dlp=lp-lastPAR[5]; // ax
    //         p1               p2          p3                 p4
    return (lastPAR[0]*dlp*dlp+lastPAR[1]+lastPAR[2]/p)/(1.+lastPAR[3]/p);
  }
  return 0.;
} // End of GetTableValues

// Returns max -t=Q2 (GeV^2) for the momentum pP(GeV) and the target nucleus (tgN,tgZ)
G4double G4QAntiBaryonElasticCrossSection::GetQ2max(G4int PDG, G4int tgZ, G4int tgN,
                                                    G4double pP)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass()*.001; // MeV to GeV
  static const G4double mProt= G4QPDGCode(2212).GetMass()*.001; // MeV to GeV
  //static const G4double mLamb= G4QPDGCode(3122).GetMass()*.001; // MeV to GeV
  //static const G4double mHe3 = G4QPDGCode(2112).GetNuclMass(2,1,0)*.001; // MeV to GeV
  //static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0)*.001; // MeV to GeV
  //static const G4double mDeut = G4QPDGCode(2112).GetNuclMass(1,1,0)*.001; // MeV to GeV
  static const G4double mNuc2= sqr((mProt+mNeut)/2);
  //static const G4double mProt2= mProt*mProt;
  //static const G4double mNeut2= mNeut*mNeut;
  //static const G4double mDeut2= mDeut*mDeut;
  G4double pP2=pP*pP;                                 // squared momentum of the projectile
  if(tgZ || tgN>-1)                                   // ---> pipA
  {
    G4double mt=G4QPDGCode(90000000+tgZ*1000+tgN).GetMass()*.001; // Target mass in GeV
    G4double dmt=mt+mt;
    G4double s=dmt*std::sqrt(pP2+mNuc2)+mNuc2+mt*mt;    // Mondelstam s (@@ other AntiBar?)
    return dmt*dmt*pP2/s;
  }
  else
  {
    G4cout<<"*Error*G4QAntiBaryonElasticCrossSection::GetQ2ma:PDG="<<PDG<<",Z="<<tgZ<<",N="
          <<tgN<<", while it is defined only for p projectiles & Z_target>0"<<G4endl;
    throw G4QException("G4QAntiBaryonElasticCrossSection::GetQ2max: only aBA implemented");
  }
}
