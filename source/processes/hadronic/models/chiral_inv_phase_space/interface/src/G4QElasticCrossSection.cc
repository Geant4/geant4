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
// $Id: G4QElasticCrossSection.cc,v 1.20 2007/01/16 14:41:40 mkossov Exp $
// GEANT4 tag $Name: geant4-08-02-patch-01 $
//
//
// G4 Physics class: G4QElasticCrossSection for N+A elastic cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Oct-06
// 
//================================================================================

//#define debug
//#define isodebug
//#define pdebug
//#define ppdebug
//#define tdebug
//#define sdebug

#include "G4QElasticCrossSection.hh"

// Initialization of the static parameters
const G4int G4QElasticCrossSection::nPoints=128;//#of points in the AMDB tables(>anyPar)(D)
const G4int G4QElasticCrossSection::nLast=nPoints-1; // the Last element in the table   (D)
G4double  G4QElasticCrossSection::lPMin=-8.;  // Min tabulated logarithmic Momentum     (D)
G4double  G4QElasticCrossSection::lPMax= 8.;  // Max tabulated logarithmic Momentum     (D)
G4double  G4QElasticCrossSection::dlnP=(lPMax-lPMin)/nLast;// Log step in the table     (D)
G4bool    G4QElasticCrossSection::onlyCS=true;// Flag to calculate only CS (not Si/Bi)		(L)
G4double  G4QElasticCrossSection::lastSIG=0.; // Last calculated cross section								  (L)
G4double  G4QElasticCrossSection::lastLP=-10.;// Last log(mom_of_the_incident_hadron)   (L)
G4double  G4QElasticCrossSection::lastTM=0.;  // Last t_maximum                         (L)
G4double  G4QElasticCrossSection::theSS=0.;   // The Last sq.slope of first difruction  (L)
G4double  G4QElasticCrossSection::theS1=0.;   // The Last mantissa of first difruction  (L)
G4double  G4QElasticCrossSection::theB1=0.;   // The Last slope of first difruction     (L)
G4double  G4QElasticCrossSection::theS2=0.;   // The Last mantissa of second difruction (L)
G4double  G4QElasticCrossSection::theB2=0.;   // The Last slope of second difruction    (L)
G4double  G4QElasticCrossSection::theS3=0.;   // The Last mantissa of third difruction  (L)
G4double  G4QElasticCrossSection::theB3=0.;   // The Last slope of third difruction     (L)
G4double  G4QElasticCrossSection::theS4=0.;   // The Last mantissa of 4-th difruction   (L)
G4double  G4QElasticCrossSection::theB4=0.;   // The Last slope of 4-th difruction      (L)
G4int     G4QElasticCrossSection::lastTZ=0;   // Last atomic number of the target				
G4int     G4QElasticCrossSection::lastTN=0;   // Last number of neutrons of the target
G4double  G4QElasticCrossSection::lastPIN=0.; // Last initialized max momentum
G4double* G4QElasticCrossSection::lastCST=0;  // Elastic cross-section table															
G4double* G4QElasticCrossSection::lastPAR=0;  // Parameters for functional calculation					
G4double* G4QElasticCrossSection::lastSST=0;  // E-dep of sq.slope of the first difruction	
G4double* G4QElasticCrossSection::lastS1T=0;  // E-dep of mantissa of the first difruction	
G4double* G4QElasticCrossSection::lastB1T=0;  // E-dep of the slope of the first difruction
G4double* G4QElasticCrossSection::lastS2T=0;  // E-dep of mantissa of the second difruction
G4double* G4QElasticCrossSection::lastB2T=0;  // E-dep of the slope of theSecond difruction
G4double* G4QElasticCrossSection::lastS3T=0;  // E-dep of mantissa of the third difruction	
G4double* G4QElasticCrossSection::lastB3T=0;  // E-dep of the slope of the third difruction
G4double* G4QElasticCrossSection::lastS4T=0;  // E-dep of mantissa of the 4-th difruction	
G4double* G4QElasticCrossSection::lastB4T=0;  // E-dep of the slope of the 4-th difruction
G4int     G4QElasticCrossSection::lastPDG=0;  // The last PDG code of the projectile
G4int     G4QElasticCrossSection::lastN=0;    // The last N of calculated nucleus
G4int     G4QElasticCrossSection::lastZ=0;    // The last Z of calculated nucleus
G4double  G4QElasticCrossSection::lastP=0.;   // Last used in cross section Momentum
G4double  G4QElasticCrossSection::lastTH=0.;  // Last threshold momentum
G4double  G4QElasticCrossSection::lastCS=0.;  // Last value of the Cross Section
G4int     G4QElasticCrossSection::lastI=0;    // The last position in the DAMDB

// Returns Pointer to the G4VQCrossSection class
G4VQCrossSection* G4QElasticCrossSection::GetPointer()
{
  static G4QElasticCrossSection theCrossSection; //**Static body of the QEl Cross Section**
  return &theCrossSection;
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4QElasticCrossSection::GetCrossSection(G4bool fCS, G4double pMom, G4int tgZ,
                                                 G4int tgN, G4int pPDG)
{
  static std::vector <G4int>    colPDG;// Vector of the projectile PDG code
  static std::vector <G4int>    colN;  // Vector of N for calculated nuclei (isotops)
  static std::vector <G4int>    colZ;  // Vector of Z for calculated nuclei (isotops)
  static std::vector <G4double> colP;  // Vector of last momenta for the reaction
  static std::vector <G4double> colTH; // Vector of energy thresholds for the reaction
  static std::vector <G4double> colCS; // Vector of last cross sections for the reaction
  // ***---*** End of the mandatory Static Definitions of the Associative Memory ***---***
  G4double pEn=pMom;
  onlyCS=fCS;
#ifdef pdebug
  G4cout<<"G4QElCS::GetCS:>>> f="<<fCS<<", p="<<pMom<<", Z="<<tgZ<<"("<<lastZ<<") ,N="<<tgN
        <<"("<<lastN<<"),PDG="<<pPDG<<"("<<lastPDG<<"), T="<<pEn<<"("<<lastTH<<")"<<",Sz="
        <<colN.size()<<G4endl;
		//CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  if(!pPDG)
  {
#ifdef pdebug
    G4cout<<"G4QElCS::GetCS: *** Found pPDG="<<pPDG<<" ====> CS=0"<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return 0.;                         // projectile PDG=0 is a mistake (?!) @@
  }
  G4bool in=false;                     // By default the isotope must be found in the AMDB
  if(tgN!=lastN || tgZ!=lastZ || pPDG!=lastPDG)// The nucleus was not the last used isotope
  {
    in = false;                        // By default the isotope haven't be found in AMDB  
    lastP   = 0.;                      // New momentum history (nothing to compare with)
    lastPDG = pPDG;                    // The last PDG of the projectile
    lastN   = tgN;                     // The last N of the calculated nucleus
    lastZ   = tgZ;                     // The last Z of the calculated nucleus
    lastI   = colN.size();             // Size of the Associative Memory DB in the heap
    if(lastI) for(G4int i=0; i<lastI; i++) // Loop over proj/tgZ/tgN lines of DB
	   {                                  // The nucleus with projPDG is found in AMDB
      if(colPDG[i]==pPDG && colN[i]==tgN && colZ[i]==tgZ)
						{
        lastI=i;
        lastTH =colTH[i];                // Last THreshold (A-dependent)
#ifdef pdebug
        G4cout<<"G4QElCS::GetCS:*Found* P="<<pMom<<",Threshold="<<lastTH<<",i="<<i<<G4endl;
        //CalculateCrossSection(fCS,-27,i,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(pEn<=lastTH)
        {
#ifdef pdebug
          G4cout<<"G4QElCS::GetCS:Found T="<<pEn<<" < Threshold="<<lastTH<<",CS=0"<<G4endl;
          //CalculateCrossSection(fCS,-27,i,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
          return 0.;                     // Energy is below the Threshold value
        }
        lastP  =colP [i];                // Last Momentum  (A-dependent)
        lastCS =colCS[i];                // Last CrossSect (A-dependent)
        if(std::fabs(lastP/pMom-1.)<tolerance)
        {
#ifdef pdebug
          G4cout<<"G4QElCS::GetCS:P="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
#endif
          CalculateCrossSection(fCS,-1,i,lastPDG,lastZ,lastN,pMom); // Update param's only
          return lastCS*millibarn;     // Use theLastCS
        }
        in = true;                       // This is the case when the isotop is found in DB
        // Momentum pMom is in IU ! @@ Units
#ifdef pdebug
        G4cout<<"G4QElCS::G:UpdateDB P="<<pMom<<",f="<<fCS<<",I="<<lastI<<",i="<<i<<G4endl;
#endif
        lastCS=CalculateCrossSection(fCS,-1,i,lastPDG,lastZ,lastN,pMom); // read & update
#ifdef pdebug
        G4cout<<"G4QElCS::GetCrosSec: *****> New (inDB) Calculated CS="<<lastCS<<G4endl;
        //CalculateCrossSection(fCS,-27,i,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(lastCS<=0. && pEn>lastTH)    // Correct the threshold
        {
#ifdef pdebug
          G4cout<<"G4QElCS::GetCS: New T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
        break;                           // Go out of the LOOP with found lastI
      }
#ifdef pdebug
      G4cout<<"---G4QElCrossSec::GetCrosSec:pPDG="<<pPDG<<",i="<<i<<",N="<<colN[i]
            <<",Z["<<i<<"]="<<colZ[i]<<",cPDG="<<colPDG[i]<<G4endl;
      //CalculateCrossSection(fCS,-27,i,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
	   }
	   if(!in)                            // This nucleus has not been calculated previously
	   {
#ifdef pdebug
      G4cout<<"G4QElCS::GetCrosSec:CalcNew P="<<pMom<<",f="<<fCS<<",lastI="<<lastI<<G4endl;
#endif
      //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
      lastCS=CalculateCrossSection(fCS,0,lastI,lastPDG,lastZ,lastN,pMom);//calculate&create
      if(lastCS<=0.)
						{
        lastTH = ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last
#ifdef pdebug
        G4cout<<"G4QElCrossSection::GetCrossSect: NewThresh="<<lastTH<<",T="<<pEn<<G4endl;
#endif
        if(pEn>lastTH)
        {
#ifdef pdebug
          G4cout<<"G4QElCS::GetCS: First T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
						}
#ifdef pdebug
      G4cout<<"G4QElCS::GetCrosSec: New CS="<<lastCS<<",lZ="<<lastN<<",lN="<<lastZ<<G4endl;
      //CalculateCrossSection(fCS,-27,lastI,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      colN.push_back(tgN);
      colZ.push_back(tgZ);
      colPDG.push_back(pPDG);
      colP.push_back(pMom);
      colTH.push_back(lastTH);
      colCS.push_back(lastCS);
#ifdef pdebug
      G4cout<<"G4QElCS::GetCS:1st,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
      //CalculateCrossSection(fCS,-27,lastI,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      return lastCS*millibarn;
	   } // End of creation of the new set of parameters
    else
				{
#ifdef pdebug
      G4cout<<"G4QElCS::GetCS: Update lastI="<<lastI<<G4endl;
#endif
      colP[lastI]=pMom;
      colPDG[lastI]=pPDG;
      colCS[lastI]=lastCS;
    }
  } // End of parameters udate
  else if(pEn<=lastTH)
  {
#ifdef pdebug
    G4cout<<"G4QElCS::GetCS: Current T="<<pEn<<" < Threshold="<<lastTH<<", CS=0"<<G4endl;
    //CalculateCrossSection(fCS,-27,lastI,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return 0.;                         // Momentum is below the Threshold Value -> CS=0
  }
  else if(std::fabs(lastP/pMom-1.)<tolerance)
  {
#ifdef pdebug
    G4cout<<"G4QElCS::GetCS:OldCur P="<<pMom<<"="<<pMom<<", CS="<<lastCS*millibarn<<G4endl;
    //CalculateCrossSection(fCS,-27,lastI,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
    G4cout<<"G4QElCS::GetCrSec:***SAME***, onlyCS="<<onlyCS<<G4endl;
#endif
    return lastCS*millibarn;     // Use theLastCS
  }
  else
  {
#ifdef pdebug
    G4cout<<"G4QElCS::GetCS:UpdatCur P="<<pMom<<",f="<<fCS<<",I="<<lastI<<",j="<<j<<G4endl;
#endif
    lastCS=CalculateCrossSection(fCS,1,lastI,lastPDG,lastZ,lastN,pMom); // Only UpdateDB
    lastP=pMom;
  }
#ifdef pdebug
  G4cout<<"G4QElCS::GetCrSec:End,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
  //CalculateCrossSection(fCS,-27,lastI,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
  G4cout<<"G4QElCS::GetCrSec:***End***, onlyCS="<<onlyCS<<G4endl;
#endif
  return lastCS*millibarn;
}

// Calculation of total elastic cross section (p in IU, CS in mb) @@ Units (?)
// F=0 - create AMDB, F=-1 - read&update AMDB, F=1 - update AMDB (sinchro with higher AMDB)
G4double G4QElasticCrossSection::CalculateCrossSection(G4bool CS,G4int F,G4int I,G4int PDG,
                                                       G4int tgZ, G4int tgN, G4double pIU)
{
  // *** Begin of Associative Memory DB for acceleration of the cross section calculations
  static std::vector <G4double>  PIN;   // Vector of max initialized log(P) in the table
  static std::vector <G4double*> PAR;   // Vector of parameters for functional calculations
  static std::vector <G4double*> CST;   // Vector of cross-section table
  static std::vector <G4double*> SST;   // Vector of the first squared slope
  static std::vector <G4double*> S1T;   // Vector of the first mantissa
  static std::vector <G4double*> B1T;   // Vector of the first slope
  static std::vector <G4double*> S2T;   // Vector of the secon mantissa
  static std::vector <G4double*> B2T;   // Vector of the second slope
  static std::vector <G4double*> S3T;   // Vector of the third mantissa
  static std::vector <G4double*> B3T;   // Vector of the third slope
  static std::vector <G4double*> S4T;   // Vector of the 4-th mantissa (gloria)
  static std::vector <G4double*> B4T;   // Vector of the 4-th slope    (gloria)
  // *** End of Static Definitions (Associative Memory Data Base) ***
  G4double pMom=pIU/GeV;                // All calculations are in GeV
  onlyCS=CS;                            // Flag to calculate only CS (not Si/Bi)
#ifdef pdebug
  G4cout<<"G4QElasticCrossSection::CalcCS:->onlyCS="<<onlyCS<<",F="<<F<<",p="<<pIU<<G4endl;
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
      G4cout<<"G4QElasticCS::CalcCS: DB is updated for I="<<I<<",*,PIN4="<<PIN[4]<<G4endl;
#endif
    }
#ifdef pdebug
    G4cout<<"G4QElasticCrossSection::CalcCS:*read*, LP="<<lastLP<<",PIN="<<lastPIN<<G4endl;
#endif
    if(lastLP>lastPIN && lastLP<lPMax)
    {
      lastPIN = GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);// Can update upper logP-Limit in tabs
#ifdef pdebug
						G4cout<<"G4QElCS::CalcCS:*updated(I)*,LP="<<lastLP<<"<IN["<<I<<"]="<<lastPIN<<G4endl;
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
    G4cout<<"G4QElasticCrossSection::CalcCS:*ini*,lastLP="<<lastLP<<",min="<<lPMin<<G4endl;
#endif
    lastPIN = GetPTables(lastLP,lPMin,PDG,tgZ,tgN); // Returns the new P-limit for tables
#ifdef pdebug
    G4cout<<"G4QElCS::CalcCS:*i*Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<",LP"<<lastPIN<<G4endl;
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
  G4cout<<"G4QElCS::CalcCS:?update?,LP="<<lastLP<<",IN="<<lastPIN<<",ML="<<lPMax<<G4endl;
#endif
  if(lastLP>lastPIN && lastLP<lPMax)
  {
    lastPIN = GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);
#ifdef pdebug
    G4cout<<"G4QElCS::CalcCS: *updated(O)*, LP="<<lastLP<<" < IN="<<lastPIN<<G4endl;
#endif
  }
#ifdef pdebug
  G4cout<<"G4QElastCS::CalcCS: lastLP="<<lastLP<<",lPM="<<lPMin<<",lPIN="<<lastPIN<<G4endl;
#endif
  if(!onlyCS) lastTM=GetQ2max(PDG, tgZ, tgN, pMom); // Calculate (-t)_max=Q2_max (GeV2)
#ifdef pdebug
  G4cout<<"G4QElasticCrosSec::CalcCS:oCS="<<onlyCS<<",-t="<<lastTM<<", p="<<lastLP<<G4endl;
#endif
  if(lastLP>lPMin && lastLP<=lastPIN)   // Linear fit is made using precalculated tables
  {
    if(std::fabs(lastLP-lastPIN)<.0001) // Just take the highest tabulated value
    {
      G4double shift=(lastLP-lPMin)/dlnP+.000001; // Log distance from lPMin
      G4int    blast=static_cast<int>(shift); // this is a bin number of the lower edge (0)
      if(blast<0 || blast>=nLast) G4cout<<"G4QEleastCS::CCS:b="<<blast<<","<<nLast<<G4endl;
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
      G4cout<<"G4QElasticCrossSection::CalculateCS:(E) S1="<<theS1<<", B1="<<theB1<<G4endl;
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
      G4cout<<"G4QElCS::CalcCrossSection: Sig="<<lastSIG<<", P="<<pMom<<", Z="<<tgZ<<", N="
            <<tgN<<", PDG="<<PDG<<", onlyCS="<<onlyCS<<G4endl;
#endif
      if(!onlyCS)                       // Skip the differential cross-section parameters
      {
        G4double SSTL=lastSST[blast];           // the low bin of the first squared slope
        theSS=SSTL+shift*(lastSST[lastL]-SSTL); // the basic value of the first sq.slope
        G4double S1TL=lastS1T[blast];           // the low bin of the first mantissa
        theS1=S1TL+shift*(lastS1T[lastL]-S1TL); // the basic value of the first mantissa
        G4double B1TL=lastB1T[blast];           // the low bin of the first slope
#ifdef pdebug
        G4cout<<"G4QElCS::CalcCrossSection:bl="<<blast<<",ls="<<lastL<<",SL="<<S1TL<<",SU="
              <<lastS1T[lastL]<<",BL="<<B1TL<<",BU="<<lastB1T[lastL]<<G4endl;
#endif
        theB1=B1TL+shift*(lastB1T[lastL]-B1TL); // the basic value of the first slope
        G4double S2TL=lastS2T[blast];           // the low bin of the second mantissa
        theS2=S2TL+shift*(lastS2T[lastL]-S2TL); // the basic value of the second mantissa
        G4double B2TL=lastB2T[blast];           // the low bin of the second slope
        theB2=B2TL+shift*(lastB2T[lastL]-B2TL); // the basic value of the second slope
        G4double S3TL=lastS3T[blast];           // the low bin of the third mantissa
        theS3=S3TL+shift*(lastS3T[lastL]-S3TL); // the basic value of the third mantissa
#ifdef pdebug
        G4cout<<"G4QElCS::CCS: s3l="<<S3TL<<",sh3="<<shift<<",s3h="<<lastS3T[lastL]<<",b="
              <<blast<<",l="<<lastL<<G4endl;
#endif
        G4double B3TL=lastB3T[blast];           // the low bin of the third slope
        theB3=B3TL+shift*(lastB3T[lastL]-B3TL); // the basic value of the third slope
        G4double S4TL=lastS4T[blast];           // the low bin of the 4-th mantissa
        theS4=S4TL+shift*(lastS4T[lastL]-S4TL); // the basic value of the 4-th mantissa
#ifdef pdebug
        G4cout<<"G4QElCS::CCS: s4l="<<S4TL<<",sh4="<<shift<<",s4h="<<lastS4T[lastL]<<",b="
              <<blast<<",l="<<lastL<<G4endl;
#endif
        G4double B4TL=lastB4T[blast];           // the low bin of the 4-th slope
        theB4=B4TL+shift*(lastB4T[lastL]-B4TL); // the basic value of the 4-th slope
      }
#ifdef pdebug
      G4cout<<"G4QElasticCrossSection::CalculateCS:(I) S1="<<theS1<<", B1="<<theB1<<G4endl;
#endif
    }
  }
  else lastSIG=GetTabValues(lastLP, PDG, tgZ, tgN); // Direct calculation beyond the table
  if(lastSIG<0.) lastSIG = 0.;                   // @@ a Warning print can be added
#ifdef pdebug
  G4cout<<"G4QElasticCrossSection::CalculateCS: END, onlyCS="<<onlyCS<<G4endl;
#endif
  return lastSIG;
}

// It has parameter sets for all tZ/tN/PDG, using them the tables can be created/updated
G4double G4QElasticCrossSection::GetPTables(G4double LP,G4double ILP, G4int PDG, G4int tgZ,
                                                                                 G4int tgN)
{
  // @@ At present all nA==pA ---------> Each neucleus can have not more than 51 parameters
  static const G4double pwd=2727;
  const G4int n_npel=24;                // #of parameters for np-elastic (<nPoints=51)
  const G4int n_ppel=32;                // #of parameters for pp-elastic (<nPoints=51)
  //const G4int n_pdel=30;                // #of parameters for pp-elastic (<nPoints=51)
  //const G4int n_phe4=32;                // #of parameters for phe4-elastic (<nPoints=51)
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
  //                      -0- -1--2- -3-  -4-  -5-  -6-    -7-  -8- -9- -10--11--12--13-
  //G4double pd_el[n_pdel]={7.5,.09,5.,.06,.0013,.1,.000006,1.287,.001,34.,.4,4.726,7.5,.1,
  //                        .05,.000017,.0004,1.15,5.5,.13,.02,.1911,4.857,46.,40.,2.,.01,
  //                        .0000007,13.,.1};
  //                       -14- -15-  -16-  -17--18- -19--20- -21- -22- -23--24--25--26-
  //                        -27-  -28- -29-
  //                      -0- -1- -2- -3-  -4-   -5-   -6-  -7-  -8- -9- -10--11--12-
  //G4double p_he4[n_phe4]={22.5,.3,5.,.037,.004,.00003,12.9,2500.,.05,22.5,.3,.04,.0065,
		//																								12.5,2500.,.053,.035,2.5e-8,29.,.6,.7,.1,1.5,.03,.00015,4.24,
  //                        .004,.005,1.5e-11,2.3,.002,23.};
  //                      -13- -14-  -15- -16- -17- -18--19--20-21--22--23-  -24- -25-
  //                      -26- -27-   -28- -29- -30- -31-
  if(PDG==2212 || PDG==2112)
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
      if     (PDG==2112&&tgZ==1&&tgN==0)
                                for(G4int ip=0;ip<n_npel;ip++)lastPAR[ip]=np_el[ip];// np
						else if(PDG==2212&&tgZ==1&&tgN==0)
                                for(G4int ip=0;ip<n_ppel;ip++)lastPAR[ip]=pp_el[ip];// pp
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
        // a17
        G4double a32=a16*a16;
        G4double r1a16=3.e16/a16;
        G4double r2a16=2.e13/a16;
        G4double r3a16=6.e13/a16;
        G4double pa10=5.e-27*a10;
        // Reaction cross-section parameters (pel=peh_fit.f)
        lastPAR[0]=.28*a;                                                     // p1
        lastPAR[1]=4.8*std::pow(a,1.14)/(1.+3.6/a3);                          // p2
        lastPAR[2]=3.3/a+.8*ssa/(1.+2.e7/a8);                                 // p3
        lastPAR[3]=1./(1.+6.e-9*a12)+.6*a/(1.+1.6e15/a16);                    // p4
        lastPAR[4]=6.e-4/(a4+3.e-6*a12)+.16/(a+9.e5/a3+r1a16*r1a16);          // p5
        lastPAR[5]=3.e-4/a2+(1.5e-3+3.e-14*a7)/(1.+r2a16*r2a16+3.e-12*a6);    // p6
        lastPAR[6]=4.e-30+(1.2e-28*a3+pa10*pa10)/(1.+r3a16*r3a16+4.e-26*a14); // p7
        lastPAR[7]=.07/(1.+1.7e-8*a16)+.5/(1.+1.7e18/a16+4.5e-7*a4);          // p8
        lastPAR[8]=(1.5e-10+2.e-18*a8)/(1.+3.e-25*a16);                       // p10
        // @@ the differential cross-section is parameterized separately for A>6 & A<7
        if(a<6.5)
								{
          // a11
          // a13
          G4double a28=a16*a12;
          // a31
          // a40
          // The main pre-exponent      (pel_sg)
          lastPAR[ 9]=4000*a;                                // p1
          lastPAR[10]=1.2e7*a8+380*a16*a;                    // p2
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
          lastPAR[24]=8.e4/(a8+2.5e12/a16);                  // p1
          lastPAR[25]=8.e7/(a12+1.e-27*a28*a28);             // p2
          lastPAR[26]=.0011*a3;                              // p3
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
          lastPAR[38]=1.05e8*std::exp(.32*asa);              // p1
          lastPAR[39]=19.5*std::exp(.45*asa);                // p2
          lastPAR[40]=7.e3+2.4e6/a5;                         // p3
          lastPAR[41]=2.5e5*std::exp(.085*a3);               // p4
          lastPAR[42]=2.5*a;                                 // p5
          // The gloria slope           (pel_ub)
          lastPAR[43]=920.+.03*a8*a3;                        // p1
          lastPAR[44]=93.+.0023*a12;                         // p2
#ifdef pdebug
         G4cout<<"G4QElCS::CalcCS:la "<<lastPAR[38]<<", "<<lastPAR[39]<<", "<<lastPAR[40]
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
         G4cout<<"G4QElCS::CalcCS:ha "<<lastPAR[41]<<", "<<lastPAR[42]<<", "<<lastPAR[43]
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
      G4cout<<"G4QElasticCrossSection::GetPTables:ip=0(init), lp="<<lp<<",S1="<<theS1
												<<",B1="<<theB1<<",S2="<<theS2<<",B2="<<theB3<<",S3="<<theS3
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
        if(fin>nPoints) fin=nLast;                // Limit of the tabular initialization
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
            G4cout<<"G4QElasticCrossSection::GetPTables:ip="<<ip<<",lp="<<lp<<",S1="<<theS1
                  <<",B1="<<theB1<<",S2="<<theS2<<",B2="<<theB2<<",S3="
                  <<theS3<<",B3="<<theB3<<",S4="<<theS4<<",B4="<<theB4<<G4endl;
#endif
          }
          return lp;
        }
        else G4cout<<"*Warning*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ
                   <<", N="<<tgN<<", i="<<ini<<" > fin="<<fin<<", LP="<<LP<<" > ILP="<<ILP
                   <<" nothing is done!"<<G4endl;
      }
      else G4cout<<"*Warning*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ
                 <<", N="<<tgN<<", i="<<ini<<">= max="<<nPoints<<", LP="<<LP<<" > ILP="
                 <<ILP<<", lPMax="<<lPMax<<" nothing is done!"<<G4endl;
    }
#ifdef pdebug
    else G4cout<<"*Warning*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ
               <<", N="<<tgN<<", LP="<<LP<<" <= ILP="<<ILP<<" nothing is done!"<<G4endl;
#endif
  }
  else
  {
    G4cout<<"*Error*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for PDG=2212(p)/2112(n)"<<G4endl;
    throw G4QException("G4QElasticCrossSection::GetPTables: only nA & pA are implemented");
  }
  return ILP;
}

// Returns Q2=-t in independent units (MeV^2) (all internal calculations are in GeV)
G4double G4QElasticCrossSection::GetExchangeT(G4int tgZ, G4int tgN, G4int PDG)
{
  static const G4double GeVSQ=gigaelectronvolt*gigaelectronvolt;
  static const G4double third=1./3.;
  static const G4double fifth=1./5.;
  static const G4double sevth=1./7.;
#ifdef tdebug
  G4cout<<"G4QElasticCS::GetExchangeT:F="<<onlyCS<<",Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<G4endl;
#endif
  if(onlyCS) G4cout<<"*Warning*G4QElasticCrossSection::GetExchangeQ2: onlyCS=true"<<G4endl;
  if(lastLP<-4.3) return lastTM*GeVSQ*G4UniformRand();// S-wave for p<14 MeV/c (kinE<.1MeV)
  G4double q2=0.;
  if(PDG==2112 && tgZ==1 && tgN==0)                // ===> n+p=n+p
  {
#ifdef tdebug
    G4cout<<"G4QElasticCS::GetExchangeT: TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1<<",S2="
          <<theS2<<",B2="<<theB2<<",GeV2="<<GeVSQ<<G4endl;
#endif
    G4double E1=lastTM*theB1;
  		G4double R1=(1.-std::exp(-E1));
#ifdef tdebug
    G4double ts1=-std::log(1.-R1)/theB1;
    G4double ds1=std::fabs(ts1-lastTM)/lastTM;
    if(ds1>.0001)
      G4cout<<"*Warning*G4QElCS::GetExT:1n "<<ts1<<"#"<<lastTM<<",d="<<ds1
            <<",R1="<<R1<<",E1="<<E1<<G4endl;
#endif
    G4double E2=lastTM*theB2;
  		G4double R2=(1.-std::exp(-E2));
#ifdef tdebug
    G4double ts2=-std::log(1.-R2)/theB2;
    G4double ds2=std::fabs(ts2-lastTM)/lastTM;
    if(ds2>.0001)
      G4cout<<"*Warning*G4QElCS::GetExT:2n "<<ts2<<"#"<<lastTM<<",d="<<ds2
            <<",R2="<<R2<<",E2="<<E2<<G4endl;
#endif
    //G4double E3=lastTM*theB3;
  		//G4double R3=(1.-std::exp(-E3));
#ifdef tdebug
    //G4double ts3=-std::log(1.-R3)/theB3;
    //G4double ds3=std::fabs(ts3-lastTM)/lastTM;
    //if(ds3>.01)G4cout<<"*Warn*G4QElCS::GetExT:3n "<<ts3<<"#"<<lastTM<<",d="<<ds3<<G4endl;
#endif
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
  else if(PDG==2212 && tgZ==1 && tgN==0)                // ===> p+p=p+p
  {
#ifdef tdebug
    G4cout<<"G4QElasticCS::GetExchangeT: TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1<<",S2="
          <<theS2<<",B2="<<theB2<<",S3="<<theS3<<",B3="<<theB3<<",GeV2="<<GeVSQ<<G4endl;
#endif
    G4double E1=lastTM*theB1;
  		G4double R1=(1.-std::exp(-E1));
#ifdef tdebug
    G4double ts1=-std::log(1.-R1)/theB1;
    G4double ds1=std::fabs(ts1-lastTM)/lastTM;
    if(ds1>.0001)
      G4cout<<"*Warning*G4QElCS::GetExT:1p "<<ts1<<"#"<<lastTM<<",d="<<ds1
            <<",R1="<<R1<<",E1="<<E1<<G4endl;
#endif
    G4double E2=lastTM*theB2;
  		G4double R2=(1.-std::exp(-E2*E2*E2));
#ifdef tdebug
    G4double ts2=std::pow(-std::log(1.-R2),.333333333)/theB2;
    G4double ds2=std::fabs(ts2-lastTM)/lastTM;
    if(ds2>.0001)
      G4cout<<"*Warning*G4QElCS::GetExT:2p "<<ts2<<"#"<<lastTM<<",d="<<ds2
            <<",R2="<<R2<<",E2="<<E2<<G4endl;
#endif
    G4double E3=lastTM*theB3;
  		G4double R3=(1.-std::exp(-E3));
#ifdef tdebug
    G4double ts3=-std::log(1.-R3)/theB3;
    G4double ds3=std::fabs(ts3-lastTM)/lastTM;
    if(ds3>.0001)
      G4cout<<"*Warning*G4QElCS::GetExT:3p "<<ts3<<"#"<<lastTM<<",d="<<ds3
            <<",R3="<<R1<<",E3="<<E3<<G4endl;
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
  else if(PDG==2212 || PDG==2112)
  {
    G4double a=tgZ+tgN;
#ifdef tdebug
    G4cout<<"G4QElCS::GetExT: a="<<a<<",t="<<lastTM<<",S1="<<theS1<<",B1="<<theB1<<",SS="
          <<theSS<<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS3<<",B3="<<theB3<<",S4="
          <<theS4<<",B4="<<theB4<<G4endl;
#endif
    G4double E1=lastTM*(theB1+lastTM*theSS);
  		G4double R1=(1.-std::exp(-E1));
    G4double tss=theSS+theSS; // for future solution of quadratic equation (imediate check)
#ifdef tdebug
    G4double ts1=-std::log(1.-R1)/theB1;
    if(std::fabs(tss)>1.e-7) ts1=(std::sqrt(theB1*(theB1+(tss+tss)*ts1))-theB1)/tss;
    G4double ds1=(ts1-lastTM)/lastTM;
    if(ds1>.0001)
      G4cout<<"*Warning*G4QElCS::GetExT:1a "<<ts1<<"#"<<lastTM<<",d="<<ds1
            <<",R1="<<R1<<",E1="<<E1<<G4endl;
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
      G4cout<<"*Warning*G4QElCS::GetExT:2a "<<ts2<<"#"<<lastTM<<",d="<<ds2
            <<",R2="<<R2<<",E2="<<E2<<G4endl;
#endif
    G4double E3=lastTM*theB3;
    if(a>6.5)E3*=tm2*tm2*tm2;                       // power 1 for lowA, 7 (2nd) for HighA
  		G4double R3=(1.-std::exp(-E3));
#ifdef tdebug
    G4double ts3=-std::log(1.-R3)/theB3;
    if(a>6.5)ts3=std::pow(ts3,sevth);
    G4double ds3=std::fabs(ts3-lastTM)/lastTM;
    if(ds3>.0001)
      G4cout<<"*Warning*G4QElCS::GetExT:3a "<<ts3<<"#"<<lastTM<<",d="<<ds3
            <<",R3="<<R3<<",E3="<<E3<<G4endl;
#endif
    G4double E4=lastTM*theB4;
  		G4double R4=(1.-std::exp(-E4));
#ifdef tdebug
    G4double ts4=-std::log(1.-R4)/theB4;
    G4double ds4=std::fabs(ts4-lastTM)/lastTM;
    if(ds4>.0001)
      G4cout<<"*Warning*G4QElCS::GetExT:4a "<<ts4<<"#"<<lastTM<<",d="<<ds4
            <<",R4="<<R4<<",E4="<<E4<<G4endl;
#endif
  		G4double I1=R1*theS1;
  		G4double I2=R2*theS2;
				G4double I3=R3*theS3;
				G4double I4=R4*theS4;
    G4double I12=I1+I2;
    G4double I13=I12+I3;
    G4double rand=(I13+I4)*G4UniformRand();
#ifdef tdebug
    G4cout<<"G4QElCS::GtExT:1="<<I1<<",2="<<I2<<",3="<<I3<<",4="<<I4<<",r="<<rand<<G4endl;
#endif
    if(rand<I1)
    {
      G4double ran=R1*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB1;
      if(std::fabs(tss)>1.e-7) q2=(std::sqrt(theB1*(theB1+(tss+tss)*q2))-theB1)/tss;
#ifdef tdebug
      G4cout<<"G4QElCS::GetExT:Q2="<<q2<<",ss="<<tss/2<<",b1="<<theB1<<",t1="<<ts1<<G4endl;
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
      G4cout<<"G4QElCS::GetExT: Q2="<<q2<<", r2="<<R2<<", b2="<<theB2<<",t2="<<ts2<<G4endl;
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
      G4cout<<"G4QElCS::GetExT:Q2="<<q2<<", r3="<<R2<<", b3="<<theB2<<",t3="<<ts2<<G4endl;
#endif
    }
    else
    {
      G4double ran=R4*G4UniformRand();
      if(ran>1.) ran=1.;
      q2=-std::log(1.-ran)/theB4;
      if(a<6.5) q2=lastTM-q2;                    // u reduced for lightA (starts from 0)
#ifdef tdebug
      G4cout<<"G4QElCS::GetExT:Q2="<<q2<<",m="<<lastTM<<",b4="<<theB3<<",t4="<<ts3<<G4endl;
#endif
    }
  }
  else
  {
    G4cout<<"*Error*G4QElasticCrossSection::GetExchangeT: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for PDG=2212, Z=1,N=(0,1) & Z=2,N=2"<<G4endl;
    throw G4QException("G4QElasticCrossSection::GetExchangeT: n-p,p-H/He are implemented");
  }
  if(q2<0.) q2=0.;
  if(!(q2>-1.||q2<1.)) G4cout<<"*NAN*G4QElasticCrossSect::GetExT: -t="<<q2<<G4endl;
  if(q2>lastTM)G4cout<<"*Warning*G4QElasticCrossSect::GetExT:-t="<<q2<<">"<<lastTM<<G4endl;
  return q2*GeVSQ;
}

// Returns half max(Q2=-t) in independent units (MeV^2)
G4double G4QElasticCrossSection::GetHMaxT()
{
  static const G4double HGeVSQ=gigaelectronvolt*gigaelectronvolt/2.;
  return lastTM*HGeVSQ;
}

// lastLP is used, so calculating tables, one need to remember and then recover lastLP
G4double G4QElasticCrossSection::GetTabValues(G4double lp, G4int PDG, G4int tgZ,G4int tgN)
{
  static const G4double GeVSQ=gigaelectronvolt*gigaelectronvolt;
  static const G4double mNeut= G4QPDGCode(2112).GetMass(); // Neutron mass in MeV
  static const G4double dNM= (mNeut+mNeut)/GeVSQ;// Doubled neutron mass in ...
  static const G4int    nZ=92;                   // #0f covered elements
  static const G4int    n01=2;                   // Z= 1 H
  static const G4int    N01[n01]={0,1};          // #of neutrons in the isotope
  static const G4double T01[n01]={0.,0.};        // Threshold in MeV
  static const G4double X01[n01]={20000.,3395.}; // XS level in mb
  static const G4int    n02=2;                   // Z= 2 He
  static const G4int    N02[n02]={1,2};
  static const G4double T02[n02]={0.,0.};
  static const G4double X02[n02]={1000.,759.};
  static const G4int    n03=2;                   // Z= 3 Li
  static const G4int    N03[n03]={3,4};
  static const G4double T03[n03]={0.,0.};
  static const G4double X03[n03]={710.,1050.};
  static const G4int    n04=1;                   // Z= 4 Be
  static const G4int    N04[n04]={5};
  static const G4double T04[n04]={0.};
  static const G4double X04[n04]={6000.};
  static const G4int    n05=2;                   // Z= 5 B
  static const G4int    N05[n05]={5,6};
  static const G4double T05[n05]={0.,0.};
  static const G4double X05[n05]={2140.,4840.};
  static const G4int    n06=2;                   // Z= 6 C
  static const G4int    N06[n06]={6,7};
  static const G4double T06[n06]={0.,0.};
  static const G4double X06[n06]={4730.,4730.};
  static const G4int    n07=2;                   // Z= 7 N
  static const G4int    N07[n07]={7,8};
  static const G4double T07[n07]={0.,0.};
  static const G4double X07[n07]={9957.,1018.};
  static const G4int    n08=3;                   // Z= 8 O
  static const G4int    N08[n08]={8,9,10};
  static const G4double T08[n08]={0.,0.,0.};
  static const G4double X08[n08]={3880.,3740.,4000.};
  static const G4int    n09=1;                   // Z= 9 F
  static const G4int    N09[n09]={10};
  static const G4double T09[n09]={0.};
  static const G4double X09[n09]={3997.};
  static const G4int    n10=3;                   // Z=10 Ne
  static const G4int    N10[n10]={10,11,12};
  static const G4double T10[n10]={0.,0.,0.};
  static const G4double X10[n10]={1000.,1000.,1000.};
  static const G4int    n11=1;                   // Z=11 Na
  static const G4int    N11[n11]={12};
  static const G4double T11[n11]={0.};
  static const G4double X11[n11]={728.};
  static const G4int    n12=3;                   // Z=12 Mg
  static const G4int    N12[n12]={12,13,14};
  static const G4double T12[n12]={.52,.22,.45};
  static const G4double X12[n12]={3500.,4000.,4000.};
  static const G4int    n13=1;                   // Z=13 Al
  static const G4int    N13[n13]={14};
  static const G4double T13[n13]={0.};
  static const G4double X13[n13]={1400.};
  static const G4int    n14=3;                   // Z=14 Si
  static const G4int    N14[n14]={14,15,16};
  static const G4double T14[n14]={1.75,.1,.5};
  static const G4double X14[n14]={1770.,4700.,5000.};
  static const G4int    n15=1;                   // Z=15 P
  static const G4int    N15[n15]={16};
  static const G4double T15[n15]={0.};
  static const G4double X15[n15]={703.};
  static const G4int    n16=4;                   // Z=16 S
  static const G4int    N16[n16]={16,17,18,20};
  static const G4double T16[n16]={1.57,.26,.48,1.};
  static const G4double X16[n16]={4230.,6090.,6750.,5000.};
  static const G4int    n17=2;                   // Z=17 Cl
  static const G4int    N17[n17]={18,20};
  static const G4double T17[n17]={.226,.202};
  static const G4double X17[n17]={2000.,2000.};
  static const G4int    n18=3;                   // Z=18 Ar
  static const G4int    N18[n18]={18,20,22};
  static const G4double T18[n18]={.00001,.3,.66};
  static const G4double X18[n18]={50.,1000.,2500.};
  static const G4int    n19=3;                   // Z=19 K
  static const G4int    N19[n19]={20,21,22};
  static const G4double T19[n19]={.2,0.,.125};
  static const G4double X19[n19]={2200.,2200.,3000.};
  static const G4int    n20=6;                   // Z=20 Ca
  static const G4int    N20[n20]={20,22,23,24,26,28};
  static const G4double T20[n20]={.5,.3,.04,.5,0.,.5};
  static const G4double X20[n20]={700.,3500.,8000.,3000.,3900.,3700.};
  static const G4int    n21=1;                   // Z=21 Sc
  static const G4int    N21[n21]={24};
  static const G4double T21[n21]={.0969};
  static const G4double X21[n21]={4600.};
  static const G4int    n22=5;                   // Z=22 Ti
  static const G4int    N22[n22]={24,25,26,27,28};
  static const G4double T22[n22]={.1,.1,.1,0.,.1};
  static const G4double X22[n22]={4700.,5000.,5500.,1000.,6000.};
  static const G4int    n23=2;                   // Z=23 V
  static const G4int    N23[n23]={27,28};
  static const G4double T23[n23]={.1,.1};
  static const G4double X23[n23]={5000.,7000.};
  static const G4int    n24=4;                   // Z=24 Cr
  static const G4int    N24[n24]={26,28,29,30};
  static const G4double T24[n24]={.595,.637,.248,.4};
  static const G4double X24[n24]={4000.,3900.,5600.,2260.};
  static const G4int    n25=1;                   // Z=25 Mn
  static const G4int    N25[n25]={30};
  static const G4double T25[n25]={.08};
  static const G4double X25[n25]={1070.};
  static const G4int    n26=4;                   // Z=26 Fe
  static const G4int    N26[n26]={28,30,31,32};
  static const G4double T26[n26]={.5,.862,0.,.37};
  static const G4double X26[n26]={4540.,2500.,600.,5500.};
  static const G4int    n27=1;                   // Z=27 Co
  static const G4int    N27[n27]={32};
  static const G4double T27[n27]={.08};
  static const G4double X27[n27]={1070.};
  static const G4int    n28=5;                   // Z=28 Ni
  static const G4int    N28[n28]={30,32,33,34,36};
  static const G4double T28[n28]={.812,.45,.07,.6,.6};
  static const G4double X28[n28]={3000.,6000.,20000.,1200.,1200.};
  static const G4int    n29=2;                   // Z=29 Cu
  static const G4int    N29[n29]={34,36};
  static const G4double T29[n29]={.0995,.0995};
  static const G4double X29[n29]={3500.,10000.};
  static const G4int    n30=5;                   // Z=30 Zn
  static const G4int    N30[n30]={34,36,37,38,40};
  static const G4double T30[n30]={.1,.1,.1,.1,.1};
  static const G4double X30[n30]={7000.,7000.,7000.,7000.,7000.};
  static const G4int    n31=2;                   // Z=31 Ga
  static const G4int    N31[n31]={38,40};
  static const G4double T31[n31]={.0059,.0056};
  static const G4double X31[n31]={11000.,10000.};
  static const G4int    n32=5;                   // Z=32 Ge
  static const G4int    N32[n32]={38,40,41,42,44};
  static const G4double T32[n32]={.015,.00553,.00059,.00627,.0175};
  static const G4double X32[n32]={10000.,1470.,10000.,13000.,10000.};
  static const G4int    n33=1;                   // Z=33 As
  static const G4int    N33[n33]={42};
  static const G4double T33[n33]={.00241};
  static const G4double X33[n33]={15000.};
  static const G4int    n34=6;                   // Z=34 Se
  static const G4int    N34[n34]={40,42,43,44,46,48};
  static const G4double T34[n34]={.0024,.00749,.00272,.0121,.00671,.0311};
  static const G4double X34[n34]={14000.,11900.,14000.,10000.,11000.,8400.};
  static const G4int    n35=2;                   // Z=35 Br
  static const G4int    N35[n35]={44,46};
  static const G4double T35[n35]={.000412,.00363};
  static const G4double X35[n35]={17000.,11200.};
  static const G4int    n36=6;                   // Z=36 Kr
  static const G4int    N36[n36]={42,44,46,47,48,50};
  static const G4double T36[n36]={.000865,.001,.0001,.000521,.002,.013};
  static const G4double X36[n36]={9000.,20000.,9000.,7200.,11180.,11000.};
  static const G4int    n37=2;                   // Z=37 Rb
  static const G4int    N37[n37]={48,50};
  static const G4double T37[n37]={.00923,.03};
  static const G4double X37[n37]={9300.,9000.};
  static const G4int    n38=4;                   // Z=38 Sr
  static const G4int    N38[n38]={46,48,49,50};
  static const G4double T38[n38]={.0035,.0199,.00289,.1655};
  static const G4double X38[n38]={18000.,8500.,10600.,8000.};
  static const G4int    n39=1;                   // Z=39 Y
  static const G4int    N39[n39]={50};
  static const G4double T39[n39]={0.0};
  static const G4double X39[n39]={583.};
  static const G4int    n40=5;                   // Z=40 Zr
  static const G4int    N40[n40]={50,51,52,54,56};
  static const G4double T40[n40]={.131,.1,.039,.0677,.1};
  static const G4double X40[n40]={8650.,11300.,8700.,9000.,8600.};
  static const G4int    n41=1;                   // Z=41 Nb
  static const G4int    N41[n41]={52};
  static const G4double T41[n41]={.1};
  static const G4double X41[n41]={9200.};
  static const G4int    n42=7;                   // Z=42 Mo
  static const G4int    N42[n42]={50,52,53,54,55,56,58};
  static const G4double T42[n42]={.000032,.00808,.00216,.00816,.00196,.05,.1};
  static const G4double X42[n42]={7870,8000.,8100.,7800.,7700.,8000.,8500.};
  static const G4int    n43=1;                   // Z=43 Te - No stable isotopes
  static const G4int    N43[n43]={-1};
  static const G4double T43[n43]={0.};
  static const G4double X43[n43]={0.};
  static const G4int    n44=7;                   // Z=44 Mo
  static const G4int    N44[n44]={52,54,55,56,57,58,60};
  static const G4double T44[n44]={0.,0.,.00101,.0005,.025,.00161,.00151};
  static const G4double X44[n44]={8000,20000.,7000.,8000.,8000.,8000.,7200.};
  static const G4int    n45=1;                   // Z=45 Rh
  static const G4int    N45[n45]={58};
  static const G4double T45[n45]={.039};
  static const G4double X45[n45]={8000.};
  static const G4int    n46=6;                   // Z=46 Pd
  static const G4int    N46[n46]={56,58,59,60,62,64};
  static const G4double T46[n46]={.0004,.1,.02,.1,.1,.1};
  static const G4double X46[n46]={12000.,8000.,7400.,7500.,7400.,7300.};
  static const G4int    n47=2;                   // Z=47 Ag
  static const G4int    N47[n47]={60,62};
  static const G4double T47[n47]={.00106,.000984};
  static const G4double X47[n47]={12000.,670.};
  static const G4int    n48=8;                   // Z=48 Cd
  static const G4int    N48[n48]={58,60,62,63,64,65,66,68};
  static const G4double T48[n48]={.00069,.1,.00718,.01,.00662,.001,.00341,.00289};
  static const G4double X48[n48]={8000.,7500,6000.,7000.,6000.,6500.,6500.,6500.};
  static const G4int    n49=2;                   // Z=49 In
  static const G4int    N49[n49]={64,66};
  static const G4double T49[n49]={.00033,.00112};
  static const G4double X49[n49]={5700.,5700.};
  static const G4int    n50=10;                  // Z=50 Sn
  static const G4int    N50[n50]={62,64,65,66,67,68,69,70,72,74};
  static const G4double T50[n50]={.15,.002,.000693,.00215,.000611,.00189,.000547,.014,
                                  .0086,.00651};
  static const G4double X50[n50]={6600.,6700,6700.,6700.,6700.,6500.,6500.,6300.,6300.,
                                  6200.};
  static const G4int    n51=2;                   // Z=51 Sb
  static const G4int    N51[n51]={70,72};
  static const G4double T51[n51]={.000922,.00151};
  static const G4double X51[n51]={6000.,7000.};
  static const G4int    n52=8;                   // Z=52 Te
  static const G4int    N52[n52]={68,70,71,72,73,74,76,78};
  static const G4double T52[n52]={0.,.00379,.000633,.0057,.00134,.00635,.0004,.034};
  static const G4double X52[n52]={8000.,6000,6000.,6000.,7000.,6200.,300.,6000.};
  static const G4int    n53=1;                   // Z=53 I
  static const G4int    N53[n53]={74};
  static const G4double T53[n53]={.06};
  static const G4double X53[n53]={5400.};
  static const G4int    n54=9;                   // Z=54 Xe
  static const G4int    N54[n54]={70,72,74,75,76,77,78,80,82};
  static const G4double T54[n54]={0.,.00065,.004,.0042,.004,.00219,.00426,.002,0.};
  static const G4double X54[n54]={4200.,2350,14000.,12000.,14000.,16000.,15000.,20000.,
                                  4300.};
  static const G4int    n55=1;                   // Z=55 Cs
  static const G4int    N55[n55]={78};
  static const G4double T55[n55]={.00351};
  static const G4double X55[n55]={3000.};
  static const G4int    n56=7;                   // Z=56 Ba
  static const G4int    N56[n56]={74,76,78,79,80,81,82};
  static const G4double T56[n56]={.1,.1,.00207,.00486,.00318,.00195,.00303};
  static const G4double X56[n56]={4000,5000.,13000.,20000.,13000.,16000.,300.};
  static const G4int    n57=2;                   // Z=57 La
  static const G4int    N57[n57]={81,82};
  static const G4double T57[n57]={.1,.25};
  static const G4double X57[n57]={5000.,5600.};
  static const G4int    n58=4;                   // Z=58 Ce
  static const G4int    N58[n58]={78,80,82,84};
  static const G4double T58[n58]={.1,.1,.02,.01};
  static const G4double X58[n58]={6000.,6000.,9000.,5290.};
  static const G4int    n59=1;                   // Z=59 Pr
  static const G4int    N59[n59]={82};
  static const G4double T59[n59]={.00577};
  static const G4double X59[n59]={600.};
  static const G4int    n60=7;                   // Z=60 Nd
  static const G4int    N60[n60]={82,83,84,85,86,88,90};
  static const G4double T60[n60]={.011,.01,.05,.05,.05,.05,.00512};
  static const G4double X60[n60]={12000,19000.,6000.,8000.,6000.,7000.,3400.};
  static const G4int    n61=1;                   // Z=61 Pm - No stable isotopes
  static const G4int    N61[n61]={-1};
  static const G4double T61[n61]={0.};
  static const G4double X61[n61]={0.};
  static const G4int    n62=7;                   // Z=62 Nd
  static const G4int    N62[n62]={82,85,86,87,88,90,92};
  static const G4double T62[n62]={0.,.07,0.,.00015,.000593,.00369,.0031};
  static const G4double X62[n62]={8000,9000.,5300.,5300.,30000.,3000.,2000.};
  static const G4int    n63=2;                   // Z=63 Eu
  static const G4int    N63[n63]={88,90};
  static const G4double T63[n63]={.01,.01};
  static const G4double X63[n63]={15000.,14000.};
  static const G4int    n64=7;                   // Z=64 Gd
  static const G4int    N64[n64]={88,90,91,92,93,94,96};
  static const G4double T64[n64]={.1,.000276,.1,.01,.000215,.00604,.00288};
  static const G4double X64[n64]={7400,30000.,6200.,11000.,20000.,20000.,40000.};
  static const G4int    n65=1;                   // Z=65 Tb
  static const G4int    N65[n65]={94};
  static const G4double T65[n65]={.01};
  static const G4double X65[n65]={13000.};
  static const G4int    n66=7;                   // Z=66 Dy
  static const G4int    N66[n66]={90,92,94,95,96,97,98};
  static const G4double T66[n66]={.1,.1,.0000024,.0000068,.000043,.000049,0.};
  static const G4double X66[n66]={7000,7000.,50000.,50000.,7000.,50000.,300000.};
  static const G4int    n67=1;                   // Z=67 Ho
  static const G4int    N67[n67]={98};
  static const G4double T67[n67]={.01};
  static const G4double X67[n67]={13000.};
  static const G4int    n68=6;                   // Z=68 Er
  static const G4int    N68[n68]={94,96,98,99,100,102};
  static const G4double T68[n68]={0.,0.,.00201,.000134,0.,0.};
  static const G4double X68[n68]={50000.,50000.,40000.,50000.,50000.,50000.};
  static const G4int    n69=1;                   // Z=69 Tm
  static const G4int    N69[n69]={100};
  static const G4double T69[n69]={0.};
  static const G4double X69[n69]={10000.};
  static const G4int    n70=7;                   // Z=70 Yb
  static const G4int    N70[n70]={98,100,101,102,103,104,106};
  static const G4double T70[n70]={0.,0.,0.,0.,0.,0.,0.};
  static const G4double X70[n70]={7000,7000.,7000.,7000.,7000.,7000.,7000.};
  static const G4int    n71=2;                   // Z=71 Lu
  static const G4int    N71[n71]={104,105};
  static const G4double T71[n71]={0.,0.};
  static const G4double X71[n71]={5280.,3160.};
  static const G4int    n72=6;                   // Z=72 Hf
  static const G4int    N72[n72]={102,104,105,106,107,108};
  static const G4double T72[n72]={.0003,.00128,.001,.0024,.001,.011};
  static const G4double X72[n72]={40000.,20000.,13000.,20000.,20000.,12000.};
  static const G4int    n73=2;                   // Z=73 Ta
  static const G4int    N73[n73]={107,108};
  static const G4double T73[n73]={.1,.005};
  static const G4double X73[n73]={10000.,10000.};
  static const G4int    n74=5;                   // Z=74 W
  static const G4int    N74[n74]={106,108,109,110,112};
  static const G4double T74[n74]={.1,.1,.045,.1,.1};
  static const G4double X74[n74]={8000.,8600.,10000.,8800.,8800.};
  static const G4int    n75=2;                   // Z=75 Re
  static const G4int    N75[n75]={110,112};
  static const G4double T75[n75]={.1,.1};
  static const G4double X75[n75]={8500.,8500.};
  static const G4int    n76=7;                   // Z=76 Os
  static const G4int    N76[n76]={108,110,111,112,113,114,116};
  static const G4double T76[n76]={.000008,.000008,.000008,.000008,.000008,.000008,.000008};
  static const G4double X76[n76]={25000,25000.,25000.,25000.,25000.,25000.,25000.};
  static const G4int    n77=2;                   // Z=77 Ir
  static const G4int    N77[n77]={114,116};
  static const G4double T77[n77]={.1,.01};
  static const G4double X77[n77]={14000.,15000.};
  static const G4int    n78=6;                   // Z=78 Pt
  static const G4int    N78[n78]={112,114,116,117,118,120};
  static const G4double T78[n78]={0.,0.,0.,0.,0.,0.};
  static const G4double X78[n78]={12000.,12000.,12000.,12000.,12000.,12000.};
  static const G4int    n79=1;                   // Z=79 Au
  static const G4int    N79[n79]={118};
  static const G4double T79[n79]={.005};
  static const G4double X79[n79]={14000.};
  static const G4int    n80=7;                   // Z=80 Hg
  static const G4int    N80[n80]={116,118,119,120,121,122,124};
  static const G4double T80[n80]={0.,0.,0.,0.,0.,0.,0.};
  static const G4double X80[n80]={30000,30000.,30000.,30000.,30000.,30000.,30000.};
  static const G4int    n81=2;                   // Z=81 Tl
  static const G4int    N81[n81]={122,124};
  static const G4double T81[n81]={0.,0.};
  static const G4double X81[n81]={10000.,10000.};
  static const G4int    n82=4;                   // Z=82 Pb
  static const G4int    N82[n82]={122,124,125,126};
  static const G4double T82[n82]={.1,.1,.475,1.};
  static const G4double X82[n82]={12000.,5000.,6000.,4700.};
  static const G4int    n83=1;                   // Z=79 Tm
  static const G4int    N83[n83]={126};
  static const G4double T83[n83]={0.};
  static const G4double X83[n83]={9000.};
  static const G4int    n84=1;                   // Z=84 Po - No stable isotopes
  static const G4int    N84[n84]={-1};
  static const G4double T84[n84]={0.};
  static const G4double X84[n84]={0.};
  static const G4int    n85=1;                   // Z=85 At - No stable isotopes
  static const G4int    N85[n85]={-1};
  static const G4double T85[n85]={0.};
  static const G4double X85[n85]={0.};
  static const G4int    n86=1;                   // Z=86 Rn - No stable isotopes
  static const G4int    N86[n86]={-1};
  static const G4double T86[n86]={0.};
  static const G4double X86[n86]={0.};
  static const G4int    n87=1;                   // Z=87 Fr - No stable isotopes
  static const G4int    N87[n87]={-1};
  static const G4double T87[n87]={0.};
  static const G4double X87[n87]={0.};
  static const G4int    n88=1;                   // Z=88 Ra - No stable isotopes
  static const G4int    N88[n88]={-1};
  static const G4double T88[n88]={0.};
  static const G4double X88[n88]={0.};
  static const G4int    n89=1;                   // Z=89 Ac - No stable isotopes
  static const G4int    N89[n89]={-1};
  static const G4double T89[n89]={0.};
  static const G4double X89[n89]={0.};
  static const G4int    n90=1;                   // Z=90 Th
  static const G4int    N90[n90]={142};
  static const G4double T90[n90]={0.};
  static const G4double X90[n90]={11800.};
  static const G4int    n91=1;                   // Z=91 Pa - No stable isotopes
  static const G4int    N91[n91]={-1};
  static const G4double T91[n91]={0.};
  static const G4double X91[n91]={0.};
  static const G4int    n92=3;                   // Z=92 U
  static const G4int    N92[n92]={142,143,146};
  static const G4double T92[n92]={.1,0.,1.};
  static const G4double X92[n92]={10400.,15100.,14700.};
  static const G4int nn[nZ]={n01,n02,n03,n04,n05,n06,n07,n08,n09,n10,n11,n12,n13,n14,n15,
   n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,n32,n33,n34,n35,n36,n37,
   n38,n39,n40,n41,n42,n43,n44,n45,n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,
   n60,n61,n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,n78,n79,n80,n81,
			n82,n83,n84,n85,n86,n87,n88,n89,n90,n91,n92};
  static const G4int* nN[nZ]={N01,N02,N03,N04,N05,N06,N07,N08,N09,N10,N11,N12,N13,N14,N15,
   N16,N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,N36,N37,
   N38,N39,N40,N41,N42,N43,N44,N45,N46,N47,N48,N49,N50,N51,N52,N53,N54,N55,N56,N57,N58,N59,
   N60,N61,N62,N63,N64,N65,N66,N67,N68,N69,N70,N71,N72,N73,N74,N75,N76,N77,N78,N79,N80,N81,
   N82,N83,N84,N85,N86,N87,N88,N89,N90,N91,N92};
  static const G4double* nT[nZ]={T01,T02,T03,T04,T05,T06,T07,T08,T09,T10,T11,T12,T13,T14,
   T15,T16,T17,T18,T19,T20,T21,T22,T23,T24,T25,T26,T27,T28,T29,T30,T31,T32,T33,T34,T35,T36,
   T37,T38,T39,T40,T41,T42,T43,T44,T45,T46,T47,T48,T49,T50,T51,T52,T53,T54,T55,T56,T57,T58,
   T59,T60,T61,T62,T63,T64,T65,T66,T67,T68,T69,T70,T71,T72,T73,T74,T75,T76,T77,T78,T79,T80,
   T81,T82,T83,T84,T85,T86,T87,T88,T89,T90,T91,T92};
  static const G4double* nX[nZ]={X01,X02,X03,X04,X05,X06,X07,X08,X09,X10,X11,X12,X13,X14,
   X15,X16,X17,X18,X19,X20,X21,X22,X23,X24,X25,X26,X27,X28,X29,X30,X31,X32,X33,X34,X35,X36,
   X37,X38,X39,X40,X41,X42,X43,X44,X45,X46,X47,X48,X49,X50,X51,X52,X53,X54,X55,X56,X57,X58,
   X59,X60,X61,X62,X63,X64,X65,X66,X67,X68,X69,X70,X71,X72,X73,X74,X75,X76,X77,X78,X79,X80,
   X81,X82,X83,X84,X85,X86,X87,X88,X89,X90,X91,X92};
  if(tgZ<1 || tgZ>92)
  {
    G4cout<<"*Warning*G4QElasticCS::GetTabValue: (1-92) No isotopes for Z="<<tgZ<<G4endl;
    return 0.;
  }
  G4int iZ=tgZ-1; // Z index
  if(nN[iZ][0] < 0)
  {
#ifdef isodebug
    G4cout<<"*Warning*G4QElasticCS::GetTabValue: No isotopes for Z="<<tgZ<<G4endl;
#endif
    return 0.;
  }
#ifdef pdebug
  G4cout<<"G4QElasticCS::GetTabVal: lp="<<lp<<",Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<G4endl;
#endif
  G4double p=std::exp(lp);              // momentum
  G4double sp=std::sqrt(p);             // sqrt(p)
  G4double p2=p*p;            
  G4double p3=p2*p;
  G4double p4=p3*p;
  if(PDG==2112 && tgZ==1 && tgN==0)       // np
  {
    G4double ssp=std::sqrt(sp);           // sqrt(sqrt(p))=p^.25
    G4double p2s=p2*sp;
		  G4double dl1=lp-lastPAR[3];
    theSS=lastPAR[21];
    theS1=(lastPAR[9]+lastPAR[10]*dl1*dl1+lastPAR[11]/p)/(1.+lastPAR[12]/p4)
          +lastPAR[13]/(p4+lastPAR[14]);
    theB1=(lastPAR[17]+lastPAR[18]/(p4*p4+lastPAR[19]*p3))/(1.+lastPAR[20]/p4);
    theS2=(lastPAR[15]+lastPAR[16]/p4/p)/p3;
    theB2=lastPAR[22]/(p*sp+lastPAR[23]); 
    theS3=0.;
    theB3=0.; 
    theS4=0.;
    theB4=0.; 
#ifdef tdebug
    G4cout<<"G4QElasticCS::GetTableValues:(np) TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1
          <<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS1<<",B3="<<theB1<<G4endl;
#endif
    // Returns the total elastic pp cross-section (to avoid spoiling lastSIG)
    return lastPAR[0]/(p2s+lastPAR[1]*p+lastPAR[2]/ssp)+lastPAR[4]/p
           +(lastPAR[5]+lastPAR[6]*dl1*dl1+lastPAR[7]/p)/(1.+lastPAR[8]/p4);
  }
  else if(PDG==2212 && tgZ==1 && tgN==0) // pp
  {
    G4double p2s=p2*sp;
		  G4double dl1=lp-lastPAR[3];
		  G4double dl2=lp-lastPAR[8];
    theSS=lastPAR[31];
    theS1=(lastPAR[9]+lastPAR[10]*dl2*dl2)/(1.+lastPAR[11]/p4/p)+
          (lastPAR[12]/p2+lastPAR[13]*p)/(p4+lastPAR[14]*sp);
    theB1=lastPAR[15]*std::pow(p,lastPAR[16])/(1.+lastPAR[17]/p3);
    theS2=lastPAR[18]+lastPAR[19]/(p4+lastPAR[20]*p);
    theB2=lastPAR[21]+lastPAR[22]/(p4+lastPAR[23]/sp); 
    theS3=lastPAR[24]+lastPAR[25]/(p4*p4+lastPAR[26]*p2+lastPAR[27]);
    theB3=lastPAR[28]+lastPAR[29]/(p4+lastPAR[30]); 
    theS4=0.;
    theB4=0.; 
#ifdef tdebug
    G4cout<<"G4QElasticCS::GetTableValues:(pp) TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1
          <<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS1<<",B3="<<theB1<<G4endl;
#endif
    // Returns the total elastic pp cross-section (to avoid spoiling lastSIG)
    return lastPAR[0]/p2s/(1.+lastPAR[7]/p2s)+(lastPAR[1]+lastPAR[2]*dl1*dl1+lastPAR[4]/p)
                                                   /(1.+lastPAR[5]*lp)/(1.+lastPAR[6]/p4);
  }
  else if(PDG==2212 || PDG==2112) // n/p+A (isotope/projectile invariant)
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
      G4cout<<"G4QElCS::GetTabV: lA, p="<<p<<",S1="<<theS1<<",B1="<<theB1<<",SS="<<theSS
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
      G4cout<<"G4QElCS::GetTabV: hA, p="<<p<<",S1="<<theS1<<",B1="<<theB1<<",SS="<<theSS
            <<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS3<<",B3="<<theB3<<",S4="<<theS4
            <<",B4="<<theB4<<G4endl;
#endif
    }
    G4double rp16=lastPAR[6]/p16;
    // Returns the total elastic (n/p)A cross-section (to avoid spoiling lastSIG)
#ifdef tdebug
    G4cout<<"G4QElCS::GetTabV: PDG="<<PDG<<",P="<<p<<",N="<<tgN<<",Z="<<tgZ<<G4endl;
#endif
    if(PDG!=2112 || !tgN || p>.1)      // No Low Energy neutron addition
     return (lastPAR[0]*dl*dl+lastPAR[1])/(1.+lastPAR[2]/p2)+lastPAR[3]/(p4+lastPAR[4]/p2)+
            lastPAR[5]/(p5+rp16*rp16)+lastPAR[7]/(p4+lastPAR[8]/p4);
    else
    {
      G4int nI=nn[iZ];
      for(G4int i=0; i<nI; i++) if(tgN==nN[iZ][i])
      {
        G4double kE=p2/dNM;
#ifdef tdebug
        G4cout<<"G4QECS::GTV:P="<<p<<",E="<<kE<<",X="<<nX[iZ][i]<<",T="<<nT[iZ][i]<<G4endl;
#endif
        if(kE<nT[iZ][i]) return 0.;    // 0 below the threshold (in MeV)
        if(p<2.) return nX[iZ][i];     // At very low momentum(<2MeV/c) -> only LECS
    				return lastPAR[3]/(p4+lastPAR[4]/p2)+lastPAR[5]/(p5+rp16*rp16)+
               lastPAR[7]/(p4+lastPAR[8]/p4)+nX[iZ][i]/(1.+lastPAR[51]*p16);
      }
      G4cerr<<"*G4QElCS::GTV:PDG="<<PDG<<",Z="<<tgZ<<",N="<<tgN<<": Sig not found"<<G4endl;
    }
  }
  else
  {
    G4cout<<"*Error*G4QElasticCrossSection::GetTabValues: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for PDG=2212, Z=1, N=0"<<G4endl;
    throw G4QException("G4QElasticCrossSection::GetTabValues: only pp is implemented");
  }
  return 0.;
} // End of GetTableValues

// Returns max -t=Q2 (GeV^2) for the momentum pP(GeV) and the target nucleus (tgN,tgZ)
G4double G4QElasticCrossSection::GetQ2max(G4int PDG, G4int tgZ, G4int tgN, G4double pP)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass()*.001; // MeV to GeV
  static const G4double mProt= G4QPDGCode(2212).GetMass()*.001; // MeV to GeV
  //static const G4double mLamb= G4QPDGCode(3122).GetMass()*.001; // MeV to GeV
  //static const G4double mHe3 = G4QPDGCode(2112).GetNuclMass(2,1,0)*.001; // MeV to GeV
  //static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0)*.001; // MeV to GeV
  //static const G4double mDeut = G4QPDGCode(2112).GetNuclMass(1,1,0)*.001; // MeV to GeV
  static const G4double mProt2= mProt*mProt;
  static const G4double mNeut2= mNeut*mNeut;
  //static const G4double mDeut2= mDeut*mDeut;
  G4double pP2=pP*pP;                                 // squared momentum of the projectile
  if(PDG==2212 && tgZ==1 && tgN==0)                   // ---> pp(identical,symmetric 90deg)
  {
    G4double tMid=std::sqrt(pP2+mProt2)*mProt-mProt2; // CMS 90deg value of -t=Q2 (GeV^2)
    return tMid+tMid;
  }
  else if(PDG==2112 && tgZ)                           // ---> nA
  {
    G4double mt=mProt;                                // Target mass in GeV
    if(tgN||tgZ>1) mt=G4QPDGCode(90000000+tgZ*1000+tgN).GetMass()*.001; // Target mass GeV
    G4double dmt=mt+mt;
    G4double s=dmt*std::sqrt(pP2+mNeut2)+mNeut2+mt*mt; // Mondelstam s (in GeV^2)
    return dmt*dmt*pP2/s;
  }
  else if(PDG==2212 && tgZ)                           // ---> pA
  {
    G4double mt=G4QPDGCode(90000000+tgZ*1000+tgN).GetMass()*.001; // Target mass in GeV
    G4double dmt=mt+mt;
    G4double s=dmt*std::sqrt(pP2+mProt2)+mProt2+mt*mt; // Mondelstam s
    return dmt*dmt*pP2/s;
  }
  else
  {
    G4cout<<"*Error*G4QElasticCrossSection::GetQ2max: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for p,n projectiles & Z_target>0"<<G4endl;
    throw G4QException("G4QElasticCrossSection::GetQ2max: only nA & pA are implemented");
  }
}
