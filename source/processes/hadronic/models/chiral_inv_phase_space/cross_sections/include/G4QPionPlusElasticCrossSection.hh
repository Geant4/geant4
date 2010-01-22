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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4QPionPlusElasticCrossSection -- header file
// M.V. Kossov, ITEP(Moscow), 21-Jan-10
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 21-Jan-2010
//
//================================================================================
// Short description: Interaction cross-sections for the G4QElastic process
// -------------------------------------------------------------------------------

#ifndef G4QPionPlusElasticCrossSection_h
#define G4QPionPlusElasticCrossSection_h 1

#include "G4QPDGCode.hh"
#include "G4QException.hh"
#include <vector>
#include "Randomize.hh"
#include "G4VQCrossSection.hh"

class G4QPionPlusElasticCrossSection : public G4VQCrossSection
{
protected:

  G4QPionPlusElasticCrossSection();               // Constructor

public:

  ~G4QPionPlusElasticCrossSection();

  static G4VQCrossSection* GetPointer(); // Gives a pointer to this singletone

  // Cross-section is mb. At present momentum (pMom) is in MeV=IU (@@ make Indep. Units)
  virtual G4double GetCrossSection(G4bool fCS, G4double pMom, G4int tgZ, G4int tgN,
                                                                          G4int pPDG=2212);

  G4double CalculateCrossSection(G4bool CS, G4int F, G4int I, G4int pPDG, G4int Z, G4int N,
                                                                              G4double pP);
  G4double GetSlope(G4int tZ, G4int tN, G4int pPDG);     // Slope of the 1st diff. maximum
  G4double GetExchangeT(G4int tZ, G4int tN, G4int pPDG); // Randomizes -t=Q2 (in IU=MeV^2)
  G4double GetHMaxT();                   // Currrent Max(-t=Q2)/2. (in IU=MeV^2)

private:
  G4double GetPTables(G4double lpP, G4double lPm, G4int PDG, G4int tZ, G4int tN); // newLP
  G4double GetTabValues(G4double lp, G4int pPDG, G4int tgZ, G4int tgN); // return CS(Si/Bi)
  G4double GetQ2max(G4int pPDG, G4int tgZ, G4int tgN, G4double pP); // return -t=Q2

// Body
private:
  // --- Data formating AMDB (define the precalculated table structure) ---
  static const G4int nPoints;// #of points in the AMDB tables     
  static const G4int nLast;  // the Last element in the table
  static G4double    lPMin;  // Min tabulated logarithmic Momentum  
  static G4double    lPMax;  // Max tabulated logarithmic Momentum  
  static G4double    dlnP;   // Log step in the table     
  // ---- Local (for particular pP, pPDG, tZ, tN) -----
  static G4bool    onlyCS;   // flag to calculate only CS (not S1/B1,S2/B2,S3/B3)
  static G4double  lastSIG;  // Last calculated cross section
  static G4double  lastLP;   // Last log(mom_of_the_incident_hadron in GeV)
  static G4double  lastTM;   // Last t_maximum                       
  static G4int     lastN;    // The last N of calculated nucleus
  static G4int     lastZ;    // The last Z of calculated nucleus
  static G4double  lastP;    // Last used in the cross section Momentum
  static G4double  lastTH;   // Last value of the Momentum Threshold
  static G4double  lastCS;   // Last value of the Cross Section
  static G4int     lastI;    // The last position in the DAMDB
  static G4double  theSS;    // The Last squared slope of first diffruction 
  static G4double  theS1;    // The Last mantissa of first diffruction 
  static G4double  theB1;    // The Last slope of first diffruction    
  static G4double  theS2;    // The Last mantissa of second diffruction
  static G4double  theB2;    // The Last slope of second diffruction   
  static G4double  theS3;    // The Last mantissa of third diffruction 
  static G4double  theB3;    // The Last slope of third diffruction    
  static G4double  theS4;    // The Last mantissa of 4-th diffruction 
  static G4double  theB4;    // The Last slope of 4-th diffruction    
  // ---- Global (AMBD of P-dependent tables for pPDG,tZ,tN) -----
  static G4int     lastTZ;   // Last atomic number of the target
  static G4int     lastTN;   // Last number of neutrons of the target
  static G4double  lastPIN;  // Last initialized max momentum
  static G4double* lastCST;  // Last cross-section table
  static G4double* lastPAR;  // Last parameters for functional calculation
  static G4double* lastSST;  // E-dep of squared slope of the first difruction 
  static G4double* lastS1T;  // E-dep of mantissa of the first difruction 
  static G4double* lastB1T;  // E-dep of the slope of the first difruction
  static G4double* lastS2T;  // E-dep of mantissa of the second difruction
  static G4double* lastB2T;  // E-dep of the slope of theSecond difruction
  static G4double* lastS3T;  // E-dep of mantissa of the third difruction 
  static G4double* lastB3T;  // E-dep of the slope of the third difruction
  static G4double* lastS4T;  // E-dep of mantissa of the 4-th difruction 
  static G4double* lastB4T;  // E-dep of the slope of the 4-th difruction

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
 };
#endif
