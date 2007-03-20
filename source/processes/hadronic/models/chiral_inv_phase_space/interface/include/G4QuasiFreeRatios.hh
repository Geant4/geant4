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
// GEANT4 physics class: G4QuasiFreeRatios -- header file
// M.V. Kossov, ITEP(Moscow), 24-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Oct-2006
//

#ifndef G4QuasiFreeRatios_h
#define G4QuasiFreeRatios_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4QPDGCode.hh"
#include "G4QException.hh"
#include <vector>
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4VQCrossSection.hh"
#include "G4QElasticCrossSection.hh"


class G4QuasiFreeRatios
{
protected:

  G4QuasiFreeRatios()  {}               // Constructor

public:

  ~G4QuasiFreeRatios() {}

  static G4QuasiFreeRatios* GetPointer(); // Gives a pointer to this singletone

  std::pair<G4double,G4double> GetRatios(G4double pMom,G4int tgZ,G4int tgN,G4int pPDG=0);

  std::pair<G4double,G4double> CalculateRatios(G4int F, G4int I, G4int pPDG,
                                               G4int Z, G4int N, G4double pP);
  std::pair<G4double,G4double> GetElTot(G4double p, G4int PDG, G4int Z, G4int N, G4bool F);
private:
  G4double GetPTables(G4double lpP, G4double lPm, G4int PDG, G4int tZ, G4int tN); // newLP
  std::pair<G4double,G4double> GetTabValues(G4double lp, G4int pPDG, G4int tgZ, G4int tgN);
  G4double QF2IN_Ratio(G4double TCSmb, G4double A);

// Body
private:
  // --- Data formating AMDB (define the precalculated table structure) ---
  static const G4int nPoints;// #of points in the AMDB tables     
  static const G4int nLast;  // the Last element in the table
  static G4double    lPMin;  // Min tabulated logarithmic Momentum  
  static G4double    lPMax;  // Max tabulated logarithmic Momentum  
  static G4double    dlnP;   // Log step in the table     
  // ---- Local (for particular pP, pPDG, tZ, tN) -----
  static std::pair<G4double,G4double> lastRAT; // Last Ratios (QE/QF,QF/IN) Table (Calc)
  static G4double  lastLP;   // Last log(mom_of_the_incident_hadron in GeV)
  static G4int     lastN;    // The last N of calculated nucleus
  static G4int     lastZ;    // The last Z of calculated nucleus
  static G4double  lastP;    // Last used for the ratio Momentum
  static std::pair<G4double,G4double> lastQIR;  // Last Ratios (QE/QF,QF/IN) Values (Get)
  static G4int     lastI;    // The last position in the DAMDB
  static G4double  theQFP;   // The Last p-dependent QF/IN function
  static G4double  theETR;   // The ast p-dependent El/Tot Ratio function
  // ---- Global (AMBD of P-dependent tables for pPDG,tZ,tN) -----
  static G4int     lastPDG;  // Last PDG code of the projectile
  static G4int     lastTZ;   // Last atomic number of the target
  static G4int     lastTN;   // Last number of neutrons of the target
  static G4double  lastPIN;  // Last initialized max momentum
  static std::pair<G4double,G4double>* lastRST;  // Last ratio table
  static G4double* lastPAR;  // Last parameters for functional calculation
  static const G4double tolerance;// relative tolerance in momentum to get old Ratio
 }; 					
#endif
