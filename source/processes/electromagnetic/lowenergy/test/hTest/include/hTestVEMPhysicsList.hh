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
#ifndef hTestVEMPhysicsList_h
#define hTestVEMPhysicsList_h 1

//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
//
// ClassName:   hTestVEMPhysicsList
//  
// Description: Virtual class to build Hadron Physics List for Geant4
//
// Authors:    08.04.01 V.Ivanchenko 
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VUserPhysicsList.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestVEMPhysicsList:  public G4VUserPhysicsList
{
public:
  hTestVEMPhysicsList() {};
  ~hTestVEMPhysicsList() {};

public:
  void ConstructEM() {ConstructProcess();};
  void ConstructParticle() {};
  void SetCuts() {};

  inline void SetVerbose(G4int val) {verbose = val;};
  inline void SetNuclearStopping(G4bool val) {nuclStop = val;};
  inline void SetBarkas(G4bool val) {barkas = val;};
  inline void SetEStoppingTable(G4String name) {table = name;};
  
private:

  // hide assignment operator 
  hTestVEMPhysicsList & operator=(const hTestVEMPhysicsList &right);
  hTestVEMPhysicsList(const hTestVEMPhysicsList&);

protected:

  G4int verbose;
  G4bool nuclStop;
  G4bool barkas;
  G4String table;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


