#ifndef HsVEMPhysicsList_h
#define HsVEMPhysicsList_h 1

//---------------------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//---------------------------------------------------------------------------
//
// ClassName:   HsVEMPhysicsList
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

class HsVEMPhysicsList:  public G4VUserPhysicsList
{
public:
  HsVEMPhysicsList() {};
  ~HsVEMPhysicsList() {};

public:
  void ConstructEM() {ConstructProcess();};
  void ConstructParticle() {};
  void SetCuts() {};

  inline void SetVerbose(G4int val) {verbose = val;};
  inline void SetMaxChargedStep(G4double val) {maxChargedStep = val;};
  inline void SetNuclearStopping(G4bool val) {nuclStop = val;};
  inline void SetBarkas(G4bool val) {barkas = val;};
  inline void SetEStoppingTable(G4String name) {table = name;};
  
private:

  // hide assignment operator 
  HsVEMPhysicsList & operator=(const HsVEMPhysicsList &right);
  HsVEMPhysicsList(const HsVEMPhysicsList&);

protected:

  G4int verbose;
  G4double maxChargedStep;
  G4bool nuclStop;
  G4bool barkas;
  G4String table;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


