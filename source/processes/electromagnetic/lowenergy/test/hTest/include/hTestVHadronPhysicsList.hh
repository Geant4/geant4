#ifndef hTestVHadronPhysicsList_h
#define hTestVHadronPhysicsList_h 1

//---------------------------------------------------------------------------
//
// ClassName:   hTestVHadronPhysicsList
//  
// Description: Virtual class to build Hadron Physics List for Geant4
//
// Authors:     V.Ivanchenko 29/03/01
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VUserPhysicsList.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestVHadronPhysicsList:  public G4VUserPhysicsList
{
public:
  hTestVHadronPhysicsList() {};
  ~hTestVHadronPhysicsList() {};

public:
  void ConstructHad() {ConstructProcess();};
  void ConstructParticle() {};
  void SetCuts() {};

  void SetVerbose(G4int val) {verbose = val;};
  
private:

  // hide assignment operator 
  hTestVHadronPhysicsList & operator=(const hTestVHadronPhysicsList &right);
  hTestVHadronPhysicsList(const hTestVHadronPhysicsList&);

protected:

  G4int verbose;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


