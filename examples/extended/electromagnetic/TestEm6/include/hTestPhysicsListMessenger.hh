// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// The messenger for physics list class.
// Class Description - end
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef hTestPhysicsListMessenger_h
#define hTestPhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class hTestPhysicsList;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestPhysicsListMessenger: public G4UImessenger
{
  public: // Without description
  
    hTestPhysicsListMessenger(hTestPhysicsList*);
   ~hTestPhysicsListMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
  
    hTestPhysicsList*            hTestList;

    G4UIcmdWithADoubleAndUnit* cutGCmd;
    G4UIcmdWithADoubleAndUnit* cutECmd;
    G4UIcmdWithADoubleAndUnit* cutPCmd;
    G4UIcmdWithADoubleAndUnit* rCmd;
    G4UIcmdWithADoubleAndUnit* eCmd;
    G4UIcmdWithADoubleAndUnit* setMaxStepCmd;     
};

#endif

