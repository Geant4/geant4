#ifndef hTestPhysicsListMessenger_h
#define hTestPhysicsListMessenger_h 1

// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestPhysicsListMessenger -------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4UImessenger.hh"

class hTestPhysicsList;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

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
    G4UIcmdWithADoubleAndUnit* eCmd;
    G4UIcmdWithADoubleAndUnit* lowLimCmd;
    G4UIcmdWithADoubleAndUnit* setMaxStepCmd;     
    G4UIcmdWithAString*        EMPhysicsCmd;
    G4UIcmdWithAString*        HadPhysicsCmd;

};

#endif

