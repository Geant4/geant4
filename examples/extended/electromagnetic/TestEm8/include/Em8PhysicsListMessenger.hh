// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8PhysicsListMessenger.hh,v 1.2 2000-06-27 13:29:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em8PhysicsListMessenger_h
#define Em8PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em8PhysicsList;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em8PhysicsListMessenger: public G4UImessenger
{
  public:
    Em8PhysicsListMessenger(Em8PhysicsList*);
   ~Em8PhysicsListMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em8PhysicsList*          Em8List;

    G4UIcmdWithADoubleAndUnit* setMaxStepCmd;

    G4UIcmdWithADoubleAndUnit* cutGCmd;
    G4UIcmdWithADoubleAndUnit* cutECmd;
    G4UIcmdWithADoubleAndUnit* cutPCmd;
    G4UIcmdWithADoubleAndUnit* rCmd;
    G4UIcmdWithADoubleAndUnit* eCmd;
 
};

#endif

