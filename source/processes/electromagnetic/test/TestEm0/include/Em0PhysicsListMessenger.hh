// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0PhysicsListMessenger.hh,v 1.1 1999-01-08 16:32:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em0PhysicsListMessenger_h
#define Em0PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em0PhysicsList;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em0PhysicsListMessenger: public G4UImessenger
{
  public:
    Em0PhysicsListMessenger(Em0PhysicsList*);
   ~Em0PhysicsListMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em0PhysicsList* physList;
    G4UIcmdWithADoubleAndUnit* cutGCmd;
    G4UIcmdWithADoubleAndUnit* cutECmd;
    G4UIcmdWithADoubleAndUnit* rCmd;
    G4UIcmdWithADoubleAndUnit* eCmd;
};

#endif
