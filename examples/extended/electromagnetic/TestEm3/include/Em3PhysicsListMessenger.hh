// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3PhysicsListMessenger.hh,v 1.1 2000-04-17 12:06:23 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em3PhysicsListMessenger_h
#define Em3PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em3PhysicsList;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em3PhysicsListMessenger: public G4UImessenger
{
  public:
  
    Em3PhysicsListMessenger(Em3PhysicsList* );
   ~Em3PhysicsListMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
  
    Em3PhysicsList* pPhysicsList;
    
    G4UIcmdWithADoubleAndUnit* gammaCutCmd;
    G4UIcmdWithADoubleAndUnit* electCutCmd;
    G4UIcmdWithADoubleAndUnit* protoCutCmd;    
};

#endif

