// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14PhysicsListMessenger.hh,v 1.1 1999-05-29 14:12:05 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst14PhysicsListMessenger_h
#define Tst14PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst14PhysicsList;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst14PhysicsListMessenger: public G4UImessenger
{
  public:
    Tst14PhysicsListMessenger(Tst14PhysicsList*);
   ~Tst14PhysicsListMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Tst14PhysicsList*          Tst14List;

    G4UIcmdWithADoubleAndUnit* setMaxStepCmd;

    G4UIcmdWithADoubleAndUnit* cutGCmd;
    G4UIcmdWithADoubleAndUnit* cutECmd;
    G4UIcmdWithADoubleAndUnit* cutPCmd;
    G4UIcmdWithADoubleAndUnit* rCmd;
    G4UIcmdWithADoubleAndUnit* eCmd;
 
};

#endif

