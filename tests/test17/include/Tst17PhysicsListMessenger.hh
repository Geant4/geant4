// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17PhysicsListMessenger.hh,v 1.1 1999-11-30 18:01:51 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst17PhysicsListMessenger_h
#define Tst17PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst17PhysicsList;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst17PhysicsListMessenger: public G4UImessenger
{
  public:
    Tst17PhysicsListMessenger(Tst17PhysicsList*);
   ~Tst17PhysicsListMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:

  Tst17PhysicsList*          Tst17List;
  G4UIdirectory* lowEnDir;
  G4UIcmdWithADoubleAndUnit* cutGLowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutELowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutGELowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutSecPhotCmd;
  G4UIcmdWithADoubleAndUnit* cutSecElecCmd;


  G4UIcmdWithADoubleAndUnit* setMaxStepCmd;
  
  G4UIcmdWithADoubleAndUnit* cutGCmd;
  G4UIcmdWithADoubleAndUnit* cutECmd;
  G4UIcmdWithADoubleAndUnit* cutPCmd;
  G4UIcmdWithADoubleAndUnit* rCmd;
  G4UIcmdWithADoubleAndUnit* eCmd;
 

};

#endif

