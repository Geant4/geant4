// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20PhysicsListMessenger.hh,v 1.1 2001-05-25 12:50:06 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst20PhysicsListMessenger_h
#define Tst20PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst20PhysicsList;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst20PhysicsListMessenger: public G4UImessenger
{
  
public:

  Tst20PhysicsListMessenger(Tst20PhysicsList*);
  ~Tst20PhysicsListMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:

  Tst20PhysicsList*          Tst20List;

  G4UIdirectory* lowEnDir;
  G4UIcmdWithADoubleAndUnit* cutGLowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutELowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutGELowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutSecPhotCmd;
  G4UIcmdWithADoubleAndUnit* cutSecElecCmd;
  G4UIcmdWithADoubleAndUnit* cutGCmd;
  G4UIcmdWithADoubleAndUnit* cutECmd;
};

#endif








