//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestPhysicsListMessenger_h
#define FluoTestPhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class FluoTestPhysicsList;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestPhysicsListMessenger: public G4UImessenger
{
  public:

  FluoTestPhysicsListMessenger(FluoTestPhysicsList*);
  ~FluoTestPhysicsListMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);

 private:

  FluoTestPhysicsList*         FluoTestList;

  G4UIdirectory* EnDir;
  G4UIcmdWithADoubleAndUnit* cutGLowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutELowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutGELowLimCmd;
};

#endif
