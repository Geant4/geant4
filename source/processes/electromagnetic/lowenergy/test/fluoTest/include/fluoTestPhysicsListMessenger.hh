//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestPhysicsListMessenger_h
#define fluoTestPhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class fluoTestPhysicsList;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestPhysicsListMessenger: public G4UImessenger
{
  public:

  fluoTestPhysicsListMessenger(fluoTestPhysicsList*);
  ~fluoTestPhysicsListMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);

 private:

  fluoTestPhysicsList*         fluoTestList;

  G4UIdirectory* EnDir;
  G4UIcmdWithADoubleAndUnit* cutGLowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutELowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutGELowLimCmd;
};

#endif
