//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef myPhysicsListMessenger_h
#define myPhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class myPhysicsList;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class myPhysicsListMessenger: public G4UImessenger
{
  public:

  myPhysicsListMessenger(myPhysicsList*);
  ~myPhysicsListMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);

 private:

  myPhysicsList*         myList;

  G4UIdirectory* EnDir;
  G4UIcmdWithADoubleAndUnit* cutGLowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutELowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutGELowLimCmd;
};

#endif
