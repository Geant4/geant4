//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoPhysicsListMessenger_h
#define  XraFluocPhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class  XrayFluoPhysicsList;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class  XrayFluoPhysicsListMessenger: public G4UImessenger
{
  public:

   XrayFluoPhysicsListMessenger(XrayFluoPhysicsList*);
  ~ XrayFluoPhysicsListMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);

 private:

   XrayFluoPhysicsList*        XrayFluoList;

  G4UIdirectory* EnDir;
  G4UIcmdWithADoubleAndUnit* cutGLowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutELowLimCmd;
  G4UIcmdWithADoubleAndUnit* cutGELowLimCmd;
};

#endif
