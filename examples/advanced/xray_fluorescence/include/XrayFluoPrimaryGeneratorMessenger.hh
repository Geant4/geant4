//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoPrimaryGeneratorMessenger_h
#define XrayFluoPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class XrayFluoPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  XrayFluoPrimaryGeneratorMessenger(XrayFluoPrimaryGeneratorAction*);
  ~XrayFluoPrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  XrayFluoPrimaryGeneratorAction* Action; 
  G4UIcmdWithAString*          RndmCmd;
  G4UIcmdWithAString*          RndmCmmd;
  G4UIcmdWithADoubleAndUnit*     SigmAngleCmd;
  G4UIcmdWithADoubleAndUnit*  SigmaMomentumCmd;
};

#endif








