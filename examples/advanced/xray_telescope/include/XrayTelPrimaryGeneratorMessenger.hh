//  XrayTelPrimaryGeneratorMessenger.hh

#ifndef XrayTelPrimaryGeneratorMessenger_h
#define XrayTelPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class XrayTelPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayTelPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  XrayTelPrimaryGeneratorMessenger(XrayTelPrimaryGeneratorAction*);
  ~XrayTelPrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  XrayTelPrimaryGeneratorAction*    XrayTelAction; 
  G4UIcmdWithAString*          RndmCmd;
  G4UIcmdWithAString*          ErndmCmd;
  G4UIcmdWithADoubleAndUnit*   SetRmin;
  G4UIcmdWithADoubleAndUnit*   SetRmax;
  G4UIcmdWithADoubleAndUnit*   SetTmax;
};

#endif

