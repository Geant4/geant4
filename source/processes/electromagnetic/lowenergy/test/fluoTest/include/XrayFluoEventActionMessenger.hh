//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoEventActionMessenger_h
#define XrayFluoEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class XrayFluoEventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoEventActionMessenger: public G4UImessenger
{
  public:
    XrayFluoEventActionMessenger(XrayFluoEventAction*);
   ~XrayFluoEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    XrayFluoEventAction*   eventAction;   
    G4UIcmdWithAString* DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;    
};

#endif
