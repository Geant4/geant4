//  XrayTelEventActionMessenger.hh

#ifndef XrayTelEventActionMessenger_h
#define XrayTelEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class XrayTelEventAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayTelEventActionMessenger: public G4UImessenger
{
  public:
    XrayTelEventActionMessenger(XrayTelEventAction*);
   ~XrayTelEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    XrayTelEventAction*   eventAction;   
    G4UIcmdWithAString* DrawCmd;
};

#endif






