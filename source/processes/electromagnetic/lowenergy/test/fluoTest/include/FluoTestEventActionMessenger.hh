//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestEventActionMessenger_h
#define FluoTestEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class FluoTestEventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestEventActionMessenger: public G4UImessenger
{
  public:
    FluoTestEventActionMessenger(FluoTestEventAction*);
   ~FluoTestEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    FluoTestEventAction*   eventAction;   
    G4UIcmdWithAString* DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;    
};

#endif
