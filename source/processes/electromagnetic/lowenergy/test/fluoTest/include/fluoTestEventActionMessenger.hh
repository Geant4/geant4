//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestEventActionMessenger_h
#define fluoTestEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class fluoTestEventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestEventActionMessenger: public G4UImessenger
{
  public:
    fluoTestEventActionMessenger(fluoTestEventAction*);
   ~fluoTestEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    fluoTestEventAction*   eventAction;   
    G4UIcmdWithAString* DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;    
};

#endif
