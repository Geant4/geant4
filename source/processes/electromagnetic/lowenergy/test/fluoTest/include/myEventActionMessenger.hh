//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef myEventActionMessenger_h
#define myEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class myEventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class myEventActionMessenger: public G4UImessenger
{
  public:
    myEventActionMessenger(myEventAction*);
   ~myEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    myEventAction*   eventAction;   
    G4UIcmdWithAString* DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;    
};

#endif
