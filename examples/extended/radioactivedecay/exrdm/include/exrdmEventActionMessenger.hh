
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 

#ifndef exrdmEventActionMessenger_h
#define exrdmEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class exrdmEventAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exrdmEventActionMessenger: public G4UImessenger
{
  public:
    exrdmEventActionMessenger(exrdmEventAction*);
   ~exrdmEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    exrdmEventAction*   eventAction;   
    G4UIcmdWithAString* DrawCmd;
};

#endif






