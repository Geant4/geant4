// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// Event action messenger.
// Class Description - end
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef hTestEventActionMessenger_h
#define hTestEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class hTestEventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestEventActionMessenger: public G4UImessenger
{
public: // Without  description
 
    hTestEventActionMessenger(hTestEventAction*);
   ~hTestEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    hTestEventAction* eventAction;   
    G4UIcmdWithAnInteger* setVerboseCmd;
    G4UIcmdWithAString*   DrawCmd;
};

#endif
