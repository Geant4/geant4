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

#ifndef Em6EventActionMessenger_h
#define Em6EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em6EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em6EventActionMessenger: public G4UImessenger
{
public: // Without  description
 
    Em6EventActionMessenger(Em6EventAction*);
   ~Em6EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em6EventAction* eventAction;   
    G4UIcmdWithAnInteger* setVerboseCmd;
    G4UIcmdWithAString*   DrawCmd;
};

#endif
