//  Em6EventActionMessenger.hh

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em6EventActionMessenger_h
#define Em6EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em6EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em6EventActionMessenger: public G4UImessenger
{
  public:
    Em6EventActionMessenger(Em6EventAction*);
   ~Em6EventActionMessenger();

    void SetNewValue(G4UIcommand*, G4String);

  private:
    Em6EventAction*   eventAction;
    G4UIcmdWithAString*   DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
