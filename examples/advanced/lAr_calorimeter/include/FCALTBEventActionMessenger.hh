// First version of  FCALTBEventActionMessenger.
// /14/11/02 P.Mendez


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef FCALTBEventActionMessenger_h
#define FCALTBEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class FCALTBEventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class FCALTBEventActionMessenger: public G4UImessenger
{
  public:
    FCALTBEventActionMessenger(FCALTBEventAction*);
   ~FCALTBEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    FCALTBEventAction* fcaltbeventAction;   
    G4UIcmdWithAString* DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


