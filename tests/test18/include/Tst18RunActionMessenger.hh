// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		Tst18RunActionMessenger.hh
//
// Author:		F Lei
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------

#ifndef Tst18RunActionMessenger_h
#define Tst18RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst18RunAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst18RunActionMessenger: public G4UImessenger
{
  public:
    Tst18RunActionMessenger(Tst18RunAction*);
   ~Tst18RunActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Tst18RunAction*   Tst18Run;   
    G4UIcmdWithAString* FileCmd;
};

#endif






