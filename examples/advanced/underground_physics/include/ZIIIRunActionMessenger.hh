// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		exrdmRunActionMessenger.hh
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

#ifndef ZIIIRunActionMessenger_h
#define ZIIIRunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ZIIIRunAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ZIIIRunActionMessenger: public G4UImessenger
{
  public:
    ZIIIRunActionMessenger(ZIIIRunAction*);
   ~ZIIIRunActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    ZIIIRunAction*   ZIIIRun;   
    G4UIcmdWithAString* FileCmd;
};

#endif






