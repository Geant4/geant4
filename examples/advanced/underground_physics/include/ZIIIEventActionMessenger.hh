// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ZIIIEventActionMessenger.hh,v 1.1 2001-06-26 11:23:17 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef ZIIIEventActionMessenger_h
#define ZIIIEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ZIIIEventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ZIIIEventActionMessenger: public G4UImessenger
{
  public:
    ZIIIEventActionMessenger(ZIIIEventAction*);
   ~ZIIIEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    ZIIIEventAction*   eventAction;   
    G4UIcmdWithAString* DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;    
};

#endif
