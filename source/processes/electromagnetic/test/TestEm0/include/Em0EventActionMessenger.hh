// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0EventActionMessenger.hh,v 1.2 1999-12-15 14:51:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em0EventActionMessenger_h
#define Em0EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em0EventAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em0EventActionMessenger: public G4UImessenger
{
  public:
    Em0EventActionMessenger(Em0EventAction*);
   ~Em0EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em0EventAction* eventAction;   
    G4UIcmdWithAString* DrawCmd;
};

#endif
