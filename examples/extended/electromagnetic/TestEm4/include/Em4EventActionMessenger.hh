// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em4EventActionMessenger.hh,v 1.2 1999-12-15 14:49:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em4EventActionMessenger_h
#define Em4EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em4EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em4EventActionMessenger: public G4UImessenger
{
  public:
    Em4EventActionMessenger(Em4EventAction*);
   ~Em4EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em4EventAction* eventAction;   
    G4UIcmdWithAString*   DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;
};

#endif
