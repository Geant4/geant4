// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2EventActionMessenger.hh,v 1.2 1999-12-15 14:48:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em2EventActionMessenger_h
#define Em2EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em2EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em2EventActionMessenger: public G4UImessenger
{
  public:
    Em2EventActionMessenger(Em2EventAction*);
   ~Em2EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em2EventAction*   eventAction;   
    G4UIcmdWithAString*   DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;
};

#endif
