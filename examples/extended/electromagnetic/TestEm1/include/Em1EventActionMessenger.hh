// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1EventActionMessenger.hh,v 1.1 1999-10-11 13:07:36 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em1EventActionMessenger_h
#define Em1EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em1EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em1EventActionMessenger: public G4UImessenger
{
  public:
    Em1EventActionMessenger(Em1EventAction*);
   ~Em1EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em1EventAction* eventAction;   
    G4UIcmdWithAString* DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;    
};

#endif
