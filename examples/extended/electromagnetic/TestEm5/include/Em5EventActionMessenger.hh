// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5EventActionMessenger.hh,v 1.1 1999-10-12 12:23:30 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em5EventActionMessenger_h
#define Em5EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em5EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em5EventActionMessenger: public G4UImessenger
{
  public:
    Em5EventActionMessenger(Em5EventAction*);
   ~Em5EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em5EventAction* eventAction;   
    G4UIcmdWithAnInteger* setVerboseCmd;
    G4UIcmdWithAString*   DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;
};

#endif
