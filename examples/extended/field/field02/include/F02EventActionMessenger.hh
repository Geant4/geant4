// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F02EventActionMessenger.hh,v 1.1 2001-03-27 16:26:18 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F02EventActionMessenger_h
#define F02EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class F02EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F02EventActionMessenger: public G4UImessenger
{
  public:
    F02EventActionMessenger(F02EventAction*);
   ~F02EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    F02EventAction* eventAction;   
    G4UIcmdWithAnInteger* setVerboseCmd;
    G4UIcmdWithAString*   DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;
};

#endif
