// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F01EventActionMessenger.hh,v 1.1 2001-03-27 16:21:26 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F01EventActionMessenger_h
#define F01EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class F01EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F01EventActionMessenger: public G4UImessenger
{
  public:
    F01EventActionMessenger(F01EventAction*);
   ~F01EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    F01EventAction* eventAction;   
    G4UIcmdWithAnInteger* setVerboseCmd;
    G4UIcmdWithAString*   DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;
};

#endif
