// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17EventActionMessenger.hh,v 1.2 1999-12-15 14:54:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst17EventActionMessenger_h
#define Tst17EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst17EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst17EventActionMessenger: public G4UImessenger
{
  public:
    Tst17EventActionMessenger(Tst17EventAction*);
   ~Tst17EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Tst17EventAction* eventAction;   
    G4UIcmdWithAnInteger* setVerboseCmd;
};

#endif
