// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14EventActionMessenger.hh,v 1.1 1999-05-29 14:12:05 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst14EventActionMessenger_h
#define Tst14EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst14EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst14EventActionMessenger: public G4UImessenger
{
  public:
    Tst14EventActionMessenger(Tst14EventAction*);
   ~Tst14EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Tst14EventAction* eventAction;   
    G4UIcmdWithAnInteger* setVerboseCmd;
};

#endif
