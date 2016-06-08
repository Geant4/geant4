// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN03EventActionMessenger.hh,v 1.2.8.1 1999/12/07 20:47:28 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef ExN03EventActionMessenger_h
#define ExN03EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ExN03EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ExN03EventActionMessenger: public G4UImessenger
{
  public:
    ExN03EventActionMessenger(ExN03EventAction*);
   ~ExN03EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    ExN03EventAction*   eventAction;   
    G4UIcmdWithAString* DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;    
};

#endif
