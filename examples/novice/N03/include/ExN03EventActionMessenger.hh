// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN03EventActionMessenger.hh,v 1.1 1999-01-07 16:05:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
};

#endif
