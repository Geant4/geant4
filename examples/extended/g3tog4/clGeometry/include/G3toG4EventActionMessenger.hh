// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4EventActionMessenger.hh,v 1.1 2000-07-24 11:23:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G3toG4EventActionMessenger_h
#define G3toG4EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G3toG4EventAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G3toG4EventActionMessenger: public G4UImessenger
{
  public:
    G3toG4EventActionMessenger(G3toG4EventAction*);
   ~G3toG4EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    G3toG4EventAction*   eventAction;   
    G4UIcmdWithAString* DrawCmd;
};

#endif
