// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: exGPSEventActionMessenger.hh,v 1.2 2004-10-28 11:21:37 flei Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef exGPSEventActionMessenger_h
#define exGPSEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class exGPSEventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exGPSEventActionMessenger: public G4UImessenger
{
  public:
    exGPSEventActionMessenger(exGPSEventAction*);
   ~exGPSEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    exGPSEventAction*   eventAction;   
    G4UIcmdWithAString* DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;    
};

#endif
