// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTEventActionMessenger.hh,v 1.1 2003-11-07 21:30:28 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTEventActionMessenger_h
#define NTSTEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class NTSTEventAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class NTSTEventActionMessenger: public G4UImessenger
{
  public:
    NTSTEventActionMessenger(NTSTEventAction*);
   ~NTSTEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    NTSTEventAction* eventAction;   
    G4UIcmdWithAString* DrawCmd;
};

#endif
