// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UnitsMessenger.hh,v 1.1 1999-01-07 16:09:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4UnitsMessenger_h
#define G4UnitsMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4UnitsMessenger: public G4UImessenger
{
  public:
    G4UnitsMessenger();
   ~G4UnitsMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:    
    G4UIdirectory*             UnitsTableDir;
    G4UIcmdWithoutParameter*   ListCmd;
};

#endif

