// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// The stepping messenger is defined
// Class Description - end
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test17SteppingMessenger_h
#define Test17SteppingMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17SteppingAction;
class G4UIdirectory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17SteppingMessenger: public G4UImessenger
{
public: // Without description

   Test17SteppingMessenger(Test17SteppingAction* );
  ~Test17SteppingMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   Test17SteppingAction* steppingAction;

   G4UIdirectory*     steppingDir;

};

#endif

