// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5SteppingMessenger.hh,v 1.2 1999-12-15 14:49:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em5SteppingMessenger_h
#define Em5SteppingMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em5SteppingAction;
class G4UIdirectory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em5SteppingMessenger: public G4UImessenger
{
  public:

   Em5SteppingMessenger(Em5SteppingAction* );
  ~Em5SteppingMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   Em5SteppingAction* steppingAction;

   G4UIdirectory*     steppingDir;

};

#endif

