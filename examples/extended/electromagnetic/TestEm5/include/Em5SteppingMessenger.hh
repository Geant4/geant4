// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5SteppingMessenger.hh,v 1.1 1999-10-12 12:23:33 maire Exp $
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

