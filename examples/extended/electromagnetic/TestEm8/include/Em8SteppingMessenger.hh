// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8SteppingMessenger.hh,v 1.1 2000-01-07 14:50:23 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em8SteppingMessenger_h
#define Em8SteppingMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em8SteppingAction;
class G4UIdirectory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em8SteppingMessenger: public G4UImessenger
{
  public:

   Em8SteppingMessenger(Em8SteppingAction* );
  ~Em8SteppingMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   Em8SteppingAction* steppingAction;

   G4UIdirectory*     steppingDir;

};

#endif

