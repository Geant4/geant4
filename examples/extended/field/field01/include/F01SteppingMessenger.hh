// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F01SteppingMessenger.hh,v 1.1 2001-03-27 16:21:30 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F01SteppingMessenger_h
#define F01SteppingMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F01SteppingAction;
class G4UIdirectory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F01SteppingMessenger: public G4UImessenger
{
  public:

   F01SteppingMessenger(F01SteppingAction* );
  ~F01SteppingMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   F01SteppingAction* steppingAction;

   G4UIdirectory*     steppingDir;

};

#endif

