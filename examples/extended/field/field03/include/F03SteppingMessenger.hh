// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F03SteppingMessenger.hh,v 1.1 2001-06-08 11:55:40 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F03SteppingMessenger_h
#define F03SteppingMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F03SteppingAction;
class G4UIdirectory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F03SteppingMessenger: public G4UImessenger
{
  public:

   F03SteppingMessenger(F03SteppingAction* );
  ~F03SteppingMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   F03SteppingAction* steppingAction;

   G4UIdirectory*     steppingDir;

};

#endif

