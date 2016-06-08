// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN03PrimaryGeneratorMessenger.hh,v 1.1.10.1 1999/12/07 20:47:29 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef ExN03PrimaryGeneratorMessenger_h
#define ExN03PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class ExN03PrimaryGeneratorAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ExN03PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    ExN03PrimaryGeneratorMessenger(ExN03PrimaryGeneratorAction*);
   ~ExN03PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    ExN03PrimaryGeneratorAction* ExN03Action; 
    G4UIcmdWithAString*          RndmCmd;
};

#endif

