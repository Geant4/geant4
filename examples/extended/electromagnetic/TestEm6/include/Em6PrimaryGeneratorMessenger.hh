// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em6PrimaryGeneratorMessenger.hh,v 1.1 2000-03-01 07:50:10 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em6PrimaryGeneratorMessenger_h
#define Em6PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class Em6PrimaryGeneratorAction;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em6PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    Em6PrimaryGeneratorMessenger(Em6PrimaryGeneratorAction*);
   ~Em6PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em6PrimaryGeneratorAction* Em6Action; 
    G4UIcmdWithoutParameter*   DefaultCmd;
};

#endif






