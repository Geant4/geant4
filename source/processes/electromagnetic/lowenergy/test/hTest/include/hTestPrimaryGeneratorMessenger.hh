// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: hTestPrimaryGeneratorMessenger.hh,v 1.1 2000-05-21 18:37:45 chauvie Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef hTestPrimaryGeneratorMessenger_h
#define hTestPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class hTestPrimaryGeneratorAction;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    hTestPrimaryGeneratorMessenger(hTestPrimaryGeneratorAction*);
   ~hTestPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    hTestPrimaryGeneratorAction* hTestAction; 
    G4UIcmdWithoutParameter*   DefaultCmd;
};

#endif






