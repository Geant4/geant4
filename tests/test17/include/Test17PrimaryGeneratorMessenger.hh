// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Test17PrimaryGeneratorMessenger.hh,v 1.1 2000-05-26 06:36:30 chauvie Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test17PrimaryGeneratorMessenger_h
#define Test17PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class Test17PrimaryGeneratorAction;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    Test17PrimaryGeneratorMessenger(Test17PrimaryGeneratorAction*);
   ~Test17PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Test17PrimaryGeneratorAction* Test17Action; 
    G4UIcmdWithoutParameter*   DefaultCmd;
};

#endif






