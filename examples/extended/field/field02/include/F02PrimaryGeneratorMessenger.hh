// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F02PrimaryGeneratorMessenger.hh,v 1.1 2001-03-27 16:26:20 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F02PrimaryGeneratorMessenger_h
#define F02PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class F02PrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F02PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    F02PrimaryGeneratorMessenger(F02PrimaryGeneratorAction*);
   ~F02PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    F02PrimaryGeneratorAction* F02Action; 
    G4UIcmdWithAString*        RndmCmd;
    G4UIcmdWithADoubleAndUnit* setxvertexCmd;
    G4UIcmdWithADoubleAndUnit* setyvertexCmd;
    G4UIcmdWithADoubleAndUnit* setzvertexCmd;
};

#endif

