// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F01PrimaryGeneratorMessenger.hh,v 1.1 2001-03-27 16:21:29 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F01PrimaryGeneratorMessenger_h
#define F01PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class F01PrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F01PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    F01PrimaryGeneratorMessenger(F01PrimaryGeneratorAction*);
   ~F01PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    F01PrimaryGeneratorAction* F01Action; 
    G4UIcmdWithAString*        RndmCmd;
    G4UIcmdWithADoubleAndUnit* setxvertexCmd;
    G4UIcmdWithADoubleAndUnit* setyvertexCmd;
    G4UIcmdWithADoubleAndUnit* setzvertexCmd;
};

#endif

