// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01PrimaryGeneratorMessenger.hh,v 1.1.1.1 2000-09-14 11:37:21 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef AnaEx01PrimaryGeneratorMessenger_h
#define AnaEx01PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class AnaEx01PrimaryGeneratorAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class AnaEx01PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    AnaEx01PrimaryGeneratorMessenger(AnaEx01PrimaryGeneratorAction*);
   ~AnaEx01PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    AnaEx01PrimaryGeneratorAction* AnaEx01Action; 
    G4UIcmdWithAString*          RndmCmd;
};

#endif

