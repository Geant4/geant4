// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FluoTestPrimaryGeneratorMessenger.hh,v 1.7 2001-10-31 12:33:43 elena Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestPrimaryGeneratorMessenger_h
#define FluoTestPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class FluoTestPrimaryGeneratorAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    FluoTestPrimaryGeneratorMessenger(FluoTestPrimaryGeneratorAction*);
   ~FluoTestPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    FluoTestPrimaryGeneratorAction* FluoTestAction; 
    G4UIcmdWithAString*          RndmCmd;
 G4UIcmdWithAString*          RndmPart;
  G4UIcmdWithAString*          RndmVert;

};

#endif

