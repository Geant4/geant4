// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3PrimaryGeneratorMessenger.hh,v 1.1 1999-10-11 16:55:49 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em3PrimaryGeneratorMessenger_h
#define Em3PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class Em3PrimaryGeneratorAction;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em3PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    Em3PrimaryGeneratorMessenger(Em3PrimaryGeneratorAction*);
   ~Em3PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em3PrimaryGeneratorAction* Em3Action; 
    G4UIcmdWithoutParameter*   DefaultCmd;
};

#endif

