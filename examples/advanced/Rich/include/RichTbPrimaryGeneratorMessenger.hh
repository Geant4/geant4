// Rich advanced example for Geant4
// RichTbPrimaryGeneratorMessenger.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbPrimaryGeneratorMessenger_h
#define RichTbPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class RichTbPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class RichTbPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    RichTbPrimaryGeneratorMessenger(RichTbPrimaryGeneratorAction*);
   ~RichTbPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    RichTbPrimaryGeneratorAction* RichTbAction; 
    G4UIcmdWithAString*        RndmCmd;
    G4UIcmdWithADoubleAndUnit* setxvertexCmd;
    G4UIcmdWithADoubleAndUnit* setyvertexCmd;
    G4UIcmdWithADoubleAndUnit* setzvertexCmd;
};

#endif

