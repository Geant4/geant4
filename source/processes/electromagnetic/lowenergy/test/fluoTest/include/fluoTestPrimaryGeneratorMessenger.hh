//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestPrimaryGeneratorMessenger_h
#define fluoTestPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class fluoTestPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    fluoTestPrimaryGeneratorMessenger(fluoTestPrimaryGeneratorAction*);
   ~fluoTestPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
   
  private:
    fluoTestPrimaryGeneratorAction* Action; 
    G4UIcmdWithAString*          RndmCmd;
    G4UIcmdWithAString*          RndmCmmd;
    G4UIcmdWithADoubleAndUnit*     SigmAngleCmd;
    G4UIcmdWithADoubleAndUnit*  SigmaMomentumCmd;
};

#endif








