//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef myPrimaryGeneratorMessenger_h
#define myPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class myPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class myPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    myPrimaryGeneratorMessenger(myPrimaryGeneratorAction*);
   ~myPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
   
  private:
    myPrimaryGeneratorAction* Action; 
    G4UIcmdWithAString*          RndmCmd;
    G4UIcmdWithAString*          RndmCmmd;
    G4UIcmdWithADoubleAndUnit*     SigmAngleCmd;
    G4UIcmdWithADoubleAndUnit*  SigmaMomentumCmd;
};

#endif








