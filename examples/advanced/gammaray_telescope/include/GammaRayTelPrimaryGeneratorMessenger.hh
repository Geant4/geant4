//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelPrimaryGeneratorMessenger_h
#define GammaRayTelPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class GammaRayTelPrimaryGeneratorAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    GammaRayTelPrimaryGeneratorMessenger(GammaRayTelPrimaryGeneratorAction*);
   ~GammaRayTelPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    GammaRayTelPrimaryGeneratorAction* GammaRayTelAction; 
    G4UIcmdWithAString*          RndmCmd;
};

#endif

