//  Em6PhysicsListMessenger.hh

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em6PhysicsListMessenger_h
#define Em6PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithADouble.hh"

class Em6PhysicsList;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em6PhysicsListMessenger: public G4UImessenger
{
  public:

    Em6PhysicsListMessenger(Em6PhysicsList* );
   ~Em6PhysicsListMessenger();

    void SetNewValue(G4UIcommand*, G4String);

  private:

    Em6PhysicsList* pPhysicsList;

    G4UIcmdWithADoubleAndUnit* gammaCutCmd;
    G4UIcmdWithADoubleAndUnit* electCutCmd;
    G4UIcmdWithADoubleAndUnit* protoCutCmd;
    G4UIcmdWithADouble* GammaToMuPairFac;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

