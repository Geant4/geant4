// Em6PrimaryGeneratorAction.hh

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em6PrimaryGeneratorAction_h
#define Em6PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4Event;
class Em6DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em6PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Em6PrimaryGeneratorAction(Em6DetectorConstruction*);
   ~Em6PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
    G4ParticleGun* GetParticleGun() {return particleGun;};

  private:
    G4ParticleGun*              particleGun;
    Em6DetectorConstruction*    Em6Detector;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


