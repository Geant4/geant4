#ifndef GB01PrimaryGeneratorAction_h
#define GB01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class GB01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  GB01PrimaryGeneratorAction();    
  virtual ~GB01PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

private:
  G4ParticleGun*           particleGun;         //pointer a to G4  class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


