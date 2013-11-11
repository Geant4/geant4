#ifndef GB02PrimaryGeneratorAction_h
#define GB02PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class GB02PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  GB02PrimaryGeneratorAction();    
  virtual ~GB02PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

private:
  G4ParticleGun*           particleGun;         //pointer a to G4  class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


