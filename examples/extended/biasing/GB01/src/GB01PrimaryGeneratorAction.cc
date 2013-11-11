#include "GB01PrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"


#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB01PrimaryGeneratorAction::GB01PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="gamma");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(100.*MeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.,0.,-50*cm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB01PrimaryGeneratorAction::~GB01PrimaryGeneratorAction()
{
  delete particleGun;
}


void GB01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->GeneratePrimaryVertex(anEvent);
}

