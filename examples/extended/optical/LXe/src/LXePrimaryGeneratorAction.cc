#include "LXePrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXePrimaryGeneratorAction::LXePrimaryGeneratorAction(){
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  
  G4String particleName;
  particleGun->SetParticleDefinition(particleTable->
				     FindParticle(particleName="gamma"));
  //Default energy,position,momentum
  particleGun->SetParticleEnergy(511.*keV);
  particleGun->SetParticlePosition(G4ThreeVector(0.0 , 0.0, -20.0*cm));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXePrimaryGeneratorAction::~LXePrimaryGeneratorAction(){
    delete particleGun;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  particleGun->GeneratePrimaryVertex(anEvent);
}


