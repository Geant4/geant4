
#include "ExN04PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4HEPEvtInterface.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "ExN04PrimaryGeneratorMessenger.hh"

ExN04PrimaryGeneratorAction::ExN04PrimaryGeneratorAction()
{
  HEPEvt = new G4HEPEvtInterface("pythia_event.data");

  G4int n_particle = 1;
  G4ParticleGun* fParticleGun = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="mu+");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.,0.));
  fParticleGun->SetParticleEnergy(100.*GeV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  particleGun = fParticleGun;

  messenger = new ExN04PrimaryGeneratorMessenger(this);
  useHEPEvt = true;
}

ExN04PrimaryGeneratorAction::~ExN04PrimaryGeneratorAction()
{
  delete HEPEvt;
  delete particleGun;
  delete messenger;
}

void ExN04PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(useHEPEvt)
  { HEPEvt->GeneratePrimaryVertex(anEvent); }
  else
  { particleGun->GeneratePrimaryVertex(anEvent); }
}


