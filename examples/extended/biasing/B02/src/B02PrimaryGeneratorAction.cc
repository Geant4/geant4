#include "B02PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4UImanager.hh"
#include "globals.hh"

B02PrimaryGeneratorAction::B02PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/gun/particle neutron");
  UI->ApplyCommand("/gun/energy 10.0 MeV");
  UI->ApplyCommand("/gun/position 0.0 0.0 -15.0005 cm");
  UI->ApplyCommand("/gun/direction 0. 0. 1.");
}

B02PrimaryGeneratorAction::~B02PrimaryGeneratorAction()
{
  delete particleGun;
}

void B02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4int i = anEvent->GetEventID() % 3;
  particleGun->GeneratePrimaryVertex(anEvent);
}


