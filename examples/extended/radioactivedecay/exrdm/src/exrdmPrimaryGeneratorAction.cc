
// **********************************************************************

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"

#include "exrdmPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmPrimaryGeneratorAction::exrdmPrimaryGeneratorAction()
{
  particleGun = new G4GeneralParticleSource ();
}

exrdmPrimaryGeneratorAction::~exrdmPrimaryGeneratorAction()
{
  delete particleGun;
}

void exrdmPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->GeneratePrimaryVertex(anEvent);
}



