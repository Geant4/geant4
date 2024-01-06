#include "generator.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
    fParticleGun = new G4ParticleGun(1);

    G4ThreeVector pos(0., 0., 0.);
    G4double energy = 0. * eV;

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleEnergy(energy);
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
    delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
    G4int Z = 9;
    G4int A = 18;

    G4double charge = 0. * eplus;
    G4double energy = 0. * keV;

    G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(Z, A, energy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(charge);

    fParticleGun->GeneratePrimaryVertex(anEvent);
}