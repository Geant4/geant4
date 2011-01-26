
#include "PrimaryGenerator.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4Electron.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "Randomize.hh"


PrimaryGenerator::PrimaryGenerator() {

  primaryKineticEnergy = 1 * MeV;
  sigmaKineticEnergy = 0 * keV;
  sigmaSpatialPlacement = 0 * mm;
  incidentAngle = 0.0 * deg;

  particleGun = new G4ParticleGun(G4Electron::Electron(),1);  
  messenger = new PrimaryGeneratorMessenger(this);
}


PrimaryGenerator::~PrimaryGenerator() {

  delete particleGun;
  delete messenger;
}


void PrimaryGenerator::GeneratePrimaries(G4Event* event) {

  G4double kineticEnergy = primaryKineticEnergy;
  G4double latVar1 = 0.0 * mm;
  G4double latVar2 = 0.0 * mm;

  if(sigmaKineticEnergy > 0.0 * keV)
     kineticEnergy = G4RandGauss::shoot(primaryKineticEnergy, 
                                        sigmaKineticEnergy);

  if(sigmaSpatialPlacement > 0.0 * mm) { 
     latVar1 = G4RandGauss::shoot(0.0 * mm, sigmaSpatialPlacement);
     latVar2 = G4RandGauss::shoot(0.0 * mm, sigmaSpatialPlacement);
  }

  particleGun -> SetParticleEnergy(kineticEnergy);
  particleGun -> SetParticlePosition(
         G4ThreeVector(latVar1, latVar2, -10.0 * mm).rotateX(incidentAngle));
  particleGun -> SetParticleMomentumDirection(
         G4ThreeVector(0.0 * mm, 0.0 * mm, 10.0 * mm).rotateX(incidentAngle));
  particleGun -> GeneratePrimaryVertex(event);
}


void PrimaryGenerator::SetPrimaryKineticEnergy(G4double kinEnergy) {

  if(kinEnergy > 0.0 * MeV) {
     primaryKineticEnergy = kinEnergy;
  }
}


void PrimaryGenerator::SetSigmaKineticEnergy(G4double sigma) {

  if(sigma >= 0.0 * keV) {
     sigmaKineticEnergy = sigma;
  }
}


void PrimaryGenerator::SetSigmaSpatialPlacement(G4double sigma) {

  if(sigma >= 0.0 * mm) {
     sigmaSpatialPlacement = sigma;
  }
}


void PrimaryGenerator::SetIncidentAngle(G4double angle) {

  if(angle >= 0.0 * deg && angle < 90.0 * deg) {
     incidentAngle = angle;
  }
}


