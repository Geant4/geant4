// $Id: A01PrimaryGeneratorAction.cc,v 1.1 2002-11-13 07:24:01 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "A01PrimaryGeneratorAction.hh"
#include "A01PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

A01PrimaryGeneratorAction::A01PrimaryGeneratorAction()
{
  momentum = 1000.*MeV;
  sigmaMomentum = 50.*MeV;
  sigmaAngle = 2.*deg;
  randomizePrimary = true;

  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);

  //create a messenger for this class
  gunMessenger = new A01PrimaryGeneratorMessenger(this);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  positron = particleTable->FindParticle(particleName="e+");
  muon = particleTable->FindParticle(particleName="mu+");
  pion = particleTable->FindParticle(particleName="pi+");
  kaon = particleTable->FindParticle(particleName="kaon+");
  proton = particleTable->FindParticle(particleName="proton");

  // default particle kinematics
  particleGun->SetParticlePosition(G4ThreeVector(0.,0.,-8.*m));
  particleGun->SetParticleDefinition(positron);
}

A01PrimaryGeneratorAction::~A01PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

void A01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* particle;

  if(randomizePrimary)
  {
    /////////////////////////////////////////////G4int i = (int)(5.*G4UniformRand());
    G4int i = (int)(2.*G4UniformRand());
    switch(i)
    {
      case 0:
        particle = positron;
        break;
      case 1:
        particle = muon;
        break;
      case 2:
        particle = pion;
        break;
      case 3:
        particle = kaon;
        break;
      default:
        particle = proton;
        break;
    }
    particleGun->SetParticleDefinition(particle);
  }
  else
  {
    particle = particleGun->GetParticleDefinition();
  }

  G4double pp = momentum + (G4UniformRand()-0.5)*sigmaMomentum;
  G4double mass = particle->GetPDGMass();
  G4double Ekin = sqrt(pp*pp+mass*mass)-mass;
  particleGun->SetParticleEnergy(Ekin);

  G4double angle = (G4UniformRand()-0.5)*sigmaAngle;
  particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(angle),0.,cos(angle)));

  particleGun->GeneratePrimaryVertex(anEvent);
}

