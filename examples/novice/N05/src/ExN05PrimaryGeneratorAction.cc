// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05PrimaryGeneratorAction.cc,v 1.2 1999-12-15 14:49:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN05PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

ExN05PrimaryGeneratorAction::ExN05PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle 
    = particleTable->FindParticle(particleName="geantino");
  particleGun->SetParticleDefinition(particle);

  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.,0.));
  particleGun->SetParticleEnergy(100.*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,-300.*cm,0.*cm));
}

ExN05PrimaryGeneratorAction::~ExN05PrimaryGeneratorAction()
{
  delete particleGun;
}

void ExN05PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->GeneratePrimaryVertex(anEvent);
}

G4ParticleGun* ExN05PrimaryGeneratorAction::GetParticleGun()
{
  return particleGun;
} 

