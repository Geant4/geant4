// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayShooter.cc,v 1.1 2000-01-29 00:44:53 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4RayShooter.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "G4Geantino.hh"

G4RayShooter::G4RayShooter()
{
  SetInitialValues();
}

void G4RayShooter::SetInitialValues()
{
  particle_definition = G4Geantino::GeantinoDefinition();
  G4ThreeVector zero;
  particle_momentum_direction = (G4ParticleMomentum)zero;
  particle_energy = 1.0*GeV;
  particle_position = zero;
  particle_time = 0.0;
  particle_polarization = zero;
}

G4RayShooter::~G4RayShooter()
{
}

void G4RayShooter::Shoot(G4Event* evt,G4ThreeVector vtx,G4ThreeVector direc)
{
  // create a new vertex
  G4PrimaryVertex* vertex = new G4PrimaryVertex(vtx,particle_time);

  // create new primaries and set them to the vertex
  G4double mass =  particle_definition->GetPDGMass();
  G4double energy = particle_energy + mass;
  G4double pmom = sqrt(energy*energy-mass*mass);
  G4double px = pmom*direc.x();
  G4double py = pmom*direc.y();
  G4double pz = pmom*direc.z();
  G4PrimaryParticle* particle =
      new G4PrimaryParticle(particle_definition,px,py,pz);
  particle->SetMass( mass );
  particle->SetPolarization(particle_polarization.x(),
                               particle_polarization.y(),
                               particle_polarization.z());
  vertex->SetPrimary( particle );

  evt->AddPrimaryVertex( vertex );
}


