// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleGun.cc,v 1.1 1999-01-07 16:06:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// G4ParticleGun
#include "G4ParticleGun.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleGunMessenger.hh"
#include "G4Event.hh"
#include "G4ios.hh"

G4ParticleGun::G4ParticleGun()
{
  SetInitialValues();
}

G4ParticleGun::G4ParticleGun(G4int numberofparticles)
{
  SetInitialValues();
  NumberOfParticlesToBeGenerated = numberofparticles;
}

G4ParticleGun::G4ParticleGun
    (G4ParticleDefinition * particleDef, G4int numberofparticles)
{
  SetInitialValues();
  NumberOfParticlesToBeGenerated = numberofparticles;
  particle_definition = particleDef;
}

void G4ParticleGun::SetInitialValues()
{
  NumberOfParticlesToBeGenerated = 1;
  particle_definition = NULL;
  G4ThreeVector zero;
  particle_momentum_direction = (G4ParticleMomentum)zero;
  particle_energy = 0.0;
  particle_position = zero;
  particle_time = 0.0;
  particle_polarization = zero;
  theMessenger = new G4ParticleGunMessenger(this);
}

G4ParticleGun::~G4ParticleGun()
{
  delete theMessenger;
}

void G4ParticleGun::SetParticleMomentum(G4ParticleMomentum aMomentum)
{
  if(particle_definition==NULL)
  {
    G4cout <<"Particle Definition not defined yet for G4ParticleGun"<< endl;
    G4cout <<"Zero Mass is assumed"<<endl;
    particle_momentum_direction =  aMomentum.unit();
    particle_energy = aMomentum.mag();
  } 
  else 
  {
    G4double mass =  particle_definition->GetPDGMass();
    G4double p = aMomentum.mag();
    particle_momentum_direction =  aMomentum.unit();
    if ((particle_energy>0.0)&&(abs(particle_energy+mass-sqrt(p*p+mass*mass))>keV))
    {
      G4cout << "G4ParticleGun::" << particle_definition->GetParticleName() << endl;
      G4cout << "  KineticEnergy and Momentum could be inconsistent" << endl;
      G4cout << " (Momentum:" << p/GeV << " GeV/c";
      G4cout << "  Mass:" << mass/GeV << " GeV/c/c)" << endl;
      G4cout << "  KineticEnergy is overwritten!! ";
      G4cout << particle_energy/GeV << "->";
      G4cout << (sqrt(p*p+mass*mass)-mass)/GeV << "GeV" << endl;
    }
    particle_energy = sqrt(p*p+mass*mass)-mass;
  }
}

void G4ParticleGun::GeneratePrimaryVertex(G4Event* evt)
{
  if(particle_definition==NULL) return;

  // create a new vertex
  G4PrimaryVertex* vertex = 
    new G4PrimaryVertex(particle_position,particle_time);

  // create new primaries and set them to the vertex
  G4double mass =  particle_definition->GetPDGMass();
  G4double energy = particle_energy + mass;
  G4double pmom = sqrt(energy*energy-mass*mass);
  G4double px = pmom*particle_momentum_direction.x();
  G4double py = pmom*particle_momentum_direction.y();
  G4double pz = pmom*particle_momentum_direction.z();
  for( int i=0; i<NumberOfParticlesToBeGenerated; i++ )
  {
    G4PrimaryParticle* particle =
      new G4PrimaryParticle(particle_definition,px,py,pz);
    particle->SetMass( mass );
    particle->SetPolarization(particle_polarization.x(),
                               particle_polarization.y(),
                               particle_polarization.z());
    vertex->SetPrimary( particle );
  }

  evt->AddPrimaryVertex( vertex );
}


