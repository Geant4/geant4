// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleGun.cc,v 1.3 2000-10-18 12:41:25 kurasige Exp $
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
  particle_charge = 0.0;
  theMessenger = new G4ParticleGunMessenger(this);
}

G4ParticleGun::~G4ParticleGun()
{
  delete theMessenger;
}

void G4ParticleGun::SetParticleDefinition
                 (G4ParticleDefinition * aParticleDefinition)
{ 
  particle_definition = aParticleDefinition; 
  particle_charge = particle_definition->GetPDGCharge();
}
 
void G4ParticleGun::SetParticleMomentum(G4ParticleMomentum aMomentum)
{
  if(particle_definition==NULL)
  {
    G4cout <<"Particle Definition not defined yet for G4ParticleGun"<< G4endl;
    G4cout <<"Zero Mass is assumed"<<G4endl;
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
      G4cout << "G4ParticleGun::" << particle_definition->GetParticleName() << G4endl;
      G4cout << "  KineticEnergy and Momentum could be inconsistent" << G4endl;
      G4cout << " (Momentum:" << p/GeV << " GeV/c";
      G4cout << "  Mass:" << mass/GeV << " GeV/c/c)" << G4endl;
      G4cout << "  KineticEnergy is overwritten!! ";
      G4cout << particle_energy/GeV << "->";
      G4cout << (sqrt(p*p+mass*mass)-mass)/GeV << "GeV" << G4endl;
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
    particle->SetCharge( particle_charge );
    particle->SetPolarization(particle_polarization.x(),
                               particle_polarization.y(),
                               particle_polarization.z());
    vertex->SetPrimary( particle );
  }

  evt->AddPrimaryVertex( vertex );
}


