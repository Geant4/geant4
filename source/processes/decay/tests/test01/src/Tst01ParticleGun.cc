//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: Tst01ParticleGun.cc,v 1.6 2008-10-24 15:24:16 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// G4ParticleGun
#include "Tst01ParticleGun.hh"
#include "Tst01ParticleGunMessenger.hh"
#include "G4Event.hh"



Tst01ParticleGun::Tst01ParticleGun()
                :G4ParticleGun()
{
  SetInitialValues();
}


Tst01ParticleGun::Tst01ParticleGun(G4int )
                 :G4ParticleGun(1)
{
  SetInitialValues();
}


Tst01ParticleGun::Tst01ParticleGun(G4ParticleDefinition * particleDef, 
				   G4int )
                 :G4ParticleGun(particleDef, 1)  
{
  SetInitialValues();
}


Tst01ParticleGun::~Tst01ParticleGun()
{
  while (!decay_products.empty()){
    G4PrimaryParticle* particle = (G4PrimaryParticle*)(decay_products.back());
    decay_products.pop_back();
    delete particle;
  }
}
    
Tst01ParticleGun::Tst01ParticleGun(const Tst01ParticleGun &right)
  : G4ParticleGun()
{
  *this = right;
}


const Tst01ParticleGun & Tst01ParticleGun::operator=(const Tst01ParticleGun &right)
{
  if ( this == &right) return *this;
  NumberOfParticlesToBeGenerated = right.NumberOfParticlesToBeGenerated;
  particle_definition = right.particle_definition;
  particle_momentum_direction = right.particle_momentum_direction; 
  particle_energy = right.particle_energy;
  particle_charge = right.particle_charge;
  particle_polarization = right.particle_polarization;
  particle_decay_time  =  right.particle_decay_time;
  decay_products.clear();
  return *this;
}

void Tst01ParticleGun::GeneratePrimaryVertex(G4Event* evt)
{
  if(particle_definition==NULL) return;

  // create a new vertex
  G4PrimaryVertex* vertex = 
    new G4PrimaryVertex(particle_position,particle_time);

  // create new primaries and set them to the vertex
  G4double mass =  particle_definition->GetPDGMass();
  G4double energy = particle_energy + mass;
  G4double pmom = std::sqrt(energy*energy-mass*mass);
  G4double px = pmom*particle_momentum_direction.x();
  G4double py = pmom*particle_momentum_direction.y();
  G4double pz = pmom*particle_momentum_direction.z();
 
  G4PrimaryParticle* particle =
    new G4PrimaryParticle(particle_definition,px,py,pz);
  particle->SetMass( mass );
  particle->SetCharge( particle_charge );
  particle->SetPolarization(particle_polarization.x(),
			    particle_polarization.y(),
			    particle_polarization.z());
  if (particle_decay_time>=0.) particle->SetProperTime(particle_decay_time);
  while (!decay_products.empty()){
    G4PrimaryParticle* daughter = (G4PrimaryParticle*)(decay_products.back());
    decay_products.pop_back();
    particle->SetDaughter(daughter);
  }
  vertex->SetPrimary( particle );

  evt->AddPrimaryVertex( vertex );
}


void Tst01ParticleGun::SetInitialValues()
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

  particle_decay_time = -1.0;

  theMessenger = new Tst01ParticleGunMessenger(this);
}











