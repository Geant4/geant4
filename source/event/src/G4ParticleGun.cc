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
// G4ParticleGun class implementation
//
// Author: Makoto Asai, 1997
// --------------------------------------------------------------------

#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
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

G4ParticleGun::G4ParticleGun(G4ParticleDefinition* particleDef,
                             G4int numberofparticles)
{
  SetInitialValues();
  NumberOfParticlesToBeGenerated = numberofparticles;
  SetParticleDefinition( particleDef );
}

void G4ParticleGun::SetInitialValues()
{
  NumberOfParticlesToBeGenerated = 1;
  particle_definition = nullptr;
  G4ThreeVector zero;
  particle_momentum_direction = (G4ParticleMomentum)zero;
  particle_energy = 0.0;
  particle_momentum = 0.0;
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

void G4ParticleGun::
SetParticleDefinition(G4ParticleDefinition* aParticleDefinition)
{ 
  if(aParticleDefinition == nullptr)
  {
    G4Exception("G4ParticleGun::SetParticleDefinition()", "Event0101",
                FatalException, "Null pointer is given.");
  }
  if(aParticleDefinition->IsShortLived())
  {
    if(aParticleDefinition->GetDecayTable() == nullptr)
    {
      G4ExceptionDescription ED;
      ED << "G4ParticleGun does not support shooting a short-lived "
         << "particle without a valid decay table." << G4endl;
      ED << "G4ParticleGun::SetParticleDefinition for "
         << aParticleDefinition->GetParticleName() << " is ignored." << G4endl;
      G4Exception("G4ParticleGun::SetParticleDefinition()", "Event0102",
                  JustWarning, ED);
      return;
    }
  }
  particle_definition = aParticleDefinition; 
  particle_charge = particle_definition->GetPDGCharge();
  if(particle_momentum>0.0)
  {
    G4double mass = particle_definition->GetPDGMass();
    particle_energy =
      std::sqrt(particle_momentum*particle_momentum+mass*mass)-mass;
  }
}

void G4ParticleGun::SetParticleEnergy(G4double aKineticEnergy)
{
  particle_energy = aKineticEnergy;
  if(particle_momentum>0.0)
  {
    if(particle_definition != nullptr)
    {
      G4cout << "G4ParticleGun::" << particle_definition->GetParticleName()
             << G4endl;
    }
    else
    {
      G4cout << "G4ParticleGun::" << " " << G4endl;
    }
    G4cout << " was defined in terms of Momentum: " 
           << particle_momentum/GeV << "GeV/c" << G4endl;
    G4cout << " is now defined in terms of KineticEnergy: " 
           << particle_energy/GeV   << "GeV"   << G4endl;
    particle_momentum = 0.0;
  }
}

void G4ParticleGun::SetParticleMomentum(G4double aMomentum)
{
  if(particle_energy>0.0)
  {
    if(particle_definition != nullptr)
    {
      G4cout << "G4ParticleGun::" << particle_definition->GetParticleName()
             << G4endl;
    }
    else
    {
      G4cout << "G4ParticleGun::" << " " << G4endl;
    }
    G4cout << " was defined in terms of KineticEnergy: "
           << particle_energy/GeV << "GeV"   << G4endl;
    G4cout << " is now defined in terms Momentum: "
           << aMomentum/GeV       << "GeV/c" << G4endl;
  }
  if(particle_definition == nullptr)
  {
    G4cout << "Particle Definition not defined yet for G4ParticleGun"
           << G4endl;
    G4cout << "Zero Mass is assumed" << G4endl;
    particle_momentum = aMomentum;
    particle_energy = aMomentum;
  }
  else
  {
    G4double mass = particle_definition->GetPDGMass();
    particle_momentum = aMomentum;
    particle_energy =
      std::sqrt(particle_momentum*particle_momentum+mass*mass)-mass;
  }
}
 
void G4ParticleGun::SetParticleMomentum(G4ParticleMomentum aMomentum)
{
  if(particle_energy>0.0)
  {
    if(particle_definition != nullptr)
    {
      G4cout << "G4ParticleGun::" << particle_definition->GetParticleName()
             << G4endl;
    }
    else
    {
      G4cout << "G4ParticleGun::" << " " << G4endl;
    }
    G4cout << " was defined in terms of KineticEnergy: "
           << particle_energy/GeV << "GeV"   << G4endl;
    G4cout << " is now defined in terms Momentum: "
           << aMomentum.mag()/GeV << "GeV/c" << G4endl;
  }
  if(particle_definition == nullptr)
  {
    G4cout << "Particle Definition not defined yet for G4ParticleGun"
           << G4endl;
    G4cout << "Zero Mass is assumed" << G4endl;
    particle_momentum_direction = aMomentum.unit();
    particle_momentum = aMomentum.mag();
    particle_energy = aMomentum.mag();
  } 
  else 
  {
    G4double mass =  particle_definition->GetPDGMass();
    particle_momentum = aMomentum.mag();
    particle_momentum_direction = aMomentum.unit();
    particle_energy = 
      std::sqrt(particle_momentum*particle_momentum+mass*mass)-mass;
  }
}

void G4ParticleGun::GeneratePrimaryVertex(G4Event* evt)
{
  if(particle_definition == nullptr) 
  {
    G4ExceptionDescription ED;
    ED << "Particle definition is not defined." << G4endl;
    ED << "G4ParticleGun::SetParticleDefinition() has to be invoked beforehand."
       << G4endl;
    G4Exception("G4ParticleGun::GeneratePrimaryVertex()", "Event0109",
                FatalException, ED);
    return;
  }

  // Create a new vertex
  //
  auto* vertex = 
    new G4PrimaryVertex(particle_position,particle_time);

  // Create new primaries and set them to the vertex
  //
  G4double mass =  particle_definition->GetPDGMass();
  for( G4int i=0; i<NumberOfParticlesToBeGenerated; ++i )
  {
    auto* particle =
      new G4PrimaryParticle(particle_definition);
    particle->SetKineticEnergy( particle_energy );
    particle->SetMass( mass );
    particle->SetMomentumDirection( particle_momentum_direction );
    particle->SetCharge( particle_charge );
    particle->SetPolarization(particle_polarization.x(),
                              particle_polarization.y(),
                              particle_polarization.z());
    vertex->SetPrimary( particle );
  }
  evt->AddPrimaryVertex( vertex );
}
