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
// G4SingleParticleSource source implementation
//
// Author: Fan Lei, QinetiQ ltd. - 05/02/2004
// Customer: ESA/ESTEC
// Revision: Andrea Dotti, SLAC
// --------------------------------------------------------------------

#include <cmath>

#include "G4SingleParticleSource.hh"

#include "G4SystemOfUnits.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4AutoLock.hh"

G4SingleParticleSource::part_prop_t::part_prop_t()
{
  momentum_direction = G4ParticleMomentum(1,0,0);
  energy = 1.*MeV;
  position = G4ThreeVector();
}

G4SingleParticleSource::G4SingleParticleSource()
{
  // Initialise all variables
  // Position distribution Variables

  NumberOfParticlesToBeGenerated = 1;
  definition = G4Geantino::GeantinoDefinition();

  charge = 0.0;
  time = 0;
  polarization = G4ThreeVector();

  biasRndm = new G4SPSRandomGenerator();
  posGenerator = new G4SPSPosDistribution();
  posGenerator->SetBiasRndm(biasRndm);
  angGenerator = new G4SPSAngDistribution();
  angGenerator->SetPosDistribution(posGenerator);
  angGenerator->SetBiasRndm(biasRndm);
  eneGenerator = new G4SPSEneDistribution();
  eneGenerator->SetBiasRndm(biasRndm);

  verbosityLevel = 0;
    
  G4MUTEXINIT(mutex);
}

G4SingleParticleSource::~G4SingleParticleSource()
{
  delete biasRndm;
  delete posGenerator;
  delete angGenerator;
  delete eneGenerator;

  G4MUTEXDESTROY(mutex);
}

void G4SingleParticleSource::SetVerbosity(G4int vL)
{
  G4AutoLock l(&mutex);
  verbosityLevel = vL;
  posGenerator->SetVerbosity(vL);
  angGenerator->SetVerbosity(vL);
  eneGenerator->SetVerbosity(vL);
}

void G4SingleParticleSource::
SetParticleDefinition(G4ParticleDefinition* aParticleDefinition)
{
  definition = aParticleDefinition;
  charge = aParticleDefinition->GetPDGCharge();
}

void G4SingleParticleSource::GeneratePrimaryVertex(G4Event* evt)
{
  if (definition == nullptr)
  {
    // TODO: Should this rise an exception???
    return;
  }

  if (verbosityLevel > 1)
  {
    G4cout << " NumberOfParticlesToBeGenerated: "
           << NumberOfParticlesToBeGenerated << G4endl;
  }

  part_prop_t& pp = ParticleProperties.Get();

  // Position stuff
  pp.position = posGenerator->GenerateOne();

  // Create a new vertex
  auto* vertex = new G4PrimaryVertex(pp.position,time);

  for (G4int i=0; i<NumberOfParticlesToBeGenerated; ++i)
  {
    // Angular stuff
    pp.momentum_direction = angGenerator->GenerateOne();

    // Energy stuff
    pp.energy = eneGenerator->GenerateOne(definition);

    if (verbosityLevel >= 2)
    {
      G4cout << "Creating primaries and assigning to vertex" << G4endl;
    }

    // Create new primaries and set them to the vertex
    //
    G4double mass = definition->GetPDGMass();
    auto* particle = new G4PrimaryParticle(definition);
    particle->SetKineticEnergy(pp.energy );
    particle->SetMass( mass );
    particle->SetMomentumDirection( pp.momentum_direction );
    particle->SetCharge( charge );
    particle->SetPolarization(polarization.x(),
                              polarization.y(),
                              polarization.z());
    if (verbosityLevel > 1)
    {
      G4cout << "Particle name: " << definition->GetParticleName() << G4endl;
      G4cout << "       Energy: " << pp.energy << G4endl;
      G4cout << "     Position: " << pp.position << G4endl;
      G4cout << "    Direction: " << pp.momentum_direction << G4endl;
    }

    // Set bweight equal to the multiple of all non-zero weights
    //
    G4double weight = eneGenerator->GetWeight()*biasRndm->GetBiasWeight();

    if(eneGenerator->IfApplyEnergyWeight())
    {
      weight *= eneGenerator->GetArbEneWeight(pp.energy);
    }

    // Pass it to primary particle
    //
    particle->SetWeight(weight);
    vertex->SetPrimary(particle);
  }

  // Now pass the weight to the primary vertex. CANNOT be used here!
  //  vertex->SetWeight(particle_weight);
  evt->AddPrimaryVertex(vertex);

  if (verbosityLevel > 1)
  {
    G4cout << " Primary Vetex generated !" << G4endl;
  }
}
