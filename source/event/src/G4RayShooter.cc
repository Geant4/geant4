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
// G4RayShooter class implementation
//
// Author: Makoto Asai, 2000
// --------------------------------------------------------------------

#include "G4RayShooter.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

void G4RayShooter::Shoot(G4Event* evt, G4ThreeVector vtx, G4ThreeVector direc)
{
  if(particle_definition == nullptr)
  {
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    particle_definition = particleTable->FindParticle(particleName="geantino");
    if(particle_definition == nullptr)
    {
      G4String msg;
      msg = "G4RayTracer uses geantino to trace the ray, but your physics list does not\n";
      msg += "define G4Geantino. Please add G4Geantino in your physics list.";
      G4Exception("G4RayShooter::Shoot()", "RayTracer001", FatalException, msg);
    }
  }

  // Create a new vertex
  //
  auto* vertex = new G4PrimaryVertex(vtx,particle_time);

  // Create new primaries and set them to the vertex
  //
  G4double mass = particle_definition->GetPDGMass();
  auto* particle = new G4PrimaryParticle(particle_definition);
  particle->SetKineticEnergy( particle_energy );
  particle->SetMass( mass );
  particle->SetMomentumDirection( direc );
  particle->SetPolarization(particle_polarization.x(),
                            particle_polarization.y(),
                            particle_polarization.z());
  vertex->SetPrimary( particle );

  evt->AddPrimaryVertex( vertex );
}
