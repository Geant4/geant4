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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the CaTS::PrimaryGeneratorAction class

// Geant4 headers:
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4HEPEvtInterface.hh"
#include "G4Version.hh"
#include "G4GenericMessenger.hh"
// project headers:
#include "PrimaryGeneratorAction.hh"
// c++ headers:
#include <cstring>

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  DefineCommands();
  G4int n_particle            = 1;
  G4ParticleGun* fParticleGun = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle =
    particleTable->FindParticle(particleName = "mu+");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  fParticleGun->SetParticleEnergy(10. * CLHEP::GeV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  gentypeMap["particleGun"] = fParticleGun;
  gentypeMap["GPS"]         = new G4GeneralParticleSource;
  fcurrentGenerator         = gentypeMap["particleGun"];
  fcurrentGeneratorName     = "particleGun";
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fMessenger;
  for( auto& gentype : gentypeMap )
  {
    delete gentype.second;
  }
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fcurrentGenerator->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::DefineCommands()
{
  fMessenger = new G4GenericMessenger(this, "/CaTS/primaryGenerator/",
                                      "select primary Generator");
  // command to select primary Generator
  fMessenger
    ->DeclareMethod("setGenerator", &PrimaryGeneratorAction::SetGenerator)
    .SetGuidance("Select the primary Generator (Available generators : "
                 "particleGun,GPS")
    .SetParameterName("generator", true)
    .SetDefaultValue("particleGun")
    .SetCandidates("particleGun GPS");
  fMessenger->DeclareMethod("Print", &PrimaryGeneratorAction::Print)
    .SetGuidance("Print name of primary generator")
    .SetStates(G4State_PreInit, G4State_Idle);
}
