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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRPrimaryGeneratorAction.cc
//   Gorad primary generator action class
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

GRPrimaryGeneratorAction::GRPrimaryGeneratorAction(
             G4bool useParticleGun, G4bool useParticleSource)
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(nullptr), fParticleSource(nullptr)     
{
  if(useParticleGun)
  {
    fParticleGun  = new G4ParticleGun(1);
  
    auto particleTable = G4ParticleTable::GetParticleTable();
    auto fPion = particleTable->FindParticle("pi+");
    fParticleGun->SetParticleDefinition(fPion);
  
    // default particle kinematics
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
    fParticleGun->SetParticleEnergy(1.*GeV);
  }

  if(useParticleSource)
  { fParticleSource = new G4GeneralParticleSource(); }
}

GRPrimaryGeneratorAction::~GRPrimaryGeneratorAction()
{
  if(fParticleGun) delete fParticleGun;
  if(fParticleSource) delete fParticleSource;
}

void GRPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  if(fParticleGun) fParticleGun->GeneratePrimaryVertex(event);
  if(fParticleSource) fParticleSource->GeneratePrimaryVertex(event);
}

