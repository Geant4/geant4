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
/// \file analysis/A01/src/A01PrimaryGeneratorAction.cc
/// \brief Implementation of the A01PrimaryGeneratorAction class
//
// $Id$
// --------------------------------------------------------------
//

#include "A01PrimaryGeneratorAction.hh"
#include "A01PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

A01PrimaryGeneratorAction::A01PrimaryGeneratorAction()
{
  fMomentum = 1000.*MeV;
  fSigmaMomentum = 50.*MeV;
  fSigmaAngle = 2.*deg;
  fRandomizePrimary = true;

  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  //create a messenger for this class
  fGunMessenger = new A01PrimaryGeneratorMessenger(this);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  fPositron = particleTable->FindParticle(particleName="e+");
  fMuon = particleTable->FindParticle(particleName="mu+");
  fPion = particleTable->FindParticle(particleName="pi+");
  fKaon = particleTable->FindParticle(particleName="kaon+");
  fProton = particleTable->FindParticle(particleName="proton");

  // default particle kinematics
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-8.*m));
  fParticleGun->SetParticleDefinition(fPositron);
}

A01PrimaryGeneratorAction::~A01PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
}

void A01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* particle;

  if(fRandomizePrimary)
  {
    /////////////////////////////////////////////G4int i = (int)(5.*G4UniformRand());
    G4int i = (int)(2.*G4UniformRand());
    switch(i)
    {
      case 0:
        particle = fPositron;
        break;
      case 1:
        particle = fMuon;
        break;
      case 2:
        particle = fPion;
        break;
      case 3:
        particle = fKaon;
        break;
      default:
        particle = fProton;
        break;
    }
    fParticleGun->SetParticleDefinition(particle);
  }
  else
  {
    particle = fParticleGun->GetParticleDefinition();
  }

  G4double pp = fMomentum + (G4UniformRand()-0.5)*fSigmaMomentum;
  G4double mass = particle->GetPDGMass();
  G4double Ekin = std::sqrt(pp*pp+mass*mass)-mass;
  fParticleGun->SetParticleEnergy(Ekin);

  G4double angle = (G4UniformRand()-0.5)*fSigmaAngle;
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(std::sin(angle),0.,std::cos(angle)));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

