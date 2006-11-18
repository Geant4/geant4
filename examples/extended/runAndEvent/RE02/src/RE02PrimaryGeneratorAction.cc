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
// $Id: RE02PrimaryGeneratorAction.cc,v 1.3 2006-11-18 01:37:24 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RE02PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"

//
RE02PrimaryGeneratorAction::RE02PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

// default particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("proton");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,1.));
  particleGun->SetParticleEnergy(150.0*MeV);
//
// default beam position
  G4double position = -200./2.*cm;
//
// Initial beam spot size in sigma.; This is not a part of ParticleGun.
  fsigmaPosition = 10.* mm;

//
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm, 0.*cm, position));
}

//
RE02PrimaryGeneratorAction::~RE02PrimaryGeneratorAction()
{
  delete particleGun;
}

//
void RE02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 

  G4ThreeVector position = particleGun->GetParticlePosition();
  G4double dx = (G4UniformRand()-0.5)*fsigmaPosition;
  G4double dy = (G4UniformRand()-0.5)*fsigmaPosition;
  position.setX(dx);
  position.setY(dy);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);
}
