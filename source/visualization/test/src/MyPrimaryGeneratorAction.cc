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
// $Id: MyPrimaryGeneratorAction.cc,v 1.6 2006-06-29 21:34:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyPrimaryGeneratorAction.hh"

#include "G4Event.hh"

#ifdef GPS
#include "G4GeneralParticleSource.hh"
#else
#include "G4ParticleGun.hh"
#endif

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

MyPrimaryGeneratorAction::MyPrimaryGeneratorAction()
{
  G4int n_particle = 1;
#ifdef GPS
  particleGun = new G4GeneralParticleSource;
#else
  particleGun = new G4ParticleGun(n_particle);
#endif
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle 
    = particleTable->FindParticle(particleName="geantino");
  particleGun->SetParticleDefinition(particle);

#ifndef GPS
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.,0.));
  particleGun->SetParticleEnergy(100.*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,-300.*cm,0.*cm));
#endif
}

MyPrimaryGeneratorAction::~MyPrimaryGeneratorAction()
{
  delete particleGun;
}

void MyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->GeneratePrimaryVertex(anEvent);
}

#ifdef GPS
G4GeneralParticleSource* MyPrimaryGeneratorAction::GetGeneralParticleSource()
{
  return particleGun;
} 
#else
G4ParticleGun* MyPrimaryGeneratorAction::GetParticleGun()
{
  return particleGun;
} 
#endif
