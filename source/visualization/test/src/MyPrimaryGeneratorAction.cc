//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: MyPrimaryGeneratorAction.cc,v 1.5 2005-02-04 16:28:49 johna Exp $
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
