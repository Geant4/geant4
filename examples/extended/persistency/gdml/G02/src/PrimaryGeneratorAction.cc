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
// $Id: PrimaryGeneratorAction.cc,v 1.2 2008-12-18 12:57:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class PrimaryGeneratorAction implementation
//
// ----------------------------------------------------------------------------

#include "PrimaryGeneratorAction.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"

// ----------------------------------------------------------------------------

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  // Particle gun and particle table 
  //
  particleGun = new G4ParticleGun();
  particleTable = G4ParticleTable::GetParticleTable();

  // Default particle
  //
  particleGun->SetParticleDefinition(particleTable->FindParticle("geantino"));
  particleGun->SetParticleEnergy( 1.0*MeV );

  G4ThreeVector err1=G4ThreeVector(-1260,-560,40); // outside
  G4ThreeVector err2=G4ThreeVector(100,-240,120);  // inside
  G4ThreeVector err2v=(err2-err1).unit();
  
  particleGun->SetParticleMomentumDirection(err2v);
  particleGun->SetParticlePosition(err1);

}

// ----------------------------------------------------------------------------

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

// ----------------------------------------------------------------------------

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->GeneratePrimaryVertex(anEvent);
}
