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
// $Id: Tst52PrimaryGeneratorAction.cc,v 1.1 2007-04-12 12:00:17 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------

#include "G4IonTable.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "Tst52PrimaryGeneratorAction.hh"
#include "Tst52PrimaryGeneratorMessenger.hh"

Tst52PrimaryGeneratorAction::Tst52PrimaryGeneratorAction()
{
  G4int particleNumber = 1;

  particleGun = new G4ParticleGun(particleNumber);
  gunMessenger = new Tst52PrimaryGeneratorMessenger(this);
  
  // default particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
 
 G4ParticleDefinition* particle = particleTable -> FindParticle("e-");

  particleGun -> SetParticleDefinition(particle);

  particleGun -> SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  particleGun -> SetParticleEnergy(1.*MeV);

  particleGun -> SetParticlePosition(G4ThreeVector(0.*m,0.*m,-1.*m));
}

Tst52PrimaryGeneratorAction::~Tst52PrimaryGeneratorAction()
{
  delete gunMessenger;
  delete particleGun;
}

void Tst52PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  particleGun -> GeneratePrimaryVertex(anEvent);
}

G4double Tst52PrimaryGeneratorAction::GetInitialEnergy()
{
  G4double primaryParticleEnergy = particleGun -> GetParticleEnergy(); 
  return primaryParticleEnergy;
}

G4String Tst52PrimaryGeneratorAction::GetParticle()
{
  G4String primaryParticleName = particleGun -> GetParticleDefinition() 
                                             ->GetParticleName();
  return primaryParticleName;
}



