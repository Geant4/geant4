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
// $Id: Tst51PrimaryGeneratorAction.cc,v 1.1 2005-07-05 11:06:27 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "Tst51PrimaryGeneratorAction.hh"
//#include "Tst51PrimaryGeneratorMessenger.hh"

Tst51PrimaryGeneratorAction::Tst51PrimaryGeneratorAction()
{
  G4int particleNumber = 1;

  particleGun = new G4ParticleGun(particleNumber);
  // gunMessenger = new Tst51PrimaryGeneratorMessenger(this);
  
  // default particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4cout << "G4ParticleTable loaded" << G4endl;
  
  G4ParticleDefinition* particle = particleTable -> FindParticle("e-");
  particleGun -> SetParticleDefinition(particle);

  // the angle beteen the e- direction and the Z axis is 45. degrees
  G4double d = 1. *m;
  G4double angle = 45. *deg;

  G4double dir = std::cos(angle);
  
  particleGun -> SetParticlePosition(G4ThreeVector(0. , d, -d));
 
  particleGun -> SetParticleMomentumDirection(G4ThreeVector(0., -dir, dir));
  particleGun -> SetParticleEnergy(70. * keV);
}

Tst51PrimaryGeneratorAction::~Tst51PrimaryGeneratorAction()
{
  //delete gunMessenger;
  delete particleGun;
}

void Tst51PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  particleGun -> GeneratePrimaryVertex(anEvent);
}

