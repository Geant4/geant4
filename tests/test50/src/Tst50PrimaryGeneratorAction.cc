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
// $Id: Tst50PrimaryGeneratorAction.cc,v 1.10 2003-05-17 18:11:53 guatelli Exp $
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
#include "G4RadioactiveDecay.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "Tst50PrimaryGeneratorAction.hh"
#include "Tst50PrimaryGeneratorMessenger.hh"

Tst50PrimaryGeneratorAction::Tst50PrimaryGeneratorAction():
  randomDirection("off")
{
  G4int particleNumber = 1;

  particleGun = new G4ParticleGun(particleNumber);
  gunMessenger = new Tst50PrimaryGeneratorMessenger(this);
  
  // default particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable -> FindParticle("gamma");
  particleGun -> SetParticleDefinition(particle);
  particleGun -> SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun -> SetParticleEnergy(1.*MeV);
  particleGun -> SetParticlePosition(G4ThreeVector(0.*m,0.*m,-1.*m));
}

Tst50PrimaryGeneratorAction::~Tst50PrimaryGeneratorAction()
{
  delete gunMessenger;
  delete particleGun;
}

void Tst50PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  if(randomDirection == "on")
    {
      G4double a,b,c;
      G4double n;
      do{
	a = (G4UniformRand()-0.5)/0.5;
	b = (G4UniformRand()-0.5)/0.5; 
	c = (G4UniformRand()-0.5)/0.5;
	n = a*a+b*b+c*c;
      }while(n > 1 || n == 0.0);
      n = sqrt(n);
      a /= n;
      b /= n;
      c /= n;
      G4ThreeVector direction(a,b,c);
      particleGun -> SetParticleMomentumDirection(direction);
    }
  
  if(randomDirection == "off")
    {
      particleGun -> SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    }

  particleGun -> GeneratePrimaryVertex(anEvent);
}

G4double Tst50PrimaryGeneratorAction::GetInitialEnergy()
{
  G4double primaryParticleEnergy = particleGun->GetParticleEnergy(); 
  return primaryParticleEnergy;
}

G4String Tst50PrimaryGeneratorAction::GetParticle()
{
  G4String primaryParticleName = particleGun->GetParticleDefinition()->GetParticleName();
  return primaryParticleName;
}





