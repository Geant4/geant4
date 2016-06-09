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
// $Id: ExN07PrimaryGeneratorAction.cc,v 1.1 2003/03/10 01:43:36 asaim Exp $
// GEANT4 tag $Name: geant4-05-01 $
//

#include "ExN07PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

ExN07PrimaryGeneratorAction::ExN07PrimaryGeneratorAction()
:serial(false)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="mu-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(100.*GeV);
}

ExN07PrimaryGeneratorAction::~ExN07PrimaryGeneratorAction()
{
  delete particleGun;
}

void ExN07PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(serial)
  {
    particleGun->SetParticlePosition(G4ThreeVector(0.,0.,-3.5*m));
    particleGun->GeneratePrimaryVertex(anEvent);
  }
  else
  {
    for(G4int i=0;i<3;i++)
    {
      particleGun->SetParticlePosition(G4ThreeVector(0.,G4double(i-1)*m,-1.5*m));
      particleGun->GeneratePrimaryVertex(anEvent);
    }
  }
}

