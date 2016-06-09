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
// $Id: ExN07PrimaryGeneratorAction.cc,v 1.3 2006-06-29 17:55:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

