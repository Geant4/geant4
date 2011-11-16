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
// $Id: Tst33PrimaryGeneratorAction.cc,v 1.7 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"

#include "Tst33PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4ThreeVector.hh"

Tst33PrimaryGeneratorAction::Tst33PrimaryGeneratorAction()
  :
  fParticleGun(new G4ParticleGun(1))
{
  if (!fParticleGun) {
    G4Exception("Tst33PrimaryGeneratorAction::Tst33PrimaryGeneratorAction()",
        "TST33-06", FatalException, " new failed to create G4ParticleGun!");
  }
  fParticleGun->SetParticleDefinition(G4Neutron::NeutronDefinition());
  //  fParticleGun->SetParticleDefinition(G4Gamma::GammaDefinition());
  //  fParticleGun->SetParticleDefinition(G4Proton::ProtonDefinition());
  fParticleGun->SetParticleEnergy(10.0*MeV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, -90.0005*cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
}

Tst33PrimaryGeneratorAction::~Tst33PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void Tst33PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
