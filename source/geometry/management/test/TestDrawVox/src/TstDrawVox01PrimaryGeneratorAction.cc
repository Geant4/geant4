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

#include "TstDrawVox01PrimaryGeneratorAction.hh"
#include "TstDrawVox01PrimaryGeneratorMessenger.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4TransportationManager.hh"
#include <CLHEP/Random/RandFlat.h>

TstDrawVox01PrimaryGeneratorAction::TstDrawVox01PrimaryGeneratorAction():
  generatorAction (standardGun),
  particleGun (0),
  messenger (0),
  worldVolume (0)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  // default particle

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle 
    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(1.*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));

  // messenger

  messenger = new TstDrawVox01PrimaryGeneratorMessenger (this);

  // world extent

  worldVolume = G4TransportationManager::GetTransportationManager ()
    -> GetNavigatorForTracking () -> GetWorldVolume ();
  if (worldVolume) worldExtent = worldVolume -> GetLogicalVolume ()
		     -> GetSolid () -> GetExtent ();

}

TstDrawVox01PrimaryGeneratorAction::~TstDrawVox01PrimaryGeneratorAction()
{
  delete particleGun;
  delete messenger;
}

void TstDrawVox01PrimaryGeneratorAction::SelectPrimaryGeneratorAction
(Action action)
{
  generatorAction = action;
}

void TstDrawVox01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4VPhysicalVolume* currentWorldVolume;

  double costheta;
  double sintheta;
  double phi;
  double cosphi;
  double sinphi;

  switch (generatorAction) {

  case standardGun:

    particleGun->GeneratePrimaryVertex(anEvent);
    break;

  case randomDirectionGun:

    costheta = CLHEP::RandFlat::shoot (-1., 1.);
    sintheta = std::sqrt (1. - costheta * costheta);
    phi      = CLHEP::RandFlat::shoot (twopi);
    cosphi   = std::cos (phi);
    sinphi   = std::sin (phi);
    particleGun->SetParticleMomentumDirection
      (G4ThreeVector (sintheta * cosphi, sintheta * sinphi, costheta));

    particleGun->GeneratePrimaryVertex(anEvent);
    break;

  case randomPositionGun:

    // Check if world is in place or has changed.
    currentWorldVolume = 
      G4TransportationManager::GetTransportationManager ()
      -> GetNavigatorForTracking () -> GetWorldVolume ();
    if (!worldVolume ||	worldVolume != currentWorldVolume) {
      worldVolume = currentWorldVolume;
      if (worldVolume) worldExtent = worldVolume -> GetLogicalVolume ()
			 -> GetSolid () -> GetExtent ();
    }

    particleGun->SetParticlePosition
      (G4ThreeVector
       (CLHEP::RandFlat::shoot (worldExtent.GetXmin (), worldExtent.GetXmax ()),
	CLHEP::RandFlat::shoot (worldExtent.GetYmin (), worldExtent.GetYmax ()),
	CLHEP::RandFlat::shoot (worldExtent.GetZmin (), worldExtent.GetZmax ())));

    particleGun->GeneratePrimaryVertex(anEvent);
    break;

  case randomPositionAndDirectionGun:

    // Check if world is in place or has changed.
    currentWorldVolume =
      G4TransportationManager::GetTransportationManager ()
      -> GetNavigatorForTracking () -> GetWorldVolume ();
    if (!worldVolume ||	worldVolume != currentWorldVolume) {
      worldVolume = currentWorldVolume;
      if (worldVolume) worldExtent = worldVolume -> GetLogicalVolume ()
			 -> GetSolid () -> GetExtent ();
    }

    particleGun->SetParticlePosition
      (G4ThreeVector
       (CLHEP::RandFlat::shoot (worldExtent.GetXmin (), worldExtent.GetXmax ()),
	CLHEP::RandFlat::shoot (worldExtent.GetYmin (), worldExtent.GetYmax ()),
	CLHEP::RandFlat::shoot (worldExtent.GetZmin (), worldExtent.GetZmax ())));

    costheta = CLHEP::RandFlat::shoot (-1., 1.);
    sintheta = std::sqrt (1. - costheta * costheta);
    phi      = CLHEP::RandFlat::shoot (twopi);
    cosphi   = std::cos (phi);
    sinphi   = std::sin (phi);
    particleGun->SetParticleMomentumDirection
      (G4ThreeVector (sintheta * cosphi, sintheta * sinphi, costheta));

    particleGun->GeneratePrimaryVertex(anEvent);
    break;

  }
}
