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

#include <CLHEP/Random/RandFlat.h>

#include "Tst01PrimaryGeneratorAction.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "Tst01PrimaryGeneratorMessenger.hh"
#include "Tst01DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4TransportationManager.hh"

////////////////////////////////////////////////////////////////////////
//
//

Tst01PrimaryGeneratorAction::Tst01PrimaryGeneratorAction()
 : generatorAction (standardGun), particleGun(0),
   messenger(0), fGunPosition(0.0, 0.0, 0.0)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  // default particle

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle 
    = particleTable->FindParticle(particleName="geantino");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(1.*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));

  // messenger

  messenger = new Tst01PrimaryGeneratorMessenger (this);

  // world extent

  Tst01DetectorConstruction detector;
  G4double wSize = detector.GetWorldSize();
  fSize = std::sqrt( wSize*wSize + wSize*wSize + wSize*wSize ) ;

  G4cout << "World Size factor = " << fSize << G4endl;
}

///////////////////////////////////////////////////////////////////////////
//
// Destructor: delets pointers

Tst01PrimaryGeneratorAction::~Tst01PrimaryGeneratorAction()
{
  delete particleGun;
  delete messenger;
}

////////////////////////////////////////////////////////////////////////
//
//

void Tst01PrimaryGeneratorAction::SelectPrimaryGeneratorAction(Action action)
{
  generatorAction = action;
}

////////////////////////////////////////////////////////////////////////
//
//

void Tst01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ThreeVector direction, position ;
  G4double costheta, sintheta, phi, cosphi, sinphi ;

  switch (generatorAction) 
  {

    case standardGun:  

    particleGun->GeneratePrimaryVertex(anEvent);
    break;
//...................................................................
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
//......................................................................
    case randomPositionGun:   

    particleGun->SetParticlePosition
      (G4ThreeVector
       (CLHEP::RandFlat::shoot(-fSize/2, fSize/2),
	CLHEP::RandFlat::shoot(-fSize/2, fSize/2),
	CLHEP::RandFlat::shoot(-fSize/2, fSize/2)));

    particleGun->GeneratePrimaryVertex(anEvent);
    break;
//.................................................................
    case randomPositionAndDirectionGun:  

    particleGun->SetParticlePosition
      (G4ThreeVector
       (CLHEP::RandFlat::shoot(-fSize/2, fSize/2),
	CLHEP::RandFlat::shoot(-fSize/2, fSize/2),
	CLHEP::RandFlat::shoot(-fSize/2, fSize/2)));

    costheta = CLHEP::RandFlat::shoot (-1., 1.);
    sintheta = std::sqrt (1. - costheta * costheta);
    phi      = CLHEP::RandFlat::shoot (twopi);
    cosphi   = std::cos (phi);
    sinphi   = std::sin (phi);
    particleGun->SetParticleMomentumDirection
      (G4ThreeVector (sintheta * cosphi, sintheta * sinphi, costheta));

    particleGun->GeneratePrimaryVertex(anEvent);
    break;
//.................................................................
    case viewerGun:

    particleGun->SetParticlePosition(fGunPosition) ;
    if(fPosition)
    {
      direction = G4ThreeVector( -fGunPosition.x()/fPosition ,
                                 -fGunPosition.y()/fPosition ,
                                 -fGunPosition.z()/fPosition    ) ;

      costheta = CLHEP::RandFlat::shoot(fPosition
                                       /std::sqrt(fPosition*fPosition +
                                                  fSize*fSize), 1.);
      sintheta = std::sqrt (1. - costheta * costheta);
      phi      = CLHEP::RandFlat::shoot (twopi);
      cosphi   = std::cos (phi) ;
      sinphi   = std::sin (phi) ;
      position = G4ThreeVector (sintheta * cosphi, 
                                sintheta * sinphi, costheta) ;

      position.rotateUz(direction) ;

      particleGun->SetParticleMomentumDirection(position);      
    }
    else
    {
      costheta = CLHEP::RandFlat::shoot (-1., 1.);
      sintheta = std::sqrt (1. - costheta * costheta);
      phi      = CLHEP::RandFlat::shoot (twopi);
      cosphi   = std::cos (phi);
      sinphi   = std::sin (phi);
      particleGun->SetParticleMomentumDirection
      (G4ThreeVector (sintheta * cosphi, sintheta * sinphi, costheta));
    }
    particleGun->GeneratePrimaryVertex(anEvent) ;
    break ;
//..................................................................
    case planeGun:

    if(fPosition)
    {
      direction = G4ThreeVector( -fGunPosition.x()/fPosition ,
                                 -fGunPosition.y()/fPosition ,
                                 -fGunPosition.z()/fPosition    ) ;
      particleGun->SetParticleMomentumDirection(direction) ;

      position = G4ThreeVector(CLHEP::RandFlat::shoot(-0.5*fSize,0.5*fSize),
                               CLHEP::RandFlat::shoot(-0.5*fSize,0.5*fSize),
                               0.0) ;
      position.rotateUz(direction) ;
      particleGun->SetParticlePosition(fGunPosition+position) ;
    }
    else
    {
      G4Exception("Tst01PrimaryGeneratorAction::GeneratePrimaries",
                  "InvalidSetup", FatalException, "Invalid setting for gun");
    }
    particleGun->GeneratePrimaryVertex(anEvent);
    break ;
  }
}


//
//
/////////////////////////////////////////////////////////////////////////
