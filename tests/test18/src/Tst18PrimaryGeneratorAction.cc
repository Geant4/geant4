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
// $Id$
//
//  File:        Tst18PrimaryGeneratorAction.cc
//  Description: Generator action for radioactive decay system test 
//  Author:      F. Lei (DERA UK)
//                  updated by Dennis Wright (SLAC)
//  Date:        14 August 2013
//

#include "Tst18PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleMomentum.hh"
#include "globals.hh"
#include "RadioactiveDecayGun.hh"


Tst18PrimaryGeneratorAction::Tst18PrimaryGeneratorAction()
{
  theParticleGun = new RadioactiveDecayGun();
}

Tst18PrimaryGeneratorAction::~Tst18PrimaryGeneratorAction()
{
  delete theParticleGun;
}

void Tst18PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleMomentum pm(1.0, 0.0, 0.0);

  G4int i = anEvent->GetEventID() % 3;
  switch(i)
  {
    case 0:
      break;
    case 1:
      pm.setX(0.0);
      pm.setY(1.0);
      break;
    case 2:
      pm.setY(0.0);
      pm.setZ(1.0);
      break;
  }

  theParticleGun->SetParticleMomentumDirection(pm);
  theParticleGun->GeneratePrimaryVertex(anEvent);
}

