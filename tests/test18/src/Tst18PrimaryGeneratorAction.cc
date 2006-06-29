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
// $Id: Tst18PrimaryGeneratorAction.cc,v 1.7 2006-06-29 21:44:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst18PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4UImanager.hh"
#include "globals.hh"

#include "RadioactiveDecayGun.hh"


Tst18PrimaryGeneratorAction::Tst18PrimaryGeneratorAction()
{
  theParticleGun = new RadioactiveDecayGun();
}

Tst18PrimaryGeneratorAction::~Tst18PrimaryGeneratorAction()
{
  delete theParticleGun;
//  delete theIonTable;
}

void Tst18PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    G4UImanager* UI = G4UImanager::GetUIpointer();
  G4int i = anEvent->GetEventID() % 3;
  switch(i)
  {
    case 0:
      UI->ApplyCommand("/gun/direction 1.0 0.0 0.0");
      break;
    case 1:
      UI->ApplyCommand("/gun/direction 0.0 1.0 0.0");
      break;
    case 2:
      UI->ApplyCommand("/gun/direction 0.0 0.0 1.0");
      break;
  }

  theParticleGun->GeneratePrimaryVertex(anEvent);
}






