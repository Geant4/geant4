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
// $Id: g4rwPrimaryGeneratorAction.cc,v 1.4 2006/06/29 17:21:52 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//

#include "g4rwPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4UImanager.hh"
#include "globals.hh"

gogdmlPrimaryGeneratorAction::gogdmlPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/gun/particle geantino");
  UI->ApplyCommand("/gun/position 0.0 0.0 0.0 cm");
}

gogdmlPrimaryGeneratorAction::~gogdmlPrimaryGeneratorAction()
{
  delete particleGun;
}

void gogdmlPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4int i = anEvent->GetEventID() % 3;
  switch(i)
  {
    case 0:
      UI->ApplyCommand("/gun/direction 0.2 0.0 1.0");
      break;
    case 1:
      UI->ApplyCommand("/gun/direction 0.0 0.3 1.0");
      break;
    case 2:
      UI->ApplyCommand("/gun/direction 0.1 0.0 1.0");
      break;
  }

  particleGun->GeneratePrimaryVertex(anEvent);
}


