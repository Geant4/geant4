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
// $Id: exrdm02PrimaryGeneratorAction.cc,v 1.2 2006-12-13 15:48:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "exrdm02PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4UImanager.hh"
#include "globals.hh"

#include "RadioactiveDecayGun.hh"


exrdm02PrimaryGeneratorAction::exrdm02PrimaryGeneratorAction()
{
  theParticleGun = new RadioactiveDecayGun();
}

exrdm02PrimaryGeneratorAction::~exrdm02PrimaryGeneratorAction()
{
  delete theParticleGun;
//  delete theIonTable;
}

void exrdm02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
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






