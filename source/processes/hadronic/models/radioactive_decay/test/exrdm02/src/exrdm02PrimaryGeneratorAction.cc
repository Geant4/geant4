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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: exrdm02PrimaryGeneratorAction.cc,v 1.1 2003-10-08 16:31:51 hpw Exp $
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






