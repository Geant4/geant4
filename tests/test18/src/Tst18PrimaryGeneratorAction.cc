// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst18PrimaryGeneratorAction.cc,v 1.5 2000-06-14 17:48:13 flei Exp $
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






