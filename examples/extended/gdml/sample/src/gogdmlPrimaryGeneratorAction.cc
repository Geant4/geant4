// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: gogdmlPrimaryGeneratorAction.cc,v 1.1.1.1 2002-05-31 00:34:43 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "gogdmlPrimaryGeneratorAction.hh"

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
  UI->ApplyCommand("/gun/energy 1.0 GeV");
  UI->ApplyCommand("/gun/position -2.0 0.0 0.0 cm");
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
      UI->ApplyCommand("/gun/direction 1.0 0.0 0.0");
      break;
    case 1:
      UI->ApplyCommand("/gun/direction 1.0 0.1 0.0");
      break;
    case 2:
      UI->ApplyCommand("/gun/direction 1.0 0.0 0.1");
      break;
  }

  particleGun->GeneratePrimaryVertex(anEvent);
}


