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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: B01PrimaryGeneratorAction.cc,v 1.2 2002-04-19 10:54:26 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"

#include "B01PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4UImanager.hh"

B01PrimaryGeneratorAction::B01PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/gun/particle neutron");
  UI->ApplyCommand("/gun/energy 10.0 MeV");
  UI->ApplyCommand("/gun/position 0.0 0.0 -15.0005 cm");
  UI->ApplyCommand("/gun/direction 0. 0. 1.");
}

B01PrimaryGeneratorAction::~B01PrimaryGeneratorAction()
{
  delete particleGun;
}

void B01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4int i = anEvent->GetEventID() % 3;
  particleGun->GeneratePrimaryVertex(anEvent);
}
