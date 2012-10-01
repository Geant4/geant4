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
#include "RadioactiveDecayGun.hh"
#include "RadioactiveDecayGunmessenger.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"


#include "G4UImanager.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
////////////////////////////////////////////////////////////////////////////////
//
RadioactiveDecayGun::RadioactiveDecayGun ()
{
  //  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  theRadioactiveDecayGunMessenger = new RadioactiveDecayGunmessenger(this);

  // G4IonTable *theIonTable = G4ParticleTable::GetParticleTable()->GetIonTable();
  //G4ParticleDefinition *aIon = NULL;

  //aIon = theIonTable->GetIon (52, 109, 0.0);

  //SetParticleDefinition(aIon);

  //G4UImanager* UI = G4UImanager::GetUIpointer();
  //      UI->ApplyCommand("/gun/particle proton");
  //UI->ApplyCommand("/gun/energy 0.00000000 MeV");
  //UI->ApplyCommand("/gun/position 0.0 0.0 0.0 m");

}
////////////////////////////////////////////////////////////////////////////////
//
RadioactiveDecayGun::~RadioactiveDecayGun()
{
   delete theRadioactiveDecayGunMessenger;
}

void RadioactiveDecayGun::SetNucleus (Nucleus theIon1)
{
  theIon = theIon1;

  //  G4IonTable *theIonTable =
  //   const_cast<G4IonTable* const>(G4ParticleTable::GetParticleTable()->GetIonTable());
  G4IonTable *theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());

  G4ParticleDefinition *aIon = NULL;

  G4int A = theIon.GetA();
  G4int Z = theIon.GetZ();
  G4double E = theIon.GetE();

  aIon = theIonTable->GetIon (Z, A, E);

  SetParticleDefinition(aIon);

  G4UImanager* UI = G4UImanager::GetUIpointer();
  //      UI->ApplyCommand("/gun/particle proton");
  UI->ApplyCommand("/gun/energy 0.000001 MeV");
  UI->ApplyCommand("/gun/position 0.0 0.0 0.0 m");

}

