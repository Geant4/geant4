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
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"

#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4EmStandardPhysics.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "PhysicsListEmSpecialized.hh"
#include "PhysicsListEmStandardTracking.hh"
#include "PhysicsListMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
{
  SetVerboseLevel(0);

  fMessenger.reset(new PhysicsListMessenger(this));
  // By default, use the standard physics list with processes.
  fEmPhysicsList.reset(new G4EmStandardPhysics(GetVerboseLevel()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  G4BaryonConstructor::ConstructParticle();
  G4BosonConstructor::ConstructParticle();
  G4IonConstructor::ConstructParticle();
  G4LeptonConstructor::ConstructParticle();
  G4MesonConstructor::ConstructParticle();
  G4ShortLivedConstructor::ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  fEmPhysicsList->ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetMode(const G4String& name)
{
  if (verboseLevel > -1) {
    G4cout << "PhysicsList::SetMode: <" << name << ">" << G4endl;
  }

  if (name == "processes") {
    fEmPhysicsList.reset(new G4EmStandardPhysics(GetVerboseLevel()));
  }
  else if (name == "tracking") {
    fEmPhysicsList.reset(new PhysicsListEmStandardTracking(GetVerboseLevel()));
  }
  else if (name == "specialized") {
    fEmPhysicsList.reset(new PhysicsListEmSpecialized(GetVerboseLevel()));
  }
  else {
    G4cout << "PhysicsList::SetMode: <" << name << ">"
           << " is not defined" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
