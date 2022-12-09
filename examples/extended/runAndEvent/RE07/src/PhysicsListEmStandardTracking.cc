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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsListEmStandardTracking.hh"

#include "G4Electron.hh"
#include "G4EmBuilder.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4Positron.hh"

#include "EmStandardPhysicsTrackingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListEmStandardTracking::PhysicsListEmStandardTracking(G4int ver)
  : G4VPhysicsConstructor("EmStandardTracking")
{
  SetVerboseLevel(ver);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListEmStandardTracking::~PhysicsListEmStandardTracking() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListEmStandardTracking::ConstructParticle()
{
  // Taken care of in PhysicsList.cc
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListEmStandardTracking::ConstructProcess()
{
  G4EmBuilder::PrepareEMPhysics();

  // Register custom tracking manager for e-/e+ and gammas. Note that tracking
  // managers take precendence over processes registered in the G4ProcessManager
  // and hides them, so additional physics constructors like G4EmExtraPhysics
  // will not work!
  auto* trackingManager = new EmStandardPhysicsTrackingManager;

  G4Electron::Definition()->SetTrackingManager(trackingManager);
  G4Positron::Definition()->SetTrackingManager(trackingManager);
  G4Gamma::Definition()->SetTrackingManager(trackingManager);
}
