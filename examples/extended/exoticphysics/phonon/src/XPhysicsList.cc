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
/// \file exoticphysics/phonon/src/XPhysicsList.cc
/// \brief Implementation of the XPhysicsList class
//
//

#include "XPhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhononDownconversion.hh"
#include "G4PhononLong.hh"
#include "G4PhononReflection.hh"
#include "G4PhononScattering.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicsList::XPhysicsList(G4int verbose) : G4VUserPhysicsList() {
  if (verbose) G4cout << "XPhysicsList::constructor" << G4endl;

  SetVerboseLevel(verbose);
  SetDefaultCutValue(100*mm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XPhysicsList::~XPhysicsList() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XPhysicsList::ConstructParticle() {
  G4PhononLong::PhononDefinition();
  G4PhononTransFast::PhononDefinition();
  G4PhononTransSlow::PhononDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XPhysicsList::ConstructProcess() {
  AddTransportation();

  // Only make processes once
  G4VProcess* phScat = new G4PhononScattering;
  G4VProcess* phRefl = new G4PhononReflection;
  G4VProcess* phDown = new G4PhononDownconversion;

  // Set process verbosity to match physics list, for diagnostics
  phScat->SetVerboseLevel(verboseLevel);
  phRefl->SetVerboseLevel(verboseLevel);
  phDown->SetVerboseLevel(verboseLevel);

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    // WARNING!  CHANGING ORDER OF REGISTRATION CAN CHANGE PHYSICS RESULTS
    if (phScat->IsApplicable(*particle)) pmanager->AddDiscreteProcess(phScat);
    if (phDown->IsApplicable(*particle)) pmanager->AddDiscreteProcess(phDown);
    if (phRefl->IsApplicable(*particle)) pmanager->AddDiscreteProcess(phRefl);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XPhysicsList::SetCuts() {
  // These values are used as the default production thresholds
  // for the world volume.
  SetCutsWithDefault();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



