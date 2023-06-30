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

#include "ActionInitialization.hh"
#include "ChemistryWorld.hh"
#include "DetectorConstruction.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAEventScheduler.hh"
#include "G4DNAScavengerMaterial.hh"
#include "G4H2O.hh"
#include "G4Molecule.hh"
#include "G4MoleculeCounter.hh"
#include "G4MoleculeGun.hh"
#include "G4Scheduler.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "StackingAction.hh"
#include "TimeStepAction.hh"
#include "PulseAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ActionInitialization::ActionInitialization(DetectorConstruction *pDetector)
    : G4VUserActionInitialization(), fpDetector(pDetector) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ActionInitialization::BuildForMaster() const {
  SetUserAction(new RunAction());
  G4DNAChemistryManager::Instance()->ResetCounterWhenRunEnds(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ActionInitialization::Build() const {
  G4MoleculeCounter::Instance()->Use(true);
  G4MoleculeCounter::Instance()->DontRegister(G4H2O::Definition());
  G4MoleculeCounter::Instance()->SetVerbose(0);
  G4MoleculeCounter::Instance()->CheckTimeForConsistency(false);
  auto pPulseAction = new PulseAction();
  SetUserAction(pPulseAction);
  SetUserAction(new PrimaryGeneratorAction(fpDetector));
  auto pRunAction = new RunAction();
  SetUserAction(pRunAction);
  SetUserAction(new StackingAction());
  auto pChemWorld = fpDetector->GetChemistryWorld();
  auto pScavenger =
      std::make_unique<G4DNAScavengerMaterial>(pChemWorld);
  // To counter Scavenger
  dynamic_cast<G4DNAScavengerMaterial *>(pScavenger.get())
      ->SetCounterAgainstTime();
  G4Scheduler::Instance()->SetScavengerMaterial(std::move(pScavenger));
  auto timeStepAction = new TimeStepAction(pChemWorld, pPulseAction);
  auto eventScheduler = timeStepAction->GetEventScheduler();
  pRunAction->SetEventScheduler(eventScheduler);
  G4Scheduler::Instance()->SetUserAction(timeStepAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
