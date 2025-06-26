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
/// \file scavenger/src/ActionInitialization.cc
/// \brief Implementation of the scavenger::ActionInitialization class

#include "ActionInitialization.hh"

#include "G4DNAChemistryManager.hh"
#include "G4H2O.hh"
#include "G4MoleculeCounter.hh"
#include "G4Scheduler.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "StackingAction.hh"
#include "TimeStepAction.hh"

namespace scavenger
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ActionInitialization::ActionInitialization() : G4VUserActionInitialization() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction());

  BuildMoleculeCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction());
  SetUserAction(new RunAction());
  SetUserAction(new EventAction);
  SetUserAction(new StackingAction());
  G4Scheduler::Instance()->SetUserAction(new TimeStepAction());

  BuildMoleculeCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildMoleculeCounter() const
{
  G4MoleculeCounterManager::Instance()->SetResetCountersBeforeEvent(true);
  G4MoleculeCounterManager::Instance()->SetResetCountersBeforeRun(true);
  G4MoleculeCounterManager::Instance()->SetAccumulateCounterIntoMaster(false);

  auto counter = std::make_unique<G4MoleculeCounter>();
  counter->SetVerbose(0);
  counter->SetTimeComparer(G4MoleculeCounterTimeComparer::CreateWithFixedPrecision(1 * ps));
  counter->IgnoreMolecule(G4H2O::Definition());
  G4MoleculeCounterManager::Instance()->RegisterCounter(std::move(counter));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

}  // namespace scavenger
