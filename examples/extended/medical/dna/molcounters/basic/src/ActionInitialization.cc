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
/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

// The `molcounters` example(s) are provided as part of Geant4-DNA
// and any report or published result obtained using it shall cite
// the respective Geant4-DNA collaboration publications.
//
// Reports or results obtained using the spatially-aware `MoleculeCounter`
// provided in this example, shall further cite:
//
// Velten & Tomé, Radiation Physics and Chemistry, 2023 (10.1016/j.radphyschem.2023.111194)
//
//
// Author: Christian Velten (2025)
//

#include "ActionInitialization.hh"

#include "ChemistrySteppingAction.hh"
#include "ChemistryTrackingManager.hh"
#include "EventAction.hh"
#include "MoleculeCounter.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "StackingAction.hh"
#include "G4SystemOfUnits.hh"

#include "G4DNAChemistryManager.hh"
#include "G4H2O.hh"
#include "G4MoleculeCounter.hh"
#include "G4MoleculeCounterManager.hh"
#include "G4MoleculeReactionCounter.hh"
#include "G4Scheduler.hh"
#include "G4UserTimeStepAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction);
  if (!fBuildMultipleAndCustomMoleculeCounters)
    BuildMoleculeCounters();
  else
    BuildMultipleAndCustomMoleculeCounters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction());

  SetUserAction(new RunAction);
  SetUserAction(new EventAction);
  SetUserAction(new StackingAction);

  if (G4DNAChemistryManager::IsActivated()) {
    G4Scheduler::Instance()->SetUserAction(new G4UserTimeStepAction);

    auto chemTrackingManager = new ChemistryTrackingManager();
    chemTrackingManager->SetUserAction(new ChemistrySteppingAction);
    G4Scheduler::Instance()->SetInteractivity(chemTrackingManager);

    if (!fBuildMultipleAndCustomMoleculeCounters)
      BuildMoleculeCounters();
    else
      BuildMultipleAndCustomMoleculeCounters();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ActionInitialization::BuildMoleculeCounters() const
{
  G4MoleculeCounterManager::Instance()->SetResetCountersBeforeEvent(false);  // defaults to false
  G4MoleculeCounterManager::Instance()->SetResetCountersBeforeRun(true);  // defaults to false

  // Register molecule counters
  {
    // Basic molecule counter
    auto counter = std::make_unique<G4MoleculeCounter>("Molecules");
    counter->IgnoreMolecule(G4H2O::Definition());
    G4MoleculeCounterManager::Instance()->RegisterCounter(std::move(counter));
  }
  {
    // Basic reaction counter
    auto counter = std::make_unique<G4MoleculeReactionCounter>("Reactions");
    G4MoleculeCounterManager::Instance()->RegisterCounter(std::move(counter));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ActionInitialization::BuildMultipleAndCustomMoleculeCounters() const
{
  G4MoleculeCounterManager::Instance()->SetResetCountersBeforeEvent(false);
  G4MoleculeCounterManager::Instance()->SetResetCountersBeforeRun(true);

  // Register molecule counters

  // Basic (built-in) Counters
  {
    // Basic molecule counter using a fixed time precision.
    // this will create many records {molecule -> {time -> conut}}
    auto counter = std::make_unique<G4MoleculeCounter>("BasicCounter");
    counter->IgnoreMolecule(G4H2O::Definition());
    counter->SetTimeComparer(G4MoleculeCounterTimeComparer::CreateWithFixedPrecision(25 * ps));
    G4MoleculeCounterManager::Instance()->RegisterCounter(std::move(counter));
  }
  {
    // Basic molecule counter using a fixed time precision [same as above].
    // However, we are activating it to only count molecules between 500 ps and 10 ns!
    auto counter = std::make_unique<G4MoleculeCounter>("BasicCounter_Restricted");
    counter->SetActiveLowerBound(500 * ps);
    counter->SetActiveUpperBound(10 * ns);  // add option to truncate after time?
    counter->IgnoreMolecule(G4H2O::Definition());
    counter->SetTimeComparer(G4MoleculeCounterTimeComparer::CreateWithFixedPrecision(25 * ps));
    G4MoleculeCounterManager::Instance()->RegisterCounter(std::move(counter));
  }
  {
    // Basic molecule counter using variable time precision without time restrition.
    // The precision is changed with respect to chemistry time.
    auto counter = std::make_unique<G4MoleculeCounter>("BasicCounter_VariablePrecision");
    counter->SetActiveLowerBound(500 * ps);
    counter->SetActiveUpperBound(10 * ns);  // add option to truncate after time?
    counter->IgnoreMolecule(G4H2O::Definition());
    counter->SetTimeComparer(G4MoleculeCounterTimeComparer::CreateWithVariablePrecision({
      {10 * ps, 5 * ps},
      {100 * ps, 50 * ps},
      {1000 * ps, 500 * ps},
      {10 * ns, 5 * ns},
      {1 * microsecond, 50 * ns},
    }));
    G4MoleculeCounterManager::Instance()->RegisterCounter(std::move(counter));
  }

  // Custom molecule counter, see: 10.1016/j.radphyschem.2023.111194
  {
    // Here we create a custom volume-aware molecule counter that uses variable time precision.
    // This counter records not just the molecules but also the encompassing volume.
    // Note:
    // - For this counter you must set `SetSensitiveToStepping(true)`.
    // - Other options (like SetNegativeCountsAreFatal) are optional but recommended
    // - SetIgnoreMoleculePosition can be set to true to override this counter's volume-awareness
    //   this would, in effect, make it the same as a BasicCounter `G4MoleculeCounter`
    auto counter = std::make_unique<MoleculeCounter>("MoleculeCounter");
    counter->SetVerbose(1);
    counter->IgnoreMolecule(G4H2O::Definition());
    counter->SetSensitiveToStepping(true);
    counter->SetIgnoreMoleculePosition(false);
    counter->SetNegativeCountsAreFatal(true);
    counter->SetTimeComparer(G4MoleculeCounterTimeComparer::CreateWithVariablePrecision({
      {10 * ps, 5 * ps},
      {100 * ps, 50 * ps},
      {1000 * ps, 500 * ps},
      {10 * ns, 5 * ns},
      {1 * microsecond, 50 * ns},
    }));
    G4MoleculeCounterManager::Instance()->RegisterCounter(std::move(counter));
  }

  // Molecule Reaction Counter with fixed time precision.
  // Set to be active from 0ps to 1us.
  {
    auto counter = std::make_unique<G4MoleculeReactionCounter>("Reactions");
    counter->SetActiveLowerBound(0 * ps);
    counter->SetActiveUpperBound(1 * microsecond);
    counter->SetTimeComparer(G4MoleculeCounterTimeComparer::CreateWithFixedPrecision(50 * ps));
    G4MoleculeCounterManager::Instance()->RegisterCounter(std::move(counter));
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
