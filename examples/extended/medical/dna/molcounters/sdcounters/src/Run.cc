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
/// \file Run.cc
/// \brief Implementation of the Run class

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

#include "Run.hh"

#include "ScoreBasicMoleculeCounts.hh"
#include "ScoreBasicReactionCounts.hh"

#include "G4MoleculeCounterManager.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Run::Run() : G4Run()
{
  auto mfdet = dynamic_cast<G4MultiFunctionalDetector*>(
    G4SDManager::GetSDMpointer()->FindSensitiveDetector("mfDetector"));

  fScorerMoleculesBasic = mfdet->GetPrimitive(
    G4SDManager::GetSDMpointer()->GetCollectionID("mfDetector/BasicMoleculeCounts"));
  fScorerMoleculesBasicVariablePrecision = mfdet->GetPrimitive(
    G4SDManager::GetSDMpointer()->GetCollectionID("mfDetector/BasicCounter_VariablePrecision"));
  fScorerReactionsBasic = mfdet->GetPrimitive(
    G4SDManager::GetSDMpointer()->GetCollectionID("mfDetector/BasicReactionCounts"));
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Run::Merge(const G4Run* aRun)
{
  if (aRun == nullptr || aRun == this) return;

  // Ideally, merging of workers -> master counters is handled by the counter manager
  // as part of its EndOfEventAction and EndOfRunAction, if and when the manager it instructed
  // to reset the counters between events or runs.
  // In the case of resets between runs, the user has two options:
  //   (1) ignore the master counter entries and merge results from sensitive detector workers; or
  //   (2) do not trigger sensitive detectors in the workers and read out from the master, only
  // If the counter manager does not accumulate, only option (1) is feasible.
  //
  // In the case of resets between events, and the counter manager does not accumulate, the user
  // must:
  //   (*) fill the scorers during the Event (or using G4VPrimitiveScorer::EndOfEvent); and
  //   (*) merge the worker scorers at the end of the run
  // This is the "current" behavior, for example, when calculating G-values.
  //
  // The code below accumulates the scorers' workers into the respective master
  // in case the manager does not accumulate, or we reset between events or runs.

  if (G4MoleculeCounterManager::Instance()->GetResetCountersBeforeEvent()
      || G4MoleculeCounterManager::Instance()->GetResetCountersBeforeRun()
      || !G4MoleculeCounterManager::Instance()->GetAccumulateCounterIntoMaster())
  {
    {
      auto masterScorer = dynamic_cast<ScoreBasicMoleculeCounts*>(this->fScorerMoleculesBasic);
      auto localScorer = dynamic_cast<ScoreBasicMoleculeCounts*>(
        static_cast<const Run*>(aRun)->fScorerMoleculesBasic);
      masterScorer->AbsorbResultsFromWorkerScorer(localScorer);
    }

    {
      auto masterScorer =
        dynamic_cast<ScoreBasicMoleculeCounts*>(this->fScorerMoleculesBasicVariablePrecision);
      auto localScorer = dynamic_cast<ScoreBasicMoleculeCounts*>(
        static_cast<const Run*>(aRun)->fScorerMoleculesBasicVariablePrecision);
      masterScorer->AbsorbResultsFromWorkerScorer(localScorer);
    }

    {
      auto masterScorer = dynamic_cast<ScoreBasicReactionCounts*>(this->fScorerReactionsBasic);
      auto localScorer = dynamic_cast<ScoreBasicReactionCounts*>(
        static_cast<const Run*>(aRun)->fScorerReactionsBasic);
      masterScorer->AbsorbResultsFromWorkerScorer(localScorer);
    }
  }

  G4Run::Merge(aRun);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......