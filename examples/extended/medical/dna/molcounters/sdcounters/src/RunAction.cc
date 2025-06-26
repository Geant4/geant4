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

#include "RunAction.hh"

#include "Run.hh"
#include "ScoreBasicMoleculeCounts.hh"
#include "ScoreBasicReactionCounts.hh"

#include "G4AnalysisManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4RunManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Run* RunAction::GenerateRun()
{
  auto run = new Run();
  return run;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::BeginOfRunAction(const G4Run* run)
{
  // ensure that the chemistry is notified!
  if (G4DNAChemistryManager::GetInstanceIfExists() != nullptr)
    G4DNAChemistryManager::GetInstanceIfExists()->BeginOfRunAction(run);

  // do your own stuff for an event here

  G4cout << "### Run " << run->GetRunID() << " starts." << G4endl;
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::EndOfRunAction(const G4Run* run)
{
  // ensure that the chemistry is notified!
  if (G4DNAChemistryManager::GetInstanceIfExists() != nullptr)
    G4DNAChemistryManager::GetInstanceIfExists()->EndOfRunAction(run);

  // do your own stuff for an event here

  auto nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  if (IsMaster()) {
    G4cout << "End of Global Run totaling nEvents = " << nofEvents << G4endl;
    auto myRun = dynamic_cast<const Run*>(run);

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    G4String fileN = "MoleculeCounters";
    analysisManager->OpenFile(fileN);

    dynamic_cast<ScoreBasicMoleculeCounts*>(myRun->GetBasicMoleculeScorer())->OutputAndClear();
    dynamic_cast<ScoreBasicMoleculeCounts*>(myRun->GetBasicMoleculeScorerWithVariablePrecision())
      ->OutputAndClear();
    dynamic_cast<ScoreBasicReactionCounts*>(myRun->GetBasicReactionScorer())->OutputAndClear();

    analysisManager->Write();
    analysisManager->CloseFile();
  }
  else {
    G4cout << "End of Local Run with nEvents = " << nofEvents << G4endl;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......