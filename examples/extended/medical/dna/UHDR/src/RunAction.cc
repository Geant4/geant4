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
#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4DNAEventScheduler.hh"
#include "G4DNAScavengerMaterial.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Scheduler.hh"
#include "G4SystemOfUnits.hh"
#include "G4VChemistryWorld.hh"
#include "Run.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

RunAction::RunAction() : G4UserRunAction() {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::SetEventScheduler(G4DNAEventScheduler *pEventScheduler) {
  fpEventScheduler = pEventScheduler;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4Run *RunAction::GenerateRun() {
  Run *run = new Run();
  return run;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::BeginOfRunAction(const G4Run *run) {
  G4cout << "### Run " << run->GetRunID() << " starts." << G4endl;
  if (G4Threading::IsMultithreadedApplication() && IsMaster()) {
    return;
  }
  auto runInfo = dynamic_cast<const Run *>(run);
  auto GvaluesScorer =
      dynamic_cast<Scorer<Gvalues> *>(runInfo->GetGvaluesScorer());

  if (fpEventScheduler != nullptr) {
    for (const auto &it: GvaluesScorer->GetpScorer()->fTimeToRecord) {
      fpEventScheduler->AddTimeToRecord(it);
    }
  }
  GvaluesScorer->SetEventScheduler(fpEventScheduler);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::EndOfRunAction(const G4Run *run) {
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) {
    return;
  }
  auto chem6Run = dynamic_cast<const Run *>(run);
  G4double sumDose = chem6Run->GetSumDose();

  if (G4Threading::IsMultithreadedApplication() && IsMaster()) {
    G4cout << G4endl
           << "--------------------------End of Global Run------------------------------"
           << G4endl << "The run has " << nofEvents << " events " << G4endl;

    auto masterGvaluesScorer =
        dynamic_cast<Scorer<Gvalues> *>(chem6Run->GetGvaluesScorer());

    auto masterDoseLimite =
        dynamic_cast<Scorer<Dose> *>(chem6Run->GetSumDoseLimit());

    G4cout << "Number of events recorded by the species scorer : "
           << masterGvaluesScorer->GetpScorer()->GetNumberOfRecordedEvents()
           << " events " << G4endl;
    G4double dose_mean = masterDoseLimite->GetpScorer()->fDosesCutOff / gray;
    auto boundingBox =
        masterDoseLimite->GetChemistryWorld()->GetChemistryBoundary();
    G4double V = boundingBox->Volume() / cm3;
    G4double DoseInGray = (sumDose / eV) / (0.001 * V * 6.242e+18);
    G4cout << "Cut-off dose for each beam line : " << dose_mean << " Gy " << G4endl;
    G4cout << "Actual dose : " << DoseInGray << " Gy for " << nofEvents
           << " events. Actual average dose : " << DoseInGray / nofEvents << " Gy" << G4endl;
    masterGvaluesScorer->OutputAndClear(std::to_string(dose_mean));

  } else {
    G4cout << G4endl
           << "--------------------------End of Local Run------------------------------"
           << G4endl << "The run has " << nofEvents << " events. Scavenger info:" << G4endl;
    auto pScavengerMaterial = dynamic_cast<G4DNAScavengerMaterial *>(
        G4Scheduler::Instance()->GetScavengerMaterial());
    pScavengerMaterial->PrintInfo();
  }

  G4cout << "Total energy deposited in the world volume : " << sumDose / keV
         << " keV" << G4endl
         << "-------------------------------------------------------------------------"
         << G4endl << G4endl;
}