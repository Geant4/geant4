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
/// \file scavenger/src/RunAction.cc
/// \brief Implementation of the scavenger::RunAction class

#include "RunAction.hh"
#include "Run.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "ScoreSpecies.hh"

namespace scavenger
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

RunAction::RunAction()
  : G4UserRunAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4Run *RunAction::GenerateRun() {
  Run *run = new Run();
  return run;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::BeginOfRunAction(const G4Run *run) {
  G4cout << "### Run " << run->GetRunID() << " starts." << G4endl;

  // informs the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void RunAction::EndOfRunAction(const G4Run *run) {
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) {
    return;
  }
  const Run *scavengerRun = dynamic_cast<const Run *>(run);
  G4double sumDose = scavengerRun->GetSumDose();
  if (IsMaster()) {
    G4cout
      << G4endl
      << "--------------------End of Global Run-----------------------"
      << G4endl
      << "  The run has " << nofEvents << " events "
      << G4endl;
    auto masterScorer = dynamic_cast<ScoreSpecies *>(scavengerRun->GetPrimitiveScorer());
    G4cout << "Number of events recorded by the species scorer = "
           << masterScorer->GetNumberOfRecordedEvents()
           << G4endl;
    masterScorer->OutputAndClear();
  } else {
    G4cout
      << G4endl
      << "--------------------End of Local Run------------------------"
      << G4endl
      << "  The run has " << nofEvents << " events"
      << G4endl;
  }
  G4cout
    << " Total energy deposited in the world volume : " << sumDose / eV << " eV"
    << G4endl
    << " ------------------------------------------------------------"
    << G4endl
    << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

}