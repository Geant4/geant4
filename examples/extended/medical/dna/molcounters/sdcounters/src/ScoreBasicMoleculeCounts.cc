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

#include "ScoreBasicMoleculeCounts.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeCounterManager.hh"
#include "G4Scheduler.hh"
#include "G4SystemOfUnits.hh"
#include "G4TScoreNtupleWriter.hh"
#include "G4UImessenger.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ScoreBasicMoleculeCounts::ScoreBasicMoleculeCounts(G4String name, G4int depth,
                                                   G4String moleculeCounterName)
  : G4VPrimitiveScorer(name, depth),
    G4UImessenger(),
    fMoleculeCounterName(moleculeCounterName),
    fRunID(),
    fNbOfScoredEvents(),
    fTimesToRecord(),
    fMoleculeCountPerIndexPerTime()
{
  fAddTimeToRecordcmd = new G4UIcmdWithADoubleAndUnit(ScorerCommand("addTimeToRecord"), this);
  fTimeBincmd = new G4UIcmdWithAnInteger(ScorerCommand("nOfTimeBins"), this);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const char* ScoreBasicMoleculeCounts::ScorerCommand(G4String command)
{
  command = "/scorer/" + GetName() + "/" + command;
  return command.c_str();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ScoreBasicMoleculeCounts::~ScoreBasicMoleculeCounts()
{
  delete fAddTimeToRecordcmd;//can be smart ptr ?
  delete fTimeBincmd;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ScoreBasicMoleculeCounts::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fAddTimeToRecordcmd) {
    G4double cmdTime = fAddTimeToRecordcmd->GetNewDoubleValue(newValue);
    AddTimeToRecord(cmdTime);
  }
  if (command == fTimeBincmd) {
    ClearTimesToRecord();
    G4int cmdBins = fTimeBincmd->GetNewIntValue(newValue);
    G4double timeMin = 1 * ps;
    G4double timeMax = G4Scheduler::Instance()->GetEndTime();
    G4double timeLogMin = std::log10(timeMin);
    G4double timeLogMax = std::log10(timeMax);
    for (G4int i = 0; i < cmdBins; i++) {
      AddTimeToRecord(std::pow(10, timeLogMin + i * (timeLogMax - timeLogMin) / (cmdBins - 1)));
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool ScoreBasicMoleculeCounts::ProcessHits(G4Step*, G4TouchableHistory*)
{
  return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ScoreBasicMoleculeCounts::EndOfEvent(G4HCofThisEvent*)
{
  if (G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted()) return;

  auto counters = G4MoleculeCounterManager::Instance()->GetMoleculeCounters(fMoleculeCounterName);

  G4MoleculeCounter* counter = nullptr;
  if (counters.rbegin() != counters.rend()) {
    counter =
      const_cast<G4MoleculeCounter*>(dynamic_cast<const G4MoleculeCounter*>(*counters.begin()));
    if (counter == nullptr)
      G4Exception("ScoreMoleculeCounter::EndOfEvent", "SCOREMOLCOUNT", FatalException,
                  "Molecule counter has wrong type!");
  }
  else {
    G4Exception("ScoreMoleculeCounter::EndOfEvent", "SCOREMOLCOUNT", FatalException,
                "No molecule counter with given name found!");
  }

  for (auto const& entry : counter->GetCounterMap()) {
    auto key = entry.first.GetMolecule()->GetName();
    for (auto const& time : fTimesToRecord) {
      G4int nbOfMoleculesAtTime = counter->GetNbMoleculesAtTime(entry.first, time);
      fMoleculeCountPerIndexPerTime[key][time] += nbOfMoleculesAtTime;
    }
  }

  ++fNbOfScoredEvents;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ScoreBasicMoleculeCounts::AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer* workerScorer)
{
  ScoreBasicMoleculeCounts* worker = dynamic_cast<ScoreBasicMoleculeCounts*>(workerScorer);
  if (worker == nullptr || worker == this) return;

  fNbOfScoredEvents += worker->fNbOfScoredEvents;

  for (auto const& workerEntry : worker->fMoleculeCountPerIndexPerTime) {
    auto emplacedPair =
      fMoleculeCountPerIndexPerTime.emplace(workerEntry.first, workerEntry.second);
    if (!emplacedPair.second) {  // workerEntry.first already exists, emplacedPair.first == it
      auto it = emplacedPair.first;

      for (auto const& workerMoleculeCountWithTimeIter : workerEntry.second) {
        auto masterMoleculeCountWithTimeIter =
          it->second.find(workerMoleculeCountWithTimeIter.first);
        masterMoleculeCountWithTimeIter->second += workerMoleculeCountWithTimeIter.second;
      }
    }
  }

  worker->clear();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ScoreBasicMoleculeCounts::clear()
{
  G4VPrimitiveScorer::clear();
  fNbOfScoredEvents = 0;
  fMoleculeCountPerIndexPerTime.clear();
}

void ScoreBasicMoleculeCounts::OutputAndClear()
{
  if (G4Threading::IsWorkerThread()) return;

  G4cout << "\n\nWriting Scorers: " << GetName() << G4endl;

  WriteWithAnalysisManager();

  ++fRunID;
  clear();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ScoreBasicMoleculeCounts::WriteWithAnalysisManager()
{
  auto analysisManager = G4AnalysisManager::Instance();

  int fNtupleID = analysisManager->CreateNtuple("BasicMoleculeCount_" + fMoleculeCounterName,
                                                "BasicMoleculeCount_" + fMoleculeCounterName);

  analysisManager->CreateNtupleDColumn("Time__ps_");
  analysisManager->CreateNtupleDColumn("Molecule_Count");
  analysisManager->CreateNtupleSColumn("Molecule_Name");
  analysisManager->FinishNtuple();

  for (auto const& entry : fMoleculeCountPerIndexPerTime) {
    G4cout << entry.first << ":" << G4endl;

    for (auto const& moleculeCountWithTime : entry.second) {
      auto moleculeCount = moleculeCountWithTime.second;
      //auto moleculeCountPerEvent = moleculeCount / fNbOfScoredEvents;//not used

      auto time = moleculeCountWithTime.first;
      auto moleculeName = entry.first;

      G4cout << "  t=" << G4BestUnit(time, "Time") << ", n=" << moleculeCount << G4endl;

      analysisManager->FillNtupleDColumn(fNtupleID, 0, time / ps);  // time
      analysisManager->FillNtupleDColumn(fNtupleID, 1, moleculeCount);  // Number
      analysisManager->FillNtupleSColumn(fNtupleID, 2, moleculeName);  // molecule name
      analysisManager->AddNtupleRow(fNtupleID);
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......