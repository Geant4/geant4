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
/// \file scavenger/src/ScoreSpecies.cc
/// \brief Implementation of the scavenger::ScoreSpecies class

#include "ScoreSpecies.hh"

#include <G4MolecularConfiguration.hh>
#include <G4MoleculeCounter.hh>
#include <G4SystemOfUnits.hh>
#include <G4EventManager.hh>
#include <globals.hh>
#include "G4UnitsTable.hh"
#include "G4Scheduler.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4UImessenger.hh"

namespace scavenger
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreSpecies::ScoreSpecies(const G4String &name, const G4int &depth)
  : G4VPrimitiveScorer(name, depth),
    G4UImessenger(),
    fpSpeciesdir(new G4UIdirectory("/scorer/species/")),
    fpTimeBincmd(new G4UIcmdWithAnInteger("/scorer/species/nOfTimeBins", this)),
    fpAddTimeToRecordcmd(
      new G4UIcmdWithADoubleAndUnit("/scorer/species/addTimeToRecord", this)),
    fpSetResultsFileNameCmd(
      new G4UIcmdWithAString("/scorer/species/setRootFileName", this)) {
  fpSpeciesdir->SetGuidance("ScoreSpecies commands");
  fpSetResultsFileNameCmd->SetGuidance("Set the root file name");
  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::SetNewValue(G4UIcommand *command, G4String newValue) {
  if (command == fpAddTimeToRecordcmd.get()) {
    G4double cmdTime = fpAddTimeToRecordcmd->GetNewDoubleValue(newValue);
    AddTimeToRecord(cmdTime);
  }
  if (command == fpTimeBincmd.get()) {
    ClearTimeToRecord();
    G4int cmdBins = fpTimeBincmd->GetNewIntValue(newValue);
    G4double timeMin = 1 * ps;
    G4double timeMax = G4Scheduler::Instance()->GetEndTime() - timeMin;
    G4double timeLogMin = std::log10(timeMin);
    G4double timeLogMax = std::log10(timeMax);
    for (G4int i = 0; i < cmdBins; i++) {
      AddTimeToRecord(std::pow(10, timeLogMin +
                                   i * (timeLogMax - timeLogMin) / (cmdBins - 1)));
    }
  }
  if (command == fpSetResultsFileNameCmd.get()) {
    fRootFileName = newValue;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool ScoreSpecies::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep == 0.) {
    return FALSE;
  }
  edep *= aStep->GetPreStepPoint()->GetWeight(); // (Particle Weight)
  auto index = GetIndex(aStep);
  fEvtMap->add(index, edep);
  fEdep += edep;
  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::Initialize(G4HCofThisEvent *HCE) {
  fEvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
                                     GetName());
  if (fHCID < 0) {
    fHCID = GetCollectionID(0);
  }
  HCE->AddHitsCollection(fHCID, (G4VHitsCollection *) fEvtMap);
  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::EndOfEvent(G4HCofThisEvent *) {
  if (G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted()) {
    fEdep = 0.;
    G4MoleculeCounter::Instance()->ResetCounter();
    return;
  }

  auto species = G4MoleculeCounter::Instance()->GetRecordedMolecules();

  if (species == nullptr || species->empty()) {
    G4cout << "No molecule recorded, energy deposited= "
           << G4BestUnit(fEdep, "Energy") << G4endl;
    ++fNEvent;
    fEdep = 0.;
    G4MoleculeCounter::Instance()->ResetCounter();
    return;
  }
  for (auto molecule: *species) {
    for (auto time_mol: fTimeToRecord) {
      G4double n_mol =
        G4MoleculeCounter::Instance()->GetNMoleculesAtTime(molecule,
                                                           time_mol);
      if (n_mol < 0) {
        G4cerr << "N molecules not valid < 0 " << G4endl;
        G4Exception("", "N<0", FatalException, "");
      }
      SpeciesInfo &molInfo = fSpeciesInfoPerTime[time_mol][molecule];
      molInfo.fNumber += n_mol;
      G4double gValue = (n_mol / (fEdep / eV)) * 100.;
      molInfo.fG += gValue;
      molInfo.fG2 += gValue * gValue;
    }
  }
  ++fNEvent;
  fEdep = 0.;
  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
ScoreSpecies::AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer *workerScorer) {
  auto right = dynamic_cast<ScoreSpecies *>(dynamic_cast<G4VPrimitiveScorer *>(workerScorer));

  if (right == nullptr) {
    return;
  }
  if (right == this) {
    return;
  }

  auto it_map1 = right->fSpeciesInfoPerTime.begin();
  auto end_map1 = right->fSpeciesInfoPerTime.end();

  for (; it_map1 != end_map1; ++it_map1) {
    InnerSpeciesMap &map2 = it_map1->second;
    auto it_map2 = map2.begin();
    auto end_map2 = map2.end();

    for (; it_map2 != end_map2; ++it_map2) {
      SpeciesInfo &molInfo =
        fSpeciesInfoPerTime[it_map1->first][it_map2->first];
      molInfo.fNumber += it_map2->second.fNumber;
      molInfo.fG += it_map2->second.fG;
      molInfo.fG2 += it_map2->second.fG2;
    }
  }
  right->fSpeciesInfoPerTime.clear();
  fNEvent += right->fNEvent;
  right->fNEvent = 0;
  right->fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::PrintAll() {
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of events " << fNEvent << G4endl;
  G4cout << " Number of energy deposition recorded "
         << fEvtMap->entries() << G4endl;

  for (auto itr : *fEvtMap->GetMap()) {
    G4cout << "  copy no.: " << itr.first
           << "  energy deposit: "
           << *(itr.second) / GetUnitValue()
           << " [" << GetUnit() << "]"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::OutputAndClear() {
  if (G4Threading::IsWorkerThread()) {
    return;
  }

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType(fOutputType);
  this->WriteWithAnalysisManager(analysisManager);
  fNEvent = 0;
  fSpeciesInfoPerTime.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
ScoreSpecies::WriteWithAnalysisManager(G4VAnalysisManager *analysisManager) {
  analysisManager->OpenFile(fRootFileName);
  G4int fNtupleID = analysisManager->CreateNtuple("species", "species");
  analysisManager->CreateNtupleIColumn(fNtupleID, "speciesID");
  analysisManager->CreateNtupleIColumn(fNtupleID, "number");
  analysisManager->CreateNtupleIColumn(fNtupleID, "nEvent");
  analysisManager->CreateNtupleSColumn(fNtupleID, "speciesName");
  analysisManager->CreateNtupleDColumn(fNtupleID, "time");
  analysisManager->CreateNtupleDColumn(fNtupleID, "sumG");
  analysisManager->CreateNtupleDColumn(fNtupleID, "sumG2");
  analysisManager->FinishNtuple(fNtupleID);

  for (const auto &it_map1: fSpeciesInfoPerTime) {
    const InnerSpeciesMap &map2 = it_map1.second;

    for (const auto &it_map2 : map2) {
      G4double time = it_map1.first;
      auto species = it_map2.first;
      const G4String &name = species->GetName();
      auto molID = it_map2.first->GetMoleculeID();
      auto number = it_map2.second.fNumber;
      auto G = it_map2.second.fG;
      auto G2 = it_map2.second.fG2;
      analysisManager->FillNtupleIColumn(fNtupleID, 0, molID);   // MolID
      analysisManager->FillNtupleIColumn(fNtupleID, 1, number);  // Number
      analysisManager->FillNtupleIColumn(fNtupleID, 2, fNEvent); // Total nb events
      analysisManager->FillNtupleSColumn(fNtupleID, 3, name);    // molName
      analysisManager->FillNtupleDColumn(fNtupleID, 4, time);    // time
      analysisManager->FillNtupleDColumn(fNtupleID, 5, G);       // G
      analysisManager->FillNtupleDColumn(fNtupleID, 6, G2);      // G2
      analysisManager->AddNtupleRow(fNtupleID);
    }
  }

  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}