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
/// \file ScoreSpecies.cc
/// \brief Implementation of the ScoreSpecies class

// This example is provided by the Geant4-DNA collaboration
// chem6 example is derived from chem4 and chem5 examples
//
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// J. Appl. Phys. 125 (2019) 104301
// Med. Phys. 45 (2018) e722-e739
// J. Comput. Phys. 274 (2014) 841-882
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157-178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: W. G. Shin and S. Incerti (CENBG, France)

#include "ScoreSpecies.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4Scheduler.hh"
#include "G4TScoreNtupleWriter.hh"
#include "G4UImessenger.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <G4EventManager.hh>
#include <G4MolecularConfiguration.hh>
#include <G4MoleculeCounter.hh>
#include <globals.hh>
/**
 \file ScoreSpecies.cc
 \class ScoreSpecies
  This is a primitive scorer class for molecular species.
  The number of species is recorded for all times (predetermined or
  user chosen). It also scores the energy deposition in order to compute the
  radiochemical yields.
*/
extern std::ofstream out;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreSpecies::ScoreSpecies(const G4String &name, const G4int depth)
  : G4VPrimitiveScorer(name, depth), G4UImessenger() {
  fSpeciesdir = std::make_unique<G4UIdirectory>("/scorer/species/");
  fSpeciesdir->SetGuidance("ScoreSpecies commands");

  fAddTimeToRecordcmd = std::make_unique<G4UIcmdWithADoubleAndUnit>
      ("/scorer/species/addTimeToRecord", this);
  fTimeBincmd = std::make_unique<G4UIcmdWithAnInteger>("/scorer/species/nOfTimeBins", this);
  fNameCmd = std::make_unique<G4UIcmdWithAString>("/scorer/species/fileName", this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::SetNewValue(G4UIcommand *command, const G4String newValue) {
  if (command == fAddTimeToRecordcmd.get()) {
    const G4double cmdTime = fAddTimeToRecordcmd->GetNewDoubleValue(newValue);
    AddTimeToRecord(cmdTime);
  }
  if (command == fTimeBincmd.get()) {
    ClearTimeToRecord();
    const G4int cmdBins = fTimeBincmd->GetNewIntValue(newValue);
    constexpr G4double timeMin = 1 * ps;
    const G4double timeMax = G4Scheduler::Instance()->GetEndTime() - 1 * ps;
    if (cmdBins <= 1) {
      // Single bin (or invalid <=0) request: record only the latest time
      AddTimeToRecord(timeMax);
    } else {
      const G4double timeLogMin = std::log10(timeMin);
      const G4double timeLogMax = std::log10(timeMax);
      for (G4int i = 0; i < cmdBins; i++) {
        AddTimeToRecord(std::pow(10, timeLogMin + i * (timeLogMax - timeLogMin) / (cmdBins - 1)));
      }
    }
  }
  if (command == fNameCmd.get()) {
    fFileName = newValue;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool ScoreSpecies::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep == 0.) { return false; }

  edep *= aStep->GetPreStepPoint()->GetWeight(); // (Particle Weight)
  const G4int index = GetIndex(aStep);
  fEvtMap->add(index, edep);
  fEdep += edep;

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::Initialize(G4HCofThisEvent *HCE) {
  fEvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(), GetName());

  if (fHCID < 0) {
    fHCID = GetCollectionID(0);
  }

  HCE->AddHitsCollection(fHCID, fEvtMap);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::EndOfEvent(G4HCofThisEvent *) {
  if (G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted()) {
    fEdep = 0.;
    return;
  }

  // get the first, and in this case only, counter
  const auto counter = G4MoleculeCounterManager::Instance()->GetMoleculeCounter<
    G4MoleculeCounter>(0);
  if (counter == nullptr) {
    G4Exception("ScoreSpecies::EndOfEvent", "BAD_REFERENCE", FatalException,
                "The molecule counter could not be received!");
  }

  const auto indices = counter->GetMapIndices();

  if (indices.empty()) {
    G4cout << "No molecule recorded, energy deposited= " << G4BestUnit(fEdep, "Energy") << G4endl;
    ++fNEvent;
    fEdep = 0.;
    return;
  }
  for (const auto &idx: indices) {
    for (const auto time_mol: fTimeToRecord) {
      const G4double n_mol = counter->GetNbMoleculesAtTime(idx, time_mol);

      if (n_mol < 0) {
        G4cerr << "N molecules not valid < 0 " << G4endl;
        G4Exception("", "N<0", FatalException, "");
      }

      SpeciesInfo &molInfo = fSpeciesInfoPerTime[time_mol][idx.Molecule];
      molInfo.fNumber += n_mol;

      if (fEdep <= 0.) {
        // Avoid division by zero; skip G-value accumulation for this event
        continue;
      }
      const G4double gValue = (n_mol / (fEdep / eV)) * 100.;
      molInfo.fG += gValue;
      molInfo.fG2 += gValue * gValue;
    }
  }

  ++fNEvent;

  fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer *workerScorer) {
  auto *right = dynamic_cast<ScoreSpecies *>(workerScorer);
  if (right == nullptr || right == this) {
    return;
  }

  for (const auto &[timeKey, innerMap]: right->fSpeciesInfoPerTime) {
    for (const auto &[speciesPtr, info]: innerMap) {
      SpeciesInfo &molInfo = fSpeciesInfoPerTime[timeKey][speciesPtr];
      molInfo.fNumber += info.fNumber;
      molInfo.fG += info.fG;
      molInfo.fG2 += info.fG2;
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
  G4cout << " Number of energy deposition recorded " << fEvtMap->entries() << G4endl;

  for (auto [fst, snd]: *fEvtMap->GetMap()) {
    G4cout << "  copy no.: " << fst << "  energy deposit: " << *(snd) / GetUnitValue()
        << " [" << GetUnit() << "]" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::OutputAndClear() {
  if (G4Threading::IsWorkerThread()) { return; }

  //---------------------------------------------------------------------------
  // Save results

  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType(fOutputType);

  if (analysisManager) {
    this->WriteWithAnalysisManager(analysisManager);
  }

  fNEvent = 0;
  fSpeciesInfoPerTime.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::WriteWithAnalysisManager(G4VAnalysisManager *analysisManager) {
  const G4String fileN = fFileName + G4UIcommand::ConvertToString(fRunID);
  analysisManager->OpenFile(fileN);
  const int fNtupleID = analysisManager->CreateNtuple("species", "species");
  analysisManager->CreateNtupleIColumn(fNtupleID, "speciesID");
  analysisManager->CreateNtupleIColumn(fNtupleID, "number");
  analysisManager->CreateNtupleIColumn(fNtupleID, "nEvent");
  analysisManager->CreateNtupleSColumn(fNtupleID, "speciesName");
  analysisManager->CreateNtupleDColumn(fNtupleID, "time");
  analysisManager->CreateNtupleDColumn(fNtupleID, "sumG");
  analysisManager->CreateNtupleDColumn(fNtupleID, "sumG2");
  analysisManager->FinishNtuple(fNtupleID);

  for (auto it_map1: fSpeciesInfoPerTime) {
    InnerSpeciesMap &map2 = it_map1.second;

    if (it_map1.first == fSpeciesInfoPerTime.begin()->first) {
      for (auto it_molname: map2) {
        const auto tmp_species = it_molname.first;
        out << std::setw(12) << tmp_species->GetName() << std::setw(12)
            << tmp_species->GetMoleculeID();
      }
      out << '\n';
    }

    for (auto it_map2: map2) {
      const G4double time = it_map1.first;
      const auto species = it_map2.first;
      const G4String &name = species->GetName();
      const G4int molID = it_map2.first->GetMoleculeID();
      const G4int number = it_map2.second.fNumber;
      const G4double G = it_map2.second.fG;
      const G4double G2 = it_map2.second.fG2;
      const G4int N = fNEvent;

      if (time == *fTimeToRecord.rbegin()) {
        if (N > 1) {
          out << std::setw(12) << G / N << std::setw(12)
              << std::sqrt(((G2 / N) - std::pow(G / N, 2)) / (N - 1));
        } else {
          out << std::setw(12) << G / N << std::setw(12)
              << std::sqrt(((G2 / N) - std::pow(G / N, 2)) / N);
        }
      }

      analysisManager->FillNtupleIColumn(fNtupleID, 0, molID); // MolID
      analysisManager->FillNtupleIColumn(fNtupleID, 1, number); // Number
      analysisManager->FillNtupleIColumn(fNtupleID, 2, fNEvent); // Total nb events
      analysisManager->FillNtupleSColumn(fNtupleID, 3, name); // molName
      analysisManager->FillNtupleDColumn(fNtupleID, 4, time); // time
      analysisManager->FillNtupleDColumn(fNtupleID, 5, G); // G
      analysisManager->FillNtupleDColumn(fNtupleID, 6, G2); // G2
      analysisManager->AddNtupleRow(fNtupleID);
    }
  }

  analysisManager->Write();
  analysisManager->CloseFile();
  fRunID++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
