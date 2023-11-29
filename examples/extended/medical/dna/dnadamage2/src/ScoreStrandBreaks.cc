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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file ScoreStrandBreaks.cc
/// \brief Implementation of the ScoreStrandBreaks class

#include "ScoreStrandBreaks.hh"

#include <G4SystemOfUnits.hh>
#include <G4EventManager.hh>
#include <globals.hh>
#include "G4UnitsTable.hh"
#include "G4Scheduler.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4UImessenger.hh"
#include "G4EventManager.hh"
#include "G4DNAChemistryManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreStrandBreaks::ScoreStrandBreaks(G4String name, 
    DetectorConstruction* detect, G4double* radius)
: G4VPrimitiveScorer(name), 
  G4UImessenger(),
  fpDetector(detect)
{
  fpOutputFileUI = new G4UIcmdWithAString("/scorer/StrandBreak/OutputFile", this);
  fpOutputFileUI->SetGuidance("Set output file name");

  fpOutputTypeUI = new G4UIcmdWithAString("/scorer/StrandBreak/OutputFormat", this);
  fpOutputTypeUI->SetGuidance("Set output file format: ASCII by default");

  fpBreakEnergyUI = 
     new G4UIcmdWithADoubleAndUnit("/scorer/StrandBreak/BreakEnergy", this);
  fpBreakEnergyUI->SetDefaultUnit("eV");
  fpBreakEnergyUI->SetGuidance("Direct SSB energy treshold");

  fRadius = radius;
  fpGun = new MoleculeInserter(true);
  //G4DNAChemistryManager::Instance()->SetGun(fpGun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreStrandBreaks::~ScoreStrandBreaks()
{
  delete fpOutputFileUI;
  delete fpBreakEnergyUI;
  delete fpOutputTypeUI;
  delete fpGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fpOutputFileUI) {
    fOutputName = newValue;
  }

  if (command == fpBreakEnergyUI) {
    fBreakEnergy = fpBreakEnergyUI->GetNewDoubleValue(newValue);
  }

  if (command == fpOutputTypeUI) {
    fOutputType = newValue;
    G4StrUtil::to_lower(fOutputType);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool ScoreStrandBreaks::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep <= 0) { return FALSE; }
  fEnergyDeposit += edep;

  G4StepPoint* aStepPoint      = aStep->GetTrack()->GetStep()->GetPreStepPoint();
  G4TouchableHandle aTouchable = aStepPoint->GetTouchableHandle();
  G4String volumeName  = aTouchable->GetVolume()->GetName();
  G4StrUtil::to_lower(volumeName);
  G4int strand = 0;

  if (G4StrUtil::contains(volumeName, "deoxyribose") || 
      G4StrUtil::contains(volumeName, "phosphate")) {
    if (G4StrUtil::contains(volumeName,"1")) 
      strand = 1;
    else if (G4StrUtil::contains(volumeName,"2"))
      strand = 2;

    G4int StrandID = strand;
    G4int BaseID = aTouchable->GetReplicaNumber(0);
    G4int PlasmidID = aTouchable->GetReplicaNumber(1);
    G4int EventID = G4EventManager::GetEventManager()->
                    GetConstCurrentEvent()->GetEventID();

    fEnergyDepositMap[EventID][PlasmidID][BaseID][StrandID]+=edep;
    return TRUE;
  }

  if (!fDNAInserted) {
    // Insert DNA molecules
    std::vector<std::vector<G4int>> DNADetails = fpDetector->GetDNADetails();
    std::vector<G4ThreeVector> DNAPositions = fpDetector->GetDNAPositions();
    std::vector<G4String> DNAMolecules = fpDetector->GetDNANames();

    fpGun->Clean();
    for(size_t i = 0; i < DNAPositions.size(); i++) {
      fpGun->AddMolecule(DNAMolecules[i],DNAPositions[i],1.0*ps);
    }
    fpGun->DefineTracks();

    fDNAInserted = true;
  }
  return FALSE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::Initialize(G4HCofThisEvent*)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::EndOfEvent(G4HCofThisEvent*)
{
  // Get EventID
  G4int EventID   = G4EventManager::GetEventManager()->
                    GetConstCurrentEvent()->GetEventID();

  // Absorbed Dose
  G4double volume  = 4.0/3 * 3.14159 * (*fRadius) * (*fRadius) * (*fRadius) / 8;
  G4double density = 1.0 * g/cm3;
  G4double mass = density * volume;
  G4double Dose = fEnergyDeposit / mass;
  fDoseArray[EventID] = Dose / gray;

  // Direct Strand Breaks
  for( auto EventAndElse : fEnergyDepositMap) {
    G4int event = EventAndElse.first;
    for ( auto PlasmidAndElse : EventAndElse.second) {
      G4int plasmid = PlasmidAndElse.first;
      for ( auto BaseAndElse : PlasmidAndElse.second) {
        G4int base = BaseAndElse.first;
        for ( auto StrandAndEdep : BaseAndElse.second) {
          G4int strand  = StrandAndEdep.first;
          G4double edep = StrandAndEdep.second;
          if (edep > fBreakEnergy)
            fDirectDamageMap[event].push_back({plasmid,base,strand});
        }
      }
    }
  }
  fEnergyDepositMap.clear();

  // Indirect Strand Breaks
  std::vector<std::vector<G4int>> DNADetails = fpDetector->GetDNADetails();
  std::vector<G4Track*> InsertedTracks = fpGun->GetInsertedTracks();
  std::vector<std::vector<G4int>> IndirectBreaks;

  for(size_t i = 0; i < InsertedTracks.size(); i++) {
    if(InsertedTracks[i]->GetTrackStatus() == fStopAndKill)
      IndirectBreaks.push_back(DNADetails[i]);
  }

  for(size_t i = 0; i < IndirectBreaks.size(); i++) {
    fIndirectDamageMap[EventID].push_back(IndirectBreaks[i]);
  }

  fEnergyDeposit = 0;
  fDNAInserted   = false;
  fnbOfEvents++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
ScoreStrandBreaks::AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer* workerScorer)
{
  ScoreStrandBreaks* right =
  dynamic_cast<ScoreStrandBreaks*>(dynamic_cast<G4VPrimitiveScorer*>(workerScorer));

  if (right == 0)
    return;
  if (right == this)
    return;

  DamageMap WorkerDirectDamageMap     = right->fDirectDamageMap;
  DamageMap WorkerIndirectDamageMap   = right->fIndirectDamageMap;
  std::map<G4int,G4double> WorkerDose = right->fDoseArray;

  for ( auto EventAndBreaks : WorkerDirectDamageMap) {
    G4int EventID = EventAndBreaks.first;
    std::vector<std::vector<G4int>> Breaks = EventAndBreaks.second;
    for (size_t i = 0; i < Breaks.size(); i++) {
      fDirectDamageMap[EventID].push_back(Breaks[i]);
    }
  }

  for ( auto EventAndBreaks : WorkerIndirectDamageMap) {
    G4int EventID = EventAndBreaks.first;
    std::vector<std::vector<G4int>> Breaks = EventAndBreaks.second;
    for (size_t i = 0; i < Breaks.size(); i++) {
      fIndirectDamageMap[EventID].push_back(Breaks[i]);
    }
  }

  for ( auto EventAndDose : WorkerDose) {
    G4int EventID = EventAndDose.first;
    G4double dose = EventAndDose.second;
    fDoseArray[EventID] = dose;
  }

  fnbOfEvents += right->fnbOfEvents;
  right->Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::Clear()
{
  if (fDirectDamageMap.size() != 0)
    fDirectDamageMap.clear();

  if (fIndirectDamageMap.size() != 0)
    fIndirectDamageMap.clear();

  if (fEnergyDepositMap.size() != 0)
    fEnergyDepositMap.clear();

  if (fDoseArray.size() != 0)
    fDoseArray.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::DrawAll()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::PrintAll()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreStrandBreaks::OutputAndClear(G4double LET, G4double LET_STD)
{
  if(G4Threading::IsWorkerThread()) return;

  if (fOutputType != "ascii") {

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    if (fOutputType == "csv")
      analysisManager->SetDefaultFileType("csv");

    else if (fOutputType == "root")
      analysisManager->SetDefaultFileType("root");

    else if (fOutputType == "xml")
      analysisManager->SetDefaultFileType("xml");

    WriteWithAnalysisManager(analysisManager);
  }

  else
    ASCII(LET, LET_STD);

  Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
ScoreStrandBreaks::WriteWithAnalysisManager(G4VAnalysisManager* analysisManager)
{
  analysisManager->OpenFile(fOutputName);
  int fNtupleID = analysisManager->CreateNtuple("SB","Direct And Indirect SBs");
  analysisManager->CreateNtupleIColumn(fNtupleID, "EventID");
  analysisManager->CreateNtupleIColumn(fNtupleID, "PlasmidID");
  analysisManager->CreateNtupleIColumn(fNtupleID, "StrandID");
  analysisManager->CreateNtupleIColumn(fNtupleID, "BaseID");
  analysisManager->CreateNtupleIColumn(fNtupleID, "DamageID");
  analysisManager->CreateNtupleIColumn(fNtupleID, "BreakID");
  analysisManager->FinishNtuple(fNtupleID);

  for (G4int i = 0; i < fnbOfEvents; i++) {
    if (fDirectDamageMap.find(i) != fDirectDamageMap.end()) {
      for(size_t j = 0; j < fDirectDamageMap[i].size(); j++) {
        std::vector<G4int> Break = fDirectDamageMap[i][j];
        analysisManager->FillNtupleIColumn(fNtupleID, 0, i);
        analysisManager->FillNtupleIColumn(fNtupleID, 1, Break[0]);
        analysisManager->FillNtupleIColumn(fNtupleID, 2, Break[1]);
        analysisManager->FillNtupleIColumn(fNtupleID, 3, Break[2]);
        analysisManager->FillNtupleIColumn(fNtupleID, 4, 1);
        analysisManager->AddNtupleRow(fNtupleID);
      }
    }

    if (fIndirectDamageMap.find(i) != fIndirectDamageMap.end()) {
      for(size_t j = 0; j < fIndirectDamageMap[i].size(); j++) {
        std::vector<G4int> Break = fIndirectDamageMap[i][j];
        analysisManager->FillNtupleIColumn(fNtupleID, 0, i);
        analysisManager->FillNtupleIColumn(fNtupleID, 1, Break[0]);
        analysisManager->FillNtupleIColumn(fNtupleID, 2, Break[1]);
        analysisManager->FillNtupleIColumn(fNtupleID, 3, Break[2]);
        analysisManager->FillNtupleIColumn(fNtupleID, 4, 2);
        analysisManager->AddNtupleRow(fNtupleID);
      }
    }
  }

  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreStrandBreaks::ASCII(G4double LET, G4double LET_STD)
{
  std::ofstream dna(fOutputName+".txt");

  dna << "# DNA SB Map File" << G4endl;
  dna << "# LET = " << LET / (keV / um) << " +- " 
      << LET_STD / (keV / um) << " keV / um" << G4endl;
  dna << "#" << std::setw(11) << "Event-ID"  
             << std::setw(12) << "Plasmid-ID" 
             << std::setw(12) << "BP-ID" 
             << std::setw(12) << "Strand-ID"
             << std::setw(10) << "Break-ID"
             << std::setw(12) << "Dose (Gy) "<< G4endl;

  for (G4int i = 0; i < fnbOfEvents; i++) {
    if (fDirectDamageMap.find(i) != fDirectDamageMap.end()) {
      for(size_t j = 0; j < fDirectDamageMap[i].size(); j++) {
        std::vector<G4int> Break = fDirectDamageMap[i][j];
        dna << std::setw(12) << i
            << std::setw(12) << Break[0]
            << std::setw(12) << Break[1]
            << std::setw(12) << Break[2]
            << std::setw(10) << 1
            << std::setw(12) << fDoseArray[i] << G4endl;
      }
    }

    if (fIndirectDamageMap.find(i) != fIndirectDamageMap.end()) {
      for(size_t j = 0; j < fIndirectDamageMap[i].size(); j++) {
        std::vector<G4int> Break = fIndirectDamageMap[i][j];
        dna << std::setw(12) << i
            << std::setw(12) << Break[0]
            << std::setw(12) << Break[1]
            << std::setw(12) << Break[2]
            << std::setw(10) << 2
            << std::setw(12) << fDoseArray[i] << G4endl;        
      }
    }
  }

  dna.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......