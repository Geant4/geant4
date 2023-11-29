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
//
// $Id$
//
/// \file ScoreSpecies.cc
/// \brief Implementation of the ScoreSpecies class

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
#include "G4TScoreNtupleWriter.hh"

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

ScoreSpecies::ScoreSpecies(G4String name, G4int depth)
: G4VPrimitiveScorer(name,depth),
  G4UImessenger(),
  fEdep(0),
  fOutputType("root"), // other options: "csv", "hdf5", "xml"
  fHCID(-1),
  fEvtMap(0)
{
  fSpeciesdir = new G4UIdirectory("/scorer/species/");
  fSpeciesdir->SetGuidance("ScoreSpecies commands");

  fAddTimeToRecordcmd =
    new G4UIcmdWithADoubleAndUnit("/scorer/species/addTimeToRecord",this);

  fTimeBincmd = new G4UIcmdWithAnInteger("/scorer/species/nOfTimeBins",this);

  fEdep = 0;
  fNEvent = 0;
  fRunID = 0;
  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreSpecies::~ScoreSpecies()
{
  delete fSpeciesdir;
  delete fAddTimeToRecordcmd;
  delete fTimeBincmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fAddTimeToRecordcmd){
    G4double cmdTime = fAddTimeToRecordcmd->GetNewDoubleValue(newValue);
    AddTimeToRecord(cmdTime);
  }
  if(command == fTimeBincmd){
    ClearTimeToRecord();
    G4int cmdBins = fTimeBincmd->GetNewIntValue(newValue);
    G4double timeMin = 1*ps;
    G4double timeMax = G4Scheduler::Instance()->GetEndTime() - 1*ps;
    G4double timeLogMin = std::log10(timeMin);
    G4double timeLogMax = std::log10(timeMax);
    for(G4int i=0;i<cmdBins;i++){
      AddTimeToRecord(std::pow(10,timeLogMin +
                               i*(timeLogMax-timeLogMin)/(cmdBins-1)));
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool ScoreSpecies::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if ( edep == 0. ) return FALSE;

  edep *= aStep->GetPreStepPoint()->GetWeight(); // (Particle Weight)
  G4int  index = GetIndex(aStep);
  fEvtMap->add(index,edep);
  fEdep+=edep;

  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::Initialize(G4HCofThisEvent* HCE)
{
  fEvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
                                    GetName());

  if(fHCID < 0)
  {
    fHCID = GetCollectionID(0);
  }

  HCE->AddHitsCollection(fHCID, (G4VHitsCollection*)fEvtMap);
  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::EndOfEvent(G4HCofThisEvent*)
{
  if(G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted())
  {
    fEdep = 0.;
    G4MoleculeCounter::Instance()->ResetCounter();
    return;
  }

  auto species = G4MoleculeCounter::Instance()->GetRecordedMolecules();

  if(species.get() == 0 || species->size() == 0)
  {
    G4cout << "No molecule recorded, energy deposited= "
           << G4BestUnit(fEdep, "Energy") << G4endl;
    ++fNEvent;
    fEdep = 0.;
    G4MoleculeCounter::Instance()->ResetCounter();
    return;
  }

//  G4cout << "ScoreSpecies::EndOfEvent"<<G4endl;
//  G4cout << "End of event, deposited energy: "
//  << G4BestUnit(fEdep, "Energy") << G4endl;

  for(auto molecule: *species)
  {
    for(auto time_mol: fTimeToRecord)
    {
      double n_mol =
          G4MoleculeCounter::Instance()->GetNMoleculesAtTime(molecule,
                                                             time_mol);

      if(n_mol < 0)
      {
        G4cerr << "N molecules not valid < 0 " << G4endl;
        G4Exception("","N<0",FatalException,"");
      }

      SpeciesInfo& molInfo = fSpeciesInfoPerTime[time_mol][molecule];
      molInfo.fNumber += n_mol;
      double gValue = (n_mol/(fEdep/eV)) * 100.;
      molInfo.fG += gValue;
      molInfo.fG2 += gValue*gValue;

      //      G4cout << "In Save molucule: fNumber " << molInfo.fNumber
      //            << " fG " << molInfo.fG
      //            << " fEdep " << fEdep/eV
      //            << G4endl;
    }
  }

  ++fNEvent;

//  G4cout << "End of event " << fNEvent
//         << ", energy deposited=" << G4BestUnit(fEdep, "Energy") << G4endl;

  fEdep = 0.;
  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
ScoreSpecies::AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer* workerScorer)
{
  ScoreSpecies* right =
  dynamic_cast<ScoreSpecies*>(dynamic_cast<G4VPrimitiveScorer*>(workerScorer));

  if(right == 0)
  {
    return;
  }
  if(right == this)
  {
    return;
  }

  // G4cout<<"ScoreSpecies::AbsorbResultsFromWorkerScorer"<<G4endl;
  SpeciesMap::iterator it_map1 = right->fSpeciesInfoPerTime.begin();
  SpeciesMap::iterator end_map1 = right->fSpeciesInfoPerTime.end();

  for(; it_map1 != end_map1; ++it_map1)
  {
    InnerSpeciesMap& map2 = it_map1->second;
    InnerSpeciesMap::iterator it_map2 = map2.begin();
    InnerSpeciesMap::iterator end_map2 = map2.end();

    for(; it_map2 != end_map2; ++it_map2)
    {
      SpeciesInfo& molInfo =
      fSpeciesInfoPerTime[it_map1->first][it_map2->first] ;
      molInfo.fNumber  += it_map2->second.fNumber;
      molInfo.fG += it_map2->second.fG;
      molInfo.fG2 += it_map2->second.fG2;

      //      G4cout << "In AbsorbeResultsFromWorkerScorer: fNumber "
      //             << molInfo.fNumber
      //             << " fG "
      //             << molInfo.fG
      //             << G4endl;
    }
  }
  right->fSpeciesInfoPerTime.clear();

  fNEvent += right->fNEvent;
  right->fNEvent = 0;
  right->fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::DrawAll()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::PrintAll()
{
   G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
   G4cout << " PrimitiveScorer " << GetName() << G4endl;
   G4cout << " Number of events " << fNEvent << G4endl;
   G4cout << " Number of energy deposition recorded "
          << fEvtMap->entries() << G4endl;

  for(auto itr : *fEvtMap->GetMap()) {
     G4cout << "  copy no.: " << itr.first
     << "  energy deposit: "
     << *(itr.second)/GetUnitValue()
     << " [" << GetUnit()<<"]"
     << G4endl;
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::OutputAndClear()
{
  if(G4Threading::IsWorkerThread()) return;

  //---------------------------------------------------------------------------
  // Save results

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType(fOutputType);

  if(analysisManager)
  {
    this->WriteWithAnalysisManager(analysisManager);
  }

  fNEvent = 0;
  fSpeciesInfoPerTime.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
ScoreSpecies::WriteWithAnalysisManager(G4VAnalysisManager* analysisManager)
{
  G4String fileN = "Species"+G4UIcommand::ConvertToString(fRunID);
  analysisManager->OpenFile(fileN);
  int fNtupleID = analysisManager->CreateNtuple("species", "species");
  analysisManager->CreateNtupleIColumn(fNtupleID, "speciesID");
  analysisManager->CreateNtupleIColumn(fNtupleID, "number");
  analysisManager->CreateNtupleIColumn(fNtupleID, "nEvent");
  analysisManager->CreateNtupleSColumn(fNtupleID, "speciesName");
  analysisManager->CreateNtupleDColumn(fNtupleID, "time");
  analysisManager->CreateNtupleDColumn(fNtupleID, "sumG");
  analysisManager->CreateNtupleDColumn(fNtupleID, "sumG2");
  analysisManager->FinishNtuple(fNtupleID);

  for(auto it_map1: fSpeciesInfoPerTime)
  {
    InnerSpeciesMap& map2 = it_map1.second;

    if(it_map1.first == fSpeciesInfoPerTime.begin()->first){
      for(auto it_molname : map2){
        auto tmp_species = it_molname.first;
        out<<std::setw(12)<<tmp_species->GetName()<<
             std::setw(12)<<tmp_species->GetMoleculeID();
      }
      out<<'\n';
    }

    for(auto it_map2 : map2)
    {
      double time = it_map1.first;
      auto species = it_map2.first;
      const G4String& name = species->GetName();
      int molID = it_map2.first->GetMoleculeID();
      int number = it_map2.second.fNumber;
      double G = it_map2.second.fG;
      double G2 = it_map2.second.fG2;
      G4int N = fNEvent;

      if(time == *fTimeToRecord.rbegin()){
        if(N > 1)
        {
          out<<std::setw(12)<<G/N<<std::setw(12)<<
            std::sqrt(((G2/N)-std::pow(G/N,2))/(N-1));
        }else{
          out<<std::setw(12)<<G/N<<std::setw(12)<<
            std::sqrt(((G2/N)-std::pow(G/N,2))/N);
        }

      }

      analysisManager->FillNtupleIColumn(fNtupleID, 0, molID);  // MolID
      analysisManager->FillNtupleIColumn(fNtupleID, 1, number); // Number
      analysisManager->FillNtupleIColumn(fNtupleID,
                                         2, fNEvent); // Total nb events
      analysisManager->FillNtupleSColumn(fNtupleID, 3, name);   // molName
      analysisManager->FillNtupleDColumn(fNtupleID, 4, time);   // time
      analysisManager->FillNtupleDColumn(fNtupleID, 5, G);      // G
      analysisManager->FillNtupleDColumn(fNtupleID, 6, G2);     // G2
      analysisManager->AddNtupleRow(fNtupleID);
    }
  }

  analysisManager->Write();
  analysisManager->CloseFile();
  fRunID++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
