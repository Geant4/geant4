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
// dnadamage3 example is derived from the chem6 example
// chem6 example authors: W. G. Shin and S. Incerti (CENBG, France)
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
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file ScoreSpecies.cc
/// \brief Implementation of the ScoreSpecies class
///
///  This is a primitive scorer class for molecular species.
///  The number of species is recorded for all times (predetermined or 
///  user chosen). It also scores the energy deposition in order to compute the 
///  radiochemical yields.

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreSpecies::ScoreSpecies(G4String name, G4int depth)
: G4VPrimitiveScorer(name,depth),
  G4UImessenger(),
  fOutputType("root"), // other options: "csv", "hdf5", "xml"
  fHCID(-1)
{
  fSpeciesdir = new G4UIdirectory("/scorer/species/");
  fSpeciesdir->SetGuidance("ScoreSpecies commands");

  fAddTimeToRecordcmd =
    new G4UIcmdWithADoubleAndUnit("/scorer/species/addTimeToRecord",this);

  fTimeBincmd = new G4UIcmdWithAnInteger("/scorer/species/nOfTimeBins",this);

  fOutputTypeUI = new G4UIcmdWithAString("/scorer/species/OutputFormat",this);
  fOutputFileUI = new G4UIcmdWithAString("/scorer/species/OutputFile",this);

  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScoreSpecies::~ScoreSpecies()
{
  delete fSpeciesdir;
  delete fAddTimeToRecordcmd;
  delete fTimeBincmd;
  delete fOutputTypeUI;
  delete fOutputFileUI;
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

  if(command == fOutputTypeUI) {
    fOutputType = newValue;
    G4StrUtil::to_lower(fOutputType);
  }

  if(command == fOutputFileUI)
    fOutputFile = newValue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool ScoreSpecies::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if ( edep == 0. ) return FALSE;

  edep *= aStep->GetPreStepPoint()->GetWeight();
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
      molInfo.fNumber  += n_mol;
      molInfo.fNumber2 += n_mol*n_mol;
      double gValue = (n_mol/(fEdep/eV)) * 100.;
      molInfo.fG  += gValue;
      molInfo.fG2 += gValue*gValue;
    }
  }

  ++fNEvent;

  fEdep = 0.;
  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
ScoreSpecies::AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer* workerScorer)
{
  ScoreSpecies* right = dynamic_cast<ScoreSpecies*>(workerScorer);

  if(right == 0)
  {
    return;
  }
  if(right == this)
  {
    return;
  }

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
      fSpeciesInfoPerTime[it_map1->first][it_map2->first];
      molInfo.fNumber  += it_map2->second.fNumber;
      molInfo.fNumber2 += it_map2->second.fNumber2;
      molInfo.fG  += it_map2->second.fG;
      molInfo.fG2 += it_map2->second.fG2;
    }
  }
  right->fSpeciesInfoPerTime.clear();

  fNEvent += right->fNEvent;
  right->fNEvent = 0;
  right->fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScoreSpecies::DrawAll()
{
}

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

  if(fOutputType != "ascii") 
  {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType(fOutputType);

    if(analysisManager)
    {
      this->WriteWithAnalysisManager(analysisManager);
    }
  }
  else
    OutputToASCII();

  fNEvent = 0;
  fSpeciesInfoPerTime.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
ScoreSpecies::WriteWithAnalysisManager(G4VAnalysisManager* analysisManager)
{
  analysisManager->OpenFile(fOutputFile);
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

    for(auto it_map2 : map2)
    {
      double time = it_map1.first;
      auto species = it_map2.first;
      const G4String& name = species->GetName();
      int molID = it_map2.first->GetMoleculeID();
      int number = it_map2.second.fNumber;
      double G = it_map2.second.fG;
      double G2 = it_map2.second.fG2;

      analysisManager->FillNtupleIColumn(fNtupleID, 0, molID);
      analysisManager->FillNtupleIColumn(fNtupleID, 1, number);
      analysisManager->FillNtupleIColumn(fNtupleID,
                                         2, fNEvent);
      analysisManager->FillNtupleSColumn(fNtupleID, 3, name);
      analysisManager->FillNtupleDColumn(fNtupleID, 4, time);
      analysisManager->FillNtupleDColumn(fNtupleID, 5, G);
      analysisManager->FillNtupleDColumn(fNtupleID, 6, G2);
      analysisManager->AddNtupleRow(fNtupleID);
    }
  }

  analysisManager->Write();
  analysisManager->CloseFile();
  fRunID++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreSpecies::OutputToASCII(){

  std::ofstream SpeciesOutput;
  SpeciesOutput.open(fOutputFile+".txt");

  std::map<G4String, std::map<G4double, std::vector<G4double>>> mol;

  for(auto it_map1: fSpeciesInfoPerTime) {
    InnerSpeciesMap& map2 = it_map1.second;
    G4double time = it_map1.first/ps;

    for(auto it_map2: map2) {
      G4double G  = it_map2.second.fG;
      G4double G2 = it_map2.second.fG2;
      G4double Nb  = it_map2.second.fNumber;
      G4double Nb2 = it_map2.second.fNumber2;
      G4double N  = fNEvent;
      if (N > 0) {
        G  /= N;
        Nb /= N;
      }

      if ( N == 1 ) {
        G2  = 0.0;
        Nb2 = 0.0;
      }

      else if (N > 0) {
        G2  = std::sqrt( N/(N-1) * ( G2/N - G*G) );
        Nb2 = std::sqrt( N/(N-1) * ( Nb2/N - Nb*Nb) );
      }

      mol[it_map2.first->GetName()][time].push_back(G);
      mol[it_map2.first->GetName()][time].push_back(G2);
      mol[it_map2.first->GetName()][time].push_back(Nb);
      mol[it_map2.first->GetName()][time].push_back(Nb2);
    }
  }

  SpeciesOutput << "# Species Scorer Output " << G4endl;
  SpeciesOutput << "# " << std::setw(10) << "Time [ps]" <<
                           std::setw(15) << "Gx(/100 eV)" <<
                           std::setw(15) << "RMS" <<
                           std::setw(15) << "Gx(#Mols)" <<
                           std::setw(20) << "Molecule" << G4endl;

  for ( auto it1 : mol )
    for ( auto it2 : it1.second )
      SpeciesOutput << std::setw(12) << it2.first      
          << "   " << std::setw(12) << it2.second[0]
          << "   " << std::setw(12) << it2.second[1]  
          << "   " << std::setw(12) << it2.second[2]
          << "   " << std::setw(17) << it1.first  << G4endl;

  SpeciesOutput.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
