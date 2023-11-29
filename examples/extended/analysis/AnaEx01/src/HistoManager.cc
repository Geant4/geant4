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
/// \file analysis/AnaEx01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  // Create or get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
    // the default file type can be overriden in run macro
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if ( ! fFactoryOn ) {
    //
    analysisManager->SetVerboseLevel(1);
    // Only merge in MT mode to avoid warning when running in Sequential mode
  #ifdef G4MULTITHREADED
    analysisManager->SetNtupleMerging(true);
  #endif

    // Create directories
    analysisManager->SetHistoDirectoryName("histo");
    analysisManager->SetNtupleDirectoryName("ntuple");
  }

  // Open an output file
  //
  G4bool fileOpen = analysisManager->OpenFile("AnaEx01");
  if (! fileOpen) {
    G4cerr << "\n---> HistoManager::Book(): cannot open "
           << analysisManager->GetFileName() << G4endl;
    return;
  }

  if ( ! fFactoryOn ) {
    // Create histograms.
    // Histogram ids are generated automatically starting from 0.
    // The start value can be changed by:
    // analysisManager->SetFirstHistoId(1);

    // id = 0
    analysisManager->CreateH1("EAbs","Edep in absorber (MeV)", 100, 0., 800*MeV);
    // id = 1
    analysisManager->CreateH1("EGap","Edep in gap (MeV)", 100, 0., 100*MeV);
    // id = 2
    analysisManager->CreateH1("LAbs","trackL in absorber (mm)", 100, 0., 1*m);
    // id = 3
    analysisManager->CreateH1("LGap","trackL in gap (mm)", 100, 0., 50*cm);

    // Create ntuples.
    // Ntuples ids are generated automatically starting from 0.
    // The start value can be changed by:
    // analysisManager->SetFirstMtupleId(1);

    // Create 1st ntuple (id = 0)
    analysisManager->CreateNtuple("Ntuple1", "Edep");
    analysisManager->CreateNtupleDColumn("Eabs"); // column Id = 0
    analysisManager->CreateNtupleDColumn("Egap"); // column Id = 1
    analysisManager->FinishNtuple();

    // Create 2nd ntuple (id = 1)
    //
    analysisManager->CreateNtuple("Ntuple2", "TrackL");
    analysisManager->CreateNtupleDColumn("Labs"); // column Id = 0
    analysisManager->CreateNtupleDColumn("Lgap"); // column Id = 1
    analysisManager->FinishNtuple();

    fFactoryOn = true;
  }

  G4cout << "\n----> Output file is open in "
         << analysisManager->GetFileName() << "."
         << analysisManager->GetFileType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Save()
{
  if (! fFactoryOn) { return; }

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(ih, xbin, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  auto h1 = analysisManager->GetH1(ih);
  if (h1 != nullptr) { h1->scale(fac);
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple(G4double energyAbs, G4double energyGap,
                              G4double trackLAbs, G4double trackLGap)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Fill 1st ntuple ( id = 0)
  analysisManager->FillNtupleDColumn(0, 0, energyAbs);
  analysisManager->FillNtupleDColumn(0, 1, energyGap);
  analysisManager->AddNtupleRow(0);
  // Fill 2nd ntuple ( id = 1)
  analysisManager->FillNtupleDColumn(1, 0, trackLAbs);
  analysisManager->FillNtupleDColumn(1, 1, trackLGap);
  analysisManager->AddNtupleRow(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
  if (! fFactoryOn) { return; }

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4cout << "\n ----> print histograms statistic \n" << G4endl;
  for ( G4int i=0; i<analysisManager->GetNofH1s(); ++i ) {
    G4String name = analysisManager->GetH1Name(i);
    auto h1 = analysisManager->GetH1(i);

    G4String unitCategory;
    if (name[0U] == 'E' ) { unitCategory = "Energy"; }
    if (name[0U] == 'L' ) { unitCategory = "Length"; }
         // we use an explicit unsigned int type for operator [] argument
         // to avoid problems with windows compiler

    G4cout << name
           << ": mean = " << G4BestUnit(h1->mean(), unitCategory)
           << " rms = " << G4BestUnit(h1->rms(), unitCategory )
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
