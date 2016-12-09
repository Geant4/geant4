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
// $Id: B5EventAction.cc 101036 2016-11-04 09:00:23Z gcosmo $
//
/// \file B5EventAction.cc
/// \brief Implementation of the B5EventAction class

#include "B5EventAction.hh"
#include "B5HodoscopeHit.hh"
#include "B5DriftChamberHit.hh"
#include "B5EmCalorimeterHit.hh"
#include "B5HadCalorimeterHit.hh"
#include "B5Constants.hh"
#include "B5Analysis.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EventAction::B5EventAction()
: G4UserEventAction(), 
  fHodHC1ID(-1),
  fHodHC2ID(-1),
  fDriftHC1ID(-1),
  fDriftHC2ID(-1),
  fEmCalHCID(-1),
  fHadCalHCID(-1),
  fEmCalEdep(kNofEmCells, 0.), 
  fHadCalEdep(kNofHadCells, 0.)
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EventAction::~B5EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5EventAction::BeginOfEventAction(const G4Event*)
{
  if (fHodHC1ID==-1) {
    auto sdManager = G4SDManager::GetSDMpointer();
    fHodHC1ID = sdManager->GetCollectionID("hodoscope1/hodoscopeColl");
    fHodHC2ID = sdManager->GetCollectionID("hodoscope2/hodoscopeColl");
    fDriftHC1ID = sdManager->GetCollectionID("chamber1/driftChamberColl");
    fDriftHC2ID = sdManager->GetCollectionID("chamber2/driftChamberColl");
    fEmCalHCID = sdManager->GetCollectionID("EMcalorimeter/EMcalorimeterColl");
    fHadCalHCID = sdManager->GetCollectionID("HadCalorimeter/HadCalorimeterColl");
  }
}     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5EventAction::EndOfEventAction(const G4Event* event)
{
  auto hce = event->GetHCofThisEvent();
  if (!hce) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl; 
      G4Exception("B5EventAction::EndOfEventAction()",
                  "B5Code001", JustWarning, msg);
      return;
  }

  // Get hits collections 
  auto hHC1 
    = static_cast<B5HodoscopeHitsCollection*>(hce->GetHC(fHodHC1ID));
    
  auto hHC2 
    = static_cast<B5HodoscopeHitsCollection*>(hce->GetHC(fHodHC2ID));
    
  auto dHC1 
    = static_cast<B5DriftChamberHitsCollection*>(hce->GetHC(fDriftHC1ID));
    
  auto dHC2 
    = static_cast<B5DriftChamberHitsCollection*>(hce->GetHC(fDriftHC2ID));
    
  auto ecHC 
    = static_cast<B5EmCalorimeterHitsCollection*>(hce->GetHC(fEmCalHCID));
    
  auto hcHC 
    = static_cast<B5HadCalorimeterHitsCollection*>(hce->GetHC(fHadCalHCID));
    
  if ( (!hHC1) || (!hHC2) || (!dHC1) || (!dHC2) || (!ecHC) || (!hcHC) ) {
      G4ExceptionDescription msg;
      msg << "Some of hits collections of this event not found." << G4endl; 
      G4Exception("B5EventAction::EndOfEventAction()",
                  "B5Code001", JustWarning, msg);
      return;
  }   
  
  //
  // Fill histograms & ntuple
  // 
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
 
  // Fill histograms
 
  auto nhit  = dHC1->entries();
  analysisManager->FillH1(0, nhit );

  for (auto i=0;i<nhit ;i++) {
    auto hit = (*dHC1)[i];
    auto localPos = hit->GetLocalPos();
    analysisManager->FillH2(0, localPos.x(), localPos.y());
  }
 
  nhit  = dHC2->entries();
  analysisManager->FillH1(1, nhit );

  for (auto i=0;i<nhit ;i++) {
    auto hit = (*dHC2)[i];
    auto localPos = hit->GetLocalPos();
    analysisManager->FillH2(1, localPos.x(), localPos.y());
  }
      
  // Fill ntuple
  
  // Dc1Hits
  analysisManager->FillNtupleIColumn(0, dHC1->entries());
  // Dc2Hits
  analysisManager->FillNtupleIColumn(1, dHC1->entries());
  
  // ECEnergy
  G4int totalEmHit = 0;
  G4double totalEmE = 0.;
  for (auto i=0;i<kNofEmCells;i++) {
    auto hit = (*ecHC)[i];
    auto edep = hit->GetEdep();
    if (edep>0.) {
      totalEmHit++;
      totalEmE += edep;
    }
    fEmCalEdep[i] = edep;
  }
  analysisManager->FillNtupleDColumn(2, totalEmE);

  // HCEnergy
  G4int totalHadHit = 0;
  G4double totalHadE = 0.;
  for (auto i=0;i<kNofHadCells;i++) {
    auto hit = (*hcHC)[i];
    auto edep = hit->GetEdep();
    if (edep>0.) {
        totalHadHit++;
        totalHadE += edep;
    }
    fHadCalEdep[i] = edep;
  }
  analysisManager->FillNtupleDColumn(3, totalHadE);

  // Time 1
  for (auto i=0;i<hHC1->entries();i++) {
    analysisManager->FillNtupleDColumn(4,(*hHC1)[i]->GetTime());
  }
    
  // Time 2
  for (auto i=0;i<hHC2->entries();i++) {
    analysisManager->FillNtupleDColumn(5,(*hHC2)[i]->GetTime());
  }
  
  analysisManager->AddNtupleRow();  
  
  //
  // Print diagnostics
  // 
  
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( printModulo==0 || event->GetEventID() % printModulo != 0) return;
  
  auto primary = event->GetPrimaryVertex(0)->GetPrimary(0);
  G4cout 
    << G4endl
    << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
    << primary->GetG4code()->GetParticleName()
    << " " << primary->GetMomentum() << G4endl;
  
  // Hodoscope 1
  nhit  = hHC1->entries();
  G4cout << "Hodoscope 1 has " << nhit  << " hits." << G4endl;
  for (auto i=0;i<nhit ;i++) {
    auto hit = (*hHC1)[i];
    hit->Print();
  }

  // Hodoscope 2
  nhit  = hHC2->entries();
  G4cout << "Hodoscope 2 has " << nhit  << " hits." << G4endl;
  for (auto i=0;i<nhit ;i++) {
    auto hit = (*hHC2)[i];
    hit->Print();
  }

  // Drift chamber 1
  nhit  = dHC1->entries();
  G4cout << "Drift Chamber 1 has " << nhit  << " hits." << G4endl;
  for (auto layer=0;layer<kNofChambers;layer++) {
    for (auto i=0;i<nhit ;i++) {
      auto hit = (*dHC1)[i];
      if (hit->GetLayerID()==layer) hit->Print();
    }
  }

  // Drift chamber 2
  nhit  = dHC2->entries();
  G4cout << "Drift Chamber 2 has " << nhit  << " hits." << G4endl;
  for (auto layer=0;layer<kNofChambers;layer++) {
    for (auto i=0;i<nhit ;i++) {
      auto hit = (*dHC2)[i];
      if (hit->GetLayerID()==layer) hit->Print();
    }
  }

  // EM calorimeter
  G4cout << "EM Calorimeter has " << totalEmHit << " hits. Total Edep is "
    << totalEmE/MeV << " (MeV)" << G4endl;

  // Had calorimeter
  G4cout << "Hadron Calorimeter has " << totalHadHit << " hits. Total Edep is "
    << totalHadE/MeV << " (MeV)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
