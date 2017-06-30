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
// $Id: B5EventAction.cc 103553 2017-04-18 09:00:54Z gcosmo $
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

using std::array;
using std::vector;


namespace {

// Utility function which finds a hit collection with the given Id
// and print warnings if not found 
G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
  auto hce = event->GetHCofThisEvent();
  if (!hce) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl; 
      G4Exception("B5EventAction::EndOfEventAction()",
                  "B5Code001", JustWarning, msg);
      return nullptr;
  }

  auto hc = hce->GetHC(collId);
  if ( ! hc) {
    G4ExceptionDescription msg;
    msg << "Hits collection " << collId << " of this event not found." << G4endl; 
    G4Exception("B5EventAction::EndOfEventAction()",
                "B5Code001", JustWarning, msg);
  }
  return hc;  
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EventAction::B5EventAction()
: G4UserEventAction(), 
  fHodHCID  {{ -1, -1 }},
  fDriftHCID{{ -1, -1 }},
  fCalHCID  {{ -1, -1 }},
  fDriftHistoID{{ {{ -1, -1 }}, {{ -1, -1 }} }},
  fCalEdep{{ vector<G4double>(kNofEmCells, 0.), vector<G4double>(kNofHadCells, 0.) }}
      // std::array<T, N> is an aggregate that contains a C array. 
      // To initialize it, we need outer braces for the class itself 
      // and inner braces for the C array
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
  // Find hit collections and histogram Ids by names (just once)
  // and save them in the data members of this class

  if (fHodHCID[0] == -1) {
    auto sdManager = G4SDManager::GetSDMpointer();
    auto analysisManager = G4AnalysisManager::Instance();

    // hits collections names
    array<G4String, kDim> hHCName 
      = {{ "hodoscope1/hodoscopeColl", "hodoscope2/hodoscopeColl" }};
    array<G4String, kDim> dHCName 
      = {{ "chamber1/driftChamberColl", "chamber2/driftChamberColl" }};
    array<G4String, kDim> cHCName 
      = {{ "EMcalorimeter/EMcalorimeterColl", "HadCalorimeter/HadCalorimeterColl" }};

    // histograms names
    array<array<G4String, kDim>, kDim> histoName 
      = {{ {{ "Chamber1", "Chamber2" }}, {{ "Chamber1 XY", "Chamber2 XY" }} }};

    for (G4int iDet = 0; iDet < kDim; ++iDet) {
      // hit collections IDs
      fHodHCID[iDet]   = sdManager->GetCollectionID(hHCName[iDet]);
      fDriftHCID[iDet] = sdManager->GetCollectionID(dHCName[iDet]);
      fCalHCID[iDet]   = sdManager->GetCollectionID(cHCName[iDet]);
      // histograms IDs
      fDriftHistoID[kH1][iDet] = analysisManager->GetH1Id(histoName[kH1][iDet]);
      fDriftHistoID[kH2][iDet] = analysisManager->GetH2Id(histoName[kH2][iDet]);
    }
  }
}     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5EventAction::EndOfEventAction(const G4Event* event)
{
  //
  // Fill histograms & ntuple
  // 

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
 
  // Drift chambers hits
  for (G4int iDet = 0; iDet < kDim; ++iDet) {
    auto hc = GetHC(event, fDriftHCID[iDet]);
    if ( ! hc ) return;

    auto nhit = hc->GetSize();
    analysisManager->FillH1(fDriftHistoID[kH1][iDet], nhit );
    // columns 0, 1
    analysisManager->FillNtupleIColumn(iDet, nhit);
  
    for (unsigned long i = 0; i < nhit; ++i) {
      auto hit = static_cast<B5DriftChamberHit*>(hc->GetHit(i));
      auto localPos = hit->GetLocalPos();
      analysisManager->FillH2(fDriftHistoID[kH2][iDet], localPos.x(), localPos.y());
    }
  }
      
  // Em/Had Calorimeters hits
  array<G4int, kDim> totalCalHit = {{ 0, 0 }}; 
  array<G4double, kDim> totalCalEdep = {{ 0., 0. }}; 

  for (G4int iDet = 0; iDet < kDim; ++iDet) {
    auto hc = GetHC(event, fCalHCID[iDet]);
    if ( ! hc ) return;

    totalCalHit[iDet] = 0;
    totalCalEdep[iDet] = 0.;
    for (unsigned long i = 0; i < hc->GetSize(); ++i) {
      G4double edep = 0.;
      // The EM and Had calorimeter hits are of different types
      if (iDet == 0) {
        auto hit = static_cast<B5EmCalorimeterHit*>(hc->GetHit(i));
        edep = hit->GetEdep();
      } else {
        auto hit = static_cast<B5HadCalorimeterHit*>(hc->GetHit(i));
        edep = hit->GetEdep();
      }
      if ( edep > 0. ) {
        totalCalHit[iDet]++;
        totalCalEdep[iDet] += edep;
      }
      fCalEdep[iDet][i] = edep;
    }
    // columns 2, 3
    analysisManager->FillNtupleDColumn(iDet + 2, totalCalEdep[iDet]);
  }

  // Hodoscopes hits
  for (G4int iDet = 0; iDet < kDim; ++iDet) {
    auto hc = GetHC(event, fHodHCID[iDet]);
    if ( ! hc ) return;

    for (unsigned int i = 0; i<hc->GetSize(); ++i) {
      auto hit = static_cast<B5HodoscopeHit*>(hc->GetHit(i));
      // columns 4, 5
      analysisManager->FillNtupleDColumn(iDet + 4, hit->GetTime());
    }
  }
    
  //
  // Print diagnostics
  // 
  
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( printModulo == 0 || event->GetEventID() % printModulo != 0) return;
  
  auto primary = event->GetPrimaryVertex(0)->GetPrimary(0);
  G4cout 
    << G4endl
    << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
    << primary->GetG4code()->GetParticleName()
    << " " << primary->GetMomentum() << G4endl;
  
  // Hodoscopes
  for (G4int iDet = 0; iDet < kDim; ++iDet) {
    auto hc = GetHC(event, fHodHCID[iDet]);
    if ( ! hc ) return;
    G4cout << "Hodoscope " << iDet + 1 << " has " << hc->GetSize()  << " hits." << G4endl;
    for (unsigned int i = 0; i<hc->GetSize(); ++i) {
      hc->GetHit(i)->Print();
    }
  }

  // Drift chambers
  for (G4int iDet = 0; iDet < kDim; ++iDet) {
    auto hc = GetHC(event, fDriftHCID[iDet]);
    if ( ! hc ) return;
    G4cout << "Drift Chamber " << iDet + 1 << " has " <<  hc->GetSize()  << " hits." << G4endl;
    for (auto layer = 0; layer < kNofChambers; ++layer) {
      for (unsigned int i = 0; i < hc->GetSize(); i++) {
        auto hit = static_cast<B5DriftChamberHit*>(hc->GetHit(i));
        if (hit->GetLayerID() == layer) hit->Print();
      }
    }
  }

  // Calorimeters
  array<G4String, kDim> calName = {{ "EM", "Hadron" }};
  for (G4int iDet = 0; iDet < kDim; ++iDet) {
    G4cout << calName[iDet] << " Calorimeter has " << totalCalHit[iDet] << " hits." 
           << " Total Edep is " << totalCalEdep[iDet]/MeV << " (MeV)" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
