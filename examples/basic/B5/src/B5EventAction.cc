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
// $Id: B5EventAction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B5EventAction.cc
/// \brief Implementation of the B5EventAction class

#include "B5EventAction.hh"
#include "B5HodoscopeHit.hh"
#include "B5DriftChamberHit.hh"
#include "B5EmCalorimeterHit.hh"
#include "B5HadCalorimeterHit.hh"
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
  fHHC1ID(-1),
  fHHC2ID(-1),
  fDHC1ID(-1),
  fDHC2ID(-1),
  fECHCID(-1),
  fHCHCID(-1)
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
    if (fHHC1ID==-1) {
      G4SDManager* sdManager = G4SDManager::GetSDMpointer();
      fHHC1ID = sdManager->GetCollectionID("hodoscope1/hodoscopeColl");
      fHHC2ID = sdManager->GetCollectionID("hodoscope2/hodoscopeColl");
      fDHC1ID = sdManager->GetCollectionID("chamber1/driftChamberColl");
      fDHC2ID = sdManager->GetCollectionID("chamber2/driftChamberColl");
      fECHCID = sdManager->GetCollectionID("EMcalorimeter/EMcalorimeterColl");
      fHCHCID = sdManager->GetCollectionID("HadCalorimeter/HadCalorimeterColl");
    }
}     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5EventAction::EndOfEventAction(const G4Event* event)
{
    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) 
    {
        G4ExceptionDescription msg;
        msg << "No hits collection of this event found." << G4endl; 
        G4Exception("B5EventAction::EndOfEventAction()",
                    "B5Code001", JustWarning, msg);
        return;
    }   


    // Get hits collections 
    B5HodoscopeHitsCollection* hHC1 
      = static_cast<B5HodoscopeHitsCollection*>(hce->GetHC(fHHC1ID));
      
    B5HodoscopeHitsCollection* hHC2 
      = static_cast<B5HodoscopeHitsCollection*>(hce->GetHC(fHHC2ID));
      
    B5DriftChamberHitsCollection* dHC1 
      = static_cast<B5DriftChamberHitsCollection*>(hce->GetHC(fDHC1ID));
      
    B5DriftChamberHitsCollection* dHC2 
      = static_cast<B5DriftChamberHitsCollection*>(hce->GetHC(fDHC2ID));
      
    B5EmCalorimeterHitsCollection* ecHC 
      = static_cast<B5EmCalorimeterHitsCollection*>(hce->GetHC(fECHCID));
      
    B5HadCalorimeterHitsCollection* hcHC 
      = static_cast<B5HadCalorimeterHitsCollection*>(hce->GetHC(fHCHCID));
      
    if ( (!hHC1) || (!hHC2) || (!dHC1) || (!dHC2) || (!ecHC) || (!hcHC) ) 
    {
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
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
    // Fill histograms
 
    G4int n_hit = dHC1->entries();
    analysisManager->FillH1(0, n_hit);

    for (G4int i=0;i<n_hit;i++)
    {
       B5DriftChamberHit* hit = (*dHC1)[i];
       G4ThreeVector localPos = hit->GetLocalPos();
       analysisManager->FillH2(0, localPos.x(), localPos.y());
    }
 
    n_hit = dHC2->entries();
    analysisManager->FillH1(1, n_hit);

    for (G4int i=0;i<n_hit;i++)
    {
       B5DriftChamberHit* hit = (*dHC2)[i];
       G4ThreeVector localPos = hit->GetLocalPos();
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
    for (G4int i=0;i<80;i++)
    {
        B5EmCalorimeterHit* hit = (*ecHC)[i];
        G4double eDep = hit->GetEdep();
        if (eDep>0.)
        {
            totalEmHit++;
            totalEmE += eDep;
        }
    }
    analysisManager->FillNtupleDColumn(2, totalEmE);

    // HCEnergy
    G4int totalHadHit = 0;
    G4double totalHadE = 0.;
    for (G4int i=0;i<20;i++)
    {
        B5HadCalorimeterHit* hit = (*hcHC)[i];
        G4double eDep = hit->GetEdep();
        if (eDep>0.)
        {
            totalHadHit++;
            totalHadE += eDep;
        }
    }
    analysisManager->FillNtupleDColumn(3, totalHadE);

    // Time 1
    for (G4int i=0;i<hHC1->entries();i++)
    {
      analysisManager->FillNtupleDColumn(4,(*hHC1)[i]->GetTime());
    }
      
    // Time 2
    for (G4int i=0;i<hHC2->entries();i++)
    {
      analysisManager->FillNtupleDColumn(5,(*hHC2)[i]->GetTime());
    }
      
    analysisManager->AddNtupleRow();  
    
    //
    // Print diagnostics
    // 
    
    G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ( printModulo==0 || event->GetEventID() % printModulo != 0) return;
    
    G4PrimaryParticle* primary = event->GetPrimaryVertex(0)->GetPrimary(0);
    G4cout << G4endl
           << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
           << primary->GetG4code()->GetParticleName()
           << " " << primary->GetMomentum() << G4endl;
    
    // Hodoscope 1
    n_hit = hHC1->entries();
    G4cout << "Hodoscope 1 has " << n_hit << " hits." << G4endl;
    for (G4int i=0;i<n_hit;i++)
    {
        B5HodoscopeHit* hit = (*hHC1)[i];
        hit->Print();
    }

    // Hodoscope 2
    n_hit = hHC2->entries();
    G4cout << "Hodoscope 2 has " << n_hit << " hits." << G4endl;
    for (G4int i=0;i<n_hit;i++)
    {
        B5HodoscopeHit* hit = (*hHC2)[i];
        hit->Print();
    }

    // Drift chamber 1
    n_hit = dHC1->entries();
    G4cout << "Drift Chamber 1 has " << n_hit << " hits." << G4endl;
    for (G4int i2=0;i2<5;i2++)
    {
        for (G4int i=0;i<n_hit;i++)
        {
            B5DriftChamberHit* hit = (*dHC1)[i];
            if (hit->GetLayerID()==i2) hit->Print();
        }
    }

    // Drift chamber 2
    n_hit = dHC2->entries();
    G4cout << "Drift Chamber 2 has " << n_hit << " hits." << G4endl;
    for (G4int i2=0;i2<5;i2++)
    {
        for (G4int i=0;i<n_hit;i++)
        {
            B5DriftChamberHit* hit = (*dHC2)[i];
            if (hit->GetLayerID()==i2) hit->Print();
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
