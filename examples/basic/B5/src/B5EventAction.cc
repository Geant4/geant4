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
// $Id$
//
/// \file B5EventAction.cc
/// \brief Implementation of the B5EventAction class

#include "B5EventAction.hh"
#include "B5EventActionMessenger.hh"
#include "B5HodoscopeHit.hh"
#include "B5DriftChamberHit.hh"
#include "B5EmCalorimeterHit.hh"
#include "B5HadCalorimeterHit.hh"
#include "B5Analysis.hh"

#include "G4Event.hh"
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
  fHCHCID(-1),
  fMessenger(0),
  fVerboseLevel(1)
{
    G4String colName;
    G4SDManager* sdMan = G4SDManager::GetSDMpointer();
    fHHC1ID = sdMan->GetCollectionID(colName="hodoscope1/hodoscopeColl");
    fHHC2ID = sdMan->GetCollectionID(colName="hodoscope2/hodoscopeColl");
    fDHC1ID = sdMan->GetCollectionID(colName="chamber1/driftChamberColl");
    fDHC2ID = sdMan->GetCollectionID(colName="chamber2/driftChamberColl");
    fECHCID = sdMan->GetCollectionID(colName="EMcalorimeter/EMcalorimeterColl");
    fHCHCID = sdMan->GetCollectionID(colName="HadCalorimeter/HadCalorimeterColl");
    fMessenger = new B5EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EventAction::~B5EventAction()
{
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5EventAction::BeginOfEventAction(const G4Event*)
{
    if ( fHHC1ID == -1 ) {
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
    G4HCofThisEvent * HCE = event->GetHCofThisEvent();
    B5HodoscopeHitsCollection* fHHC1 = 0;
    B5HodoscopeHitsCollection* fHHC2 = 0;
    B5DriftChamberHitsCollection* fDHC1 = 0;
    B5DriftChamberHitsCollection* fDHC2 = 0;
    B5EmCalorimeterHitsCollection* ECHC = 0;
    B5HadCalorimeterHitsCollection* HCHC = 0;
    if(HCE)
    {
        fHHC1 = (B5HodoscopeHitsCollection*)(HCE->GetHC(fHHC1ID));
        fHHC2 = (B5HodoscopeHitsCollection*)(HCE->GetHC(fHHC2ID));
        fDHC1 = (B5DriftChamberHitsCollection*)(HCE->GetHC(fDHC1ID));
        fDHC2 = (B5DriftChamberHitsCollection*)(HCE->GetHC(fDHC2ID));
        ECHC = (B5EmCalorimeterHitsCollection*)(HCE->GetHC(fECHCID));
        HCHC = (B5HadCalorimeterHitsCollection*)(HCE->GetHC(fHCHCID));
    }
    
    // Fill histograms
    
    // get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
    if (fDHC1)
    {
        int n_hit = fDHC1->entries();
        analysisManager->FillH1(0, n_hit);

        for(int i1=0;i1<n_hit;i1++)
        {
            B5DriftChamberHit* aHit = (*fDHC1)[i1];
            G4ThreeVector localPos = aHit->GetLocalPos();
            analysisManager->FillH2(0, localPos.y(), localPos.x());
        }
    }
    if (fDHC2)
    {
        int n_hit = fDHC2->entries();
        analysisManager->FillH1(1, n_hit);

        for(int i1=0;i1<n_hit;i1++)
        {
            B5DriftChamberHit* aHit = (*fDHC2)[i1];
            G4ThreeVector localPos = aHit->GetLocalPos();
            analysisManager->FillH2(1, localPos.y(), localPos.x());
        }
    }
    
    // Fill the tuple
    // could we fill on condition ?
    
    if (fDHC1) analysisManager->FillNtupleIColumn(0, fDHC1->entries());
    if (fDHC2) analysisManager->FillNtupleIColumn(1, fDHC1->entries());
    
    if (ECHC)
    {
        int iHit = 0;
        double totalE = 0.;
        for(int i1=0;i1<80;i1++)
        {
            B5EmCalorimeterHit* aHit = (*ECHC)[i1];
            double eDep = aHit->GetEdep();
            if(eDep>0.)
            {
                iHit++;
                totalE += eDep;
            }
        }
        analysisManager->FillNtupleDColumn(2, totalE);
        
        if (fHHC1 && fHHC2 && fHHC1->entries()==1 && fHHC2->entries()==1)
        {
            double tof = (*fHHC2)[0]->GetTime() - (*fHHC1)[0]->GetTime();
            analysisManager->FillH2(2, totalE, tof);
        }
    }
    if (HCHC)
    {
        int iHit = 0;
        double totalE = 0.;
        for(int i1=0;i1<20;i1++)
        {
            B5HadCalorimeterHit* aHit = (*HCHC)[i1];
            double eDep = aHit->GetEdep();
            if(eDep>0.)
            {
                iHit++;
                totalE += eDep;
            }
        }
        analysisManager->FillNtupleDColumn(3, totalE);
    }
    if (fHHC1 && fHHC1->entries()==1) 
    {
      analysisManager->FillNtupleDColumn(4,(*fHHC1)[0]->GetTime());
    }
      
    if (fHHC2 && fHHC2->entries()==1) 
    {
      analysisManager->FillNtupleDColumn(5,(*fHHC2)[0]->GetTime());
    }
      
    analysisManager->AddNtupleRow();  
    
    // Diagnostics
    
    if (fVerboseLevel==0 || event->GetEventID() % fVerboseLevel != 0) return;
    
    G4PrimaryParticle* primary = event->GetPrimaryVertex(0)->GetPrimary(0);
    G4cout << G4endl
    << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
    << primary->GetG4code()->GetParticleName()
    << " " << primary->GetMomentum() << G4endl;
    
    if(fHHC1)
    {
        int n_hit = fHHC1->entries();
        G4cout << "Hodoscope 1 has " << n_hit << " hits." << G4endl;
        for(int i1=0;i1<n_hit;i1++)
        {
            B5HodoscopeHit* aHit = (*fHHC1)[i1];
            aHit->Print();
        }
    }
    if(fHHC2)
    {
        int n_hit = fHHC2->entries();
        G4cout << "Hodoscope 2 has " << n_hit << " hits." << G4endl;
        for(int i1=0;i1<n_hit;i1++)
        {
            B5HodoscopeHit* aHit = (*fHHC2)[i1];
            aHit->Print();
        }
    }
    if(fDHC1)
    {
        int n_hit = fDHC1->entries();
        G4cout << "Drift Chamber 1 has " << n_hit << " hits." << G4endl;
        for(int i2=0;i2<5;i2++)
        {
            for(int i1=0;i1<n_hit;i1++)
            {
                B5DriftChamberHit* aHit = (*fDHC1)[i1];
                if(aHit->GetLayerID()==i2) aHit->Print();
            }
        }
    }
    if(fDHC2)
    {
        int n_hit = fDHC2->entries();
        G4cout << "Drift Chamber 2 has " << n_hit << " hits." << G4endl;
        for(int i2=0;i2<5;i2++)
        {
            for(int i1=0;i1<n_hit;i1++)
            {
                B5DriftChamberHit* aHit = (*fDHC2)[i1];
                if(aHit->GetLayerID()==i2) aHit->Print();
            }
        }
    }
    if(ECHC)
    {
        int iHit = 0;
        double totalE = 0.;
        for(int i1=0;i1<80;i1++)
        {
            B5EmCalorimeterHit* aHit = (*ECHC)[i1];
            double eDep = aHit->GetEdep();
            if(eDep>0.)
            {
                iHit++;
                totalE += eDep;
            }
        }
        G4cout << "EM Calorimeter has " << iHit << " hits. Total Edep is "
        << totalE/MeV << " (MeV)" << G4endl;
    }
    if(HCHC)
    {
        int iHit = 0;
        double totalE = 0.;
        for(int i1=0;i1<20;i1++)
        {
            B5HadCalorimeterHit* aHit = (*HCHC)[i1];
            double eDep = aHit->GetEdep();
            if(eDep>0.)
            {
                iHit++;
                totalE += eDep;
            }
        }
        G4cout << "Hadron Calorimeter has " << iHit << " hits. Total Edep is "
        << totalE/MeV << " (MeV)" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
