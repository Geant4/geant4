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

#include "EventAction.hh"

#include "G4RunManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "SensitiveDetectorHit.hh"

#include "Analysis.hh"

EventAction::EventAction():
sdht_ID(-1){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::BeginOfEventAction(const G4Event*){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event* evt){
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    
    G4ThreeVector ssd[3];
    ssd[0]= G4ThreeVector(0.,0.,0.);
    ssd[1]= G4ThreeVector(0.,0.,0.);
    ssd[2]= G4ThreeVector(0.,0.,0.);

    if(sdht_ID == -1) {
        G4String sdName;
        if(SDman->FindSensitiveDetector(sdName="telescope",0)){
            sdht_ID = SDman->GetCollectionID(sdName="telescope/collection");
        }
    }
    
    SensitiveDetectorHitsCollection* sdht = 0;
    G4HCofThisEvent *hce = evt->GetHCofThisEvent();
    
    if(hce){
        if(sdht_ID != -1){
            G4VHitsCollection* aHCSD = hce->GetHC(sdht_ID);
            sdht = (SensitiveDetectorHitsCollection*)(aHCSD);
        }
    }
    
    int bTotalHits = 0;
    if(sdht){
        
        int n_hit_sd = sdht->entries();
        for(int i2=0;i2<3;i2++){
            for(int i1=0;i1<n_hit_sd;i1++)
            {
                SensitiveDetectorHit* aHit = (*sdht)[i1];
                if(aHit->GetLayerID()==i2) {
                    ssd[i2] = aHit->GetWorldPos();
                    bTotalHits++;
                }
            }
        }
    }

    if(bTotalHits > 2){
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        G4double angXin  = (ssd[1].x() - ssd[0].x()) / (ssd[1].z() - ssd[0].z());
        G4double angYin  = (ssd[1].y() - ssd[0].y()) / (ssd[1].z() - ssd[0].z());

        analysisManager->FillNtupleDColumn(0, angXin * 1.E6 * CLHEP::rad);
        analysisManager->FillNtupleDColumn(1, angYin * 1.E6 * CLHEP::rad);

        double posXin = ssd[1].x() - angXin * ssd[1].z();
        double posYin = ssd[1].y() - angYin * ssd[1].z();

        analysisManager->FillNtupleDColumn(2, posXin / CLHEP::mm);
        analysisManager->FillNtupleDColumn(3, posYin / CLHEP::mm);

        G4double angXout = (ssd[2].x() - posXin) / (ssd[2].z());
        G4double angYout = (ssd[2].y() - posYin) / (ssd[2].z());
        analysisManager->FillNtupleDColumn(4, angXout * 1.E6 * CLHEP::rad);
        analysisManager->FillNtupleDColumn(5, angYout * 1.E6 * CLHEP::rad);
        
        analysisManager->AddNtupleRow();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
