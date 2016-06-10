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

#include "ExExChEventAction.hh"

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

#include "ExExChSensitiveDetectorHit.hh"

#include "ExExChTrackingAction.hh"
#include "ExExChAnalysis.hh"

ExExChEventAction::ExExChEventAction()
{
    fSD_ID = -1;
    fVerboseLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExExChEventAction::~ExExChEventAction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChEventAction::BeginOfEventAction(const G4Event*){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChEventAction::EndOfEventAction(const G4Event* evt)
{
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    if(fSD_ID == -1) {
        G4String sdName;
        if(SDman->FindSensitiveDetector(sdName="telescope",0)){
            fSD_ID = SDman->GetCollectionID(sdName="telescope/collection");
        }
    }
    
    ExExChSensitiveDetectorHitsCollection* fSD = 0;
    
    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    
    if(HCE)
    {
        if(fSD_ID != -1){
            G4VHitsCollection* aHCSD = HCE->GetHC(fSD_ID);
            fSD = (ExExChSensitiveDetectorHitsCollection*)(aHCSD);
        }
    }
    
    G4ThreeVector SSDposition[3];
    SSDposition[0]= G4ThreeVector(0.,0.,0.);
    SSDposition[1]= G4ThreeVector(0.,0.,0.);
    SSDposition[2]= G4ThreeVector(0.,0.,0.);
    
    if(fSD)
    {
        int bTotalHits = 0;
        
        int n_hit_sd = fSD->entries();
        for(int i2=0;i2<3;i2++){
            for(int i1=0;i1<n_hit_sd;i1++)
            {
                ExExChSensitiveDetectorHit* aHit = (*fSD)[i1];
                if(aHit->GetLayerID()==i2) {
                    SSDposition[i2] = aHit->GetWorldPos();
                    bTotalHits++;
                }
            }
        }
        
        if(bTotalHits > 2){
            double fAngXin = (SSDposition[1].x() - SSDposition[0].x());
            fAngXin /= (SSDposition[1].z() - SSDposition[0].z());
            double fAngYin = (SSDposition[1].y() - SSDposition[0].y());
            fAngYin /= (SSDposition[1].z() - SSDposition[0].z());
            double fPosXin = SSDposition[1].x();
            double fPosYin = SSDposition[1].y();
            double fAngXout = (SSDposition[2].x() - SSDposition[1].x());
            fAngXout /= (SSDposition[2].z() - SSDposition[1].z());
            double fAngYout = (SSDposition[2].y() - SSDposition[1].y());
            fAngYout /= (SSDposition[2].z() - SSDposition[1].z());
            
            G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
            
            analysisManager->FillNtupleDColumn(0, fAngXin * 1.E6 * CLHEP::rad);
            analysisManager->FillNtupleDColumn(1, fAngYin * 1.E6 * CLHEP::rad);
            analysisManager->FillNtupleDColumn(2, fPosXin / CLHEP::mm);
            analysisManager->FillNtupleDColumn(3, fPosYin / CLHEP::mm);
            analysisManager->FillNtupleDColumn(4, fAngXout * 1.E6 * CLHEP::rad);
            analysisManager->FillNtupleDColumn(5, fAngYout * 1.E6 * CLHEP::rad);
            analysisManager->AddNtupleRow();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
