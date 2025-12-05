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
// gpaterno, October 2025
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "EventAction.hh"
#include "RunAction.hh"
#include "SensitiveDetectorHit.hh"
#include "DetectorConstruction.hh"

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
#include "globals.hh"
#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
{
    //An instance of the DetectorConstruction
    const DetectorConstruction* detectorConstruction = 
        static_cast<const DetectorConstruction*>
            (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    fNSpheres = detectorConstruction->GetNSpheres();

    //Inizialize the Edep map in the Spheres (scored through SteppingAction)
    for (int i = 0; i < fNSpheres; i++) {    
        fEdepSpheres[i] = 0.;     
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
    //reset Edep in the Radiator and Converter Crystals
    fEdepRad = 0.;
    fEdepConv = 0.;

    //reset the Edep map in the Spheres (scored through SteppingAction)
    if (fVerboseLevel > 0) {
        G4int eventID = GetEventID();
        G4cout << "EventAction::BeginOfEventAction(), "
               << "EventID: " << eventID << G4endl;
    }

    for (int i = 0; i < fNSpheres; i++) {
        fEdepSpheres[i] = 0.;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* aEvent)
{   
    //instantiating The Sensitive Detector Manager
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    //Hit Detection System
    if (fSensitiveDetector_ID == -1) {
        G4String SensitiveDetectorName;
        if (SDman->FindSensitiveDetector(SensitiveDetectorName="det",0)) {
            fSensitiveDetector_ID = 
                SDman->GetCollectionID(SensitiveDetectorName="det/collection");
        }
    }

    SensitiveDetectorHitsCollection* sensitiveDetectorHC = 0;
    G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();

    if (HCE) {
        if (fSensitiveDetector_ID != -1) {
            G4VHitsCollection* aHC = HCE->GetHC(fSensitiveDetector_ID);
            sensitiveDetectorHC = (SensitiveDetectorHitsCollection*)(aHC);
        }
    }
    
    //instantiating The Analysis Manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   
    //get the Event Number
    G4int eventID = GetEventID();

    //filling the SD scoring ntuple
    if (sensitiveDetectorHC) {
        int vNumberOfHit = sensitiveDetectorHC->entries();
        for (int i=0; i<vNumberOfHit; i++) {
            SensitiveDetectorHit* aHit = (*sensitiveDetectorHC)[i];
            analysisManager->FillNtupleIColumn(0,0,aHit->GetDetID());
            analysisManager->FillNtupleSColumn(0,1,aHit->GetParticle());
            analysisManager->FillNtupleDColumn(0,2,aHit->GetPos().x()/CLHEP::mm);
            analysisManager->FillNtupleDColumn(0,3,aHit->GetPos().y()/CLHEP::mm);
            analysisManager->FillNtupleDColumn(0,4,aHit->GetMom().x()/CLHEP::MeV);
            analysisManager->FillNtupleDColumn(0,5,aHit->GetMom().y()/CLHEP::MeV);
            analysisManager->FillNtupleDColumn(0,6,aHit->GetMom().z()/CLHEP::MeV);
            analysisManager->FillNtupleDColumn(0,7,aHit->GetTime()/CLHEP::ns);
            analysisManager->FillNtupleIColumn(0,8,eventID);
            analysisManager->AddNtupleRow(0);
        }
    }

    //filling the Edep (in the Radiator Crystal) ntuple
    if (fEdepRad > 0) {
        analysisManager->FillNtupleDColumn(1,0,fEdepRad/MeV);
        analysisManager->FillNtupleIColumn(1,1,eventID);
        analysisManager->AddNtupleRow(1);
    }

    //filling the Edep (in the Converter Crystal) ntuple
    if (fEdepConv > 0) {
        analysisManager->FillNtupleDColumn(2,0,fEdepConv/MeV);
        analysisManager->FillNtupleIColumn(2,1,eventID);
        analysisManager->AddNtupleRow(2);
    }

    //Get the results accumulated through the SteppingAction class
    //and fill the local (one for thread) ntuple for Edep map in 
    //the Spheres of thew granular target.
    for (int i = 0; i < fNSpheres; i++) {    
        analysisManager->FillNtupleIColumn(3,0,i);
        analysisManager->FillNtupleDColumn(3,1,fEdepSpheres[i]/MeV);
        analysisManager->FillNtupleIColumn(3,2,eventID);
        analysisManager->AddNtupleRow(3);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEdepInSpheres(G4int volID, G4double edep)
{    
    G4double temp = fEdepSpheres.find(volID)->second;
    fEdepSpheres[volID] = temp + edep;

    if (fVerboseLevel > 1) {
        G4int eventID = GetEventID();
        G4cout << "Event: " << eventID
               << ", Edep[" << volID << "]: "
               << edep/MeV << " MeV" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

