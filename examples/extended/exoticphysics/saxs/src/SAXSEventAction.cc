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
/// \file SAXSEventAction.cc
/// \brief Implementation of the SAXSEventAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

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

#include "SAXSAnalysis.hh"
#include "SAXSEventAction.hh"
#include "SAXSSensitiveDetectorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SAXSEventAction::SAXSEventAction():
  G4UserEventAction(),
  fSensitiveDetector_ID(-1),
  fVerboseLevel(0),
  fNRi(0),
  fNCi(0),
  fNDi(0),
  fEventWeight(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SAXSEventAction::~SAXSEventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SAXSEventAction::BeginOfEventAction(const G4Event*)
{
  fNRi = 0;
  fNCi = 0;
  fNDi = 0;
  fEventWeight = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SAXSEventAction::EndOfEventAction(const G4Event* aEvent)
{   
  //instantiating The Sensitive Detector Manager
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  //Hit Detection System
  if (fSensitiveDetector_ID==-1) {
    G4String SensitiveDetectorName;
    if (SDman->FindSensitiveDetector(SensitiveDetectorName="det",0)) {
      fSensitiveDetector_ID =
        SDman->GetCollectionID(SensitiveDetectorName="det/collection");
    }
  }
  
  SensitiveDetectorHitsCollection* fSensitiveDetectorHC = 0;
  G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
  
  if (HCE) {
    if (fSensitiveDetector_ID != -1) {
      G4VHitsCollection* aHC = HCE->GetHC(fSensitiveDetector_ID);
      fSensitiveDetectorHC = (SensitiveDetectorHitsCollection*)(aHC);
    }
  }
  
  //instantiating The Analysis Manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  //get the Event Number
  G4int eventNumber =
    G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  
  //filling the SD scoring ntuple
  if (fSensitiveDetectorHC) {
    size_t vNumberOfHit = fSensitiveDetectorHC->entries();
    for (size_t i=0; i<vNumberOfHit; i++) {
      SAXSSensitiveDetectorHit* aHit = (*fSensitiveDetectorHC)[i];
      analysisManager->FillNtupleDColumn(0,0,aHit->GetEnergy()/CLHEP::keV);
      analysisManager->FillNtupleDColumn(0,1,aHit->GetPos().x()/CLHEP::mm);
      analysisManager->FillNtupleDColumn(0,2,aHit->GetPos().y()/CLHEP::mm);
      analysisManager->FillNtupleDColumn(0,3,aHit->GetPos().z()/CLHEP::mm);
      analysisManager->FillNtupleDColumn(0,4,aHit->GetMom().x());
      analysisManager->FillNtupleDColumn(0,5,aHit->GetMom().y());
      analysisManager->FillNtupleDColumn(0,6,aHit->GetMom().z());
      analysisManager->FillNtupleDColumn(0,7,aHit->GetTime()/CLHEP::ns);
      analysisManager->FillNtupleIColumn(0,8,aHit->GetType());
      analysisManager->FillNtupleIColumn(0,9,aHit->GetTrackID());
      analysisManager->FillNtupleIColumn(0,10,fNRi);
      analysisManager->FillNtupleIColumn(0,11,fNCi);
      analysisManager->FillNtupleIColumn(0,12,fNDi);
      analysisManager->FillNtupleIColumn(0,13,eventNumber);
      analysisManager->FillNtupleDColumn(0,14,aHit->GetWeight());
      analysisManager->AddNtupleRow(0);  
    }
  }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

