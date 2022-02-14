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
/// \file SAXSSensitiveDetector.cc
/// \brief Implementation of the SAXSSensitiveDetector class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "SAXSSensitiveDetector.hh"
#include "SAXSSensitiveDetectorHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SAXSSensitiveDetector::SAXSSensitiveDetector(G4String name):
  G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="collection");
  fHCID = -1;
  
  fVarStopAndKill = true;
  
  fSDMessenger = new SAXSSensitiveDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SAXSSensitiveDetector::~SAXSSensitiveDetector() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SAXSSensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
  fHitsCollection = 
          new SensitiveDetectorHitsCollection(SensitiveDetectorName,collectionName[0]);
  if(fHCID<0) 
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);  
  HCE->AddHitsCollection(fHCID,fHitsCollection);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

bool SAXSSensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{   
  G4Track* vTrack = aStep->GetTrack();
  
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  
  G4ThreeVector worldPosition = preStepPoint->GetPosition();        
  G4ThreeVector localPosition = preStepPoint->GetTouchableHandle()->
                GetHistory()->GetTopTransform().TransformPoint(worldPosition);
  
  G4ThreeVector mom = preStepPoint->GetMomentumDirection();
  
  //if the partcile is not on the boundary of the sensitive detector, return
  if (!(preStepPoint->GetStepStatus() == fGeomBoundary)) {return false;}        
  
  G4int vType = -1;
  if(preStepPoint->GetMass() == 0 && preStepPoint->GetCharge() == 0) {
    vType = 0;
  }
  else if(preStepPoint->GetMass() != 0 && preStepPoint->GetCharge() == 0) {
    vType = 1;
  }
  else if(preStepPoint->GetMass() != 0 && preStepPoint->GetCharge() > 0) {
    vType = 2;
  }
  else if(preStepPoint->GetMass() != 0 && preStepPoint->GetCharge() < 0) {
    vType = 3;
  }

  //insert a hit        
  SAXSSensitiveDetectorHit* aHit = new SAXSSensitiveDetectorHit();
  aHit->SetPos(localPosition);
  aHit->SetMom(mom);
  aHit->SetType(vType);        
  aHit->SetEnergy(preStepPoint->GetKineticEnergy()); 
  aHit->SetTime(preStepPoint->GetGlobalTime());
  aHit->SetTrackID(vTrack->GetTrackID());
  aHit->SetWeight(vTrack->GetWeight());
  fHitsCollection->insert(aHit);
  
  if (fVarStopAndKill) {vTrack->SetTrackStatus(fStopAndKill);}
  
  return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SAXSSensitiveDetector::EndOfEvent(G4HCofThisEvent*) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

