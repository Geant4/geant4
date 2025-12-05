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
/// \file SensitiveDetector.cc
/// \brief Implementation of the SensitiveDetector class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "SensitiveDetector.hh"
#include "SensitiveDetectorHit.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SensitiveDetector::SensitiveDetector(G4String name):G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname="collection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
    fHitsCollection = 
        new SensitiveDetectorHitsCollection(SensitiveDetectorName,collectionName[0]);
    if (fHCID < 0)
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
    HCE->AddHitsCollection(fHCID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

bool SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{   
    G4Track* vTrack = aStep->GetTrack();
    
    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
                   
    if (!(preStepPoint->GetStepStatus() == fGeomBoundary)) {return true;}    

    G4String particle = 
        vTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
    
    G4ThreeVector position = preStepPoint->GetPosition();
       
    G4int detID = preStepPoint->GetTouchableHandle()->GetVolume()->GetCopyNo(); 

    G4ThreeVector momentum = preStepPoint->GetMomentum();            
    
    if (momentum.getZ() > 0) {    
        SensitiveDetectorHit* aHit = new SensitiveDetectorHit();
        aHit->SetParticle(particle);
        aHit->SetDetID(detID);
        aHit->SetPos(position);
        aHit->SetMom(momentum);    
        aHit->SetEnergy(preStepPoint->GetKineticEnergy());
        aHit->SetTime(preStepPoint->GetGlobalTime());
        aHit->SetWeight(vTrack->GetWeight());
        aHit->SetTrackID(vTrack->GetTrackID());
        aHit->SetTrackIDP(vTrack->GetParentID());
        fHitsCollection->insert(aHit);
    }
    
    return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

