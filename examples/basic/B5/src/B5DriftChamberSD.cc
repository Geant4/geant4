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
/// \file B5DriftChamber.cc
/// \brief Implementation of the B5DriftChamber class

#include "B5DriftChamberSD.hh"
#include "B5DriftChamberHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DriftChamberSD::B5DriftChamberSD(G4String name)
:G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
    G4String HCname = "driftChamberColl";
    collectionName.insert(HCname);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DriftChamberSD::~B5DriftChamberSD(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DriftChamberSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new B5DriftChamberHitsCollection
    (SensitiveDetectorName,collectionName[0]);
    if(fHCID<0)
    { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }
    hce->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B5DriftChamberSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    G4double charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
    if(charge==0.) return true;
    
    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4VPhysicalVolume* theMotherPhysical = theTouchable->GetVolume(1); // mother
    G4int copyNo = theMotherPhysical->GetCopyNo();
    G4ThreeVector worldPos = preStepPoint->GetPosition();
    G4ThreeVector localPos
    = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
    
    B5DriftChamberHit* hit = new B5DriftChamberHit(copyNo);
    hit->SetWorldPos(worldPos);
    hit->SetLocalPos(localPos);
    hit->SetTime(preStepPoint->GetGlobalTime());
    
    fHitsCollection->insert(hit);
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
