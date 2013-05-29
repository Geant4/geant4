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
/// \file analysis/A01/src/A01DriftChamber.cc
/// \brief Implementation of the A01DriftChamber class

#include "A01DriftChamber.hh"
#include "A01DriftChamberHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

A01DriftChamber::A01DriftChamber(G4String name)
:G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
    G4String HCname = "driftChamberColl";
    collectionName.insert(HCname);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

A01DriftChamber::~A01DriftChamber(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void A01DriftChamber::Initialize(G4HCofThisEvent*HCE)
{
    fHitsCollection = new A01DriftChamberHitsCollection
    (SensitiveDetectorName,collectionName[0]);
    if(fHCID<0)
    { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }
    HCE->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool A01DriftChamber::ProcessHits(G4Step*aStep,G4TouchableHistory* /*ROhist*/)
{
    G4double charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
    if(charge==0.) return true;
    
    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4VPhysicalVolume* theMotherPhysical = theTouchable->GetVolume(1); // mother
    G4int copyNo = theMotherPhysical->GetCopyNo();
    G4ThreeVector worldPos = preStepPoint->GetPosition();
    G4ThreeVector localPos
    = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
    
    A01DriftChamberHit* aHit = new A01DriftChamberHit(copyNo);
    aHit->SetWorldPos(worldPos);
    aHit->SetLocalPos(localPos);
    aHit->SetTime(preStepPoint->GetGlobalTime());
    
    fHitsCollection->insert(aHit);
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
