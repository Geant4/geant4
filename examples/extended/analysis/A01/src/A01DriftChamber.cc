//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: A01DriftChamber.cc,v 1.3 2002-12-13 11:34:33 gunter Exp $
// --------------------------------------------------------------
//
#include "A01DriftChamber.hh"
#include "A01DriftChamberHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"

A01DriftChamber::A01DriftChamber(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="driftChamberColl");
  HCID = -1;
}

A01DriftChamber::~A01DriftChamber(){;}

void A01DriftChamber::Initialize(G4HCofThisEvent*HCE)
{
  hitsCollection = new A01DriftChamberHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); }
  HCE->AddHitsCollection(HCID,hitsCollection);
}

G4bool A01DriftChamber::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
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

  hitsCollection->insert(aHit);

  return true;
}

void A01DriftChamber::EndOfEvent(G4HCofThisEvent*HCE)
{;}

