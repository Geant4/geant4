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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#include "A01EmCalorimeter.hh"
#include "A01EmCalorimeterHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

A01EmCalorimeter::A01EmCalorimeter(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="EMcalorimeterColl");
  HCID = -1;
}

A01EmCalorimeter::~A01EmCalorimeter(){;}

void A01EmCalorimeter::Initialize(G4HCofThisEvent*HCE)
{
  hitsCollection = new A01EmCalorimeterHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); }
  HCE->AddHitsCollection(HCID,hitsCollection);

  // fill calorimeter hits with zero energy deposition
  for(int i=0;i<80;i++)
  {
    A01EmCalorimeterHit* aHit = new A01EmCalorimeterHit(i);
    hitsCollection->insert( aHit );
  }
}

G4bool A01EmCalorimeter::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return true;

  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
  G4int copyNo = thePhysical->GetCopyNo();

  A01EmCalorimeterHit* aHit = (*hitsCollection)[copyNo];
  // check if it is first touch
  if(!(aHit->GetLogV()))
  {
    // fill volume information
    aHit->SetLogV(thePhysical->GetLogicalVolume());
    G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
    aTrans.Invert();
    aHit->SetRot(aTrans.NetRotation());
    aHit->SetPos(aTrans.NetTranslation());
  }
  // add energy deposition
  aHit->AddEdep(edep);

  return true;
}

void A01EmCalorimeter::EndOfEvent(G4HCofThisEvent*HCE)
{;}

