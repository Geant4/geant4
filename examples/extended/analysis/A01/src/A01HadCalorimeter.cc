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
// $Id: A01HadCalorimeter.cc,v 1.3 2002-12-13 11:34:34 gunter Exp $
// --------------------------------------------------------------
//
#include "A01HadCalorimeter.hh"
#include "A01HadCalorimeterHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

A01HadCalorimeter::A01HadCalorimeter(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="HadCalorimeterColl");
  HCID = -1;
}

A01HadCalorimeter::~A01HadCalorimeter(){;}

void A01HadCalorimeter::Initialize(G4HCofThisEvent*HCE)
{
  hitsCollection = new A01HadCalorimeterHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); }
  HCE->AddHitsCollection(HCID,hitsCollection);

  // fill calorimeter hits with zero energy deposition
  for(int iColumn=0;iColumn<10;iColumn++)
  for(int iRow=0;iRow<2;iRow++)
  {
    A01HadCalorimeterHit* aHit = new A01HadCalorimeterHit();
    hitsCollection->insert( aHit );
  }
}

G4bool A01HadCalorimeter::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return true;

  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4VPhysicalVolume* theCellPhysical = theTouchable->GetVolume(2);
  G4int rowNo = theCellPhysical->GetCopyNo();
  G4VPhysicalVolume* theColumnPhysical = theTouchable->GetVolume(3);
  G4int columnNo = theColumnPhysical->GetCopyNo();
  G4int hitID = 2*columnNo+rowNo;
  A01HadCalorimeterHit* aHit = (*hitsCollection)[hitID];

  // check if it is first touch
  if(aHit->GetColumnID()<0)
  {
    aHit->SetColumnID(columnNo);
    aHit->SetRowID(rowNo);
    G4int depth = theTouchable->GetHistory()->GetDepth();
    G4AffineTransform aTrans = theTouchable->GetHistory()->GetTransform(depth-2);
    aTrans.Invert();
    aHit->SetRot(aTrans.NetRotation());
    aHit->SetPos(aTrans.NetTranslation());
  }
  // add energy deposition
  aHit->AddEdep(edep);

  return true;
}

void A01HadCalorimeter::EndOfEvent(G4HCofThisEvent*HCE)
{;}

