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
// $Id: A01Hodoscope.cc,v 1.3 2002-12-13 11:34:34 gunter Exp $
// --------------------------------------------------------------
//
#include "A01Hodoscope.hh"
#include "A01HodoscopeHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

A01Hodoscope::A01Hodoscope(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="hodoscopeColl");
  HCID = -1;
}

A01Hodoscope::~A01Hodoscope(){;}

void A01Hodoscope::Initialize(G4HCofThisEvent*HCE)
{
  hitsCollection = new A01HodoscopeHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); }
  HCE->AddHitsCollection(HCID,hitsCollection);
}

G4bool A01Hodoscope::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return true;

  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4int copyNo = theTouchable->GetVolume()->GetCopyNo();
  G4double hitTime = preStepPoint->GetGlobalTime();

  // check if this finger already has a hit
  G4int ix = -1;
  for(G4int i=0;i<hitsCollection->entries();i++)
  {
    if((*hitsCollection)[i]->GetID()==copyNo)
    {
      ix = i;
      break;
    }
  }
  // if it has, then take the earlier time
  if(ix>=0)
  {
    if((*hitsCollection)[ix]->GetTime()>hitTime)
    { (*hitsCollection)[ix]->SetTime(hitTime); }
  }
  else
  // if not, create a new hit and set it to the collection
  {
    A01HodoscopeHit* aHit = new A01HodoscopeHit(copyNo,hitTime);
    G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
    aHit->SetLogV(thePhysical->GetLogicalVolume());
    G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
    aTrans.Invert();
    aHit->SetRot(aTrans.NetRotation());
    aHit->SetPos(aTrans.NetTranslation());
    hitsCollection->insert( aHit );
  }

  return true;
}

void A01Hodoscope::EndOfEvent(G4HCofThisEvent*HCE)
{;}

