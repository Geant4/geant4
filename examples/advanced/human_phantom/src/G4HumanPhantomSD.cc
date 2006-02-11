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

#include "G4HumanPhantomSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"


G4HumanPhantomSD::G4HumanPhantomSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="G4HumanPhantomHitsCollection");
}


G4HumanPhantomSD::~G4HumanPhantomSD(){ }


void G4HumanPhantomSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new G4HumanPhantomHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, trackerCollection ); 
}


G4bool G4HumanPhantomSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep==0.) return false;

  G4HumanPhantomHit* newHit = new G4HumanPhantomHit();
  newHit->SetTrackID   (aStep->GetTrack()->GetTrackID());
  newHit->SetBodyPartID(aStep->GetPreStepPoint()->GetTouchable()
		                                ->GetReplicaNumber());
  newHit->SetBodyPartName(aStep->GetPreStepPoint()->GetTouchable()
  			                          ->GetVolume()->GetLogicalVolume()
			                          ->GetName());
  newHit->SetEdep     (edep);
  newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());

  trackerCollection->insert( newHit );
  
  //  newHit->Print();
  //  newHit->Draw();

  return true;
}


void G4HumanPhantomSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel>0) { 
     G4int NbHits = trackerCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
            << " hits in the: " << G4endl;
     for (G4int i=0;i<NbHits;i++) (*trackerCollection)[i]->Print();
    } 
}
