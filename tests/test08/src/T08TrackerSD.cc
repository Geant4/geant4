// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08TrackerSD.cc,v 1.1 1999-01-08 16:35:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "T08TrackerSD.hh"
#include "T08TrackerHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

T08TrackerSD::T08TrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
}

T08TrackerSD::~T08TrackerSD(){;}

void T08TrackerSD::Initialize(G4HCofThisEvent*HCE)
{
  trackerCollection = new T08TrackerHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, trackerCollection ); 
}

G4bool T08TrackerSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4Track* aTrack = aStep->GetTrack();
  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep==0.) return true;

  T08TrackerHit* newHit = new T08TrackerHit();
  newHit->SetEdep( edep );
  newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  trackerCollection->insert( newHit );

  return true;
}

void T08TrackerSD::EndOfEvent(G4HCofThisEvent*HCE)
{
}

void T08TrackerSD::clear()
{
} 

void T08TrackerSD::DrawAll()
{
} 

void T08TrackerSD::PrintAll()
{
} 

