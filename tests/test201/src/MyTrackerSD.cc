// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyTrackerSD.cc,v 1.3 2000-05-26 13:11:42 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"

#include "MyTrackerSD.hh"
#include "MyTrackerHit.hh"

MyTrackerSD::MyTrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  collectionName.insert("TrackerCollection");
}

MyTrackerSD::~MyTrackerSD(){;}

void MyTrackerSD::Initialize(G4HCofThisEvent*)
{
  TrackerCollection = new MyTrackerHitsCollection(SensitiveDetectorName,
						  collectionName[0]); 
}

G4bool MyTrackerSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4ThreeVector hitPoint
   = ( aStep->GetPreStepPoint()->GetPosition()
     + aStep->GetPostStepPoint()->GetPosition() ) * 0.5;
  G4double edep = 1.0;

  if(verboseLevel>0)
  { G4cout << " New Tracker Hit at " << hitPoint << G4endl; }

  MyTrackerHit* newHit = new MyTrackerHit;
  newHit->SetEdep( edep );
  newHit->SetPos( hitPoint);
  TrackerCollection->insert(newHit);

  return true;
}

void MyTrackerSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID, TrackerCollection );
}

void MyTrackerSD::clear()
{
} 

void MyTrackerSD::DrawAll()
{
} 

void MyTrackerSD::PrintAll()
{
} 

