// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02TrackerSD.cc,v 1.1 1999-01-07 16:05:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN02TrackerSD.hh"
#include "ExN02TrackerHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

ExN02TrackerSD::ExN02TrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
}

ExN02TrackerSD::~ExN02TrackerSD(){;}

void ExN02TrackerSD::Initialize(G4HCofThisEvent*HCE)
{
  trackerCollection = new ExN02TrackerHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, trackerCollection ); 
}

G4bool ExN02TrackerSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4Track* aTrack = aStep->GetTrack();
  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep==0.) return true;

  ExN02TrackerHit* newHit = new ExN02TrackerHit();
  newHit->SetEdep( edep );
  newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  trackerCollection->insert( newHit );

  return true;
}

void ExN02TrackerSD::EndOfEvent(G4HCofThisEvent*HCE)
{
}

void ExN02TrackerSD::clear()
{
} 

void ExN02TrackerSD::DrawAll()
{
} 

void ExN02TrackerSD::PrintAll()
{
} 

