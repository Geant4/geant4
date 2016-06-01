
#include "ExN04TrackerSD.hh"
#include "ExN04TrackerHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"

ExN04TrackerSD::ExN04TrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
}

ExN04TrackerSD::~ExN04TrackerSD(){;}

void ExN04TrackerSD::Initialize(G4HCofThisEvent*HCE)
{
  static int HCID = -1;
  trackerCollection = new ExN04TrackerHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID,trackerCollection);
}

G4bool ExN04TrackerSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;

  ExN04TrackerHit* newHit = new ExN04TrackerHit();
  newHit->SetEdep( edep );
  newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  trackerCollection->insert( newHit );

  return true;
}

void ExN04TrackerSD::EndOfEvent(G4HCofThisEvent*HCE)
{;}

void ExN04TrackerSD::clear()
{
} 

void ExN04TrackerSD::DrawAll()
{
} 

void ExN04TrackerSD::PrintAll()
{
} 

