
#include "MyTrackerSD.hh"
#include "MyTrackerHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

MyTrackerSD::MyTrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="EvenCollection");
  collectionName.insert(HCname="OddCollection");
}

MyTrackerSD::~MyTrackerSD(){;}

void MyTrackerSD::Initialize(G4HCofThisEvent*HCE)
{
  EvenCollection = new 
        MyTrackerHitsCollection(SensitiveDetectorName,collectionName[0]); 
  OddCollection = new 
        MyTrackerHitsCollection(SensitiveDetectorName,collectionName[1]); 
}

G4bool MyTrackerSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  const G4VPhysicalVolume* physVol 
    = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4double edep = aStep->GetTotalEnergyDeposit();

  MyTrackerHit* newHit = new MyTrackerHit();
  newHit->SetEdep( edep );
  newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  if( physVol->GetName() == "TPhys1" )
  { EvenCollection->insert( newHit ); }
  else if( physVol->GetName() == "TPhys2" )
  { OddCollection->insert( newHit ); }
  return true;
}

void MyTrackerSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int IDe = -1;
  static G4int IDo = -1;
  if(IDe<0)
  { IDe = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  if(IDo<0)
  { IDo = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[1]); }
  HCE->AddHitsCollection( IDe, EvenCollection );
  HCE->AddHitsCollection( IDo, OddCollection ); 
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

