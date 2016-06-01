
#include "ExE02TrackerSD.hh"
#include "ExE02TrackerHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

ExE02TrackerSD::ExE02TrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="EvenCollection");
  collectionName.insert(HCname="OddCollection");
}

ExE02TrackerSD::~ExE02TrackerSD(){;}

void ExE02TrackerSD::Initialize(G4HCofThisEvent*HCE)
{
  EvenCollection = new 
        ExE02TrackerHitsCollection(SensitiveDetectorName,collectionName[0]); 
  OddCollection = new 
        ExE02TrackerHitsCollection(SensitiveDetectorName,collectionName[1]); 
}

G4bool ExE02TrackerSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  const G4VPhysicalVolume* physVol 
    = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4double edep = aStep->GetTotalEnergyDeposit();

  ExE02TrackerHit* newHit = new ExE02TrackerHit();
  newHit->SetEdep( edep );
  newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  if( physVol->GetName() == "TPhys1" )
  { EvenCollection->insert( newHit ); }
  else if( physVol->GetName() == "TPhys2" )
  { OddCollection->insert( newHit ); }
  return true;
}

void ExE02TrackerSD::EndOfEvent(G4HCofThisEvent*HCE)
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

void ExE02TrackerSD::clear()
{
} 

void ExE02TrackerSD::DrawAll()
{
} 

void ExE02TrackerSD::PrintAll()
{
} 

