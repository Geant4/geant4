#include "MySensitiveCalorimeter.hh"
#include "MyCalorimeterHit.hh"
#include "G4SDManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"


MySensitiveCalorimeter::MySensitiveCalorimeter(G4String name)
  : G4VSensitiveDetector(name), 
  calCollection(0), hitExists(false) {
  G4String HCname;
  collectionName.insert(HCname="calCollection");
}


MySensitiveCalorimeter::~MySensitiveCalorimeter() {}


void MySensitiveCalorimeter::Initialize(G4HCofThisEvent* HCE) {
  calCollection = new MyCalorimeterHitsCollection
    (SensitiveDetectorName,collectionName[0]);
  hitExists = false;
}


G4bool MySensitiveCalorimeter::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist) {
  G4double edep = aStep->GetTotalEnergyDeposit();
  if( edep < 0.001*eV ) return false;

  if ( ! hitExists ) {
    MyCalorimeterHit* calHit = new MyCalorimeterHit();
    calHit->SetEdep( edep );
    calCollection->insert( calHit );
    hitExists = true;
    // G4cout << " MySensitiveCalorimeter::ProcessHits  >>> New Calorimeter Hit <<< " 
    //	      << G4endl; //***DEBUG*** 
  } else { 
    (*calCollection)[0]->AddEdep(edep);
    // G4cout << " MySensitiveCalorimeter::ProcessHits      >>> Added energy = " 
    //	      << edep / MeV << " (MeV) <<< " << G4endl; //***DEBUG*** 
  }

  return true;
}


void MySensitiveCalorimeter::EndOfEvent(G4HCofThisEvent* HCE) {
  static G4int HCID = -1;
  if ( HCID < 0 ) { 
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection( HCID, calCollection );
  if ( ! calCollection->entries() ) {
    // G4cout << " MySensitiveCalorimeter::EndOfEvent  >>> No calorimeter hit <<< " 
    //        << G4endl; //***DEBUG***
  }   
}


void MySensitiveCalorimeter::clear() {} 


void MySensitiveCalorimeter::DrawAll() {} 


void MySensitiveCalorimeter::PrintAll() {} 

