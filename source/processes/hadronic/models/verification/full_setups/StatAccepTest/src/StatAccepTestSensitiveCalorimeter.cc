#include "StatAccepTestSensitiveCalorimeter.hh"
#include "StatAccepTestCalorimeterHit.hh"
#include "G4SDManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "StatAccepTestAnalysis.hh"


StatAccepTestSensitiveCalorimeter::
StatAccepTestSensitiveCalorimeter( const G4String name, const G4int numberOfReplicasIn )
  : G4VSensitiveDetector(name), calCollection(0), 
    numberOfReplicas( numberOfReplicasIn )
{
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    positionID.push_back( -1 );
  }
  G4String HCname;
  collectionName.insert(HCname="calCollection");
}


StatAccepTestSensitiveCalorimeter::~StatAccepTestSensitiveCalorimeter() {}


void StatAccepTestSensitiveCalorimeter::Initialize(G4HCofThisEvent* HCE) {
  calCollection = new StatAccepTestCalorimeterHitsCollection
    (SensitiveDetectorName,collectionName[0]);
  positionID.clear();
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    positionID.push_back( -1 );
  }
}


G4bool StatAccepTestSensitiveCalorimeter::
ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist) {

  G4double edep = aStep->GetTotalEnergyDeposit() * aStep->GetTrack()->GetWeight();
  // Multiply the energy deposit with the weight of the track,
  // to allow the use of biasing.

  if( edep < 0.001*eV ) return false;

  G4int replica = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(1);

  G4double radius = sqrt( aStep->GetPreStepPoint()->GetPosition().y() *
                          aStep->GetPreStepPoint()->GetPosition().y() +
                          aStep->GetPreStepPoint()->GetPosition().z() *
                          aStep->GetPreStepPoint()->GetPosition().z() );  

  // G4cout << " StatAccepTestSensitiveCalorimeter::ProcessHits : DEBUG Info" << G4endl
  //        << "\t x = " << aStep->GetPreStepPoint()->GetPosition().x()/mm 
  //        << " mm ;  replica = " << replica << " ;  radius = " << radius/mm 
  //        << " mm" << G4endl; //***DEBUG***

  if ( replica >= numberOfReplicas ) {
    G4cout << " StatAccepTestSensitiveCalorimeter::ProcessHits: ***ERROR*** "
           << " \t replica=" << replica 
	   << "  >=  numberOfReplicas=" << numberOfReplicas << G4endl;
    replica = numberOfReplicas - 1;  // To avoid a crash
  }

  if ( positionID[ replica ] == -1 ) {
    StatAccepTestCalorimeterHit* calHit = new StatAccepTestCalorimeterHit();
    calHit->SetEdep( edep );
    calHit->SetLayer( replica );
    G4int iPos = calCollection->insert( calHit );
    positionID[ replica ] = iPos - 1;
    // G4cout << " StatAccepTestSensitiveCalorimeter::ProcessHits: DEBUG Info " << G4endl
    //        << "\t New Calorimeter Hit : Energy = " 
    //        << edep / MeV << " MeV " << G4endl; //***DEBUG*** 
  } else { 
    (*calCollection)[ positionID[ replica ] ]->AddEdep(edep);
    // G4cout << " StatAccepTestSensitiveCalorimeter::ProcessHits : DEBUG Info " 
    //        << G4endl << "\t Old Calorimter Hit : Added energy = " 
    //        << " \t Added energy = " << edep / MeV << " MeV " << G4endl; //***DEBUG*** 
  }

  StatAccepTestAnalysis* analysis = StatAccepTestAnalysis::getInstance();
  if ( analysis ) analysis->fillShowerProfile( replica, radius, edep );

  return true;
}


void StatAccepTestSensitiveCalorimeter::EndOfEvent(G4HCofThisEvent* HCE) {
  static G4int HCID = -1;
  if ( HCID < 0 ) { 
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection( HCID, calCollection);
  if ( ! calCollection->entries() ) {
    // G4cout << " StatAccepTestSensitiveCalorimeter::EndOfEvent : DEBUG Info " << G4endl
    //        << "\t NO calorimeter hit! " << G4endl; //***DEBUG***
  }   
}


void StatAccepTestSensitiveCalorimeter::clear() {} 


void StatAccepTestSensitiveCalorimeter::DrawAll() {} 


void StatAccepTestSensitiveCalorimeter::PrintAll() {} 

