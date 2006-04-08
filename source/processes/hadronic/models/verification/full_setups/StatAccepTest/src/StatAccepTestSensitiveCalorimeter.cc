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


void StatAccepTestSensitiveCalorimeter::Initialize(G4HCofThisEvent* ) {
  calCollection = new StatAccepTestCalorimeterHitsCollection
    (SensitiveDetectorName,collectionName[0]);
  positionID.clear();
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    positionID.push_back( -1 );
  }
}


G4bool StatAccepTestSensitiveCalorimeter::
ProcessHits(G4Step* aStep, G4TouchableHistory* ) {

  G4double edep = aStep->GetTotalEnergyDeposit() * aStep->GetTrack()->GetWeight();
  // Multiply the energy deposit with the weight of the track,
  // to allow the use of biasing.

  if ( edep < 0.001*eV ) return false;

  // -----------------------------------------------------------------
  // 6-Apr-2006 : In the case of Scintillator as active medium,
  //              implement by hand Birks' law, using the expression
  //              and the coefficients taken from the paper 
  //              NIM 80 (1970) 239-244 for the organic scintillator
  //              NE-102:
  //                                     S*dE/dr 
  //                    dL/dr = -----------------------------------
  //                               1 + C1*(dE/dr) + C2*(dE/dr)^2 
  //              with:
  //                    C1 = 1.29 x 10^-2  g*cm^-2*MeV^-1 
  //                    C2 = 9.59 x 10^-6  g^2*cm^-4*MeV^-2 
  //              These are the same values used by ATLAS TileCal 
  //              and CMS HCAL (and also the default in Geant3).
  //              To get the "dE/dr" that appears in the formula,
  //              which has the dimensions
  //                   [ dE/dr ] = MeV * cm^2 / g
  //              we have to divide the energy deposit in MeV by the
  //              product of the step length (in cm) and the density
  //              of the scintillator, which is 1.032*g/cm3 .
  //              Of course, in our case we use only the denominator
  //              of the Birks' formula above, because we do not have
  //              digitization, i.e. we only have energy deposit but
  //              not conversion to photons.
  //              If you do not want to apply Birks' law, simply set
  //              isBirksOn = false.
  bool isBirksOn = false;  //***LOOKHERE***
  if ( isBirksOn ) {
    double C1 = 1.29e-2;                                  // [g*cm^-2*MeV^-1]  
    double C2 = 9.59e-6;                                  // [g^2*cm^-4*MeV^-2] 
    double rho = 1.032;                                   // [g/cm^3]
    double stepLength_cm = aStep->GetStepLength() / cm ;  // [cm] 
    if ( stepLength_cm > 1.0e-8 ) {
      double dedx = edep / ( rho * stepLength_cm );       // [MeV*cm^2/g]
      double birksFactor = 1.0 / ( 1.0 + C1*dedx + C2*dedx*dedx );
      //std::cout << " birksFactor=" << birksFactor << std::endl; //***DEBUG***
      edep = edep * birksFactor;
    }
  }
  // -----------------------------------------------------------------

  G4int replica = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(1);

  G4double radius = sqrt( aStep->GetPreStepPoint()->GetPosition().x() *
                          aStep->GetPreStepPoint()->GetPosition().x() +
                          aStep->GetPreStepPoint()->GetPosition().y() *
                          aStep->GetPreStepPoint()->GetPosition().y() );  

  // G4cout << " StatAccepTestSensitiveCalorimeter::ProcessHits : DEBUG Info" << G4endl
  //        << "\t z = " << aStep->GetPreStepPoint()->GetPosition().z()/mm 
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

  // We need to fill the shower profile information at this point,
  // because in the EventActiononly the hit collection is available
  // and the hits do not have information about the radius of each
  // energy deposition that contributed.
  StatAccepTestAnalysis* analysis = StatAccepTestAnalysis::getInstance();
  if ( analysis ) {
    analysis->fillShowerProfile( replica, radius, edep );
  }

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

