#include "StatAccepTestCalorimeterHit.hh"
#include "G4ios.hh"


StatAccepTestCalorimeterHit::StatAccepTestCalorimeterHit() : 
  edep( 0.0 ), layer( -1 ) {}


StatAccepTestCalorimeterHit::~StatAccepTestCalorimeterHit() {}


StatAccepTestCalorimeterHit::
StatAccepTestCalorimeterHit( const StatAccepTestCalorimeterHit &right ) :
  G4VHit( right ), edep( right.edep ), layer( right.layer ) {}


const StatAccepTestCalorimeterHit& 
StatAccepTestCalorimeterHit::operator=( const StatAccepTestCalorimeterHit &right ) {
  edep = right.edep;
  layer = right.layer;
  return *this;
}


void StatAccepTestCalorimeterHit::Draw() {}


void StatAccepTestCalorimeterHit::Print() {}


