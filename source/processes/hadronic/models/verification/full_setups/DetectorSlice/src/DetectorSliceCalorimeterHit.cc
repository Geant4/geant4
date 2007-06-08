#include "DetectorSliceCalorimeterHit.hh"
#include "G4ios.hh"


DetectorSliceCalorimeterHit::DetectorSliceCalorimeterHit() : 
  edep( 0.0 ) {}


DetectorSliceCalorimeterHit::~DetectorSliceCalorimeterHit() {}


DetectorSliceCalorimeterHit::
DetectorSliceCalorimeterHit( const DetectorSliceCalorimeterHit &right ) :
  G4VHit( right ), edep( right.edep ) {}


const DetectorSliceCalorimeterHit& 
DetectorSliceCalorimeterHit::operator=( const DetectorSliceCalorimeterHit &right ) {
  edep = right.edep;
  return *this;
}


void DetectorSliceCalorimeterHit::Draw() {}


void DetectorSliceCalorimeterHit::Print() {}


