#include "MyCalorimeterHit.hh"
#include "G4ios.hh"


MyCalorimeterHit::MyCalorimeterHit() : edep(0.0) {}


MyCalorimeterHit::~MyCalorimeterHit() {}


MyCalorimeterHit::MyCalorimeterHit(const MyCalorimeterHit &right) :
  edep(right.edep) 
{}


const MyCalorimeterHit& 
MyCalorimeterHit::operator=(const MyCalorimeterHit &right) {
  edep = right.edep;
  return *this;
}


void MyCalorimeterHit::Draw() {}


void MyCalorimeterHit::Print() {}


