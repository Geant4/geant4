#include "SimpleGeometry.hh"


SimpleGeometry::SimpleGeometry(G4VPhysicalVolume *vol) :
  fVolume(vol)
{}

SimpleGeometry::~SimpleGeometry()
{}

G4VPhysicalVolume* SimpleGeometry::Construct(){
  return fVolume;
}

