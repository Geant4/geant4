#include "G4Axis1PlacementCreator.hh"

G4Axis1PlacementCreator G4Axis1PlacementCreator::csc;

G4Axis1PlacementCreator::G4Axis1PlacementCreator(){
  G4GeometryTable::RegisterObject(this);
}

G4Axis1PlacementCreator::~G4Axis1PlacementCreator(){}

void G4Axis1PlacementCreator::CreateG4Geometry(STEPentity& Ent)
{

}

void G4Axis1PlacementCreator::CreateSTEPGeometry(void* G4obj)
{

}
