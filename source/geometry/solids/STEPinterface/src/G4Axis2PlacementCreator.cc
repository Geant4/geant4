#include "G4Axis2PlacementCreator.hh"

G4Axis2PlacementCreator G4Axis2PlacementCreator::csc;

G4Axis2PlacementCreator::G4Axis2PlacementCreator(){G4GeometryTable::RegisterObject(this);}

G4Axis2PlacementCreator::~G4Axis2PlacementCreator(){}

void G4Axis2PlacementCreator::CreateG4Geometry(STEPentity& Ent)
{

}

void G4Axis2PlacementCreator::CreateSTEPGeometry(void* G4obj)
{

}
