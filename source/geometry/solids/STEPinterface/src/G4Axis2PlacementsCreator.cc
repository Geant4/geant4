#include "G4Axis2PlacementsCreator.hh"

G4Axis2PlacementsCreator G4Axis2PlacementsCreator::csc;

G4Axis2PlacementsCreator::G4Axis2PlacementsCreator(){
  // G4cerr << " G4Axis2PlacementsCreator default constructor " 
  // << " calling G4GeometryTable::RegisterObject(this) " << G4endl;
  G4GeometryTable::RegisterObject(this);
}

G4Axis2PlacementsCreator::~G4Axis2PlacementsCreator(){}

void G4Axis2PlacementsCreator::CreateG4Geometry(STEPentity& Ent)
{

}

void G4Axis2PlacementsCreator::CreateSTEPGeometry(void* G4obj)
{

}
