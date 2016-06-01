#include "G4PointCreator.hh"

G4PointCreator G4PointCreator::csc;

G4PointCreator::G4PointCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4PointCreator::~G4PointCreator(){}

void G4PointCreator::CreateG4Geometry(STEPentity& Ent){}

void G4PointCreator::CreateSTEPGeometry(void* G4obj){}
