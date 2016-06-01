#include "G4ConnectedFaceSetCreator.hh"

G4ConnectedFaceSetCreator G4ConnectedFaceSetCreator::csc;

G4ConnectedFaceSetCreator::G4ConnectedFaceSetCreator(){G4GeometryTable::RegisterObject(this);}

G4ConnectedFaceSetCreator::~G4ConnectedFaceSetCreator(){}

void G4ConnectedFaceSetCreator::CreateG4Geometry(STEPentity& Ent)
{

}

void G4ConnectedFaceSetCreator::CreateSTEPGeometry(void* G4obj)
{

}
