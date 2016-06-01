#include "G4OpenShellCreator.hh"

G4OpenShellCreator G4OpenShellCreator::csc;

G4OpenShellCreator::G4OpenShellCreator(){G4GeometryTable::RegisterObject(this);}

G4OpenShellCreator::~G4OpenShellCreator(){}

void G4OpenShellCreator::CreateG4Geometry(STEPentity& Ent)
{

}

void G4OpenShellCreator::CreateSTEPGeometry(void* G4obj)
{

}
