#include "G4PointReplicaCreator.hh"

G4PointReplicaCreator G4PointReplicaCreator::csc;

G4PointReplicaCreator::G4PointReplicaCreator(){G4GeometryTable::RegisterObject(this);}

G4PointReplicaCreator::~G4PointReplicaCreator(){}

void G4PointReplicaCreator::CreateG4Geometry(STEPentity& Ent){}

void G4PointReplicaCreator::CreateSTEPGeometry(void* G4obj){}
