#include "G4CurveCreator.hh"

G4CurveCreator G4CurveCreator::csc;

G4CurveCreator::G4CurveCreator(){G4GeometryTable::RegisterObject(this);}

G4CurveCreator::~G4CurveCreator(){}

void G4CurveCreator::CreateG4Geometry(STEPentity& Ent){}

void G4CurveCreator::CreateSTEPGeometry(void* G4obj){}
