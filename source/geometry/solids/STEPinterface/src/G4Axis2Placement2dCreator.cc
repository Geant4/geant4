#include "G4Axis2Placement2dCreator.hh"

G4Axis2Placement2dCreator G4Axis2Placement2dCreator::csc;

G4Axis2Placement2dCreator::G4Axis2Placement2dCreator(){G4GeometryTable::RegisterObject(this);}

G4Axis2Placement2dCreator::~G4Axis2Placement2dCreator(){}

void G4Axis2Placement2dCreator::CreateG4Geometry(STEPentity& Ent)
{

}

void G4Axis2Placement2dCreator::CreateSTEPGeometry(void* G4obj)
{

}
