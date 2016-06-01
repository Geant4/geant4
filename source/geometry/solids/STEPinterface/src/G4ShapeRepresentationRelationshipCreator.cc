#include "G4ShapeRepresentationRelationshipCreator.hh"

G4ShapeRepresentationRelationshipCreator G4ShapeRepresentationRelationshipCreator::csc;

G4ShapeRepresentationRelationshipCreator::G4ShapeRepresentationRelationshipCreator(){G4GeometryTable::RegisterObject(this);}

G4ShapeRepresentationRelationshipCreator::~G4ShapeRepresentationRelationshipCreator(){}

void G4ShapeRepresentationRelationshipCreator::CreateG4Geometry(STEPentity& Ent)
{

}

void G4ShapeRepresentationRelationshipCreator::CreateSTEPGeometry(void* G4obj)
{

}
