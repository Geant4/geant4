#include "G4FaceOuterBoundCreator.hh"

G4FaceOuterBoundCreator G4FaceOuterBoundCreator::csc;

G4FaceOuterBoundCreator::G4FaceOuterBoundCreator(){G4GeometryTable::RegisterObject(this);}

G4FaceOuterBoundCreator::~G4FaceOuterBoundCreator(){}

void G4FaceOuterBoundCreator::CreateSTEPGeometry(void* G4obj)
{

}
