#include "G4PlaneCreator.hh"

G4PlaneCreator G4PlaneCreator::csc;

G4PlaneCreator::G4PlaneCreator(){G4GeometryTable::RegisterObject(this);}

G4PlaneCreator::~G4PlaneCreator(){}

void G4PlaneCreator::CreateG4Geometry(STEPentity& Ent)
{
  // L. Broglia
  // G4Placement* place;
  G4Axis2Placement3D* place;
  
  Ent.ResetAttributes();
  STEPattribute *Attr;
  Attr = Ent.NextAttribute();    
  while (Attr->NonRefType() == STRING_TYPE)
    Attr = Ent.NextAttribute();

  // Get placement
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  
  place = (G4Axis2Placement3D*)tmp;
  createdObject = new G4FPlane( place->GetRefDirection(),
				place->GetAxis(),
				place->GetLocation());
}

void G4PlaneCreator::CreateSTEPGeometry(void* G4obj)
{

}

