#include "G4FaceSurfaceCreator.hh"

G4FaceSurfaceCreator G4FaceSurfaceCreator::csc;

G4FaceSurfaceCreator::G4FaceSurfaceCreator(){G4GeometryTable::RegisterObject(this);}

G4FaceSurfaceCreator::~G4FaceSurfaceCreator(){}

void G4FaceSurfaceCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4Surface* srf=0;
  G4bool sense;
  
  G4String attrName("face_geometry");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  srf = (G4Surface*)tmp;
  
  Attr = Ent.NextAttribute();
  //  sense = *Attr->ptr.b;
  createdObject = srf;
}

void G4FaceSurfaceCreator::CreateSTEPGeometry(void* G4obj)
{

}



