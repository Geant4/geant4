#include "G4VertexPointCreator.hh"

G4VertexPointCreator G4VertexPointCreator::csc;

G4VertexPointCreator::G4VertexPointCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4VertexPointCreator::~G4VertexPointCreator(){}

void G4VertexPointCreator::CreateG4Geometry(STEPentity& Ent)
{
  // G4Point3d *point;
  G4Point3D *point;

  G4String attrName("vertex_geometry");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  // Get the geometric point  
  STEPentity* TmpEnt= *Attr->ptr.c;
  void * tmp =G4GeometryTable::CreateObject(*TmpEnt);
  
  // point = (G4Point3d*)tmp;
  point = (G4Point3D*)tmp;
  createdObject = point;
}

void G4VertexPointCreator::CreateSTEPGeometry(void* G4obj){}



