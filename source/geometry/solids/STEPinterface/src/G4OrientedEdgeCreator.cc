#include "G4OrientedEdgeCreator.hh"

G4OrientedEdgeCreator G4OrientedEdgeCreator::csc;

G4OrientedEdgeCreator::G4OrientedEdgeCreator(){G4GeometryTable::RegisterObject(this);}

G4OrientedEdgeCreator::~G4OrientedEdgeCreator(){}

void G4OrientedEdgeCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4int orientation;
  G4Curve* crv;
  
  G4String attrName("edge_element");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  // Get curve
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  crv = (G4Curve*)tmp; 

  // Get orientation
  Attr = Ent.NextAttribute();
  orientation = *Attr->ptr.i;// INTEGER_TYPE

  crv->SetSameSense(orientation);

  createdObject = crv;
}

void G4OrientedEdgeCreator::CreateSTEPGeometry(void* G4obj)
{

}
