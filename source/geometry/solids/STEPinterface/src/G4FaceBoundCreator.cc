#include "G4FaceBoundCreator.hh"

G4FaceBoundCreator G4FaceBoundCreator::csc;

G4FaceBoundCreator::G4FaceBoundCreator(){G4GeometryTable::RegisterObject(this);}

G4FaceBoundCreator::~G4FaceBoundCreator(){}

void G4FaceBoundCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4int orientation;
  G4String attrName("bound");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  // Get curve
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);

  // L. Broglia
  // Mistake : the created object tmp is a G4CurveVector
  // G4Curve* crv = (G4Curve*)tmp; 
  G4CurveVector* crv = (G4CurveVector*)tmp; 

  Attr = Ent.NextAttribute();
  //orientation = *Attr->ptr.i;
  //crv->SetSameSense(orientation);
  createdObject = crv;

}

void G4FaceBoundCreator::CreateSTEPGeometry(void* G4obj)
{

}
