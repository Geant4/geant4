#include "G4GeometricRepresentationContextCreator.hh"

G4GeometricRepresentationContextCreator G4GeometricRepresentationContextCreator
::csc;

G4GeometricRepresentationContextCreator::
G4GeometricRepresentationContextCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4GeometricRepresentationContextCreator::
~G4GeometricRepresentationContextCreator(){}

void G4GeometricRepresentationContextCreator::CreateG4Geometry(STEPentity& Ent)
{
  // Made by L. Broglia
  STEPattribute *Attr;

  G4String attrName("coordinate_space_dimension");
  Attr = GetNamedAttribute(attrName, Ent);

  // Get coordinate space dimension
  SdaiReal Tmpdim = *Attr->ptr.r; // ptr.r --> REAL_TYPE

  // void *tmp = G4GeometryTable::CreateObject(*TmpEnt);
  // place = (G4Axis2Placement3D*)tmp...;

}

void G4GeometricRepresentationContextCreator::CreateSTEPGeometry(void* G4obj)
{

}
