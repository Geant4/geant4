#include "G4ProductDefinitionShapeCreator.hh"

G4ProductDefinitionShapeCreator G4ProductDefinitionShapeCreator::csc;

G4ProductDefinitionShapeCreator::G4ProductDefinitionShapeCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4ProductDefinitionShapeCreator::~G4ProductDefinitionShapeCreator(){}

void G4ProductDefinitionShapeCreator::CreateG4Geometry(STEPentity& Ent)
{
  // Made by L. Broglia
  STEPattribute *Attr;

  G4String attrName("description");
  Attr = GetNamedAttribute(attrName, Ent);

  // Get description
  SdaiString Tmpstring = *Attr->ptr.S; // ptr.S --> STRING_TYPE

  attrName = "definition";
  Attr = GetNamedAttribute(attrName, Ent);

  // Get definition
  SCLundefined Tmpdef = *Attr->ptr.u; // ptr.u --> UNKNOWN_TYPE

  // void *tmp = G4GeometryTable::CreateObject(*TmpEnt);
  // place = (G4Axis2Placement3D*)tmp;
}

void G4ProductDefinitionShapeCreator::CreateSTEPGeometry(void* G4obj)
{

}
