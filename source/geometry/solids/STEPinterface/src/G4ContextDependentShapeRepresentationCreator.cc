#include "G4ContextDependentShapeRepresentationCreator.hh"

G4ContextDependentShapeRepresentationCreator G4ContextDependentShapeRepresentationCreator::csc;

G4ContextDependentShapeRepresentationCreator::G4ContextDependentShapeRepresentationCreator(){G4GeometryTable::RegisterObject(this);}

G4ContextDependentShapeRepresentationCreator::~G4ContextDependentShapeRepresentationCreator(){}

void G4ContextDependentShapeRepresentationCreator::CreateG4Geometry(STEPentity& Ent)
{
  
  G4String attrName("representation_relation");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  G4PlacedSolidVector* psV = (G4PlacedSolidVector*)tmp;
  
  
  // L. Broglia
  // For the moment, the represented_product_relation is not needed
  // So don`t create it
  attrName="represented_product_relation";
  Attr = GetNamedAttribute(attrName, Ent);
  TmpEnt = *Attr->ptr.c;
  // tmp =G4GeometryTable::CreateObject(*TmpEnt);

  createdObject = psV;
}

void G4ContextDependentShapeRepresentationCreator::CreateSTEPGeometry(void* obj)
{

}
