#include "G4ItemDefinedTransformationCreator.hh"

G4ItemDefinedTransformationCreator G4ItemDefinedTransformationCreator::csc;

G4ItemDefinedTransformationCreator::G4ItemDefinedTransformationCreator(){G4GeometryTable::RegisterObject(this);}

G4ItemDefinedTransformationCreator::~G4ItemDefinedTransformationCreator(){}

void G4ItemDefinedTransformationCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4PlacementVector* places = new G4PlacementVector();

  // L. Broglia
  // G4Placement* place=0;
  G4Axis2Placement3D* place=0;
  
  G4String attrName("transform_item_1");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);

  // L. Broglia
  // place = (G4Placement*)tmp;
  place = (G4Axis2Placement3D*)tmp;

  places->append(place);
  attrName="transform_item_2";
    Attr = GetNamedAttribute(attrName, Ent);
  TmpEnt = *Attr->ptr.c;
  tmp =G4GeometryTable::CreateObject(*TmpEnt);

  // L. Broglia
  // place = (G4Placement*)tmp;
  place = (G4Axis2Placement3D*)tmp;

  places->append(place);

  createdObject = places;
}

void G4ItemDefinedTransformationCreator::CreateSTEPGeometry(void* G4obj)
{

}
