#include "G4ConicCreator.hh"

G4ConicCreator G4ConicCreator::csc;

G4ConicCreator::G4ConicCreator(){G4GeometryTable::RegisterObject(this);}

G4ConicCreator::~G4ConicCreator(){}

void G4ConicCreator::CreateG4Geometry(STEPentity& Ent)
{
//  G4Placement *place;
  G4String attrName("edge_start");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  // Get the placement
    STEPentity* TmpEnt= *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
//  place = (G4Placement*)place;
  
}


void G4ConicCreator::CreateSTEPGeometry(void* G4obj)
{

}
