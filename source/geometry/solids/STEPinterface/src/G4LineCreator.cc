#include "G4LineCreator.hh"

G4LineCreator G4LineCreator::csc;

G4LineCreator::G4LineCreator(){G4GeometryTable::RegisterObject(this);}

G4LineCreator::~G4LineCreator(){}

void G4LineCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4ThreeVector origin;
  G4ThreeVector dir;
  
  Ent.ResetAttributes();
  STEPattribute *Attr;
  Attr = Ent.NextAttribute();    
  while (Attr->NonRefType() == STRING_TYPE)
    Attr = Ent.NextAttribute();

  // Get the line origin  
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  origin = *(G4ThreeVector*)tmp;

  // Get the line direction  
  Attr = Ent.NextAttribute();  
  TmpEnt = *Attr->ptr.c;
  tmp =G4GeometryTable::CreateObject(*TmpEnt);
  dir = *(G4ThreeVector*)tmp;


  G4Line* line = new G4Line();
  line->Init(origin, dir);
  createdObject = line;
}

void G4LineCreator::CreateSTEPGeometry(void* G4obj)
{

}
