#include "G4EllipseCreator.hh"

G4EllipseCreator G4EllipseCreator::csc;

G4EllipseCreator::G4EllipseCreator(){G4GeometryTable::RegisterObject(this);}

G4EllipseCreator::~G4EllipseCreator(){}

void G4EllipseCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4double semi1,semi2;
  
  // L. Broglia
  // G4Placement* place;
  G4Axis2Placement3D* place;

  Ent.ResetAttributes();
  STEPattribute* Attr = Ent.NextAttribute();
  while(Attr->NonRefType() == STRING_TYPE ||
	Attr->NonRefType() == sdaiSTRING )
    Attr = Ent.NextAttribute();	

  // Get the placement
    STEPentity* TmpEnt= *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);

  // L. Broglia
  // place = (G4Placement*)tmp;
  place = (G4Axis2Placement3D*)tmp;

  // get semi axises
  Attr = Ent.NextAttribute();	
  semi1 = *Attr->ptr.r;

  Attr = Ent.NextAttribute();	
  semi2 = *Attr->ptr.r;  
}

void G4EllipseCreator::CreateSTEPGeometry(void* G4obj)
{

}
