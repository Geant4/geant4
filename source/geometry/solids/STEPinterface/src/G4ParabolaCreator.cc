#include "G4ParabolaCreator.hh"

G4ParabolaCreator G4ParabolaCreator::csc;

G4ParabolaCreator::G4ParabolaCreator(){G4GeometryTable::RegisterObject(this);}

G4ParabolaCreator::~G4ParabolaCreator(){}

void G4ParabolaCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4double focal_dist;
  
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
  focal_dist = *Attr->ptr.r;

}

void G4ParabolaCreator::CreateSTEPGeometry(void* G4obj)
{

}
