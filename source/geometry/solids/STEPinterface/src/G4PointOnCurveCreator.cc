#include "G4PointOnCurveCreator.hh"

G4PointOnCurveCreator G4PointOnCurveCreator::csc;

G4PointOnCurveCreator::G4PointOnCurveCreator(){G4GeometryTable::RegisterObject(this);}

G4PointOnCurveCreator::~G4PointOnCurveCreator(){}

void G4PointOnCurveCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4Curve* crv=0;
  G4double pval = 0;
  
  Ent.ResetAttributes();
  STEPattribute* Attr = Ent.NextAttribute();
  while(Attr->NonRefType() == STRING_TYPE ||
	Attr->NonRefType() == sdaiSTRING )
    Attr = Ent.NextAttribute();	

  // Get basis curve
    STEPentity* TmpEnt= *Attr->ptr.c;
  void* tmp =G4GeometryTable::CreateObject(*TmpEnt);
  crv = (G4Curve*)tmp;  
  
  // Get parameter value
    Attr = Ent.NextAttribute();	
  pval = *Attr->ptr.r;
  
  createdObject = new G4double(pval);
}

void G4PointOnCurveCreator::CreateSTEPGeometry(void* G4obj){}





