#include "G4PointOnSurfaceCreator.hh"

G4PointOnSurfaceCreator G4PointOnSurfaceCreator::csc;

G4PointOnSurfaceCreator::G4PointOnSurfaceCreator(){G4GeometryTable::RegisterObject(this);}

G4PointOnSurfaceCreator::~G4PointOnSurfaceCreator(){}

void G4PointOnSurfaceCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4Surface* srf=0;
  G4double pval1 = 0;
  G4double pval2 = 0;
  
  Ent.ResetAttributes();
  STEPattribute* Attr = Ent.NextAttribute();
  while(Attr->NonRefType() == STRING_TYPE ||
	Attr->NonRefType() == sdaiSTRING )
    Attr = Ent.NextAttribute();	

  // Get basis surface
    STEPentity* TmpEnt= *Attr->ptr.c;
  void * tmp =G4GeometryTable::CreateObject(*TmpEnt);
  srf = (G4Surface*)tmp;

  // Get parameter values
    Attr = Ent.NextAttribute();	
  pval1 = *Attr->ptr.r;
  Attr = Ent.NextAttribute();	
  pval2 = *Attr->ptr.r;

  G4double* dtmp = new G4double[2];
  dtmp[0] = pval1;
  dtmp[1] = pval2;  
  
  createdObject = dtmp;  
}

void G4PointOnSurfaceCreator::CreateSTEPGeometry(void* G4obj){}
