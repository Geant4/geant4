#include "G4ManifoldSolidBrepCreator.hh"

G4ManifoldSolidBrepCreator G4ManifoldSolidBrepCreator::csc;

G4ManifoldSolidBrepCreator::G4ManifoldSolidBrepCreator(){G4GeometryTable::RegisterObject(this);}

G4ManifoldSolidBrepCreator::~G4ManifoldSolidBrepCreator(){}

void G4ManifoldSolidBrepCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4BREPSolid* sld;
  
  G4String attrName("outer");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  // Get solid
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  sld = (G4BREPSolid*)tmp; 
  createdObject = sld;
}

void G4ManifoldSolidBrepCreator::CreateSTEPGeometry(void* G4obj)
{

}
