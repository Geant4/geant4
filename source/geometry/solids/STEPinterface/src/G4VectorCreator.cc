#include "G4VectorCreator.hh"

G4VectorCreator G4VectorCreator::csc;

G4VectorCreator::G4VectorCreator(){G4GeometryTable::RegisterObject(this);}

G4VectorCreator::~G4VectorCreator(){}

void G4VectorCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4double Magnitude;
  G4ThreeVector* place;

  Ent.ResetAttributes();
  STEPattribute* Attr = Ent.NextAttribute();
  Attr = Ent.NextAttribute();

  // get the orientation
  STEPentity* TmpEnt= *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  place = (G4ThreeVector*)tmp;

  // get Magnitude
  Attr = Ent.NextAttribute();
  Magnitude = *Attr->ptr.r;

  createdObject = new G4ThreeVector(Magnitude * place->x(),
				    Magnitude * place->y(),
				    Magnitude * place->z()
				    );
}

void G4VectorCreator::CreateSTEPGeometry(void* G4obj)
{

}
