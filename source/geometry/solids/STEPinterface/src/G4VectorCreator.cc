// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VectorCreator.cc,v 1.3 2000-02-25 16:36:20 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4VectorCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4VectorCreator.hh"
#include "G4GeometryTable.hh"

G4VectorCreator G4VectorCreator::csc;

G4VectorCreator::G4VectorCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4VectorCreator::~G4VectorCreator() {}

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

  G4ThreeVector* vec = 0;
  if (place)
    vec = new G4ThreeVector(Magnitude * place->x(),
			    Magnitude * place->y(),
			    Magnitude * place->z() );
  else
    G4cerr << "WARNING - G4VectorCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL orientation (G4ThreeVector) !" << G4endl
	   << "\tEntity NOT created." << G4endl;

  createdObject = vec;
}

void G4VectorCreator::CreateSTEPGeometry(void* G4obj)
{
}
