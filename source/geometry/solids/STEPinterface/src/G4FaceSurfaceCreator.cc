// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FaceSurfaceCreator.cc,v 1.2 2000-01-21 13:46:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4FaceSurfaceCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4FaceSurfaceCreator.hh"
#include "G4GeometryTable.hh"

G4FaceSurfaceCreator G4FaceSurfaceCreator::csc;

G4FaceSurfaceCreator::G4FaceSurfaceCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4FaceSurfaceCreator::~G4FaceSurfaceCreator() {}

void G4FaceSurfaceCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4Surface* srf=0;
  G4bool sense;
  
  G4String attrName("face_geometry");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  srf = (G4Surface*)tmp;
  
  Attr = Ent.NextAttribute();
  //  sense = *Attr->ptr.b;
  createdObject = srf;
}

void G4FaceSurfaceCreator::CreateSTEPGeometry(void* G4obj)
{
}
