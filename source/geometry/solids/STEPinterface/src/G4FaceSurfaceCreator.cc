// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FaceSurfaceCreator.cc,v 1.3 2000-02-25 16:36:19 gcosmo Exp $
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
  if (!tmp)
    G4cerr << "WARNING - G4FaceSurfaceCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL pointer to G4Surface !" << G4endl
	   << "\tFace Surface NOT created." << G4endl;
  
  Attr = Ent.NextAttribute();
  //  sense = *Attr->ptr.b;
  createdObject = srf;
}

void G4FaceSurfaceCreator::CreateSTEPGeometry(void* G4obj)
{
}
