// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PlaneCreator.cc,v 1.2 2000-01-21 13:46:04 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4PlaneCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4PlaneCreator.hh"
#include "G4GeometryTable.hh"
#include "G4FPlane.hh"

G4PlaneCreator G4PlaneCreator::csc;

G4PlaneCreator::G4PlaneCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4PlaneCreator::~G4PlaneCreator() {}

void G4PlaneCreator::CreateG4Geometry(STEPentity& Ent)
{
  // L. Broglia
  // G4Placement* place;
  G4Axis2Placement3D* place;
  
  Ent.ResetAttributes();
  STEPattribute *Attr;
  Attr = Ent.NextAttribute();    
  while (Attr->NonRefType() == STRING_TYPE)
    Attr = Ent.NextAttribute();

  // Get placement
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  
  place = (G4Axis2Placement3D*)tmp;
  createdObject = new G4FPlane( place->GetRefDirection(),
				place->GetAxis(),
				place->GetLocation());
}

void G4PlaneCreator::CreateSTEPGeometry(void* G4obj)
{
}

