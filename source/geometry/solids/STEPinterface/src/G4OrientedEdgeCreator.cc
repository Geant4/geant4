// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OrientedEdgeCreator.cc,v 1.2 2000-01-21 13:46:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4OrientedEdgeCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4OrientedEdgeCreator.hh"
#include "G4GeometryTable.hh"

G4OrientedEdgeCreator G4OrientedEdgeCreator::csc;

G4OrientedEdgeCreator::G4OrientedEdgeCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4OrientedEdgeCreator::~G4OrientedEdgeCreator() {}

void G4OrientedEdgeCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4int orientation;
  G4Curve* crv;
  
  G4String attrName("edge_element");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  // Get curve
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  crv = (G4Curve*)tmp; 

  // Get orientation
  Attr = Ent.NextAttribute();
  orientation = *Attr->ptr.i;// INTEGER_TYPE

  crv->SetSameSense(orientation);

  createdObject = crv;
}

void G4OrientedEdgeCreator::CreateSTEPGeometry(void* G4obj)
{
}
