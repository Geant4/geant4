// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FaceBoundCreator.cc,v 1.4 2000-11-20 18:17:30 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4FaceBoundCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4FaceBoundCreator.hh"
#include "G4GeometryTable.hh"

G4FaceBoundCreator G4FaceBoundCreator::csc;

G4FaceBoundCreator::G4FaceBoundCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4FaceBoundCreator::~G4FaceBoundCreator() {}

void G4FaceBoundCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4String attrName("bound");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  // Get curve
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  if (!tmp)
    G4cerr << "WARNING - G4FaceBoundCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL curve vector (G4CurveVector) !" << G4endl
	   << "\tEntity NOT created." << G4endl;

  // L. Broglia
  // Mistake : the created object tmp is a G4CurveVector
  // G4Curve* crv = (G4Curve*)tmp; 
  G4CurveVector* crv = (G4CurveVector*)tmp; 

  Attr = Ent.NextAttribute();
  //orientation = *Attr->ptr.i;
  //crv->SetSameSense(orientation);
  createdObject = crv;
}

void G4FaceBoundCreator::CreateSTEPGeometry(void* G4obj)
{
}
