// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ContextDependentShapeRepresentationCreator.cc,v 1.2 2000-01-21 13:46:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ContextDependentShapeRepresentationCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4ContextDependentShapeRepresentationCreator.hh"
#include "G4GeometryTable.hh"

G4ContextDependentShapeRepresentationCreator
  G4ContextDependentShapeRepresentationCreator::csc;

G4ContextDependentShapeRepresentationCreator::
  G4ContextDependentShapeRepresentationCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4ContextDependentShapeRepresentationCreator::
  ~G4ContextDependentShapeRepresentationCreator() {}

void G4ContextDependentShapeRepresentationCreator::
  CreateG4Geometry(STEPentity& Ent)
{
  G4String attrName("representation_relation");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  G4PlacedSolidVector* psV = (G4PlacedSolidVector*)tmp;
    
  // L. Broglia
  // For the moment, the represented_product_relation is not needed
  // So don`t create it
  attrName="represented_product_relation";
  Attr = GetNamedAttribute(attrName, Ent);
  TmpEnt = *Attr->ptr.c;
  // tmp =G4GeometryTable::CreateObject(*TmpEnt);

  createdObject = psV;
}

void G4ContextDependentShapeRepresentationCreator::
  CreateSTEPGeometry(void* obj)
{
}
