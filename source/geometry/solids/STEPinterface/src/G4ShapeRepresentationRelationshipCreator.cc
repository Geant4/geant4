// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ShapeRepresentationRelationshipCreator.cc,v 1.2 2000-01-21 13:46:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4ShapeDefinitionRepresentationCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4ShapeRepresentationRelationshipCreator.hh"
#include "G4GeometryTable.hh"

G4ShapeRepresentationRelationshipCreator
  G4ShapeRepresentationRelationshipCreator::csc;

G4ShapeRepresentationRelationshipCreator::
  G4ShapeRepresentationRelationshipCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4ShapeRepresentationRelationshipCreator::
  ~G4ShapeRepresentationRelationshipCreator() {}

void G4ShapeRepresentationRelationshipCreator::CreateG4Geometry(STEPentity& Ent)
{
}

void G4ShapeRepresentationRelationshipCreator::CreateSTEPGeometry(void* G4obj)
{
}
