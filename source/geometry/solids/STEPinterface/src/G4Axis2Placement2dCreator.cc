// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Axis2Placement2dCreator.cc,v 1.2 2000-01-21 13:45:57 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Axis2Placement2dCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4Axis2Placement2dCreator.hh"
#include "G4GeometryTable.hh"

G4Axis2Placement2dCreator G4Axis2Placement2dCreator::csc;

G4Axis2Placement2dCreator::G4Axis2Placement2dCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4Axis2Placement2dCreator::~G4Axis2Placement2dCreator() {}

void G4Axis2Placement2dCreator::CreateG4Geometry(STEPentity& Ent)
{
}

void G4Axis2Placement2dCreator::CreateSTEPGeometry(void* G4obj)
{
}
