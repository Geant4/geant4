// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Axis1PlacementCreator.cc,v 1.2 2000-01-21 13:45:57 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Axis1PlacementCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4Axis1PlacementCreator.hh"
#include "G4GeometryTable.hh"

G4Axis1PlacementCreator G4Axis1PlacementCreator::csc;

G4Axis1PlacementCreator::G4Axis1PlacementCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4Axis1PlacementCreator::~G4Axis1PlacementCreator() {}

void G4Axis1PlacementCreator::CreateG4Geometry(STEPentity& Ent)
{
}

void G4Axis1PlacementCreator::CreateSTEPGeometry(void* G4obj)
{
}
