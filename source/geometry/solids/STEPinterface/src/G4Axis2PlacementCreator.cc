// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Axis2PlacementCreator.cc,v 1.2 2000-01-21 13:45:58 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Axis2PlacementCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4Axis2PlacementCreator.hh"
#include "G4GeometryTable.hh"

G4Axis2PlacementCreator G4Axis2PlacementCreator::csc;

G4Axis2PlacementCreator::G4Axis2PlacementCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4Axis2PlacementCreator::~G4Axis2PlacementCreator() {}

void G4Axis2PlacementCreator::CreateG4Geometry(STEPentity& Ent)
{
}

void G4Axis2PlacementCreator::CreateSTEPGeometry(void* G4obj)
{
}
