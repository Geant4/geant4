// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PointCreator.cc,v 1.2 2000-01-21 13:46:04 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4PointCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4PointCreator.hh"
#include "G4GeometryTable.hh"

G4PointCreator G4PointCreator::csc;

G4PointCreator::G4PointCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4PointCreator::~G4PointCreator() {}

void G4PointCreator::CreateG4Geometry(STEPentity& Ent) {}

void G4PointCreator::CreateSTEPGeometry(void* G4obj) {}
