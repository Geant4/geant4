// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FaceOuterBoundCreator.cc,v 1.2 2000-01-21 13:46:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4FaceOuterBoundCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4FaceOuterBoundCreator.hh"
#include "G4GeometryTable.hh"

G4FaceOuterBoundCreator G4FaceOuterBoundCreator::csc;

G4FaceOuterBoundCreator::G4FaceOuterBoundCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4FaceOuterBoundCreator::~G4FaceOuterBoundCreator() {}

void G4FaceOuterBoundCreator::CreateSTEPGeometry(void* G4obj)
{
}
