// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Axis2PlacementsCreator.cc,v 1.3 2000-01-21 13:45:58 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Axis2PlacementsCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4Axis2PlacementsCreator.hh"
#include "G4GeometryTable.hh"

G4Axis2PlacementsCreator G4Axis2PlacementsCreator::csc;

G4Axis2PlacementsCreator::G4Axis2PlacementsCreator()
{
  // G4cerr << " G4Axis2PlacementsCreator default constructor " 
  // << " calling G4GeometryTable::RegisterObject(this) " << G4endl;
  G4GeometryTable::RegisterObject(this);
}

G4Axis2PlacementsCreator::~G4Axis2PlacementsCreator() {}

void G4Axis2PlacementsCreator::CreateG4Geometry(STEPentity& Ent)
{
}

void G4Axis2PlacementsCreator::CreateSTEPGeometry(void* G4obj)
{
}
