// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenShellCreator.cc,v 1.2 2000-01-21 13:46:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4OpenShellCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4OpenShellCreator.hh"
#include "G4GeometryTable.hh"

G4OpenShellCreator G4OpenShellCreator::csc;

G4OpenShellCreator::G4OpenShellCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4OpenShellCreator::~G4OpenShellCreator() {}

void G4OpenShellCreator::CreateG4Geometry(STEPentity& Ent)
{
}

void G4OpenShellCreator::CreateSTEPGeometry(void* G4obj)
{
}
