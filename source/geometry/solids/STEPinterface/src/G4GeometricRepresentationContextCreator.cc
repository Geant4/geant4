// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GeometricRepresentationContextCreator.cc,v 1.4 2000-11-20 18:17:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4GeometricRepresentationContextCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4GeometricRepresentationContextCreator.hh"
#include "G4GeometryTable.hh"

G4GeometricRepresentationContextCreator
  G4GeometricRepresentationContextCreator::csc;

G4GeometricRepresentationContextCreator::
  G4GeometricRepresentationContextCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4GeometricRepresentationContextCreator::
  ~G4GeometricRepresentationContextCreator() {}

void G4GeometricRepresentationContextCreator::CreateG4Geometry(STEPentity& Ent)
{
  // Made by L. Broglia
  STEPattribute *Attr;

  G4String attrName("coordinate_space_dimension");
  Attr = GetNamedAttribute(attrName, Ent);

  // Get coordinate space dimension (ptr.r --> REAL_TYPE)
  // SdaiReal Tmpdim = *Attr->ptr.r;

  // void *tmp = G4GeometryTable::CreateObject(*TmpEnt);
  // place = (G4Axis2Placement3D*)tmp...;

}

void G4GeometricRepresentationContextCreator::CreateSTEPGeometry(void* G4obj)
{
}
