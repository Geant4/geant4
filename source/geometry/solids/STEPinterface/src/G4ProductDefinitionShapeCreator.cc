// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProductDefinitionShapeCreator.cc,v 1.3 2000-01-21 13:46:05 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4ProductDefinitionShapeCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4ProductDefinitionShapeCreator.hh"
#include "G4GeometryTable.hh"

G4ProductDefinitionShapeCreator G4ProductDefinitionShapeCreator::csc;

G4ProductDefinitionShapeCreator::G4ProductDefinitionShapeCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4ProductDefinitionShapeCreator::~G4ProductDefinitionShapeCreator() {}

void G4ProductDefinitionShapeCreator::CreateG4Geometry(STEPentity& Ent)
{
  // Made by L. Broglia
  STEPattribute *Attr;

  G4String attrName("description");
  Attr = GetNamedAttribute(attrName, Ent);

  // Get description
  SdaiString Tmpstring = *Attr->ptr.S; // ptr.S --> STRING_TYPE

  attrName = "definition";
  Attr = GetNamedAttribute(attrName, Ent);

  // Get definition
  SCLundefined Tmpdef = *Attr->ptr.u; // ptr.u --> UNKNOWN_TYPE

  // void *tmp = G4GeometryTable::CreateObject(*TmpEnt);
  // place = (G4Axis2Placement3D*)tmp;
}

void G4ProductDefinitionShapeCreator::CreateSTEPGeometry(void* G4obj)
{
}
