// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VertexPointCreator.cc,v 1.2 2000-01-21 13:46:07 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4VertexPointCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4VertexPointCreator.hh"
#include "G4GeometryTable.hh"

G4VertexPointCreator G4VertexPointCreator::csc;

G4VertexPointCreator::G4VertexPointCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4VertexPointCreator::~G4VertexPointCreator() {}

void G4VertexPointCreator::CreateG4Geometry(STEPentity& Ent)
{
  // G4Point3d *point;
  G4Point3D *point;

  G4String attrName("vertex_geometry");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);

  // Get the geometric point  
  STEPentity* TmpEnt= *Attr->ptr.c;
  void * tmp =G4GeometryTable::CreateObject(*TmpEnt);
  
  // point = (G4Point3d*)tmp;
  point = (G4Point3D*)tmp;
  createdObject = point;
}

void G4VertexPointCreator::CreateSTEPGeometry(void* G4obj) {}
