// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EdgeCurveCreator.cc,v 1.2 2000-01-21 13:46:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4EdgeCurveCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4EdgeCurveCreator.hh"
#include "G4GeometryTable.hh"

G4EdgeCurveCreator G4EdgeCurveCreator::csc;

G4EdgeCurveCreator::G4EdgeCurveCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4EdgeCurveCreator::~G4EdgeCurveCreator() {}

void G4EdgeCurveCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4Point3D *pt1,*pt2;
  G4Curve* crv;
  G4bool sameSense;

  // Get start point
  G4String attrName("edge_start");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  pt1 = (G4Point3D*)tmp;

  // Get end point
  attrName = "edge_end";
  Attr = GetNamedAttribute(attrName, Ent);
  TmpEnt = *Attr->ptr.c;
  tmp =G4GeometryTable::CreateObject(*TmpEnt);
  pt2 = (G4Point3D*)tmp;
  
  // Get curve
  attrName = "edge_geometry";
  Attr = GetNamedAttribute(attrName, Ent);  
  TmpEnt = *Attr->ptr.c;
  tmp =G4GeometryTable::CreateObject(*TmpEnt);
  crv = (G4Curve*)tmp;
  
  // set the bounds
  G4Point3D p1(pt1->x(), pt1->y(), pt1->z());
  G4Point3D p2(pt2->x(), pt2->y(), pt1->z());

  crv->SetBounds(p1, p2);

  //get sense info
  Attr = Ent.NextAttribute();  
  //  sameSense = *Attr->ptr.b;

  createdObject = crv;
}


void G4EdgeCurveCreator::CreateSTEPGeometry(void* G4obj) {}
