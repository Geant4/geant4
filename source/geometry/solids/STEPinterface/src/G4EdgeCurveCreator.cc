//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4EdgeCurveCreator.cc,v 1.6 2002-11-21 16:49:48 gcosmo Exp $
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

G4EdgeCurveCreator G4EdgeCurveCreator::GetInstance()
{
  return csc;
}

void G4EdgeCurveCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4Point3D *pt1,*pt2;
  G4Curve* crv=0;

  // Get start point
  G4String attrName("edge_start");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPentity* TmpEnt = *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  pt1 = (G4Point3D*)tmp;
  if (!tmp)
    G4cerr << "WARNING - G4EdgeCurveCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL edge-start (G4Point3D) !" << G4endl;

  // Get end point
  attrName = "edge_end";
  Attr = GetNamedAttribute(attrName, Ent);
  TmpEnt = *Attr->ptr.c;
  tmp =G4GeometryTable::CreateObject(*TmpEnt);
  pt2 = (G4Point3D*)tmp;
  if (!tmp)
    G4cerr << "WARNING - G4EdgeCurveCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL edge-end (G4Point3D) !" << G4endl;
  
  // Get curve
  attrName = "edge_geometry";
  Attr = GetNamedAttribute(attrName, Ent);  
  TmpEnt = *Attr->ptr.c;
  tmp =G4GeometryTable::CreateObject(*TmpEnt);
  
  if ((!pt1) || (!pt2))
  {
    G4cerr << "\tEdge Curve NOT created." << G4endl;
    createdObject = 0;
    return;
  }

  if (tmp)
  {
    crv = (G4Curve*)tmp;
  
    // set the bounds
    G4Point3D p1(pt1->x(), pt1->y(), pt1->z());
    G4Point3D p2(pt2->x(), pt2->y(), pt1->z());

    crv->SetBounds(p1, p2);
  }
  else
    G4cerr << "WARNING - G4EdgeCurveCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL curve (G4Curve) !" << G4endl
	   << "\tEdge Curve NOT created." << G4endl;
  
  //get sense info
  Attr = Ent.NextAttribute();  
  //  sameSense = *Attr->ptr.b;

  createdObject = crv;
}


void G4EdgeCurveCreator::CreateSTEPGeometry(void* G4obj) {}
