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
// $Id: G4CartesianPointCreator.cc,v 1.4 2002-11-21 16:49:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4CartesianPointCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4CartesianPointCreator.hh"
#include "G4GeometryTable.hh"

G4CartesianPointCreator G4CartesianPointCreator::csc;

G4CartesianPointCreator::G4CartesianPointCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4CartesianPointCreator::~G4CartesianPointCreator() {}

G4CartesianPointCreator G4CartesianPointCreator::GetInstance()
{
  return csc;
}

void G4CartesianPointCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4double x,y,z;
  Ent.ResetAttributes();
  
  STEPattribute* Attr = Ent.NextAttribute();
  Attr = Ent.NextAttribute();
  STEPaggregate& XYZAggr = *Attr->ptr.a;
  SingleLinkNode* RNode = XYZAggr.GetHead();
  
  x = ((RealNode*)RNode)->value;
  RNode = RNode->NextNode();
  y = ((RealNode*)RNode)->value;	
  RNode = RNode->NextNode();
  z = ((RealNode*)RNode)->value;
  // createdObject = new G4Point3d(x,y,z);
  createdObject = new G4Point3D(x,y,z);
}

void G4CartesianPointCreator::CreateSTEPGeometry(void* G4obj)
{
  //G4Point3d* pt =  (G4Point3d*)G4obj;
  G4Point3D* pt =  (G4Point3D*)G4obj;
  RealNode *xa,*ya,*za;

  RealAggregate pointAggr;  
  xa = (RealNode*)pointAggr.NewNode();
  // xa->value = pt->X();
  xa->value = pt->x();

  ya = (RealNode*)pointAggr.NewNode();
  // ya->value = pt->Y();
  ya->value = pt->y();
  
  za = (RealNode*)pointAggr.NewNode();  
  // za->value = pt->Z();
  za->value = pt->z();

  pointAggr.AddNode(xa);
  pointAggr.AddNode(ya);
  pointAggr.AddNode(za);

  SdaiCartesian_point *sPoint = new SdaiCartesian_point();
  sPoint->coordinates_(&pointAggr);
  sPoint->SetFileId(GetNextId());
  sPoint->name_("");

  sPoint->STEPwrite(G4cout);
  createdObject = sPoint;
}
