// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DirectionCreator.cc,v 1.2 2000-01-21 13:46:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4DirectionCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4DirectionCreator.hh"
#include "G4GeometryTable.hh"

G4DirectionCreator G4DirectionCreator::csc;

G4DirectionCreator::G4DirectionCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4DirectionCreator::~G4DirectionCreator() {}

void G4DirectionCreator::CreateG4Geometry(STEPentity& Ent)
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
  createdObject = new G4ThreeVector(x,y,z);
}


void G4DirectionCreator::CreateSTEPGeometry(void* G4obj)
{
  //G4ThreeVec* dir =  (G4ThreeVec*)G4obj;
  G4Vector3D* dir =  (G4Vector3D*)G4obj;
  RealNode *xa,*ya,*za;

  RealAggregate dirAggr;  
  xa = (RealNode*)dirAggr.NewNode();
  // xa->value = dir->X();
  xa->value = dir->x();

  ya = (RealNode*)dirAggr.NewNode();
  // ya->value = dir->Y();
  ya->value = dir->y();

  za = (RealNode*)dirAggr.NewNode();  
  // za->value = dir->Z();
  za->value = dir->z();

  dirAggr.AddNode(xa);
  dirAggr.AddNode(ya);
  dirAggr.AddNode(za);

  SdaiDirection *sDir = new SdaiDirection();
  sDir->direction_ratios_(&dirAggr);
  sDir->SetFileId(GetNextId());
  sDir->name_("");
  sDir->STEPwrite(G4cout);
  createdObject = sDir;
}
