// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EdgeLoopCreator.cc,v 1.2 2000-01-21 13:46:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4EdgeLoopCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4EdgeLoopCreator.hh"
#include "G4GeometryTable.hh"

G4EdgeLoopCreator G4EdgeLoopCreator::csc;

G4EdgeLoopCreator::G4EdgeLoopCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4EdgeLoopCreator::~G4EdgeLoopCreator() {}

void G4EdgeLoopCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4String attrName("edge_list");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  
  STEPaggregate *Aggr = Attr->ptr.a;
  SingleLinkNode *Node = Aggr->GetHead();
  STEPentity* TmpEnt=0;

  G4int EdgeCount = Aggr->EntryCount();
  G4CurveVector* CurveVec = new G4CurveVector[EdgeCount];

  for(G4int a=0; a<EdgeCount;a++)
    {
      TmpEnt = ((EntityNode*)Node)->node;
      void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
      CurveVec->append((G4Curve*)tmp);
      Node = Node->NextNode();
    }      

  createdObject = CurveVec;
}

void G4EdgeLoopCreator::CreateSTEPGeometry(void* G4obj)
{
}
