// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ClosedShellCreator.cc,v 1.5 2000-11-20 18:17:29 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ClosedShellCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4ClosedShellCreator.hh"
#include "G4GeometryTable.hh"
#include "G4BREPSolid.hh"

G4ClosedShellCreator G4ClosedShellCreator::csc;

G4ClosedShellCreator::G4ClosedShellCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4ClosedShellCreator::~G4ClosedShellCreator() {}

void G4ClosedShellCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4String attrName("cfs_faces");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPaggregate *Aggr = Attr->ptr.a;
  SingleLinkNode *Node = Aggr->GetHead();
  STEPentity* TmpEnt=0;
  
  G4int FaceCount = Aggr->EntryCount();
 
  G4SurfaceVector SurfaceVec;  

  for(G4int a=0; a<FaceCount;a++)
  {
    TmpEnt = ((EntityNode*)Node)->node;
    void *tmp = G4GeometryTable::CreateObject(*TmpEnt);
    if (tmp) SurfaceVec.append((G4Surface*)tmp);
    Node = Node->NextNode();
  }      

  // create G4solid
  G4int SurfNum = SurfaceVec.entries();
  G4Surface** srfVec =  new G4Surface*[SurfNum];
  if (SurfNum != FaceCount)
    G4cerr << "WARNING - G4ClosedShellCreator::CreateG4Geometry" << G4endl
           << "\tTotal of " << SurfNum << " G4Surface components created, out of "
	   << FaceCount << " expected !" << G4endl
	   << "\tBREP Solid not correctly instantiated." << G4endl;
  for(G4int b=0;b<SurfNum;b++)
    srfVec[b] = (SurfaceVec[b]);    
    
  G4BREPSolid* sld = new G4BREPSolid(" ", srfVec, SurfNum);

  createdObject = sld;
}

void G4ClosedShellCreator::CreateSTEPGeometry(void* G4obj)
{
  G4BREPSolid* bSld = (G4BREPSolid*)G4obj;
  SdaiClosed_shell* sld = new SdaiClosed_shell();

  G4int surfaceCount  = bSld->GetNumberOfFaces();

  EntityAggregate eAggr;
  EntityNode* eNode=0;
  void* tmp=0;
  G4Surface* srf=0;
  
  for(G4int a=0;a<surfaceCount;a++)
    {
      srf = bSld->GetSurface(a);
      G4String srfType(srf->GetEntityType());
      tmp = G4GeometryTable::CreateSTEPObject(srf,srfType);
      
      eNode = new EntityNode((STEPentity*)tmp); 
      eAggr.AddNode(eNode);
    }

  sld->cfs_faces_(&eAggr);  
  sld->SetFileId(GetNextId());
  sld->name_("");
  sld->STEPwrite(G4cout);

  createdObject = sld;
}
