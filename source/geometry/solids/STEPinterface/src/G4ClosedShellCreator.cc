#include "G4ClosedShellCreator.hh"

G4ClosedShellCreator G4ClosedShellCreator::csc;

G4ClosedShellCreator::G4ClosedShellCreator(){G4GeometryTable::RegisterObject(this);}

G4ClosedShellCreator::~G4ClosedShellCreator(){}

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
    void *tmp =::G4GeometryTable::CreateObject(*TmpEnt);
    SurfaceVec.append((G4Surface*)tmp);     
    Node = Node->NextNode();
  }      

  // create G4solid
  G4Surface** srfVec =  new G4Surface*[FaceCount];
  
  for(G4int b=0;b<FaceCount;b++)
    srfVec[b] = (SurfaceVec[b]);    
    
  G4BREPSolid* sld = new G4BREPSolid(" ", srfVec, FaceCount);

  createdObject = sld;

  // delete [] srfVec;
  // delete sld;
}

void G4ClosedShellCreator::CreateSTEPGeometry(void* G4obj)
{
  G4BREPSolid* bSld = (G4BREPSolid*)G4obj;
  SdaiClosed_shell* sld = new SdaiClosed_shell();

  G4int surfaceCount  = bSld->NumberOfFaces();

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
  


  sld->Cfs_faces(&eAggr);

  
  sld->SetFileId(GetNextId());
  sld->Name("");
  sld->STEPwrite(G4cout);

  createdObject = sld;
}
