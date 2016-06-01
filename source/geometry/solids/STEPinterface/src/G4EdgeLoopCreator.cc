#include "G4EdgeLoopCreator.hh"

G4EdgeLoopCreator G4EdgeLoopCreator::csc;

G4EdgeLoopCreator::G4EdgeLoopCreator(){G4GeometryTable::RegisterObject(this);}

G4EdgeLoopCreator::~G4EdgeLoopCreator(){}

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
