#include "G4AdvancedFaceCreator.hh"
#include "G4CurveVector.hh"
#include "G4CompositeCurve.hh"
#include "G4Surface.hh"

G4AdvancedFaceCreator G4AdvancedFaceCreator::csc;

G4AdvancedFaceCreator::G4AdvancedFaceCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4AdvancedFaceCreator::~G4AdvancedFaceCreator(){}

void G4AdvancedFaceCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4Surface* srf=0;
  G4int sense;

  STEPentity* TmpEnt=0;  
  G4String attrName("bounds");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPaggregate *Aggr = Attr->ptr.a;
  SingleLinkNode *Node = Aggr->GetHead();

  G4int boundaryCount = Aggr->EntryCount();
  G4CurveVector* boundaryVec = new G4CurveVector;
  void *tmp;
 
  for(G4int a=0; a<boundaryCount;a++)
  {
    TmpEnt = ((EntityNode*)Node)->node;
    tmp =G4GeometryTable::CreateObject(*TmpEnt);
    G4CompositeCurve* pcc= new G4CompositeCurve;
    pcc->Init(*((G4CurveVector*)tmp));
    boundaryVec->insert(pcc);
    Node = Node->NextNode();
  }      
  
  attrName = "face_geometry";
  Attr = GetNamedAttribute(attrName, Ent);
  
  TmpEnt = *Attr->ptr.c;
  tmp =G4GeometryTable::CreateObject(*TmpEnt);
  srf = (G4Surface*)tmp;
  
  Attr = Ent.NextAttribute();
  // if(*Attr->ptr.b) // BINARY_TYPE
  //   sense = 1;
  // else
  //   sense = 0;

  // Utilize the INTEGER_TYPE attribute which is 0 or 1
  sense = *Attr->ptr.i;
      
  // L. Broglia
  srf->SetSameSense(sense);
  
  // Now, re-calculate the normal of the plane
  srf->CalcNormal();

  srf->SetBoundaries(boundaryVec);
  
  createdObject= srf;
}

void G4AdvancedFaceCreator::CreateSTEPGeometry(void* G4obj)
{

}


