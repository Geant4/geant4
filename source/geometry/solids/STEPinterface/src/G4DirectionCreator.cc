#include "G4DirectionCreator.hh"

G4DirectionCreator G4DirectionCreator::csc;

G4DirectionCreator::G4DirectionCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4DirectionCreator::~G4DirectionCreator(){}

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
  sDir->Direction_ratios(&dirAggr);
  sDir->SetFileId(GetNextId());
  sDir->Name("");
  sDir->STEPwrite(G4cout);
  createdObject = sDir;
}



