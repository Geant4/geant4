#include "G4CartesianPointCreator.hh"

G4CartesianPointCreator G4CartesianPointCreator::csc;

G4CartesianPointCreator::G4CartesianPointCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4CartesianPointCreator::~G4CartesianPointCreator(){}

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
  sPoint->Coordinates(&pointAggr);
  sPoint->SetFileId(GetNextId());
  sPoint->Name("");

  sPoint->STEPwrite(G4cout);
  createdObject = sPoint;
  
}





