#include "G4CircleCreator.hh"
#include "G4CircularCurve.hh"


G4CircleCreator G4CircleCreator::csc;

G4CircleCreator::G4CircleCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4CircleCreator::~G4CircleCreator(){}

void G4CircleCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4Axis2Placement3D place;

  G4double radius;
  
  Ent.ResetAttributes();
  STEPattribute* Attr = Ent.NextAttribute();
  while(Attr->NonRefType() == STRING_TYPE ||
	Attr->NonRefType() == sdaiSTRING )
    Attr = Ent.NextAttribute();
  
  // Get the placement
  SelectNode* SelectN = (SelectNode*)Attr->ptr.sh;
  SdaiSelect* Select = (SdaiSelect*)SelectN;
  SdaiAxis2_placement* Place = (SdaiAxis2_placement*)Select;
  STEPentity* TmpEnt = Place->operator SdaiAxis2_placement_3dH();
  void * tmp =G4GeometryTable::CreateObject(*TmpEnt);
 
  place = *(G4Axis2Placement3D*)tmp;
  
  Attr = Ent.NextAttribute();	
  radius = *Attr->ptr.r;

  G4CircularCurve* circle = new G4CircularCurve();
  circle->Init(place , radius);
  createdObject = circle;
}

void G4CircleCreator::CreateSTEPGeometry(void* G4obj)
{

}

