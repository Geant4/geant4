#include "G4Axis2Placement3dCreator.hh"
G4Axis2Placement3dCreator G4Axis2Placement3dCreator::csc;

G4Axis2Placement3dCreator::G4Axis2Placement3dCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4Axis2Placement3dCreator::~G4Axis2Placement3dCreator(){}

void G4Axis2Placement3dCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4Axis2Placement3D place;
  G4Vector3D         dir;
  G4Vector3D         axis;  
  G4Point3D          srfPoint;

  Ent.ResetAttributes();
  STEPattribute* Attr = Ent.NextAttribute();
  if(Attr->NonRefType() != ENTITY_TYPE) Attr = Ent.NextAttribute();

  // Get point on surface
  STEPentity* TmpEnt = *Attr->ptr.c;
  void * tmp = G4GeometryTable::CreateObject(*TmpEnt);
  srfPoint = *(G4Point3D*)tmp;
    
  // Get axis
  Attr = Ent.NextAttribute();
  TmpEnt = *Attr->ptr.c;
  tmp =G4GeometryTable::CreateObject(*TmpEnt);
  axis = *(G4Vector3D*)tmp;

  // Get direction
  Attr = Ent.NextAttribute();
  TmpEnt = *Attr->ptr.c;
  tmp =G4GeometryTable::CreateObject(*TmpEnt);
  dir = *(G4Vector3D*)tmp;  
  
  createdObject = new G4Axis2Placement3D(dir, axis, srfPoint);
}


void G4Axis2Placement3dCreator::CreateSTEPGeometry(void* G4obj)
{
  G4Axis2Placement3D *plc = (G4Axis2Placement3D*)G4obj;

  SdaiAxis2_placement_3d *place = new SdaiAxis2_placement_3d();  

  G4Vector3D axis = plc->GetAxis();
  G4Vector3D dir  = plc->GetRefDirection();
  G4Vector3D pt   = plc->GetLocation();
  
  G4String pointName("Cartesian_Point");
  G4String directionName("Direction");  

  // Get location
  void * tmp =G4GeometryTable::CreateSTEPObject(&pt, pointName);
  place->Location((SdaiCartesian_point*)tmp);

  // Get axis
  tmp = G4GeometryTable::CreateSTEPObject(&axis, directionName);  
  place->Axis((SdaiDirection*)tmp);

  // Get dir
  tmp = G4GeometryTable::CreateSTEPObject(&dir, directionName);
  place->Ref_direction((SdaiDirection*)tmp);  

  // Set STEP info
  place->SetFileId(GetNextId());
  place->Name("");

  // Write out object & subobjects
  place->STEPwrite(G4cout);

  createdObject = place;
}





