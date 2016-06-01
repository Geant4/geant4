#include "G4ConicalSurfaceCreator.hh"
#include "G4ConicalSurface.hh"

G4ConicalSurfaceCreator G4ConicalSurfaceCreator::csc;

G4ConicalSurfaceCreator::G4ConicalSurfaceCreator(){G4GeometryTable::RegisterObject(this);}

G4ConicalSurfaceCreator::~G4ConicalSurfaceCreator(){}

void G4ConicalSurfaceCreator::CreateG4Geometry(STEPentity& Ent)
{
  // Made by L. Broglia

  STEPattribute *Attr;
  G4Axis2Placement3D* place; 
  G4double angle;

  G4String attrName("position");
  Attr = GetNamedAttribute(attrName, Ent);

  // Get placement
  STEPentity* TmpEnt = *Attr->ptr.c; // ptr.c --> ENTITY_TYPE
  void *tmp = G4GeometryTable::CreateObject(*TmpEnt);
  place = (G4Axis2Placement3D*)tmp;
  
  // Get angle
  G4String attrNameA("semi_angle");
  Attr = GetNamedAttribute(attrNameA, Ent);
  angle = *Attr->ptr.r; // Be careful, angle is in degree
   

  G4ConicalSurface* aG4cone = 
    new G4ConicalSurface( (*place).GetLocation() ,
			  (*place).GetAxis()     ,
			  angle * 2*M_PI/360      );

  createdObject = aG4cone;

 
}

void G4ConicalSurfaceCreator::CreateSTEPGeometry(void* G4obj)
{
  G4FConicalSurface* fCon = (G4FConicalSurface*)G4obj;
  SdaiConical_surface* srf= new SdaiConical_surface();
  // Get placement
  G4String placementName("Axis2Placement3d");
  void * tmp =G4GeometryTable::CreateSTEPObject(&fCon, placementName);
  SdaiAxis2_placement_3d *place = (SdaiAxis2_placement_3d*)tmp;
  srf->Position(place);
  // radius
    srf->Radius(fCon->GetSmallRadius());
  // semi-angle
  //srf->Semi_angle(fCon->GetAngle());
  // Set STEP info
  srf->SetFileId(GetNextId());
  srf->Name("");
  
  // Write out object 
  srf->STEPwrite(G4cout);

  createdObject = srf;
}



