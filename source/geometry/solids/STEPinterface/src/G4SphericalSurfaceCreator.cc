#include "G4SphericalSurfaceCreator.hh"

G4SphericalSurfaceCreator G4SphericalSurfaceCreator::csc;

G4SphericalSurfaceCreator::G4SphericalSurfaceCreator(){G4GeometryTable::RegisterObject(this);}

G4SphericalSurfaceCreator::~G4SphericalSurfaceCreator(){}

void G4SphericalSurfaceCreator::CreateG4Geometry(STEPentity& Ent)
{
  // Made by L. Broglia

  STEPattribute *Attr;
  G4Axis2Placement3D* place; 
  G4double radius;
  
  G4String attrName("position");
  Attr = GetNamedAttribute(attrName, Ent);

  // Get placement
  STEPentity* TmpEnt = *Attr->ptr.c; // ptr.c --> ENTITY_TYPE
  void *tmp = G4GeometryTable::CreateObject(*TmpEnt);
  place = (G4Axis2Placement3D*)tmp;
  
  // Get radius
  G4String attrNameR("radius");
  Attr = GetNamedAttribute(attrNameR, Ent);
  radius = *Attr->ptr.r;

  
  G4SphericalSurface* aG4sphere = 
    new G4SphericalSurface( (*place).GetLocation() ,
			    (*place).GetPX()       ,
			    (*place).GetPZ()       , 
			    radius                 ,
			    0, 2*M_PI              ,
			    0, M_PI                 );

  createdObject = aG4sphere;

}

void G4SphericalSurfaceCreator::CreateSTEPGeometry(void* G4obj)
{
  G4SphericalSurface* fSph = (G4SphericalSurface*)G4obj;
  SdaiSpherical_surface* srf= new SdaiSpherical_surface();
  // Get placement
  G4String placementName("Axis2Placement3d");
  void * tmp =G4GeometryTable::CreateSTEPObject(&fSph, placementName);
  SdaiAxis2_placement_3d *place = (SdaiAxis2_placement_3d*)tmp;
  srf->Position(place);
  // radius
    srf->Radius(fSph->GetRadius());

  // Set STEP info
  srf->SetFileId(GetNextId());
  srf->Name("");
  
  // Write out object 
  srf->STEPwrite(G4cout);

  createdObject = srf;
}



