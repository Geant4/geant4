#include "G4ToroidalSurfaceCreator.hh"

G4ToroidalSurfaceCreator G4ToroidalSurfaceCreator::csc;

G4ToroidalSurfaceCreator::G4ToroidalSurfaceCreator(){G4GeometryTable::RegisterObject(this);}

G4ToroidalSurfaceCreator::~G4ToroidalSurfaceCreator(){}

void G4ToroidalSurfaceCreator::CreateG4Geometry(STEPentity& Ent)
{
  G4Axis2Placement3D* place;
  G4double majorRadius, minorRadius;

  // Get the placement
  G4String attrName("position");
  STEPattribute *Attr = GetNamedAttribute(attrName, Ent);
  STEPentity* TmpEnt= *Attr->ptr.c;
  void *tmp =G4GeometryTable::CreateObject(*TmpEnt);
  place = (G4Axis2Placement3D*)tmp;
  
  // get radii
  attrName = "major_radius";
  Attr = GetNamedAttribute(attrName, Ent);
  majorRadius = *Attr->ptr.r;
  attrName = "minor_radius";
  Attr = GetNamedAttribute(attrName, Ent);  
  minorRadius = *Attr->ptr.r;

  
  G4ToroidalSurface* aTorus = new G4ToroidalSurface( (*place).GetLocation() ,
						     (*place).GetPZ()       ,
						     (*place).GetPX()       ,
						     minorRadius          ,
						     majorRadius            );
  createdObject = aTorus;

}


void G4ToroidalSurfaceCreator::CreateSTEPGeometry(void* G4obj)
{
  G4ToroidalSurface* tor = (G4ToroidalSurface*)G4obj;
  SdaiToroidal_surface* srf= new SdaiToroidal_surface();
  // Get placement
  G4String placementName("Axis2Placement3d");
  void * tmp =G4GeometryTable::CreateSTEPObject(&tor, placementName);
  SdaiAxis2_placement_3d *place = (SdaiAxis2_placement_3d*)tmp;
  srf->Position(place);
  // radiis
    srf->Minor_radius(tor->GetMinRadius());
  srf->Major_radius(tor->GetMaxRadius());  

  // Set STEP info
  srf->SetFileId(GetNextId());
  srf->Name("");
  
  // Write out object 
  srf->STEPwrite(G4cout);

  createdObject = srf;
}




