// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ConicalSurfaceCreator.cc,v 1.2 2000-01-21 13:45:59 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ConicalSurfaceCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4ConicalSurfaceCreator.hh"
#include "G4GeometryTable.hh"
#include "G4ConicalSurface.hh"
#include "G4FConicalSurface.hh"

G4ConicalSurfaceCreator G4ConicalSurfaceCreator::csc;

G4ConicalSurfaceCreator::G4ConicalSurfaceCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4ConicalSurfaceCreator::~G4ConicalSurfaceCreator() {}

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
  srf->position_(place);
  // radius
  srf->radius_(fCon->GetSmallRadius());
  // semi-angle
  //srf->Semi_angle(fCon->GetAngle());
  // Set STEP info
  srf->SetFileId(GetNextId());
  srf->name_("");
  
  // Write out object 
  srf->STEPwrite(G4cout);

  createdObject = srf;
}
