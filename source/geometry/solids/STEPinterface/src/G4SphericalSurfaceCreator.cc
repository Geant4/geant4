//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4SphericalSurfaceCreator.cc,v 1.4 2001-07-11 10:00:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class G4SphericalSurfaceCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4SphericalSurfaceCreator.hh"
#include "G4GeometryTable.hh"
#include "G4SphericalSurface.hh"

G4SphericalSurfaceCreator G4SphericalSurfaceCreator::csc;

G4SphericalSurfaceCreator::G4SphericalSurfaceCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4SphericalSurfaceCreator::~G4SphericalSurfaceCreator() {}

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

  
  G4SphericalSurface* aG4sphere = 0;
  if (place)
    aG4sphere = new G4SphericalSurface( (*place).GetLocation() ,
					(*place).GetPX()       ,
					(*place).GetPZ()       , 
					radius                 ,
					0, 2*M_PI              ,
					0, M_PI                 );
  else
    G4cerr << "WARNING - G4SphericalSurfaceCreator::CreateG4Geometry" << G4endl
           << "\tUnexpected NULL axis placement (G4Axis2Placement3D) !" << G4endl
	   << "\tSpherical Surface NOT created." << G4endl;

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
  srf->position_(place);
  // radius
  srf->radius_(fSph->GetRadius());

  // Set STEP info
  srf->SetFileId(GetNextId());
  srf->name_("");
  
  // Write out object 
  srf->STEPwrite(G4cout);

  createdObject = srf;
}
