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
// $Id: G4Axis2Placement3dCreator.cc,v 1.7 2002-11-21 16:49:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Axis2Placement3dCreator
//
// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------

#include "G4Axis2Placement3dCreator.hh"
#include "G4GeometryTable.hh"

G4Axis2Placement3dCreator G4Axis2Placement3dCreator::csc;

G4Axis2Placement3dCreator::G4Axis2Placement3dCreator()
{
  G4GeometryTable::RegisterObject(this);
}

G4Axis2Placement3dCreator::~G4Axis2Placement3dCreator() {}

G4Axis2Placement3dCreator G4Axis2Placement3dCreator::GetInstance()
{
  return csc;
}

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
  if (tmp)
    srfPoint = *(G4Point3D*)tmp;
  else
    G4cerr << "WARNING - G4Axis2Placement3dCreator::CreateG4Geometry" << G4endl
           << "\tPoint on surface (G4Point3D) NULL !" << G4endl
	   << "\tThis will result as an incorrect placement." << G4endl;
    
  // Get axis
  Attr = Ent.NextAttribute();
  TmpEnt = *Attr->ptr.c;
  tmp = G4GeometryTable::CreateObject(*TmpEnt);
  if (tmp)
    axis = *(G4Vector3D*)tmp;
  else
    G4cerr << "WARNING - G4Axis2Placement3dCreator::CreateG4Geometry" << G4endl
           << "\tPlacement axis (G4Vector3D) NULL !" << G4endl
	   << "\tThis will result as an incorrect placement." << G4endl;

  // Get direction
  Attr = Ent.NextAttribute();
  TmpEnt = *Attr->ptr.c;
  if (Attr->is_null())  // No direction specified!  Set it to (0,0,0) and warn.
  {
    tmp = 0;
    dir = G4Vector3D(0.,0.,0.);
  }
  else
    tmp = G4GeometryTable::CreateObject(*TmpEnt);
  if (tmp)
    dir = *(G4Vector3D*)tmp;  
  else
    G4cerr << "WARNING - G4Axis2Placement3dCreator::CreateG4Geometry" << G4endl
           << "\tPlacement direction (G4Vector3D) NULL !" << G4endl
	   << "\tThis will result as an incorrect placement." << G4endl;
  
  createdObject = new G4Axis2Placement3D(dir, axis, srfPoint);
}


void G4Axis2Placement3dCreator::CreateSTEPGeometry(void* G4obj)
{
  G4Axis2Placement3D *plc = (G4Axis2Placement3D*)G4obj;

  SdaiAxis2_placement_3d *place = new SdaiAxis2_placement_3d();  

  G4Vector3D axis = plc->GetAxis();
  G4Vector3D dir  = plc->GetRefDirection();
  G4Point3D pt    = plc->GetLocation();
  
  G4String pointName("Cartesian_Point");
  G4String directionName("Direction");  

  // Get location
  void * tmp =G4GeometryTable::CreateSTEPObject(&pt, pointName);
  place->location_((SdaiCartesian_point*)tmp);

  // Get axis
  tmp = G4GeometryTable::CreateSTEPObject(&axis, directionName);  
  place->axis_((SdaiDirection*)tmp);

  // Get dir
  tmp = G4GeometryTable::CreateSTEPObject(&dir, directionName);
  place->ref_direction_((SdaiDirection*)tmp);  

  // Set STEP info
  place->SetFileId(GetNextId());
  place->name_("");

  // Write out object & subobjects
  place->STEPwrite(G4cout);

  createdObject = place;
}
