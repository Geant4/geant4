//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4EzVolume.cc 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   G4EzVolume.cc
//
//                                         2005 Q
// ====================================================================
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4VisAttributes.hh"

#include "G4EzWorld.hh"
#include "G4EzVolume.hh"
#include "G4EzVoxelParameterization.hh"

// ====================================================================
//
// class description
//
// ====================================================================

////////////////////////
G4EzVolume::G4EzVolume()
  : name("MyVolume"),
    solid(0), 
    lv(0), lvsub(0),
    nplacement(0)
////////////////////////
{
}


/////////////////////////////////////////////
G4EzVolume::G4EzVolume(const G4String& aname)
  : name(aname), solid(0), lv(0), lvsub(0),
    nplacement(0)
/////////////////////////////////////////////
{
}


/////////////////////////
G4EzVolume::~G4EzVolume()
/////////////////////////
{
}


///////////////////////////////////////////////////////////////////////
void G4EzVolume::CreateBoxVolume(G4Material* amaterial, 
				 G4double dx, G4double dy, G4double dz)
///////////////////////////////////////////////////////////////////////
{
  if(lv !=0 ) {
    G4cout << "%%% Warning (G4EzVolume): volume is already created."
	   << G4endl;
    return;
  }

  solid= new G4Box(name, dx/2., dy/2., dz/2.);
  lv= new G4LogicalVolume(solid, amaterial, name);

  // vis. attributes
  va= new G4VisAttributes();
  lv-> SetVisAttributes(va);
}


////////////////////////////////////////////////////////////////////////////
void G4EzVolume::CreateTubeVolume(G4Material* amaterial,
				  G4double rmin, G4double rmax, G4double dz,
				  G4double phi0, G4double dphi)
////////////////////////////////////////////////////////////////////////////
{
  if(lv !=0 ) {
    G4cout << "%%% Warning (G4EzVolume): volume is already created."
	   << G4endl;
    return;
  }

  solid= new G4Tubs(name, rmin, rmax, dz, phi0, dphi);
  lv= new G4LogicalVolume(solid, amaterial, name);

  // vis. attributes
  va= new G4VisAttributes();
  lv-> SetVisAttributes(va);
}


/////////////////////////////////////////////////////////////////
void G4EzVolume::CreateConeVolume(G4Material* amaterial,
				  G4double rmin1, G4double rmax1,
				  G4double rmin2, G4double rmax2,
				  G4double dz,
				  G4double phi0, G4double dphi)
/////////////////////////////////////////////////////////////////
{
  if(lv !=0 ) {
    G4cout << "%%% Warning (G4EzVolume): volume is already created."
	   << G4endl;
    return;
  }

  solid= new G4Cons(name, rmin1, rmax1, rmin2, rmax2, 
		    dz, phi0, dphi);
  lv= new G4LogicalVolume(solid, amaterial, name);

  // vis. attributes
  va= new G4VisAttributes();
  lv-> SetVisAttributes(va);
}


/////////////////////////////////////////////////////////////////////
void G4EzVolume::CreateSphereVolume(G4Material* amaterial,
				    G4double rmin, G4double rmax,
				    G4double phi0, G4double dphi,
				    G4double theta0, G4double dtheta)
/////////////////////////////////////////////////////////////////////
{
  if(lv !=0 ) {
    G4cout << "%%% Warning (G4EzVolume): volume is already created."
	   << G4endl;
    return;
  }

  solid= new G4Sphere(name, rmin, rmax, phi0, dphi, theta0, dtheta);
  lv= new G4LogicalVolume(solid, amaterial, name);

  // vis. attributes
  va= new G4VisAttributes();
  lv-> SetVisAttributes(va);
}


//////////////////////////////////////////////////////////////////////
void G4EzVolume::CreateOrbVolume(G4Material* amaterial, G4double rmax)
//////////////////////////////////////////////////////////////////////
{
  if(lv !=0 ) {
    G4cout << "%%% Warning (G4EzVolume): volume is already created."
	   << G4endl;
    return;
  }

  solid= new G4Orb(name, rmax);
  lv= new G4LogicalVolume(solid, amaterial, name);

  // vis. attributes
  va= new G4VisAttributes();
  lv-> SetVisAttributes(va);
}


////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G4EzVolume::PlaceIt(const G4ThreeVector& pos, 
				       G4int ncopy, 
				       G4EzVolume* parent)
////////////////////////////////////////////////////////////////
{
  if(lv==0) {
    G4cout << "%%% Warning (G4EzVolume): volume is not yet created."
	   << G4endl;
    return 0;
  }
  
  G4PVPlacement* pv;
  if(parent==0) { // place it in the world
    G4VPhysicalVolume* world= G4EzWorld::GetWorldVolume();
    pv= new G4PVPlacement(0, pos, name, lv, world, false, ncopy);
  } else {
    pv= new G4PVPlacement(0, pos, lv, name, parent->lv, false, ncopy);
  }

  nplacement++;
  return pv;
}


//////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G4EzVolume::PlaceIt(const G4Transform3D& transform, 
				       G4int ncopy,
				       G4EzVolume* parent)
//////////////////////////////////////////////////////////////////////
{
  if(lv==0) {
    G4cout << "%%% Warning (G4EzVolume): volume is not yet created."
	   << G4endl;
    return 0;
  }

  G4PVPlacement* pv;
  if(parent==0) { // place it in the world
    G4VPhysicalVolume* world= G4EzWorld::GetWorldVolume();
    pv= new G4PVPlacement(transform, name, lv, world, false, ncopy);
  } else {
    pv= new G4PVPlacement(transform, lv, name, parent->lv, false, ncopy);
  }

  nplacement++;
  return pv;
}


///////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G4EzVolume::ReplicateIt(G4EzVolume* parent,
					   EAxis pAxis, G4int nReplicas,
					   G4double width, G4double offset)
///////////////////////////////////////////////////////////////////////////
{
  if(lv==0) {
    G4cout << "%%% Warning (G4EzVolume): volume is not yet created."
	   << G4endl;
    return 0;
  }

  G4PVReplica* pv= 
    new G4PVReplica(name, lv, parent->lv, pAxis, nReplicas, width, offset);

  nplacement += nReplicas;
  return pv;
}


//////////////////////////////////////////////////////////////////
G4ThreeVector G4EzVolume::VoxelizeIt(G4int nx, G4int ny, G4int nz)
//////////////////////////////////////////////////////////////////
{
  // creating voxel volume...
  G4Box* avolume= dynamic_cast<G4Box*>(solid);
  if(avolume ==0 ) {
    G4cout << "%%% Error (G4EzVolume): voxelization is valid " 
	   << "only for Box geometry." << G4endl;
    return G4ThreeVector();
  }

  if(lvsub !=0) {
    G4cout << "%%% Error (G4EzVolume): already voxelized." << G4endl;
    return G4ThreeVector();
  }

  G4double dx= (avolume-> GetXHalfLength())*2.;
  G4double dy= (avolume-> GetYHalfLength())*2.;
  G4double dz= (avolume-> GetZHalfLength())*2.;

  // voxel size
  G4double ddx= dx/nx;
  G4double ddy= dy/ny;
  G4double ddz= dz/nz;

  G4Box* voxel= new G4Box("voxel", ddx/2., ddy/2., ddz/2.);
  G4Material* voxelMaterial= lv-> GetMaterial();
  lvsub= new G4LogicalVolume(voxel, voxelMaterial, "voxel");

  G4VisAttributes* vavoxel= new G4VisAttributes(G4Color(1.,0.,0.));
  lvsub-> SetVisAttributes(vavoxel);
  
  G4EzVoxelParameterization* voxelParam=
    new G4EzVoxelParameterization(ddx, ddy, ddz, nx, ny, nz);
  G4int nvoxel= nx*ny*nz;
  new G4PVParameterised(name+"_voxel", lvsub, lv, kXAxis, 
			nvoxel, voxelParam);

  return G4ThreeVector(ddx, ddy, ddz);
}



////////////////////////////////////////////////////////////////
void G4EzVolume::SetSensitiveDetector(G4VSensitiveDetector* asd)
////////////////////////////////////////////////////////////////
{
  if(lvsub!=0) {
    lvsub-> SetSensitiveDetector(asd);
    return;
  }

  if(lv!=0) lv-> SetSensitiveDetector(asd);

}

