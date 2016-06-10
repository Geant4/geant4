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
// $Id: G4EzVolume.hh 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   G4EzVolume.hh
//
//   Utility class for supporting creating user geometry
//
//   * wrapper  of logical volume
//
//                                         2005 Q
// ====================================================================
#ifndef G4_EZ_VOLUME_H
#define G4_EZ_VOLUME_H

#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Transform3D.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
class G4Material;
class G4VSolid;
class G4VPhysicalVolume;

class G4EzVolume {
protected:
  // volume information
  G4String name;
  G4VSolid* solid;
  G4LogicalVolume* lv;
  G4LogicalVolume* lvsub; // logical volume for voxels
  G4VisAttributes* va;
  G4int nplacement;  // # of placements

public:
  G4EzVolume();
  G4EzVolume(const G4String& aname);
  ~G4EzVolume();

  // createing volume
  void CreateBoxVolume(G4Material* amaterial,
		       G4double dx, G4double dy, G4double dz);

  void CreateTubeVolume(G4Material* amaterial,
			G4double rmin, G4double rmax, 
			G4double dz,
			G4double phi0=0., G4double dphi=360*deg);

  void CreateConeVolume(G4Material* amaterial,
			G4double rmin1, G4double rmax1,
			G4double rmin2, G4double rmax2,
			G4double dz,
			G4double phi0=0., G4double dphi=360.*deg);
			
  void CreateSphereVolume(G4Material* amaterial,
			  G4double rmin, G4double rmax,
			  G4double phi0=0., G4double dphi=360.*deg,
			  G4double theta0=0., G4double dtheta=180.*deg);

  void CreateOrbVolume(G4Material* amaterial, G4double rmax);


  // placement
  G4VPhysicalVolume* PlaceIt(const G4ThreeVector& pos, G4int ncopy=0, 
			     G4EzVolume* parent=0);

  G4VPhysicalVolume* PlaceIt(const G4Transform3D& transform, G4int ncopy=0,
			     G4EzVolume* parent=0);
  // replica
  // parent volume should be exactly filled with replicated volumes.
  G4VPhysicalVolume* ReplicateIt(G4EzVolume* parent,
				 EAxis pAxis, G4int nReplicas,
				 G4double width, G4double offset=0);
  // voxelize
  // volume should be "BOX". otherwise, error.
  // returning voxel size in each dimension
  G4ThreeVector VoxelizeIt(G4int nx, G4int ny, G4int nz);
  
  // sensitivity
  // in the case of voxelized, SD will be set to a voxel.
  void SetSensitiveDetector(G4VSensitiveDetector* asd);

  // direct access to properties
  const G4String& GetName() const;

  void SetSolid(G4VSolid* asolid);
  const G4VSolid* GetSolid() const;

  void SetMaterial(G4Material* amaterial);
  G4Material* GetMaterial() const;

  G4int GetNofPlacements() const;
  
  void SetVisibility(G4bool qvisible);
  void SetColor(const G4Color& color);
  void SetColor(G4double red, G4double green, G4double blue);

};

// ====================================================================
//   inline functions
// ====================================================================

inline const G4String& G4EzVolume::GetName() const { return name; }

inline void G4EzVolume::SetSolid(G4VSolid* asolid) { solid= asolid; }

inline const G4VSolid* G4EzVolume::GetSolid() const { return solid; }

inline void G4EzVolume::SetMaterial(G4Material* amaterial) 
{
  if(lv!= 0) lv-> SetMaterial(amaterial);
  if(lvsub!= 0) lvsub-> SetMaterial(amaterial);
}

inline G4Material* G4EzVolume::GetMaterial() const 
{ 
  if(lv!=0) return lv-> GetMaterial();
  else return 0;
}

inline G4int G4EzVolume::GetNofPlacements() const { return nplacement; }

inline void G4EzVolume::SetVisibility(G4bool qvisible)
{
  if(va!=0) va-> SetVisibility(qvisible);
}

inline void G4EzVolume::SetColor(const G4Color& color)
{
  if(va!=0) va-> SetColor(color);
}

inline void G4EzVolume::SetColor(G4double red, G4double green, G4double blue)
{
  if(va!=0) va-> SetColor(red, green, blue);  
}

#endif
