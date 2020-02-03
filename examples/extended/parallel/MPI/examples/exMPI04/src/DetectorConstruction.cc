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
//
/// @file DetectorConstruction.hh
/// @brief Define geometry

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "DetectorConstruction.hh"
#include "VoxelParam.hh"
#include "VoxelSD.hh"

typedef G4LogicalVolume G4LV;
typedef G4PVPlacement G4PVP;
typedef G4VisAttributes G4VA;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
  : flv_voxel(NULL)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nistManager = G4NistManager::Instance();
  G4VisAttributes* va;

  // world volume
  const G4double DXYZ_WORLD = 200.*cm;
  G4Box* sld_world = new G4Box("world",
                               DXYZ_WORLD/2., DXYZ_WORLD/2., DXYZ_WORLD/2.);

  G4Material* vacuum = nistManager-> FindOrBuildMaterial("G4_Galactic");
  G4LV* lv_world = new G4LogicalVolume(sld_world, vacuum, "world");
  G4PVP* world = new G4PVPlacement(0, G4ThreeVector(), "AREA",
                                   lv_world, 0, false, 0);
  // vis. attributes
  va = new G4VA(G4Color(1.,1.,1.));
  va-> SetVisibility(false);
  lv_world-> SetVisAttributes(va);

  // water phantom
  const G4double DXY_PHANTOM = 20.*cm;
  const G4double DZ_PHANTOM = 50.*cm;

  G4Box* sld_phantom = new G4Box("phantom",
                                 DXY_PHANTOM/2., DXY_PHANTOM/2., DZ_PHANTOM/2.);

  G4Material* water = nistManager-> FindOrBuildMaterial("G4_WATER");
  G4LV* lv_phantom = new G4LogicalVolume(sld_phantom,
                                         water, "phantom");

  va = new G4VA(G4Color(0.,1.,1.));
  lv_phantom-> SetVisAttributes(va);

  new G4PVP(0, G4ThreeVector(), lv_phantom, "phantom", lv_world, false, 0);

  // voxel planes
  const G4double DXY_VXP = 20.*cm;
  const G4double DZ_VXP = 1.*mm;

  G4Box* sld_vxp = new G4Box("vxplane", DXY_VXP/2., DXY_VXP/2., DZ_VXP/2.);
  G4LV* lv_vxp = new G4LV(sld_vxp, water, "vxplane");

  va = new G4VA(G4Color(1.,0.,0.));
  va-> SetVisibility(false);
  lv_vxp-> SetVisAttributes(va);

  for (G4int iz =0; iz < 500; iz++) {
    G4double z0 = -DZ_PHANTOM/2. + (iz+0.5)*DZ_VXP;
    new G4PVP(0, G4ThreeVector(0.,0.,z0),
              lv_vxp, "vxplane", lv_phantom, false, 1000+iz);
  }

  // voxel parameterized
  G4Box* sld_voxel = new G4Box("voxel",1.,1.,1.); // dummy
  flv_voxel = new G4LV(sld_voxel, water, "voxel");

  va = new G4VA(G4Color(0.,1.,1.));
  va-> SetVisibility(false);
  flv_voxel-> SetVisAttributes(va);

  const G4int nvoxels = 100*100;
  new G4PVParameterised("voxle", flv_voxel, lv_vxp, kUndefined, nvoxels,
                        new VoxelParam());

  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
  flv_voxel-> SetSensitiveDetector(new VoxelSD("voxel"));
}
