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
/// @file DetectorConstruction.cc
/// @brief Define geometry

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "DetectorConstruction.hh"
#include "Materials.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // material definition
  Materials* materialConstruction = new Materials;
  materialConstruction-> Construct();

  G4Material* mate;
  G4VisAttributes* va;

  // world volume
  const G4double DXYZ_AREA = 32.*cm;
  G4Box* areaSolid = new G4Box("AREA", DXYZ_AREA/2., DXYZ_AREA/2.,
                                       DXYZ_AREA/2.);

  G4Material* vacuum = G4Material::GetMaterial("Vacuum");
  G4LogicalVolume* areaLV = new G4LogicalVolume(areaSolid, vacuum, "AREA_LV");
  G4PVPlacement* area = new G4PVPlacement(0, G4ThreeVector(), "AREA_PV",
                                          areaLV, 0, false, 0);
  // vis. attributes
  va = new G4VisAttributes(G4Color(1.,1.,1.));
  va-> SetVisibility(false);
  areaLV-> SetVisAttributes(va);

  // detectors
  // voxel
  const G4double dvoxel = 10.*mm;
  const G4double dl = 10.*cm;

  G4Box* svoxel = new G4Box("voxel", dvoxel, dl, dvoxel);
  mate = G4Material::GetMaterial("Vacuum");
  G4LogicalVolume* lvoxel = new G4LogicalVolume(svoxel, mate, "voxel");
  va = new G4VisAttributes(G4Color(0.,0.8,0.8));
  va-> SetVisibility(false);
  lvoxel-> SetVisAttributes(va);

  G4int ix, iz;
  G4int index = 0;
  for ( iz = 0; iz < 5; iz++ ) {
    for ( ix = -7; ix <= 7; ix++ ) {
      G4double x0 = (2.*ix) * cm;
      G4double z0 = (-13.+2.*iz) * cm;
      new G4PVPlacement(0, G4ThreeVector(x0, 0., z0),
                        lvoxel, "voxel", areaLV, false, index);
      index++;
    }
  }

  // tube
  G4Tubs* stube = new G4Tubs("tube", 0.*mm, 19./2.*mm, dl, 0., 360.*deg);
  mate = G4Material::GetMaterial("Al");
  G4LogicalVolume* ltube = new G4LogicalVolume(stube, mate, "tube");
  va = new G4VisAttributes(G4Color(0.,0.8,0.8));
  ltube-> SetVisAttributes(va);

  G4RotationMatrix* rmtube = new G4RotationMatrix;
  rmtube-> rotateX(-90.*deg);
  new G4PVPlacement(rmtube, G4ThreeVector(),
                    ltube, "tube", lvoxel, false, 0);

  // cal
  const G4double dxycal = 25.*mm;
  const G4double dzcal = 3.*cm;

  G4Box* scal = new G4Box("cal", dxycal, dxycal, dzcal);
  mate = G4Material::GetMaterial("CsI");
  G4LogicalVolume* lcal = new G4LogicalVolume(scal, mate, "cal");
  va = new G4VisAttributes(G4Color(0.5,0.5,0.));
  lcal-> SetVisAttributes(va);

  index = 0;
  for ( ix = -2; ix <= 2; ix++ ) {
    G4double x0 = (5.*ix)*cm;
    new G4PVPlacement(0, G4ThreeVector(x0, 0., 2.*cm),
                      lcal, "cal", areaLV, false, index);
    index++;
  }

  return area;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
}
