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
// $Id: BuildHouse.cc,v 1.1 2004-07-01 15:52:08 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "BuildHouse.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

G4VPhysicalVolume* BuildHouse()
{
  // Let's have x = east, y = up and z = south.
  const double worldx = 5 * m;         // World half-length east-west.
  const double worldy = 5 * m;         // World half-length vertical.
  const double worldz = 5 * m;         // World half-length north-south.
  const double roomx = 2 * m;          // Inner room half-length east-west.
  const double roomy = 1 * m;          // Inner room half-length vertical.
  const double roomz = 2 * m;          // Inner room half-length north-south.
  const double wall = 10 * cm;         // Wall thickness.

  // For other things, lengths, widths and/or depths are horizontal,
  // heights are vertical.
  const double wx = 1 * m;          // Window half-length.
  const double wy = 50 * cm;        // Window half-height.
  const double wz = wall + 1 * mm;  // Window half-depth, larger than
				    // wall for safe subtraction.
  const double ledgeProtrusion = 5 * cm;
  const double lx = 1.1 * m;        // Window ledge half-length.
  const double ly = 2.5 * cm;       // Window ledge half-height.
  const double lz = wall + ledgeProtrusion;  // Window ledge
				             // half-depth (protrudes
				             // into room and out into
				             // world).
  //const double ttx;  // Table top half-length.
  //const double tty;  // Table top half-height.
  //const double ttz;  // Table top half-width.
  //const double tlx;  // Table leg half-width.
  //const double tly;  // Table leg half-height (normally "length" of leg).
  //const double tlz;  // Table leg half-width.

  const G4ThreeVector nullVector;

  G4Material* air =
    new G4Material
    ("Air", 7., 14.4 * g/mole, 0.001 * g/cm3);
  G4Material* wallMaterial =
    new G4Material
    ("Brick", 20., 40. * g/mole, 3 * g/cm3);
  G4Material* tableMaterial =
    new G4Material
    ("Wood", 7., 14. * g/mole, 0.9 * g/cm3);
  
  // World...
  G4Box * world =
    new G4Box
    ("World", worldx, worldy, worldz);
  G4LogicalVolume * world_log =
    new G4LogicalVolume
    (world, air, "World", 0, 0, 0);
  world_log -> SetVisAttributes (G4VisAttributes::Invisible);
  G4VPhysicalVolume * world_phys =
    new G4PVPlacement
    (0, G4ThreeVector(), "worldP", world_log, 0, false, 0);

  // Window hole...
  G4VSolid* windowHole =
    new G4Box
    ("Window hole", wx, wy, wz);

  // Make a house out of walls...

  // North/south wall...
  G4VSolid* wallNS =
    new G4Box
    ("North/south wall", roomx + wall, roomy, wall / 2.);
  G4VSolid* wallNSWithWindowHole =
    new G4SubtractionSolid
    ("North/south wall with window hole", wallNS, windowHole,
     0, G4ThreeVector());

  // East/west wall...
  G4VSolid* wallEW =
    new G4Box
    ("North/south wall", roomx, roomy, wall / 2.);
  G4VSolid* wallEWWithWindowHole =
    new G4SubtractionSolid
    ("North/south wall with window hole", wallEW, windowHole,
     0, G4ThreeVector());

  // North wall...
  G4LogicalVolume * northWall_log =
    new G4LogicalVolume
    (wallNSWithWindowHole,
     wallMaterial, "North wall", 0, 0, 0);
  // G4VPhysicalVolume* northWall_phys =
    new G4PVPlacement
    (0, G4ThreeVector(0., 0., -(roomz + wall / 2.)),
     northWall_log, "North wall",
     world_log, false, 0);

  // South wall...
  G4LogicalVolume * southWall_log =
    new G4LogicalVolume
    (wallNSWithWindowHole,
     wallMaterial, "South wall", 0, 0, 0);
  G4RotationMatrix* rmS = new G4RotationMatrix;
  rmS->rotateY(180*deg);
  // G4VPhysicalVolume* southWall_phys =
    new G4PVPlacement
    (rmS, G4ThreeVector(0., 0., roomz + wall / 2.),
     southWall_log, "South wall",
     world_log, false, 0);

  // East wall...
  G4LogicalVolume * eastWall_log =
    new G4LogicalVolume
    (wallEWWithWindowHole,
     wallMaterial, "East wall", 0, 0, 0);
  G4RotationMatrix* rmE = new G4RotationMatrix;
  rmE->rotateY(-90*deg);
  // G4VPhysicalVolume* eastWall_phys =
    new G4PVPlacement
    (rmE, G4ThreeVector(roomz + wall / 2., 0., 0.),
     eastWall_log, "East wall",
     world_log, false, 0);

  // West wall...
  G4LogicalVolume * westWall_log =
    new G4LogicalVolume
    (wallEWWithWindowHole,
     wallMaterial, "West wall", 0, 0, 0);
  G4RotationMatrix* rmW = new G4RotationMatrix;
  rmW->rotateY(90*deg);
  // G4VPhysicalVolume* westWall_phys =
    new G4PVPlacement
    (rmW, G4ThreeVector(-(roomz + wall / 2.), 0., 0.),
     westWall_log, "West wall",
     world_log, false, 0);

  /*
  // House with holes for windows
  G4VSolid* houseOuter =
    new G4Box
    ("House outer", roomx + wall, roomy + wall, roomz + wall);
  G4VSolid* windowHole =
    new G4Box
    ("Window hole", wx, wy, wz);
  // East window...
  G4RotationMatrix* rm1 = new G4RotationMatrix;
  rm1->rotateY(-90*deg);
  G4VSolid* houseOuterWithWindowHoleE =
    new G4SubtractionSolid
    ("House with window E", houseOuter, windowHole,
     rm1, G4ThreeVector(roomx + wall / 2., 0., 0.));
  // North window...
  G4VSolid* houseOuterWithWindowHoleEN =
    new G4SubtractionSolid
    ("House with windows E and N", houseOuterWithWindowHoleE, windowHole,
     0, G4ThreeVector(0., 0., -(roomz + wall / 2.)));
  // Other windows...
  // Solid house...
  G4LogicalVolume * houseOuterWithWindowHoles_log =
    new G4LogicalVolume
    (houseOuterWithWindowHoleEN,
     wallMaterial, "House with window holes", 0, 0, 0);
  // G4VPhysicalVolume*  houseOuterWithWindowHoles_phys =
    new G4PVPlacement
    (0, G4ThreeVector(), houseOuterWithWindowHoles_log, "House",
     world_log, false, 0);
  // Inner...
  G4VSolid* houseInner =
    new G4Box
    ("House inner", roomx, roomy, roomz);
  // Room...
  G4LogicalVolume * room_log =
    new G4LogicalVolume
    (houseInner, air, "Room", 0, 0, 0);
  room_log -> SetVisAttributes (G4VisAttributes::Invisible);
  // G4VPhysicalVolume*  room_phys =
    new G4PVPlacement
    (0, G4ThreeVector(), room_log, "Room",
     houseOuterWithWindowHoles_log, false, 0);
  */

  /*
  G4VSolid* houseInner =
    new G4Box
    ("House inner", roomx, roomy, roomz);
  
  G4VSolid* walls
    new G4SubtractionSolid
    ("Walls", houseOuter, houseInner, 0, nullVector);
  
  G4VSolid* windowHole =
    new G4Box
    ("Window hole", wx, wy, wz);
  
  G4RotationMatrix* rm1 = new G4RotationMatrix;
  rm1->rotateZ(90*deg);
  G4VSolid* wallsWithWindowHole1 =
    new G4SubtractionSolid
    ("Walls with windows", walls, windowHole, rm1, roomx + wall / 2.);

  G4LogicalVolume * wallsWithWindowHole1_log =
    new G4LogicalVolume
    //(wallsWithWindowHole1, wallMaterial, "Walls with windows", 0, 0, 0);
    (walls, wallMaterial, "Walls", 0, 0, 0);

  G4VPhysicalVolume*  house =
    new G4PVPlacement
    (0, G4ThreeVector(), "House", wallsWithWindowHole1_log, world, false, 0);
  */

  return world_phys;
}
