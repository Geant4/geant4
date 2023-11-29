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
/// \file DetectorConstruction.cc
/// \brief DetectorConstruction class
//
//
#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4ProductionCuts.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"

const bool check_intersections = true; // to control geometry for errors

// new unit
const G4double ug = 1.e-6 * g;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct() {
  auto man = G4NistManager::Instance();

  // registering new useful unit
  new G4UnitDefinition("microgram", "ug", "Mass", ug);

  // geometric parameters
  G4double worldSize = 8 * cm;
  G4double worldDensity = 0.15 * mg / m3; // imperfect vacuum as in a real-world experiment

  G4double IV_diameter = 10 * mm;
  G4double IV_height = 20 * mm;
  G4double IV_vertical_offset =
      -5 * mm; // instead of moving beam axis, the IV is moved
  G4double IV_density = 0.45 * ug / cm3;

  G4double wall_tot_mass_thickness =
      0.98 * mg / cm2; // Mylar wall thickness -> 7 µm
  G4double wall_inner_layer_mass_thickness =
      0.02 * mg /
      cm2; // inner layer of wall made of water material at Mylar density

  fCollDiameter =
      3 * mm; // inner diameter of the collimator, the outer is twice that
  fCollLength = 23 * mm;
  fCollExitPosition =
      -5.5 * mm; // exit of the collimator, where the beam is already formed

  G4double enDet_diameter = 10 * mm;
  G4double enDet_thickness = 150 * um;
  G4double enDet_distance = 10 * mm;

  // materials
  auto worldMaterial = man->BuildMaterialWithNewDensity(
      "WORLD_WATER_VACUUM", "G4_WATER", worldDensity);
  auto targetMaterial =
      man->BuildMaterialWithNewDensity("TARGET_WATER", "G4_WATER", IV_density);
  auto wallMaterial = man->FindOrBuildMaterial("G4_MYLAR");
  auto collMaterial = man->FindOrBuildMaterial("G4_BRASS");
  auto enDetMaterial = man->FindOrBuildMaterial("G4_Si");

  // auxiliary variables:
  G4ThreeVector origin(0, 0, 0);
  G4ThreeVector IV_centre(0, 0, IV_vertical_offset);
  G4ThreeVector enDet_centre(enDet_distance + IV_diameter / 2., 0, 0);

  G4double mylarDensity = wallMaterial->GetDensity();
  G4double wall_total_thickness = wall_tot_mass_thickness / mylarDensity;
  G4double wall_inner_layer_thickness =
      wall_inner_layer_mass_thickness / mylarDensity;

  // materials continued
  auto innerWallMaterial =
      man->BuildMaterialWithNewDensity("WALL_WATER", "G4_WATER", mylarDensity);

  // SAVING DETECTOR SETTINGS TO FILE:
  auto filename = "SIM_GEOMETRY_SETTINGS.txt";
  std::ofstream TextFile;
  TextFile.open(filename, std::fstream::out);
  TextFile << "worldSize = " << worldSize / cm << " cm\n";
  TextFile << "worldDensity = " << worldDensity / (mg / m3) << " mg/m3\n";
  TextFile << "IV_diameter = " << IV_diameter / mm << " mm\n";
  TextFile << "IV_height = " << IV_height / mm << " mm\n";
  TextFile << "IV_vertical_offset = " << IV_vertical_offset / mm << " mm\n";
  TextFile << "IV_density = " << IV_density / (ug / cm3) << " ug/cm3\n";
  TextFile << "wall_tot_mass_thickness = "
           << wall_tot_mass_thickness / (mg / cm2) << " mg/cm2\n";
  TextFile << "wall_inner_layer_mass_thickness = "
           << wall_inner_layer_mass_thickness / (mg / cm2) << " mg/cm2\n";
  TextFile << "fCollDiameter = " << fCollDiameter / mm << " mm\n";
  TextFile << "fCollLength = " << fCollLength / mm << " mm\n";
  TextFile << "fCollExitPosition = " << fCollExitPosition / mm << " mm\n";
  TextFile << "enDet_diameter = " << enDet_diameter / mm << " mm\n";
  TextFile << "enDet_thickness = " << enDet_thickness / mm << " mm\n";
  TextFile << "enDet_distance = " << enDet_distance / mm << " mm\n";
  TextFile << "worldMaterial: " << worldMaterial->GetName() << "\n";
  TextFile << "targetMaterial: " << targetMaterial->GetName() << "\n";
  TextFile << "wallMaterial: " << wallMaterial->GetName() << "\n";
  TextFile << "collMaterial: " << collMaterial->GetName() << "\n";
  TextFile << "enDetMaterial: " << enDetMaterial->GetName() << "\n";
  TextFile << "innerWallMaterial: " << innerWallMaterial->GetName() << "\n";

  // WORLD VOLUME

  auto worldSolid = new G4Box("worldSolid", // its name
                              worldSize / 2, worldSize / 2,
                              worldSize / 2); // its size

  auto worldLogic = new G4LogicalVolume(worldSolid,    // its solid
                                        worldMaterial, // its material
                                        "worldLogic"); // its name

  auto worldPhys =
      new G4PVPlacement(nullptr,                // no rotation
                        G4ThreeVector(0, 0, 0), // placement
                        worldLogic,             // its logical volume
                        "worldPhys",            // its name
                        nullptr,                // its mother volume
                        false,                  // no boolean operation
                        0,                      // copy number
                        check_intersections);   // check intersections

  // INTERACTION VOLUME (IV, also called chamber)
  auto chamberSolid = new G4Tubs("chamberSolid",   // its name
                                 0,                // rMin
                                 IV_diameter / 2., // rMax
                                 IV_height / 2.,   // height/2
                                 0 * deg,          // phiMin
                                 360 * deg);       // phiMax

  auto chamberLogic = new G4LogicalVolume(chamberSolid,    // its solid
                                          targetMaterial,  // its material
                                          "chamberLogic"); // its name

  new G4PVPlacement(nullptr, IV_centre, chamberLogic, "chamberPhys", worldLogic,
                    false, 0, check_intersections);

  // SENSITIVE VOLUME (SV, also called target) - in general a sub-volume of the
  // IV, but in this case entire IV is the SV only ionisations that occur in the
  // SV will be counted
  auto targetSolid = new G4Tubs("targetSolid",    // name
                                0,                // rMin
                                IV_diameter / 2., // rMax
                                IV_height / 2.,   // height/2
                                0 * deg,          // phiMin
                                360 * deg);       // phiMax

  auto targetLogic =
      new G4LogicalVolume(targetSolid, targetMaterial, "targetLogic");

  new G4PVPlacement(nullptr, origin, targetLogic, "targetPhys", chamberLogic,
                    false, 0, check_intersections);

  // WALLS:
  //

  auto wallSolid = new G4Tubs("wallSolid",                             // name
                              IV_diameter / 2.,                        // rMin
                              IV_diameter / 2. + wall_total_thickness, // rMax
                              IV_height / 2., // height/2
                              0 * deg,        // phiMin
                              360 * deg);     // phiMax

  auto wallLogic = new G4LogicalVolume(wallSolid, wallMaterial, "wallSolid");
  new G4PVPlacement(nullptr, IV_centre, wallLogic, "wallPhys", worldLogic,
                    false, 0, check_intersections);

  // wall inner layer made of water
  auto innerWallSolid =
      new G4Tubs("innerWallSolid",                              // name
                 IV_diameter / 2.,                              // rMin
                 IV_diameter / 2. + wall_inner_layer_thickness, // rMax
                 IV_height / 2.,                                // height/2
                 0 * deg,                                       // phiMin
                 360 * deg);                                    // phiMax

  auto innerWallLogic =
      new G4LogicalVolume(innerWallSolid, innerWallMaterial, "innerWallLogic");
  new G4PVPlacement(nullptr, origin, innerWallLogic, "innerWallPhys", wallLogic,
                    false, 0, check_intersections);

  // COLLIMATOR:
  auto collSolid = new G4Tubs("collSolid",        // name
                              fCollDiameter / 2., // rMin
                              fCollDiameter,      // rMax
                              fCollLength / 2.,   // height/2
                              0 * deg,            // phiMin
                              360 * deg);         // phiMax

  auto collLogic = new G4LogicalVolume(collSolid, collMaterial, "collSolid");

  auto rot = new G4RotationMatrix();
  rot->rotateY(90 * deg);
  new G4PVPlacement(
      rot, G4ThreeVector(-fCollLength / 2 + fCollExitPosition, 0, 0), collLogic,
      "collPhys", worldLogic, false, 0, check_intersections);

  // SILICON DETECTOR - present only in macrometric geometry

  auto enDetSolid = new G4Tubs("enDetSolid",         // name
                               0,                    // rMin
                               enDet_diameter / 2.,  // rMax
                               enDet_thickness / 2., // height/2
                               0 * deg, 360 * deg);

  auto enDetLogic =
      new G4LogicalVolume(enDetSolid, enDetMaterial, "enDetLogic");
  new G4PVPlacement(rot, enDet_centre, enDetLogic, "enDetPhys", worldLogic,
                    false, 0, check_intersections);

  // VISUALISATION ATTRIBUTES

  auto worldVisAtt = new G4VisAttributes(G4Colour(1, 1, 1, 0.1));
  worldLogic->SetVisAttributes(worldVisAtt);

  auto targetVisAtt = new G4VisAttributes(G4Colour(0.1, 0.5, 1, 0.7));
  targetLogic->SetVisAttributes(targetVisAtt);

  auto innerWallVisAtt = new G4VisAttributes(G4Colour(0, 1, 1, 0.7));
  innerWallLogic->SetVisAttributes(innerWallVisAtt);

  auto collVisAtt = new G4VisAttributes(G4Colour(0.7, 0.65, .25, 0.5));
  collLogic->SetVisAttributes(collVisAtt);

  auto enDetVisAtt = new G4VisAttributes(G4Colour(0.5, 0.7, .5, 1.));
  enDetLogic->SetVisAttributes(enDetVisAtt);

  // Create Target G4Region and add logical volume

  auto region = new G4Region("Target");

  auto cuts = new G4ProductionCuts();

  G4double defCut = 1 * nanometer;
  cuts->SetProductionCut(defCut, "gamma");
  cuts->SetProductionCut(defCut, "e-");
  cuts->SetProductionCut(defCut, "e+");
  cuts->SetProductionCut(defCut, "proton");

  region->SetProductionCuts(cuts);
  region->AddRootLogicalVolume(chamberLogic);
  region->AddRootLogicalVolume(targetLogic);
  region->AddRootLogicalVolume(innerWallLogic);

  return worldPhys;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
