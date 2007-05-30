//******************************************************************************
// DetectorConstruction.cc
//
// 1.00 JMV, LLNL, MAR-2002:  First version.
//******************************************************************************
//
#include "DetectorConstruction.hh"

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4MaterialTable.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"

#include "globals.hh"

DetectorConstruction::DetectorConstruction()
{;}

DetectorConstruction::~DetectorConstruction()
{;}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  //------------------------------------------------------ materials

  // Air, from PhysicalConstants.h
  //
  G4NistManager* man = G4NistManager::Instance();
  G4Material* matAir = man->FindOrBuildMaterial("G4_AIR");

  // Polyethylene
  //
  G4Material* matPE  = man->FindOrBuildMaterial("G4_POLYETHYLENE");

  // He-3 detector materials
  //
  G4Material* matHe3  = new G4Material("He3",  2., 3.*g/mole, 0.00049*g/cm3, kStateGas);

  // Uranium ball material
  //
  G4Element* elHEU = new G4Element("Highly-enriched Uranium", "HEU", 2);
  elHEU->AddIsotope(new G4Isotope("U-235", 92, 235, 235.*g/mole), .93);
  elHEU->AddIsotope(new G4Isotope("U-238", 92, 238, 238.*g/mole), .07);
  G4Material* matHEU = new G4Material("HEU", 19.1*g/cm3, 1, kStateSolid);
  matHEU->AddElement(elHEU, 1.00);

  //------------------------------------------------------ volumes

  // ---------------------------------------------
  // World volume, filled with air
  // 
  G4Sphere* solidWorld = new G4Sphere("World", 0.0*cm, 80.0*cm, 
                                               0.0*deg, 360.0*deg,
                                               0.0*deg, 180.0*deg);
  G4LogicalVolume* logicWorld
    = new G4LogicalVolume(solidWorld, matAir, "World");
  G4VPhysicalVolume* physWorld 
    = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicWorld,
			"World ", 0, false, 0);

  (*logicWorld).SetVisAttributes(G4VisAttributes::Invisible);

  // ---------------------------------------------
  // 3.97-cm-radius highly-enriched Uranium sphere
  //
  G4Sphere* solidHEUball = new G4Sphere("HEU ball", 0.0*cm, 3.97*cm, 
                                               0.0*deg, 360.0*deg,
                                               0.0*deg, 180.0*deg);
  G4LogicalVolume* logicHEUball
    = new G4LogicalVolume(solidHEUball, matHEU, "HEU ball");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicHEUball,
 		    "HEU ball", logicWorld, false, 0);

  // ---------------------------------------------
  // 3" PE layer surrounding the HEU ball
  // 
  G4double PEThick = 3.*2.54*cm;
  G4Sphere* solidPELayer = new G4Sphere("PE layer", 3.97*cm, 3.97*cm+PEThick, 
                                               0.0*deg, 360.0*deg,
                                               0.0*deg, 180.0*deg);
  G4LogicalVolume* logicPELayer
    = new G4LogicalVolume(solidPELayer, matPE, "PE layer");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicPELayer,
 		    "PE layer", logicWorld, false, 0);

  // ---------------------------------------------
  // He3 detectors
  // 
  G4double He3Thick = 10.*cm;
  G4Sphere* solidHe3Detector = new G4Sphere("He-3 detector layer", 
                                               50.0*cm, 50.0*cm+He3Thick, 
                                               0.0*deg, 360.0*deg,
                                               0.0*deg, 180.0*deg);
  G4LogicalVolume* logicHe3Layer
    = new G4LogicalVolume(solidHe3Detector, matHe3, "He-3 detector");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicHe3Layer,
 		    "He3 detector", logicWorld, false, 0);

  //------------------------------------------------------------------
  // Must return pointer to the master physical volume
  //
  return physWorld;
}

