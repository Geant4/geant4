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
#include "Par03DetectorConstruction.hh"
#include "Par03DetectorMessenger.hh"
#include "Par03SensitiveDetector.hh"
#include "Par03EMShowerModel.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"

#include "G4SDManager.hh"

#include "G4UnitsTable.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03DetectorConstruction::Par03DetectorConstruction()
  : G4VUserDetectorConstruction()
{
  fDetectorMessenger = new Par03DetectorMessenger(this);

  G4NistManager* nistManager = G4NistManager::Instance();
  fDetectorMaterial          = nistManager->FindOrBuildMaterial("G4_Fe");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03DetectorConstruction::~Par03DetectorConstruction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Par03DetectorConstruction::Construct()
{
  //--------- Material definition ---------
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* air            = nistManager->FindOrBuildMaterial("G4_AIR");

  //--------- Derived dimensions ---------
  G4double full2Pi        = 2. * CLHEP::pi * rad;
  G4double layerThickness = fDetectorLength / fNbOfLayers;
  G4double cellPhi        = full2Pi / fNbOfPhiCells;
  G4double cellDR         = fDetectorRadius / fNbOfRhoCells;

  //--------- World ---------
  auto fSolidWorld  = new G4Box("World",                  // name
                               fWorldSize / 2.,          // half-width in X
                               fWorldSize / 2.,          // half-width in Y
                               fWorldSize / 2.);         // half-width in Z
  auto fLogicWorld  = new G4LogicalVolume(fSolidWorld,    // solid
                                         air,            // material
                                         "World");       // name
  auto fPhysicWorld = new G4PVPlacement(0,                // no rotation
                                        G4ThreeVector(),  // at (0,0,0)
                                        fLogicWorld,      // logical volume
                                        "World",          // name
                                        0,                // mother volume
                                        false,            // not used
                                        999,              // copy number
                                        true);            // copy number

  //--------- Detector envelope ---------
  auto fSolidDetector = new G4Tubs("Detector",               // name
                                   0,                        // inner radius
                                   fDetectorRadius,          // outer radius
                                   fDetectorLength / 2.,     // half-width in Z
                                   0,                        // start angle
                                   full2Pi);                 // delta angle
  auto fLogicDetector = new G4LogicalVolume(fSolidDetector,  // solid
                                            fDetectorMaterial,  // material
                                            "Detector");        // name
  new G4PVPlacement(
    0,  // no rotation
    G4ThreeVector(0, 0,
                  fDetectorLength / 2),  // detector face starts at (0,0,0)
    fLogicDetector,                      // logical volume
    "Detector",                          // name
    fLogicWorld,                         // mother volume
    false,                               // not used
    99,                                  // copy number
    true);                               // check overlaps

  // Region for fast simulation
  auto detectorRegion = new G4Region("DetectorRegion");
  detectorRegion->AddRootLogicalVolume(fLogicDetector);

  //--------- Readout geometry ---------
  // Layers (along z)
  auto fSolidLayer = new G4Tubs("Layer",               // name
                                0,                     // inner radius
                                fDetectorRadius,       // outer radius
                                layerThickness / 2.,   // half-width in Z
                                0,                     // start angle
                                full2Pi);              // delta angle
  auto fLogicLayer = new G4LogicalVolume(fSolidLayer,  // solid
                                         air,          // material
                                         "Layer");     // name
  if(fNbOfLayers > 1)
    new G4PVReplica("Layer",          // name
                    fLogicLayer,      // logical volume
                    fLogicDetector,   // mother volume
                    kZAxis,           // axis of replication
                    fNbOfLayers,      // number of replicas
                    layerThickness);  // width of single replica
  else
    new G4PVPlacement(0,                // no rotation
                      G4ThreeVector(),  // place at centre of mother volume
                      fLogicLayer,      // logical volume
                      "Layer",          // name
                      fLogicDetector,   // mother volume
                      false,            // not used
                      0,                // copy number
                      true);            // check overlaps

  // Layer segment (division in phi)
  auto fSolidRow = new G4Tubs("Row",                // name
                              0,                    // inner radius
                              fDetectorRadius,      // outer radius
                              layerThickness / 2.,  // half-width in Z
                              0,                    // start angle
                              cellPhi);             // delta angle

  auto fLogicRow = new G4LogicalVolume(fSolidRow,   // solid
                                       air,         // material
                                       "Segment");  // name
  if(fNbOfPhiCells > 1)
    new G4PVReplica("Segment",      // name
                    fLogicRow,      // logical volume
                    fLogicLayer,    // mother volume
                    kPhi,           // axis of replication
                    fNbOfPhiCells,  // number of replicas
                    cellPhi);       // width of single replica
  else
    new G4PVPlacement(0,                // no rotation
                      G4ThreeVector(),  // place at centre of mother volume
                      fLogicRow,        // logical volume
                      "Row",            // name
                      fLogicLayer,      // mother volume
                      false,            // not used
                      0,                // copy number
                      true);            // check overlaps

  // Final cells (segment slices in radius)
  // No volume can be placed inside a radial replication
  auto fSolidCell = new G4Tubs("Cell",               // name
                               0,                    // inner radius
                               cellDR,               // outer radius
                               layerThickness / 2.,  // half-width in Z
                               0,                    // start angle
                               cellPhi);             // delta angle

  fLogicCell = new G4LogicalVolume(fSolidCell,         // solid
                                   fDetectorMaterial,  // material
                                   "Cell");            // name
  if(fNbOfRhoCells > 1)
    new G4PVReplica("Cell",         // name
                    fLogicCell,     // logical volume
                    fLogicRow,      // mother volume
                    kRho,           // axis of replication
                    fNbOfRhoCells,  // number of replicas
                    cellDR);        // width of single replica
  else
    new G4PVPlacement(0,                // no rotation
                      G4ThreeVector(),  // place at centre of mother volume
                      fLogicCell,       // logical volume
                      "Cell",           // name
                      fLogicRow,        // mother volume
                      false,            // not used
                      0,                // copy number
                      true);            // check overlaps

  //--------- Visualisation settings ---------
  fLogicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  fLogicLayer->SetVisAttributes(G4VisAttributes::GetInvisible());
  fLogicRow->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VisAttributes attribs;
  attribs.SetColour(G4Colour(0, 0, 1, 0.3));
  attribs.SetForceSolid(true);
  fLogicCell->SetVisAttributes(attribs);

  Print();
  return fPhysicWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03DetectorConstruction::ConstructSDandField()
{
  Par03SensitiveDetector* caloSD = new Par03SensitiveDetector(
    "sensitiveDetector", fNbOfLayers, fNbOfPhiCells, fNbOfRhoCells);
  G4SDManager::GetSDMpointer()->AddNewDetector(caloSD);
  SetSensitiveDetector(fLogicCell, caloSD);

  auto detectorRegion =
    G4RegionStore::GetInstance()->GetRegion("DetectorRegion");
  new Par03EMShowerModel("model", detectorRegion);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03DetectorConstruction::Print() const
{
  G4cout << "\n------------------------------------------------------"
         << "\n--- Detector material:\t" << fDetectorMaterial->GetName()
         << "\n--- Detector length:\t" << G4BestUnit(fDetectorLength, "Length")
         << "\n--- Detector radius:\t" << G4BestUnit(fDetectorRadius, "Length")
         << "\n--- Number of layers:\t" << fNbOfLayers
         << "\n--- Number of R-cells:\t" << fNbOfRhoCells
         << "\n--- Number of phi-cells:\t" << fNbOfPhiCells << G4endl;
  G4cout << "-----------------------------------------------------" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03DetectorConstruction::SetMaterial(const G4String& aName)
{
  // search material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(aName);
  if(material)
    fDetectorMaterial = material;
  else
    G4Exception("Par03DetectorConstruction::SetMaterial()", "InvalidSetup",
                FatalException, ("Unknown material name: " + aName).c_str());
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03DetectorConstruction::SetRadius(G4double aRadius)
{
  // check if fits within world volume
  if(aRadius >= fWorldSize / 2.)
    G4Exception("Par03DetectorConstruction::SetRadius()", "InvalidSetup",
                FatalException,
                ("Detector radius cannot be larger than the world size (" +
                 G4String(G4BestUnit(fWorldSize / 2., "Length")) + ")")
                  .c_str());
  fDetectorRadius = aRadius;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03DetectorConstruction::SetLength(G4double aLength)
{
  // check if fits within world volume
  if(aLength >= fWorldSize / 2.)
    G4Exception("Par03DetectorConstruction::SetLength()", "InvalidSetup",
                FatalException,
                ("Detector length cannot be larger than the world size (" +
                 G4String(G4BestUnit(fWorldSize / 2., "Length")) + ")")
                  .c_str());
  fDetectorLength = aLength;
}