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

// User Classes
#include "Par04ParallelFastWorld.hh"
#include "Par04ParallelFullWorld.hh"
#include "Par04ParallelFastSensitiveDetector.hh"

// G4 Classes
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoDelete.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04ParallelFastWorld::Par04ParallelFastWorld(G4String aWorldName,
                                               const Par04DetectorConstruction* aMassDetector,
                                               const Par04ParallelFullWorld* aParallelFull)
  :G4VUserParallelWorld(aWorldName), fMassDetector(aMassDetector), fParallelFull(aParallelFull) {
  fNbOfLayers = fMassDetector->GetNbOfLayers();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04ParallelFastWorld::~Par04ParallelFastWorld() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04ParallelFastWorld::Construct()
{
  // In parallel world material does not matter
  G4Material* dummy   = nullptr;

  // Build parallel geometry:
  auto parallelLogicalVolume = GetWorld()->GetLogicalVolume();

  G4double detectorInnerRadius = fMassDetector->GetInnerRadius();
  G4double detectorLength = fMassDetector->GetLength();
  G4double fullLayerThickness = fMassDetector->GetAbsorberThickness(0)
    + fMassDetector->GetAbsorberThickness(1);
  fLayerThickness = fullLayerThickness;
  // Get an updated value
  fNbOfLayers = fMassDetector->GetNbOfLayers();
  fNbOfSlices = fParallelFull->GetNbOfSlices();
  fNbOfRows = fParallelFull->GetNbOfRows();
  G4double detectorRadius = fNbOfLayers * fullLayerThickness;
  G4double detectorOuterRadius = detectorInnerRadius + detectorRadius;
  G4double rowThickness = detectorLength / fNbOfRows;
  G4double full2Pi = 2.* CLHEP::pi * rad;
  Print();

  // Insert cells to create a readout structure that contains both passive and active materials
  // Mostly a copy from the detector construction
  auto solidDetector = new G4Tubs("Detector",                // name
                                   detectorInnerRadius,      // inner radius
                                   detectorOuterRadius,      // outer radius
                                   detectorLength / 2.,      // half-width in Z
                                   0,                        // start angle
                                   full2Pi);                 // delta angle
  auto logicDetector = new G4LogicalVolume(solidDetector,    // solid
                                            dummy,           // material
                                            "Detector");     // name
  new G4PVPlacement(0,                                       // no rotation
                    G4ThreeVector(0, 0, 0),                  // detector centre at (0,0,0)
                    logicDetector,                           // logical volume
                    "Detector",                              // name
                    parallelLogicalVolume,                   // mother volume
                    false,                                   // not used
                    9999,                                    // copy number
                    true);                                   // check overlaps


  //--------- Detector cylinder (division along z axis)  ---------
  auto solidRow = new G4Tubs("Row", detectorInnerRadius, detectorOuterRadius, rowThickness / 2.,
                             0,  full2Pi);

  auto logicRow = new G4LogicalVolume(solidRow, dummy, "Row");
  if (fNbOfRows > 1)
    new G4PVReplica("Row",
                    logicRow,
                    logicDetector,
                    kZAxis,
                    fNbOfRows,
                    rowThickness);
  else
    new G4PVPlacement(0,
                      G4ThreeVector(),
                      logicRow,
                      "Row",
                      logicDetector,
                      false,
                      0);

  //--------- Detector slices (division in azimuthal angle)  ---------
  G4double cellPhi = full2Pi / fNbOfSlices;
  auto solidSlice = new G4Tubs("Slice", detectorInnerRadius, detectorOuterRadius, rowThickness/2,
                               0, cellPhi);
  auto logicSlice = new G4LogicalVolume(solidSlice,
                                        dummy,
                                        "Slice");
  if(fNbOfSlices>1) {
    new G4PVReplica("Slice",
                    logicSlice,
                    logicRow,
                    kPhi,
                    fNbOfSlices,
                    cellPhi,
                    -cellPhi);
  } else {
      new G4PVPlacement(0,
                        G4ThreeVector(),
                        logicSlice,
                        "Slice",
                        logicRow,
                        false,
                        0);
  }

  //--------- Detector cells (division along radial axis)  ---------
  G4VisAttributes attribs;
  attribs.SetColour(G4Colour(0, 1, 0, 0.1));
  attribs.SetForceSolid(true);
  if(fNbOfLayers>1) {
    auto solidCell = new G4Tubs("Cell", detectorInnerRadius,
                                detectorInnerRadius + fLayerThickness,
                                rowThickness/2, 0, cellPhi);
    fLogicalCell.push_back(new G4LogicalVolume(solidCell, dummy, "Cell_0"));
    new G4PVReplica("Cell",
                    fLogicalCell.back(),
                    logicSlice,
                    kRho,
                    fNbOfLayers,
                    fLayerThickness,
                    detectorInnerRadius);
  } else {
    auto solidCell = new G4Tubs("Cell", detectorInnerRadius,
                                  detectorInnerRadius + fLayerThickness,
                                  rowThickness/2, 0, cellPhi);
      fLogicalCell.push_back(new G4LogicalVolume(solidCell, dummy, "Cell"));
      fLogicalCell.back()->SetVisAttributes(attribs);
      new G4PVPlacement(0,
                        G4ThreeVector(),
                        fLogicalCell.back(),
                        "Cell",
                        logicSlice,
                        false,
                        0);
  }
  Print();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04ParallelFastWorld::ConstructSD()
{
  // -- sensitive detectors:
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  Par04ParallelFastSensitiveDetector* caloSD =
    new Par04ParallelFastSensitiveDetector("parallelFastSD", fNbOfLayers, fNbOfSlices);
  SDman->AddNewDetector(caloSD);
  for(const auto& logicalCell: fLogicalCell)
   logicalCell->SetSensitiveDetector(caloSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04ParallelFastWorld::Print() {
  G4cout << "\n------------------------------------------------------"
         << "\n Readout geometry with physics layout is set in parallel geometry:\t"
         << "\n Cylindrical detector is divided along radius (layers), phi (slices), and z (rows)."
         << "\n Number of layers is determined by number of layers set in detector construction. "
         << "\n- Number of layers: " << fNbOfLayers << "\n------- Number of slices: " << fNbOfSlices
         << "\n- Number of rows: " << fNbOfRows;
  G4cout << "\n Readout will collect energy from fast simulation.\n------- Thickness is "
         <<  "a sum of all absorbers" << " = " << G4BestUnit(fLayerThickness, "Length")
         << "\n-----------------------------------------------------" << G4endl;
}
