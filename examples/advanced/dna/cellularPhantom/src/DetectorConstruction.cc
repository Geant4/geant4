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
// -----------------------------------------------------------------------------
//       MONTE CARLO SIMULATION OF REALISTIC GEOMETRY FROM MICROSCOPES IMAGES
//
// Authors and contributors:
// P. Barberet (a), S. Incerti (a), N. H. Tran (a), L. Morelli (a,b)
//
// a) University of Bordeaux, CNRS, LP2i, UMR5797, Gradignan, France
// b) Politecnico di Milano, Italy
//
// If you use this code, please cite the following publication:
// P. Barberet et al.,
// "Monte-Carlo dosimetry on a realistic cell monolayer
// geometry exposed to alpha particles."
// Ph. Barberet et al 2012 Phys. Med. Biol. 57 2189
// doi: 110.1088/0031-9155/57/8/2189
// -----------------------------------------------------------------------------

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"
#include "G4ProductionCuts.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction()
{
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructLine();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::DefineMaterials()
{
  G4String name, symbol;

  // Water and air are defined from NIST material database
  G4NistManager *man = G4NistManager::Instance();

  G4Material *H2O = man->FindOrBuildMaterial("G4_WATER");
  G4Material *Air = man->FindOrBuildMaterial("G4_AIR");

  fDefaultMaterial = Air;
  fPhantomMaterial = H2O;  // material is not relevant
                           // it will be changed by the ComputeMaterial
                           // method of the CellParameterisation

  // Default materials
  if (fMediumMaterial == nullptr) {fMediumMaterial = H2O;}
  if (fRedMaterial == nullptr) {fRedMaterial = H2O;}
  if (fGreenMaterial == nullptr) {fGreenMaterial = H2O;}
  if (fBlueMaterial == nullptr) {fBlueMaterial = H2O;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume *DetectorConstruction::ConstructLine() {

  //*************
  // World volume
  //*************

  fSolidWorld = new G4Box("World",                         //its name
                          fWorldSizeXY / 2, fWorldSizeXY / 2, fWorldSizeZ / 2);  //its size

  fLogicWorld = new G4LogicalVolume(fSolidWorld,           //its solid
                                    fDefaultMaterial,      //its material
                                    "World");              //its name

  fPhysiWorld = new G4PVPlacement(nullptr,                 //no rotation
                                  G4ThreeVector(),         //at (0,0,0)
                                  "World",                 //its name
                                  fLogicWorld,             //its logical volume
                                  nullptr,                 //its mother  volume
                                  false,                   //no boolean operation
                                  0);                      //copy number

  //********************
  // Cell culture medium
  //********************

  fSolidMedium = new G4Box("Medium", fMediumSizeXY / 2, fMediumSizeXY / 2, fMediumSizeZ / 2);

  fLogicMedium = new G4LogicalVolume(fSolidMedium, fMediumMaterial, "Medium");

  fPhysiMedium = new G4PVPlacement(nullptr,
                                   G4ThreeVector(0, 0, 0),
                                   "Medium",
                                   fLogicMedium,
                                   fPhysiWorld,
                                   false,
                                   0);

  // ************
  // Cell phantom
  // ************

  // The cell phantom is placed in the middle of the parent volume (fLogicMedium here)

  fPhantomParam = new CellParameterisation
    (fPhantomFileName, fRedMaterial, fGreenMaterial, fBlueMaterial, fShiftX, fShiftY, fShiftZ);

  fSolidPhantom = new G4Box("Phantom",
                            fPhantomParam->GetPixelSizeX() / 2,
                            fPhantomParam->GetPixelSizeY() / 2,
                            fPhantomParam->GetPixelSizeZ() / 2);

  fLogicPhantom = new G4LogicalVolume(fSolidPhantom,
                                      fPhantomMaterial, // material is not relevant,
                                                        // it will be changed by the
                                                        // ComputeMaterial method
                                                        // of the CellParameterisation
                                      "Phantom",
                                      nullptr,
                                      nullptr,
                                      nullptr);

  fPhysiPhantom = new G4PVParameterised(
      "Phantom",         // name
      fLogicPhantom,     // logical volume
      fLogicMedium,      // mother logical volume
      kUndefined,        // kUndefined: three-dimensional optimization
      fPhantomParam->GetPhantomTotalPixels(), // number of voxels
      fPhantomParam,     // the parametrisation
      false);

  G4cout << " #########################################################################" << G4endl;
  G4cout << "                               Phantom information                        " << G4endl;
  G4cout << " #########################################################################" << G4endl;
  G4cout << G4endl;

  G4cout << " ==========> The phantom contains " << fPhantomParam->GetPhantomTotalPixels()
         << " voxels " << G4endl;
  G4cout << " ==========> Voxel size X (um) = " << fPhantomParam->GetPixelSizeX()/um << G4endl;
  G4cout << " ==========> Voxel size Y (um) = " << fPhantomParam->GetPixelSizeY()/um << G4endl;
  G4cout << " ==========> Voxel size Z (um) = " << fPhantomParam->GetPixelSizeZ()/um << G4endl;
  G4cout << G4endl;

  G4cout << " ==========> Number of red voxels = "
         << fPhantomParam->GetRedTotalPixels() << G4endl;
  G4cout << " ==========> Number of green voxels = "
         << fPhantomParam->GetGreenTotalPixels() << G4endl;
  G4cout << " ==========> Number of blue voxels = "
         << fPhantomParam->GetBlueTotalPixels() << G4endl;
  G4cout << G4endl;

  G4cout << " ==========> Tolal mass of red voxels (kg) = "
         << fPhantomParam->GetRedMass() / kg << G4endl;
  G4cout << " ==========> Tolal mass of green voxels (kg) = "
         << fPhantomParam->GetGreenMass() / kg << G4endl;
  G4cout << " ==========> Tolal mass of blue voxels (kg) = "
        << fPhantomParam->GetBlueMass() / kg << G4endl;
  G4cout << G4endl;
  G4cout << " #########################################################################" << G4endl;
  G4cout << G4endl;

  // USER LIMITS ON STEP LENGTH

  // fLogicWorld->SetUserLimits(new G4UserLimits(100 * mm));
  // fLogicPhantom->SetUserLimits(new G4UserLimits(0.5 * micrometer));
  // fLogicMedium->SetUserLimits(new G4UserLimits(1 * micrometer));

  // Create a phantom G4Region and add logical volume

  fPhantomRegion = new G4Region("phantomRegion");

  G4ProductionCuts* cuts = new G4ProductionCuts();

  G4double defCut = 1*nanometer;
  cuts->SetProductionCut(defCut,"gamma");
  cuts->SetProductionCut(defCut,"e-");
  cuts->SetProductionCut(defCut,"e+");
  cuts->SetProductionCut(defCut,"proton");

  fPhantomRegion->SetProductionCuts(cuts);
  fPhantomRegion->AddRootLogicalVolume(fLogicMedium);

  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(const G4String& mat)
{
  if (G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat))
  {
    if (material && mat != "G4_WATER")
    {
      fMediumMaterial = material;
      G4cout << " #########################################################################"
             << G4endl;
      G4cout << "                               Cell culture medium material               "
             << G4endl;
      G4cout << fMediumMaterial << G4endl;
      G4cout << " #########################################################################"
             << G4endl;
      G4cout << G4endl;
    }
  }
  else
  {
    G4cout << G4endl;
    G4cout << "WARNING: material \"" << mat << "\" doesn't exist in NIST elements/materials"
           << G4endl;
    G4cout << " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]"
           << G4endl;
    G4cout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetRedDensity(const G4double& value)
{
  fDensityRed = value;
  if (fDensityRed != 1.0)
  {
    G4NistManager *man = G4NistManager::Instance();
    G4Material * H2O_red = man->BuildMaterialWithNewDensity("G4_WATER_red","G4_WATER",
                                                            fDensityRed);
    fRedMaterial = H2O_red;
  }
  else
  {
    G4NistManager *man = G4NistManager::Instance();
    fRedMaterial = man->FindOrBuildMaterial("G4_WATER");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGreenDensity(const G4double& value)
{
  fDensityGreen = value;
  if (fDensityGreen != 1.0)
  {
    G4NistManager *man = G4NistManager::Instance();
    G4Material * H2O_green = man->BuildMaterialWithNewDensity("G4_WATER_green","G4_WATER",
                                                              fDensityGreen);
    fGreenMaterial = H2O_green;
  }
  else
  {
    G4NistManager *man = G4NistManager::Instance();
    fGreenMaterial = man->FindOrBuildMaterial("G4_WATER");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetBlueDensity(const G4double& value)
{
  fDensityBlue = value;
  if (fDensityBlue != 1.0)
  {
    G4NistManager *man = G4NistManager::Instance();
    G4Material * H2O_blue = man->BuildMaterialWithNewDensity("G4_WATER_blue","G4_WATER",
                                                             fDensityBlue);
    fBlueMaterial = H2O_blue;
  }
  else
  {
    G4NistManager *man = G4NistManager::Instance();
    fBlueMaterial = man->FindOrBuildMaterial("G4_WATER");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetShiftX(const G4double& value)
{
  fShiftX = value;
  G4cout << "... setting phantom shift: X = " << fShiftX/um << " um" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetShiftY(const G4double& value)
{
  fShiftY = value;
  G4cout << "... setting phantom shift: Y = " << fShiftY/um << " um" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetShiftZ(const G4double& value)
{
  fShiftZ = value;
  G4cout << "... setting phantom shift: Y = " << fShiftZ/um << " um" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMediumSizeXY(const G4double& value)
{
  fMediumSizeXY = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMediumSizeZ(const G4double& value)
{
  fMediumSizeZ = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeXY(const G4double& value)
{
  fWorldSizeXY = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeZ(const G4double& value)
{
  fWorldSizeZ = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPhantomFileName(const G4String& phantomName)
{
  fPhantomFileName = phantomName;
  G4cout << " #########################################################################"
         << G4endl;
  G4cout << "                          Loading cell phantom from file:  "
         << fPhantomFileName << G4endl;
  G4cout << " #########################################################################"
         << G4endl;
  G4cout << G4endl;
}
