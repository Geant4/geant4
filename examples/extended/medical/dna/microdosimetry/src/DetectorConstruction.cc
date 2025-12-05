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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 45, (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//

#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"
#include "PhysicsList.hh"

#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction(PhysicsList* ptr)
  : G4VUserDetectorConstruction(),
    fpWaterMaterial(nullptr),
    fLogicWorld(nullptr),
    fPhysiWorld(nullptr)
{
  fDetectorMessenger = new DetectorMessenger(this, ptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::DefineMaterials()
{
  // Water is defined from NIST material database
  G4NistManager* man = G4NistManager::Instance();

  G4Material* H2O = man->FindOrBuildMaterial("G4_WATER");

  /*
   If one wishes to test other density value for water material,
   one should use instead:
   G4Material * H2O = man->BuildMaterialWithNewDensity("G4_WATER_MODIFIED",
   "G4_WATER",1.100*g/cm3);

   Note: any string for "G4_WATER_MODIFIED" parameter is accepted
   and "G4_WATER" parameter should not be changed
   Both materials are created and can be selected from dna.mac
  */
  fpWaterMaterial = H2O;

  // G4cout << "-> Density of water material (g/cm3)="
  //  << fpWaterMaterial->GetDensity()/(g/cm/cm/cm) << G4endl;

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if (fPhysiWorld) {
    return fPhysiWorld;
  }
  DefineMaterials();

  // World volume

  G4double worldSizeX = 10 * micrometer;
  G4double worldSizeY = worldSizeX;
  G4double worldSizeZ = worldSizeX;

  G4Box* solidWorld = new G4Box("World",  // its name
                                worldSizeX / 2, worldSizeY / 2, worldSizeZ / 2);  // its size

  fLogicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                    fpWaterMaterial,  // its material
                                    "World");  // its name

  fPhysiWorld = new G4PVPlacement(0,  // no rotation
                                  G4ThreeVector(),  // at (0,0,0)
                                  "World",  // its name
                                  fLogicWorld,  // its logical volume
                                  0,  // its mother  volume
                                  false,  // no boolean operation
                                  0);  // copy number

  // Target volume

  G4Box* solidTarget = new G4Box("volumeTarget",  // its name
                                 worldSizeX / 10, worldSizeY / 2, worldSizeZ / 2);  // its size

  fLogicTarget = new G4LogicalVolume(solidTarget,  // its solid
                                     fpWaterMaterial,  // its material
                                     "volumeTarget");  // its name

  fPhysiTarget = new G4PVPlacement(0,  // no rotation
                                   G4ThreeVector(),  // at (0,0,0)
                                   "volumeTarget",  // its name
                                   fLogicTarget,  // its logical volume
                                   fPhysiWorld,  // its mother  volume
                                   false,  // no boolean operation
                                   0);  // copy number

  // Target region

  G4Region* targetRegion = new G4Region("regionTarget");

  // Cuts needed to allow combination of Physics
  // Default EM settings for fast simulation for e+, e-, gammas
  // 1 um proton cut avoids creation of very slow recoil ions
  G4ProductionCuts* cuts = new G4ProductionCuts();
  cuts->SetProductionCut(0.7 * mm, "gamma");
  cuts->SetProductionCut(0.7 * mm, "e-");
  cuts->SetProductionCut(0.7 * mm, "e+");
  cuts->SetProductionCut(1 * micrometer, "proton");
  targetRegion->SetProductionCuts(cuts);

  targetRegion->AddRootLogicalVolume(fLogicTarget);

  // Visualization attributes - white
  G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  worldVisAtt->SetVisibility(true);
  fLogicWorld->SetVisAttributes(worldVisAtt);

  G4VisAttributes* targetVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  targetVisAtt->SetVisibility(true);
  fLogicTarget->SetVisAttributes(targetVisAtt);

  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& materialChoice)
{
  // Search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fpWaterMaterial = pttoMaterial;
    if (fLogicWorld) {
      fLogicWorld->SetMaterial(fpWaterMaterial);
    }
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}
