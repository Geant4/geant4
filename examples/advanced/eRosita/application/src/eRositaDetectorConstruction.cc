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

#include "eRositaDetectorConstruction.hh"
#include "eRositaTrackerSD.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

eRositaDetectorConstruction::eRositaDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

eRositaDetectorConstruction::~eRositaDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto eRositaDetectorConstruction::Construct() -> G4VPhysicalVolume*
{
    // -----------------------
    // - Material definition -
    // -----------------------

/*
    // Nitrogen
    constexpr auto NITROGEN_ATOMIC_NUMBER{7.};
    constexpr auto NITROGEN_MOLAR_MASS{14.01 * g / mole};

    auto *nitrogen = new G4Element("Nitrogen", "N", NITROGEN_ATOMIC_NUMBER, NITROGEN_MOLAR_MASS);
    
    // Oxygen
    constexpr auto OXYGEN_ATOMIC_NUMBER{8.};
    constexpr auto OXYGEN_MOLAR_MASS{16.00 * g / mole};

    auto *oxygen = new G4Element("Oxygen", "O", OXYGEN_ATOMIC_NUMBER, OXYGEN_MOLAR_MASS);
    
    // Air
    constexpr auto AIR_DENSITY{1.29 * mg / cm3};
    constexpr auto AIR_NUMBER_OF_ELEMENTS{2};
    constexpr auto NITROGEN_PERCENTAGE{70};
    constexpr auto OXYGEN_PERCENTAGE{30};

    auto *air = new G4Material("Air", AIR_DENSITY, AIR_NUMBER_OF_ELEMENTS);
    air->AddElement(nitrogen, NITROGEN_PERCENTAGE * perCent);
    air->AddElement(oxygen, OXYGEN_PERCENTAGE * perCent);
*/

    // Copper
    constexpr auto COPPER_ATOMIC_NUMBER{29.};
    constexpr auto COPPER_DENSITY{8.92 * g / cm3};
    constexpr auto COPPER_MOLAR_MASS{63.55 * g / mole};

    auto *copper = new G4Material("Copper", COPPER_ATOMIC_NUMBER, COPPER_MOLAR_MASS, COPPER_DENSITY);

    // // Aluminium (for testing)
    // constexpr auto COPPER_ATOMIC_NUMBER{13.};
    // constexpr auto COPPER_DENSITY{2.7 * g / cm3};
    // constexpr auto COPPER_MOLAR_MASS{26.98 * g / mole};

    // auto *copper = new G4Material("Aluminium", COPPER_ATOMIC_NUMBER, COPPER_MOLAR_MASS, COPPER_DENSITY);

    // Silicon
    constexpr auto SILICON_ATOMIC_NUMBER{14.};
    constexpr auto SILICON_DENSITY{2.33 * g / cm3};
    constexpr auto SILICON_MOLAR_MASS{28.09 * g / mole};

    auto *silicon = new G4Material("Silicon", SILICON_ATOMIC_NUMBER, SILICON_MOLAR_MASS, SILICON_DENSITY);

    // Vacuum

    constexpr auto VACUUM_ATOMIC_NUMBER{1};
    constexpr auto VACUUM_DENSITY{universe_mean_density}; // from PhysicalConstants.h
    constexpr auto VACUUM_MOLAR_MASS{1.01 * g / mole};
    constexpr auto VACUUM_PRESSURE{3.e-18 * pascal};
    constexpr auto VACUUM_TEMPERATURE{2.73 * kelvin};

    vacuum = new G4Material("Galactic", VACUUM_ATOMIC_NUMBER, VACUUM_MOLAR_MASS, VACUUM_DENSITY, kStateGas, VACUUM_TEMPERATURE, VACUUM_PRESSURE);

    // Print all the materials defined.
    G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

    // --------- Sizes of the principal geometrical components (solids) ---------

    // world volume
    constexpr auto WORLD_HALF_LENGTH{50.0 * mm};

    worldHalfLength = WORLD_HALF_LENGTH;

    // target: 0.5 cm * 0.5 cm * 3 cm
    constexpr auto TARGET_HALF_LENGTH{2.5 * mm};
    constexpr auto TARGET_HALF_DEPTH{15.0 * mm};

    targetHalfLength = TARGET_HALF_LENGTH;
    targetHalfDepth = TARGET_HALF_DEPTH;

    // tracker (CCD): 4 cm * 450 mu_m * 4 cm
    constexpr auto TRACKER_HALF_LENGTH{20.0 * mm};
    constexpr auto TRACKER_HALF_DEPTH{0.225 * mm};

    trackerHalfLength = TRACKER_HALF_LENGTH;
    trackerHalfDepth = TRACKER_HALF_DEPTH;

    //--------- positions of the principal geometrical components --------

    // target position: (0 cm, 0 cm, 0 cm)
    constexpr auto TARGET_INITIAL_POSITION_X{0.0 * cm};
    constexpr auto TARGET_INITIAL_POSITION_Y{0.0 * cm};
    constexpr auto TARGET_INITIAL_POSITION_Z{0.0 * cm};

    targetPositionX = TARGET_INITIAL_POSITION_X;
    targetPositionY = TARGET_INITIAL_POSITION_Y;
    targetPositionZ = TARGET_INITIAL_POSITION_Z;

    // tracker position: (0 cm, 2 cm, 0 cm)
    constexpr auto TRACKER_INITIAL_POSITION_X{0.0 * cm};
    constexpr auto TRACKER_INITIAL_POSITION_Y{2.0 * cm};
    constexpr auto TRACKER_INITIAL_POSITION_Z{0.0 * cm};

    trackerPositionX = TRACKER_INITIAL_POSITION_X;
    trackerPositionY = TRACKER_INITIAL_POSITION_Y;
    trackerPositionZ = TRACKER_INITIAL_POSITION_Z;
   
    // ------------
    // - Material -
    // ------------

    // worldMaterial = air;
    worldMaterial = vacuum;
    targetMaterial = copper;
    trackerMaterial = silicon;

    //--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

    // ---------
    // - World -
    // ---------

    worldSolid = new G4Box("world", worldHalfLength, worldHalfLength, worldHalfLength); // cube
    worldLogicalVolume = new G4LogicalVolume(worldSolid, worldMaterial, "World", nullptr, nullptr, nullptr);

    // Place the world physical volume unrotated at (0, 0, 0).
    worldPhysicalVolume = new G4PVPlacement(nullptr, // no rotation
        G4ThreeVector(), // at (0, 0, 0)
        worldLogicalVolume, // its logical volume
        "World", // its name
        nullptr, // its mother volume
        false, // no boolean operations
        0); // copy number

    // ----------
    // - Target -
    // ----------

    auto targetPosition = G4ThreeVector(targetPositionX, targetPositionY, targetPositionZ);

    targetSolid = new G4Box("target", targetHalfLength, targetHalfLength, targetHalfDepth);
    targetLogicalVolume = new G4LogicalVolume(targetSolid, targetMaterial, "Target", nullptr, nullptr, nullptr);
    targetPhysicalVolume = new G4PVPlacement(nullptr, // no rotation
        targetPosition, // at (x, y, z)
        targetLogicalVolume, // its logical volume
        "Target", // its name
        worldLogicalVolume, // its mother volume
        false, // no boolean operations
        0); // copy number

    // -----------
    // - Tracker -
    // -----------

    auto trackerPosition = G4ThreeVector(trackerPositionX, trackerPositionY, trackerPositionZ);

    trackerSolid = new G4Box("tracker", trackerHalfLength, trackerHalfDepth, trackerHalfLength);
    trackerLogicalVolume = new G4LogicalVolume(trackerSolid, trackerMaterial, "Tracker", nullptr, nullptr, nullptr);
    trackerPhysicalVolume = new G4PVPlacement(nullptr, // no rotation
        trackerPosition, // at (x, y, z)
        trackerLogicalVolume, // its logical volume
        "Tracker", // its name
        worldLogicalVolume, // its mother volume
        false, // no boolean operations
        0); // copy number

    //-----------------------
    // - Sensitive detector -
    //-----------------------

    constexpr auto TRACKER_SENSITIVE_DETECTOR_NAME{"eRosita/TrackerChamberSD"};
    auto *trackerSensitiveDetector = new eRositaTrackerSD(TRACKER_SENSITIVE_DETECTOR_NAME);
    trackerLogicalVolume->SetSensitiveDetector(trackerSensitiveDetector);

    //------------------
    // - Visualization -
    //------------------

    // make world volume invisible
    worldVisualizationStyle = new G4VisAttributes();
    worldVisualizationStyle->SetVisibility(false);
    worldLogicalVolume->SetVisAttributes(worldVisualizationStyle);

    // render the target in reddish color
    targetVisualizationStyle = new G4VisAttributes();
    // targetVisualizationStyle->SetColor(G4Color(1.0, 0.3, 0.3)); // reddish
    targetVisualizationStyle->SetColor(G4Color(1.0, 1.0, 1.0)); // black
    targetLogicalVolume->SetVisAttributes(targetVisualizationStyle);

    // render the tracker in blueish color
    trackerVisualizationStyle = new G4VisAttributes();
    // trackerVisualizationStyle->SetColor(G4Color(0.3, 0.3, 1.0)); // blueish
    trackerVisualizationStyle->SetColor(G4Color(1.0, 1.0, 1.0)); // black
    trackerLogicalVolume->SetVisAttributes(trackerVisualizationStyle);

    return worldPhysicalVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaDetectorConstruction::ConstructSDandField()
{
    G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
    auto *sensitiveDetectorManager = G4SDManager::GetSDMpointer();

    constexpr auto TRACKER_SENSITIVE_DETECTOR_NAME{"eRosita/TrackerChamberSD"};
    auto *trackerSensitiveDetector = new eRositaTrackerSD(TRACKER_SENSITIVE_DETECTOR_NAME);
    sensitiveDetectorManager->AddNewDetector(trackerSensitiveDetector);
    trackerLogicalVolume->SetSensitiveDetector(trackerSensitiveDetector);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/*
void eRositaDetectorConstruction::SetTargetMaterial(G4String materialName)
{
    // search the material by name
    auto *material = G4Material::GetMaterial(materialName);
    if (material != nullptr) {
        targetMaterial = material;
        targetLogicalVolume->SetMaterial(material);
        G4cout << "\n----> The target is " << targetFullLength / cm << " cm long and made of " << materialName << G4endl;
    }
}
*/
