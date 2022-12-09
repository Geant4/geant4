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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDetectorConstruction  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#include "GammaRayTelAnticoincidenceSD.hh"
#include "GammaRayTelCalorimeterSD.hh"
#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelDetectorMessenger.hh"
#include "GammaRayTelTrackerSD.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4UImanager.hh"
#include "G4UniformMagField.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreadLocal G4GlobalMagFieldMessenger *GammaRayTelDetectorConstruction::fMagFieldMessenger = nullptr;

GammaRayTelDetectorConstruction::GammaRayTelDetectorConstruction()
{
	// Initialize thread-local sensitive detectors
	trackerSD.Put(nullptr);
	calorimeterSD.Put(nullptr);
	anticoincidenceSD.Put(nullptr);

	ComputePayloadParameters();

	// create commands for interactive definition of the payload
	detectorMessenger = new GammaRayTelDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDetectorConstruction::~GammaRayTelDetectorConstruction() {
	delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto GammaRayTelDetectorConstruction::Construct() -> G4VPhysicalVolume* {
	DefineMaterials();
	return ConstructPayload();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::DefineMaterials() {
	G4String name;
	G4String symbol;

	G4int numberOfAtoms;
	G4int numberOfComponents;

	//
	// define Elements
	//

	constexpr auto HYDROGEN_ATOMIC_NUMBER{1.};
	constexpr auto HYDROGEN_MOLAR_MASS{1.01 * g / mole};
	auto *hydrogen = new G4Element(name = "Hydrogen", symbol = "H", HYDROGEN_ATOMIC_NUMBER, HYDROGEN_MOLAR_MASS);

	constexpr auto CARBON_ATOMIC_NUMBER{6.};
    constexpr auto CARBON_MOLAR_MASS{12.01 * g / mole};
	auto *carbon = new G4Element(name = "Carbon", symbol = "C", CARBON_ATOMIC_NUMBER, CARBON_MOLAR_MASS);

	constexpr auto NITROGEN_ATOMIC_NUMBER{7.};
    constexpr auto NITROGEN_MOLAR_MASS{14.006 * g / mole};
	auto *nitrogen = new G4Element(name = "Nitrogen", symbol = "N", NITROGEN_ATOMIC_NUMBER, NITROGEN_MOLAR_MASS);

	constexpr auto OXYGEN_ATOMIC_NUMBER{8.};
	constexpr auto OXYGEN_MOLAR_MASS{15.99 * g / mole};
	auto *oxygen = new G4Element(name = "Oxygen", symbol = "O", OXYGEN_ATOMIC_NUMBER, OXYGEN_MOLAR_MASS);

	constexpr auto ALUMINIUM_ATOMIC_NUMBER{13.};
    constexpr auto ALUMINIUM_MOLAR_MASS{26.98 * g / mole};
	auto *aluminium = new G4Element(name = "Aluminum", symbol = "Al", ALUMINIUM_ATOMIC_NUMBER, ALUMINIUM_MOLAR_MASS);

    constexpr auto SILICON_ATOMIC_NUMBER{14.};
    constexpr auto SILICON_MOLAR_MASS{28.09 * g / mole};
	auto *silicon = new G4Element(name = "Silicon", symbol = "Si", SILICON_ATOMIC_NUMBER, SILICON_MOLAR_MASS);

    constexpr auto IRON_ATOMIC_NUMBER{26.};
    constexpr auto IRON_MOLAR_MASS{55.845 * g / mole};
	auto *iron = new G4Element(name = "Iron", symbol = "Fe", IRON_ATOMIC_NUMBER, IRON_MOLAR_MASS);

    constexpr auto IODINE_ATOMIC_NUMBER{53.};
    constexpr auto IODINE_MOLAR_MASS{126.904 * g / mole};
	auto *iodine = new G4Element(name = "Iodine", symbol = "I", IODINE_ATOMIC_NUMBER, IODINE_MOLAR_MASS);

    constexpr auto CESIUM_ATOMIC_NUMBER{55};
    constexpr auto CESIUM_MOLAR_MASS{132.905 * g / mole};
	auto *cesium = new G4Element(name = "Cesium", symbol = "Cs", CESIUM_ATOMIC_NUMBER, CESIUM_MOLAR_MASS);

    constexpr auto LEAD_ATOMIC_NUMBER{82};
    constexpr auto LEAD_MOLAR_MASS{207.19 * g / mole};
	auto *lead = new G4Element(name = "Lead", symbol = "Pb", LEAD_ATOMIC_NUMBER, LEAD_MOLAR_MASS);

	//
	// define simple materials
	//

    constexpr auto TUNGSTEN_ATOMIC_NUMBER{74.};
    constexpr auto TUNGSTEN_DENSITY{19.3 * g / cm3};
    constexpr auto TUNGSTEN_MOLAR_MASS{183.84 * g / mole};
	auto *tungsten = new G4Material(name = "Tungsten", TUNGSTEN_ATOMIC_NUMBER, TUNGSTEN_MOLAR_MASS, TUNGSTEN_DENSITY);

	//
	// Define a material from elements.
	// Case 1: chemical molecule
	//

	constexpr auto SCINTILLATOR_MATERIAL_DENSITY{1.032 * g / cm3};
	auto *scintillatorMaterial = new G4Material(name = "Scintillator", SCINTILLATOR_MATERIAL_DENSITY, numberOfComponents = 2);
	scintillatorMaterial->AddElement(carbon, numberOfAtoms = 9);
	scintillatorMaterial->AddElement(hydrogen, numberOfAtoms = 10);

	constexpr auto CESIUM_IODIDE_DENSITY{4.53 * g / cm3};
	auto *cesiumIodide = new G4Material(name = "CesiumIodide", CESIUM_IODIDE_DENSITY, numberOfComponents = 2);
	cesiumIodide->AddElement(cesium, numberOfAtoms = 5);
	cesiumIodide->AddElement(iodine, numberOfAtoms = 5);

	//
	// Define a material from elements.
	// Case 2: mixture by fractional mass
	//

	constexpr auto AIR_DENSITY{1.290 * mg / cm3};
	constexpr auto AIR_NITROGEN_MASS_FRACTION{0.7};
	constexpr auto AIR_OXYGEN_MASS_FRACTION{0.3};

	auto *air = new G4Material(name = "Air", AIR_DENSITY, numberOfComponents = 2);
	air->AddElement(nitrogen, AIR_NITROGEN_MASS_FRACTION);
	air->AddElement(oxygen, AIR_OXYGEN_MASS_FRACTION);

	constexpr auto ALUMINIUM_DENSITY{2.700 * g / cm3};
	constexpr auto ALUMINIUM_MASS_FRACTION{1.};
	auto *Al = new G4Material(name = "Aluminum", ALUMINIUM_DENSITY, numberOfComponents = 1);
	Al->AddElement(aluminium, ALUMINIUM_MASS_FRACTION);

	constexpr auto SILICON_DENSITY{2.333 * g / cm3};
	constexpr auto SILICON_MASS_FRACTION{1.};
	auto *Si = new G4Material(name = "Silicon", SILICON_DENSITY, numberOfComponents = 1);
	Si->AddElement(silicon, SILICON_MASS_FRACTION);

	constexpr auto IRON_DENSITY{7.87 * g / cm3};
	constexpr auto IRON_MASS_FRACTION{1.};
	auto *Fe = new G4Material(name = "Iron", IRON_DENSITY, numberOfComponents = 1);
	Fe->AddElement(iron, IRON_MASS_FRACTION);

	constexpr auto LEAD_DENSITY{11.35 * g / cm3};
	constexpr auto LEAD_MASS_FRACTION{1.};
	auto *Pb = new G4Material(name = "Lead", LEAD_DENSITY, numberOfComponents = 1);
	Pb->AddElement(lead, LEAD_MASS_FRACTION);

	//
	// examples of vacuum
	//
	constexpr auto VACUUM_ATOMIC_NUMBER{1.};
	constexpr auto VACUUM_DENSITY{universe_mean_density}; // from PhysicalConstants.h
	constexpr auto VACUUM_MOLAR_MASS{1.01 * g / mole};
    constexpr auto VACUUM_PRESSURE{3.e-18 * pascal};
    constexpr auto VACUUM_TEMPERATURE{2.73 * kelvin};
	auto *vacuum = new G4Material(name = "Galactic", VACUUM_ATOMIC_NUMBER, VACUUM_MOLAR_MASS, VACUUM_DENSITY, kStateGas, VACUUM_TEMPERATURE, VACUUM_PRESSURE);

	constexpr auto BEAM_DENSITY{1.e-5 * g / cm3};
	constexpr auto BEAM_MASS_FRACTION{1.};
	constexpr auto BEAM_PRESSURE{2.e-2 * bar};
	constexpr auto BEAM_TEMPERATURE{STP_Temperature}; // from PhysicalConstants.h
	auto *beam = new G4Material(name = "Beam", BEAM_DENSITY, numberOfComponents = 1, kStateGas, BEAM_TEMPERATURE, BEAM_PRESSURE);
	beam->AddMaterial(air, BEAM_MASS_FRACTION);

	G4cout << *(G4Material::GetMaterialTable()) << G4endl;

	// default materials of the payload

    defaultMaterial = vacuum;

	converterMaterial = tungsten;
	acdMaterial = scintillatorMaterial; // anticoincidence (ACD)
	calMaterial = cesiumIodide; // calorimeter (CAL)
	tkrMaterial = Si; // tracker (TKR)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto GammaRayTelDetectorConstruction::ConstructPayload() -> G4VPhysicalVolume* {
	// complete the payload parameters definition
	ComputePayloadParameters();

	//
	// World
	//

	solidWorld = new G4Box("World", worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2);
	logicWorld = new G4LogicalVolume(solidWorld, defaultMaterial, "World");
	physiWorld = new G4PVPlacement(nullptr, G4ThreeVector(), "World", logicWorld, nullptr, false, 0);

	//
	// Payload
	//

	solidPayload = new G4Box("Payload", payloadSizeXY / 2, payloadSizeXY / 2, payloadSizeZ / 2);
	logicPayload = new G4LogicalVolume(solidPayload, defaultMaterial, "Payload");
	physiPayload = new G4PVPlacement(nullptr, G4ThreeVector(), "Payload", logicPayload, physiWorld, false, 0);

	//                                 
	// Calorimeter (CAL)
	//

	solidCAL = new G4Box("CAL", calSizeXY / 2, calSizeXY / 2, calSizeZ / 2);
	logicCAL = new G4LogicalVolume(solidCAL, defaultMaterial, "CAL");
	physiCAL = new G4PVPlacement(nullptr,
	        G4ThreeVector(0, 0,
	                -payloadSizeZ / 2 + calSizeZ / 2),
            "CAL", logicCAL, physiPayload, false, 0);

	//                                 
	// Tracker (TKR)
	//

	solidTKR = new G4Box("TKR", tkrSizeXY / 2, tkrSizeXY / 2, tkrSizeZ / 2);
	logicTKR = new G4LogicalVolume(solidTKR, defaultMaterial, "TKR");
	physiTKR = new G4PVPlacement(nullptr,
	        G4ThreeVector(0, 0,
	                -payloadSizeZ / 2 + calSizeZ + calTKRDistance
	                    + tkrSizeZ / 2),
            "TKR", logicTKR, physiPayload, false, 0);

	//                               
	// Anticoincidence, Top Side (ACT)
	//

	solidACT = new G4Box("ACT", actSizeXY / 2, actSizeXY / 2, actSizeZ / 2);
	logicACT = new G4LogicalVolume(solidACT, acdMaterial, "ACT");
	physiACT = new G4PVPlacement(nullptr,
			G4ThreeVector(0, 0,
					-payloadSizeZ / 2 + calSizeZ + calTKRDistance + tkrSizeZ
							+ acdTKRDistance + actSizeZ / 2),
            "ACT", logicACT, physiPayload, false, 0);

	//                               
	// Anticoincidence, Lateral Side (ACL)
	//

	solidACL1 = new G4Box("ACL1", acl1SizeX / 2, acl1SizeY / 2, acl1SizeZ / 2);
	logicACL1 = new G4LogicalVolume(solidACL1, acdMaterial, "ACL");

	physiACL1 = new G4PVPlacement(nullptr,
			G4ThreeVector(-payloadSizeXY / 2 + acl1SizeX / 2,
					-payloadSizeXY / 2 + acl1SizeY / 2,
					-payloadSizeZ / 2 + acl1SizeZ / 2),
            "ACL1", logicACL1, physiPayload, false, 0);

	physiACL1 = new G4PVPlacement(nullptr,
			G4ThreeVector(payloadSizeXY / 2 - acl1SizeX / 2,
					payloadSizeXY / 2 - acl1SizeY / 2,
					-payloadSizeZ / 2 + acl1SizeZ / 2),
            "ACL1", logicACL1, physiPayload, false, 1);

	solidACL2 = new G4Box("ACL2", acl2SizeX / 2, acl2SizeY / 2, acl2SizeZ / 2);
	logicACL2 = new G4LogicalVolume(solidACL2, acdMaterial, "ACL2");

	physiACL2 = new G4PVPlacement(nullptr,
			G4ThreeVector(-payloadSizeXY / 2 + acl2SizeX / 2,
					payloadSizeXY / 2 - acl2SizeY / 2,
					-payloadSizeZ / 2 + acl2SizeZ / 2),
            "ACL2", logicACL2, physiPayload, false, 0);

	physiACL2 = new G4PVPlacement(nullptr,
			G4ThreeVector(payloadSizeXY / 2 - acl2SizeX / 2,
					-payloadSizeXY / 2 + acl2SizeY / 2,
					-payloadSizeZ / 2 + acl2SizeZ / 2),
            "ACL2", logicACL2, physiPayload, false, 1);

	//
	// Tracker Structure (Plane + Converter + TKRDetectorX + TKRDetectorY)
	//

	solidPlane = new G4Box("Plane", tkrSizeXY / 2, tkrSizeXY / 2, tkrSupportThickness / 2);
	logicPlane = new G4LogicalVolume(solidPlane, defaultMaterial, "Plane");

	solidTKRDetectorY = new G4Box("TKRDetectorY", tkrSizeXY / 2, tkrSizeXY / 2, tkrSiliconThickness / 2);
	logicTKRDetectorY = new G4LogicalVolume(solidTKRDetectorY, tkrMaterial, "TKRDetector Y");

	solidTKRDetectorX = new G4Box("TKRDetectorX", tkrSizeXY / 2, tkrSizeXY / 2, tkrSiliconThickness / 2);
	logicTKRDetectorX = new G4LogicalVolume(solidTKRDetectorX, tkrMaterial, "TKRDetector X");

	solidConverter = new G4Box("Converter", tkrSizeXY / 2, tkrSizeXY / 2, converterThickness / 2);
	logicConverter = new G4LogicalVolume(solidConverter, converterMaterial, "Converter");

	G4int i = 0;

	for (i = 0; i < numberOfTKRLayers; i++) {
		physiTKRDetectorY = new G4PVPlacement(nullptr,
            G4ThreeVector(0., 0.,
                -tkrSizeZ / 2 + tkrSiliconThickness / 2
                    + (i) * tkrLayerDistance),
            "TKRDetectorY", logicTKRDetectorY, physiTKR, false, i);

		physiTKRDetectorX = new G4PVPlacement(nullptr,
            G4ThreeVector(0., 0.,
                -tkrSizeZ / 2 + tkrSiliconThickness / 2
                    + tkrViewsDistance + tkrSiliconThickness
                    + (i) * tkrLayerDistance),
            "TKRDetectorX", logicTKRDetectorX, physiTKR, false, i);

		physiConverter = new G4PVPlacement(nullptr,
            G4ThreeVector(0., 0.,
                -tkrSizeZ / 2 + 2 * tkrSiliconThickness
                    + tkrViewsDistance + converterThickness / 2
                    + (i) * tkrLayerDistance),
            "Converter", logicConverter, physiTKR, false, i);

		physiPlane = new G4PVPlacement(nullptr,
            G4ThreeVector(0., 0.,
                -tkrSizeZ / 2 + 2 * tkrSiliconThickness
                    + tkrViewsDistance + converterThickness
                    + tkrSupportThickness / 2
                    + (i) * tkrLayerDistance),
            "Plane", logicPlane, physiTKR, false, i);
	}

	auto *solidTKRActiveTileX = new G4Box("Active Tile X", tkrActiveTileXY / 2, tkrActiveTileXY / 2, tkrActiveTileZ / 2);
	auto *solidTKRActiveTileY = new G4Box("Active Tile Y", tkrActiveTileXY / 2, tkrActiveTileXY / 2, tkrActiveTileZ / 2);

	auto *logicTKRActiveTileX = new G4LogicalVolume(solidTKRActiveTileX, tkrMaterial, "Active Tile X", nullptr, nullptr, nullptr);
	auto *logicTKRActiveTileY = new G4LogicalVolume(solidTKRActiveTileY, tkrMaterial, "Active Tile Y", nullptr, nullptr, nullptr);

	G4int j = 0;
	G4int k = 0;

	G4double x = 0.;
	G4double y = 0.;
	G4double z = 0.;

	for (i = 0; i < numberOfTKRTiles; i++) {
		for (j = 0; j < numberOfTKRTiles; j++) {
			k = i * numberOfTKRTiles + j;

			x = -tkrSizeXY / 2 + tilesSeparation + siliconGuardRing
                    + tkrActiveTileXY / 2
                    + (i)
                        * ((2 * siliconGuardRing) + tilesSeparation
                            + tkrActiveTileXY);

			y = -tkrSizeXY / 2 + tilesSeparation + siliconGuardRing
					+ tkrActiveTileXY / 2
					+ (j)
                        * ((2 * siliconGuardRing) + tilesSeparation
                            + tkrActiveTileXY);
			z = 0.;

			new G4PVPlacement(nullptr, G4ThreeVector(x, y, z), logicTKRActiveTileY, "Active Tile Y", logicTKRDetectorY, false, k);

			x = -tkrSizeXY / 2 + tilesSeparation + siliconGuardRing
					+ tkrActiveTileXY / 2
					+ (j)
                        * ((2 * siliconGuardRing) + tilesSeparation
                            + tkrActiveTileXY);

			y = -tkrSizeXY / 2 + tilesSeparation + siliconGuardRing
					+ tkrActiveTileXY / 2
					+ (i)
                        * ((2 * siliconGuardRing) + tilesSeparation
                            + tkrActiveTileXY);
			z = 0.;

			new G4PVPlacement(nullptr, G4ThreeVector(x, y, z), logicTKRActiveTileX, "Active Tile X", logicTKRDetectorX, false, k);
		}
	}

	// Strips

	// Silicon Strips 

	/*
    G4double tkrXStripX{0.};
    G4double tkrYStripY{0.};
    G4double tkrYStripX{0.};
    G4double tkrXStripY{0.};
	*/

	tkrXStripX = tkrYStripY = tkrSiliconPitch;
	tkrYStripX = tkrXStripY = tkrActiveTileXY;
	tkrZStrip = tkrSiliconThickness;

	auto *solidTKRStripX = new G4Box("Strip X", tkrXStripX / 2, tkrYStripX / 2, tkrZStrip / 2);
	logicTKRStripX = new G4LogicalVolume(solidTKRStripX, tkrMaterial, "Strip X", nullptr, nullptr, nullptr);

	auto *solidTKRStripY = new G4Box("Strip Y", tkrXStripY / 2, tkrYStripY / 2, tkrZStrip / 2);
	logicTKRStripY = new G4LogicalVolume(solidTKRStripY, tkrMaterial, "Strip Y", nullptr, nullptr, nullptr);

	for (i = 0; i < numberOfTKRStrips; i++) {
		new G4PVPlacement(nullptr,
            G4ThreeVector(
                -tkrActiveTileXY / 2 + tkrSiliconPitch / 2
                    + (i) * tkrSiliconPitch, 0., 0.),
            logicTKRStripX, "Strip X", logicTKRActiveTileX, false, i);

		new G4PVPlacement(nullptr,
            G4ThreeVector(0.,
                -tkrActiveTileXY / 2 + tkrSiliconPitch / 2
                    + (i) * tkrSiliconPitch, 0.),
            logicTKRStripY, "Strip Y", logicTKRActiveTileY, false, i);
	}

	//
	// Calorimeter Structure (CALLayerX + CALLayerY)
	//

	solidCALLayerX = new G4Box("CALLayerX", calSizeXY / 2, calSizeXY / 2, calBarThickness / 2);
	logicCALLayerX = new G4LogicalVolume(solidCALLayerX, calMaterial, "CALLayerX");

	solidCALLayerY = new G4Box("CALLayerY", calSizeXY / 2, calSizeXY / 2, calBarThickness / 2);
	logicCALLayerY = new G4LogicalVolume(solidCALLayerY, calMaterial, "CALLayerY");

	for (i = 0; i < numberOfCALLayers; i++) {
		physiCALLayerY = new G4PVPlacement(nullptr,
            G4ThreeVector(0, 0,
                -calSizeZ / 2 + calBarThickness / 2
                    + (i) * 2 * calBarThickness),
            "CALLayerY", logicCALLayerY, physiCAL, false, i);

		physiCALLayerX = new G4PVPlacement(nullptr,
            G4ThreeVector(0, 0,
                -calSizeZ / 2 + calBarThickness / 2 + calBarThickness
                    + (i) * 2 * calBarThickness),
            "CALLayerX", logicCALLayerX, physiCAL, false, i);
	}

	//
	// Calorimeter Structure (CALDetectorX + CALDetectorY)
	//

	solidCALDetectorX = new G4Box("CALDetectorX", calBarX / 2, calBarY / 2, calBarThickness / 2);
	logicCALDetectorX = new G4LogicalVolume(solidCALDetectorX, calMaterial, "CALDetectorX");

	solidCALDetectorY = new G4Box("CALDetectorY", calBarY / 2, calBarX / 2, calBarThickness / 2);
	logicCALDetectorY = new G4LogicalVolume(solidCALDetectorY, calMaterial, "CALDetectorY");

	for (i = 0; i < numberOfCALBars; i++) {
		physiCALDetectorY = new G4PVPlacement(nullptr,
            G4ThreeVector(-calSizeXY / 2 + calBarY / 2 + (i) * calBarY, 0, 0),
            logicCALDetectorY, "CALDetectorY", logicCALLayerY, false, i);

		physiCALDetectorX = new G4PVPlacement(nullptr,
            G4ThreeVector(0, -calSizeXY / 2 + calBarY / 2 + (i) * calBarY, 0),
            logicCALDetectorX, "CALDetectorX", logicCALLayerX, false, i);
	}

/*
     // Cuts by Region

    G4String regionName[] = {"Calorimeter", "Tracker"};

    if (calorimeterCutRegion) {
        delete calorimeterCutRegion;
    }
    calorimeterCutRegion = new G4Region(regionName[0]);
    logicCAL->SetRegion(calorimeterCutRegion);
    calorimeterCutRegion->AddRootLogicalVolume(logicCAL);

    if (trackerCutRegion != nullptr) {
        delete trackerCutRegion;
    }
    trackerCutRegion = new G4Region(regionName[1]);
    logicTKR->SetRegion(trackerCutRegion);
    trackerCutRegion->AddRootLogicalVolume(logicTKR);
*/

	//                                        
	// Visualization attributes
	//
	// Invisible Volume
	logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicPayload->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicTKR->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicTKRActiveTileX->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicTKRActiveTileY->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicPlane->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicConverter->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicCAL->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicCALLayerX->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicCALLayerY->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicTKRStripX->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicTKRStripY->SetVisAttributes(G4VisAttributes::GetInvisible());

	// Some visualization styles

	auto *visualizationStyle1 = new G4VisAttributes(G4Colour(0.3, 0.8, 0.1));
	visualizationStyle1->SetVisibility(true);
	visualizationStyle1->SetForceSolid(TRUE);

	auto *visualizationStyle2 = new G4VisAttributes(G4Colour(0.2, 0.3, 0.8));
	visualizationStyle2->SetVisibility(true);
	visualizationStyle2->SetForceSolid(FALSE);

	auto *visualizationStyle3 = new G4VisAttributes(G4Colour(0.8, 0.2, 0.3));
	visualizationStyle3->SetVisibility(true);
	visualizationStyle3->SetForceWireframe(TRUE);

	// Visible Volumes

	logicCALDetectorX->SetVisAttributes(visualizationStyle1);
	logicCALDetectorY->SetVisAttributes(visualizationStyle1);
	logicTKRDetectorX->SetVisAttributes(visualizationStyle2);
	logicTKRDetectorY->SetVisAttributes(visualizationStyle2);
	logicACT->SetVisAttributes(visualizationStyle3);
	logicACL1->SetVisAttributes(visualizationStyle3);
	logicACL2->SetVisAttributes(visualizationStyle3);

	//
	// always return the physical World
	//
	PrintPayloadParameters();

	return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::ConstructSDandField() {
	//
	// Sensitive detector: Tracker
	//                                
	if (trackerSD.Get() == nullptr) {
	    constexpr auto TRACKER_SENSITIVE_DETECTOR_NAME = "TrackerSD";
		auto *sensitiveDetector = new GammaRayTelTrackerSD(TRACKER_SENSITIVE_DETECTOR_NAME);
		trackerSD.Put(sensitiveDetector);
	}

	G4SDManager::GetSDMpointer()->AddNewDetector(trackerSD.Get());

	// Flags the strips as sensitive .
	if (logicTKRStripX != nullptr) {
		SetSensitiveDetector(logicTKRStripX, trackerSD.Get()); // ActiveStripX
	}
	if (logicTKRStripY != nullptr) {
		SetSensitiveDetector(logicTKRStripY, trackerSD.Get()); // ActiveStripY
	}

	//
	// Sensitive detector: Calorimeter
	// 
	if (calorimeterSD.Get() == nullptr) {
	    constexpr auto CALORIMETER_SENSITIVE_DETECTOR_NAME = "CalorimeterSD";
		auto *sensitiveDetector = new GammaRayTelCalorimeterSD(CALORIMETER_SENSITIVE_DETECTOR_NAME);
		calorimeterSD.Put(sensitiveDetector);
	}

	G4SDManager::GetSDMpointer()->AddNewDetector(calorimeterSD.Get());
	if (logicCALDetectorX != nullptr) {
		SetSensitiveDetector(logicCALDetectorX, calorimeterSD.Get()); // CAL BarX
	}
	if (logicCALDetectorY != nullptr) {
		SetSensitiveDetector(logicCALDetectorY, calorimeterSD.Get()); // CAL BarY
	}

	//
	// Sensitive detector: Anticoincidence
	//
	if (anticoincidenceSD.Get() == nullptr) {
	    constexpr auto ANTICOINCIDENCE_SENSITIVE_DETECTOR_NAME = "AnticoincidenceSD";
		auto *sensitiveDetector = new GammaRayTelAnticoincidenceSD(ANTICOINCIDENCE_SENSITIVE_DETECTOR_NAME);
		anticoincidenceSD.Put(sensitiveDetector);
	}

	G4SDManager::GetSDMpointer()->AddNewDetector(anticoincidenceSD.Get());
	if (logicACT != nullptr) {
		SetSensitiveDetector(logicACT, anticoincidenceSD.Get()); // ACD top
	}
	if (logicACL1 != nullptr) {
		SetSensitiveDetector(logicACL1, anticoincidenceSD.Get()); // ACD lateral side
	}
	if (logicACL2 != nullptr) {
		SetSensitiveDetector(logicACL2, anticoincidenceSD.Get()); // ACD lateral side
	}

	// Create global magnetic field messenger.
	// Uniform magnetic field is then created automatically if the field value is not zero.
	auto fieldValue = G4ThreeVector();
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	fMagFieldMessenger->SetVerboseLevel(1);

	// Register the field messenger for deleting
	G4AutoDelete::Register (fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::PrintPayloadParameters() {
	G4cout
	    << "\n---------------------------------------------------------------------------------\n"
		<< "---> The tracker is composed by " << numberOfTKRLayers << " layers"
		<< ", each made of " << converterMaterial->GetName()
		<< " and " << converterThickness / mm << " mm long."
		<< "\n---------------------------------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetConverterMaterial(G4String materialChoice) {
	// search the material by its name   
	G4Material *material = G4Material::GetMaterial(materialChoice);
	if (material != nullptr) {
		converterMaterial = material;
		logicConverter->SetMaterial(material);
		PrintPayloadParameters();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetConverterThickness(G4double value) {
	converterThickness = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetTKRSiliconThickness(G4double value) {
	tkrSiliconThickness = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetTKRSiliconPitch(G4double value) {
	tkrSiliconPitch = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetTKRTileSizeXY(G4double value) {
	tkrSiliconTileXY = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetNbOfTKRLayers(G4int value) {
	numberOfTKRLayers = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetNbOfTKRTiles(G4int value) {
	numberOfTKRTiles = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetTKRLayerDistance(G4double value) {
    tkrLayerDistance = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetTKRViewsDistance(G4double value) {
    tkrViewsDistance = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetNbOfCALLayers(G4int value) {
    numberOfCALLayers = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetNbOfCALBars(G4int value) {
    numberOfCALBars = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetCALBarThickness(G4double value) {
    calBarThickness = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetACDThickness(G4double value) {
    acdThickness = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetMagField(G4double fieldValue) {
    // Just invoke manually the MT-safe command /globalField/setValue instantiated by the GlobalFieldMessenger
    std::stringstream stream;
    stream << "/globalField/setValue 0 0 " << fieldValue / tesla << " tesla";

    G4String command = stream.str();
    G4cout << "Going to execute: " << command << G4endl;

    auto *uiManager = G4UImanager::GetUIpointer();
    uiManager->ApplyCommand(command);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::UpdateGeometry() {
    // delete payloadSD;
    G4RunManager::GetRunManager()->DefineWorldVolume(ConstructPayload());
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4RegionStore::GetInstance()->UpdateMaterialList(physiWorld);
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}
