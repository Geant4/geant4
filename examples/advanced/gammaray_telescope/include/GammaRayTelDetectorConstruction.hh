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
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDetectorConstruction  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#ifndef GammaRayTelDetectorConstruction_h
#define GammaRayTelDetectorConstruction_h 1

#include "G4Cache.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class GammaRayTelAnticoincidenceSD;
class GammaRayTelCalorimeterSD;
class GammaRayTelDetectorMessenger;
class GammaRayTelTrackerSD;

class G4Box;
class G4GlobalMagFieldMessenger;
class G4LogicalVolume;
class G4Material;
class G4Region;
class G4UniformMagField;
class G4VPhysicalVolume;

class GammaRayTelDetectorConstruction: public G4VUserDetectorConstruction {
public:
    explicit GammaRayTelDetectorConstruction();

	~GammaRayTelDetectorConstruction() override;

	void SetNbOfTKRLayers(G4int value); // tracker (TKR) number of layers

	void SetTKRTileSizeXY(G4double value);

	void SetNbOfTKRTiles(G4int value); // tracker (TKR) number of tiles

	void SetTKRSiliconThickness(G4double value);

	void SetTKRSiliconPitch(G4double value);

	void SetTKRLayerDistance(G4double value);

	void SetTKRViewsDistance(G4double value);

	void SetConverterMaterial(G4String materialChoice); // tracker (TKR) converter material

	void SetConverterThickness(G4double value); // tracker (TKR) converter thickness

	void SetNbOfCALLayers(G4int value); // calorimeter (CAL) material, length, thickness

	void SetNbOfCALBars(G4int value);

	void SetCALBarThickness(G4double value);

	void SetACDThickness(G4double value); // anticoincidence (ACD) thickness

	void SetMagField(G4double fieldValue); // magnetic field

	auto Construct() -> G4VPhysicalVolume* override;

	void ConstructSDandField() override;

	void PrintPayloadParameters();

	void UpdateGeometry();

	[[nodiscard]]
	auto GetWorldSizeZ() const -> G4double {
		return worldSizeZ;
	}

	[[nodiscard]]
	auto GetWorldSizeXY() const -> G4double {
		return worldSizeXY;
	}

	[[nodiscard]]
	auto GetPayloadSizeZ() const -> G4double {
		return payloadSizeZ;
	}

	[[nodiscard]]
	auto GetPayloadSizeXY() const -> G4double {
		return payloadSizeXY;
	}

	[[nodiscard]]
	auto GetTKRSizeZ() const -> G4double {
		return tkrSizeZ;
	}

	[[nodiscard]]
	auto GetTKRSizeXY() const -> G4double {
		return tkrSizeXY;
	}

	[[nodiscard]]
	auto GetCALSizeZ() const -> G4double {
		return calSizeZ;
	}

	[[nodiscard]]
	auto GetCALTKRDistance() const -> G4double {
		return calTKRDistance;
	}

	[[nodiscard]]
	auto GetTKRSiliconThickness() const -> G4double {
		return tkrSiliconThickness;
	}

	[[nodiscard]]
	auto GetTKRSiliconTileXY() const -> G4double {
		return tkrSiliconTileXY;
	}

	[[nodiscard]]
	auto GetTKRSiliconPitch() const -> G4double {
		return tkrSiliconPitch;
	}

	[[nodiscard]]
	auto GetNbOfTKRLayers() const -> G4int {
		return numberOfTKRLayers;
	}

	auto GetNbOfTKRTiles() const -> G4int {
		return numberOfTKRTiles;
	}

	[[nodiscard]]
	auto GetNbOfTKRStrips() const -> G4int {
		return numberOfTKRStrips;
	}

	[[nodiscard]]
	auto GetTKRLayerDistance() const -> G4double {
		return tkrLayerDistance;
	}

	[[nodiscard]]
	auto GetTKRViewsDistance() const -> G4double {
		return tkrViewsDistance;
	}

	[[nodiscard]]
	auto GetTKRActiveTileXY() const -> G4double {
		return tkrActiveTileXY;
	}

	[[nodiscard]]
	auto GetTKRActiveTileZ() const -> G4double {
		return tkrActiveTileZ;
	}

	[[nodiscard]]
	auto GetSiliconGuardRing() const -> G4double {
		return siliconGuardRing;
	}

	[[nodiscard]]
	auto GetTilesSeparation() const -> G4double {
		return tilesSeparation;
	}

	[[nodiscard]]
	auto GetConverterMaterial() const -> G4Material* {
		return converterMaterial;
	}

	[[nodiscard]]
	auto GetConverterThickness() const -> G4double {
		return converterThickness;
	}

	[[nodiscard]]
	auto GetCALBarThickness() const -> G4double {
		return calBarThickness;
	}

	[[nodiscard]]
	auto GetNbOfCALLayers() const -> G4int {
		return numberOfCALLayers;
	}

	[[nodiscard]]
	auto GetNbOfCALBars() const -> G4int {
		return numberOfCALBars;
	}

	[[nodiscard]]
	auto GetACDThickness() const -> G4double {
		return acdThickness;
	}

	[[nodiscard]]
	auto GetNbOfACDTopTiles() const -> G4int {
		return numberOfACDTopTiles;
	}

	[[nodiscard]]
	auto GetNbOfACDLateralTiles() const -> G4int {
		return numberOfACDLateralTiles;
	}

private:
	G4Material *converterMaterial;
	G4double converterThickness{300. * micrometer};

	// Tracker (TKR)

	G4double tkrSiliconThickness{400. * micrometer};
	G4double tkrSiliconTileXY{9. * cm};
	G4double tkrSiliconPitch{200. * micrometer};

	G4double tkrSizeXY;
	G4double tkrSizeZ;
	G4double tkrLayerDistance{3. * cm};
	G4double tkrViewsDistance{1. * mm};
	G4double tkrSupportThickness;

	G4int numberOfTKRLayers{15};
	G4int numberOfTKRTiles{4};

	// Calorimeter (CAL)

	G4double calBarThickness{1.5 * cm};

	G4int numberOfCALLayers{5};
	G4int numberOfCALBars{12};

	G4double calSizeXY;
	G4double calSizeZ;
	G4double calBarX;
	G4double calBarY;
	G4double calBarZ;

	// Anticoincidence (ACD)

	G4double acdThickness{1. * cm};

	G4double actSizeXY;
	G4double actSizeZ;

	G4double acl1SizeX;
	G4double acl1SizeY;
	G4double acl1SizeZ;

	G4double acl2SizeX;
	G4double acl2SizeY;
	G4double acl2SizeZ;

	G4int numberOfACDTopTiles{1};
	G4int numberOfACDLateralTiles{2};

	G4double tilesSeparation{100. * micrometer};
	G4double acdTKRDistance{5. * cm};
	G4double calTKRDistance{1.5 * cm};
	G4double tkrActiveTileXY;
	G4double tkrActiveTileZ;

	G4double siliconGuardRing{1.5 * mm};
	G4int numberOfTKRStrips;

	G4double tkrXStripX;
	G4double tkrYStripX;
	G4double tkrXStripY;
	G4double tkrYStripY;
	G4double tkrZStrip;

	G4double payloadSizeZ;
	G4double payloadSizeXY;

	G4double worldSizeXY;
	G4double worldSizeZ;

	// Material
    G4Material *defaultMaterial;
    G4Material *calMaterial;
    G4Material *tkrMaterial;
    G4Material *acdMaterial;

	// World
	G4Box *solidWorld{nullptr};
	G4LogicalVolume *logicWorld{nullptr};
	G4VPhysicalVolume *physiWorld{nullptr};

	// Payload
	G4Box *solidPayload{nullptr};
	G4LogicalVolume *logicPayload{nullptr};
	G4VPhysicalVolume *physiPayload{nullptr};

	// Tracker
	G4Box *solidTKR{nullptr};
	G4LogicalVolume *logicTKR{nullptr};
	G4VPhysicalVolume *physiTKR{nullptr};

	// Calorimeter
	G4Box *solidCAL{nullptr};
	G4LogicalVolume *logicCAL{nullptr};
	G4VPhysicalVolume *physiCAL{nullptr};

	// Top Anticoincidence
	G4Box *solidACT{nullptr};
	G4LogicalVolume *logicACT{nullptr};
	G4VPhysicalVolume *physiACT{nullptr};

	// Lateral Anticoincidence
	G4Box *solidACL1{nullptr};
	G4LogicalVolume *logicACL1{nullptr};
	G4VPhysicalVolume *physiACL1{nullptr};

	G4Box *solidACL2{nullptr};
	G4LogicalVolume *logicACL2{nullptr};
	G4VPhysicalVolume *physiACL2{nullptr};

	// Tracker PLANE X
	G4Box *solidTKRDetectorX{nullptr};
	G4LogicalVolume *logicTKRDetectorX{nullptr};
	G4VPhysicalVolume *physiTKRDetectorX{nullptr};

	// Tracker PLANE Y
	G4Box *solidTKRDetectorY{nullptr};
	G4LogicalVolume *logicTKRDetectorY{nullptr};
	G4VPhysicalVolume *physiTKRDetectorY{nullptr};

	// Calorimeter PLANE X
	G4Box *solidCALLayerX{nullptr};
	G4LogicalVolume *logicCALLayerX{nullptr};
	G4VPhysicalVolume *physiCALLayerX{nullptr};

	// Calorimeter PLANE Y
	G4Box *solidCALLayerY{nullptr};
	G4LogicalVolume *logicCALLayerY{nullptr};
	G4VPhysicalVolume *physiCALLayerY{nullptr};

	// Calorimeter DETECTOR X
	G4Box *solidCALDetectorX{nullptr};
	G4LogicalVolume *logicCALDetectorX{nullptr};
	G4VPhysicalVolume *physiCALDetectorX{nullptr};

	// Calorimeter DETECTOR Y
	G4Box *solidCALDetectorY{nullptr};
	G4LogicalVolume *logicCALDetectorY{nullptr};
	G4VPhysicalVolume *physiCALDetectorY{nullptr};

	// Support Plane
	G4Box *solidPlane{nullptr};
	G4LogicalVolume *logicPlane{nullptr};
	G4VPhysicalVolume *physiPlane{nullptr};

	// Converter
	G4Box *solidConverter{nullptr};
	G4LogicalVolume *logicConverter{nullptr};
	G4VPhysicalVolume *physiConverter{nullptr};

	G4LogicalVolume *logicTKRStripX{nullptr};
	G4LogicalVolume *logicTKRStripY{nullptr};

	// magnetic field messenger
	static G4ThreadLocal G4GlobalMagFieldMessenger* fMagFieldMessenger;

	GammaRayTelDetectorMessenger *detectorMessenger;  // pointer to the messenger

	G4Cache<GammaRayTelTrackerSD*> trackerSD; // pointer to the sensitive detector, tracker (TRK)
	G4Cache<GammaRayTelCalorimeterSD*> calorimeterSD; // pointer to the sensitive detector, calorimeter (CAL)
	G4Cache<GammaRayTelAnticoincidenceSD*> anticoincidenceSD; // pointer to the sensitive detector, anticoincidence (ACD)

//    G4Region* trackerCutRegion; // tracker (TKR) cut region
//    G4Region* calorimeterCutRegion; // calorimeter (CAL) cut region

	void ComputePayloadParameters();

	auto ConstructPayload() -> G4VPhysicalVolume*;

	void DefineMaterials();
};

inline auto GammaRayTelDetectorConstruction::ComputePayloadParameters() -> void {
	// Compute derived parameters of the payload

	tkrSupportThickness = tkrLayerDistance - 2 * tkrSiliconThickness - tkrViewsDistance - converterThickness;
	tkrSizeXY = numberOfTKRTiles * tkrSiliconTileXY + (numberOfTKRTiles + 1) * tilesSeparation;
	tkrSizeZ = numberOfTKRLayers * tkrLayerDistance;

	tkrActiveTileXY = tkrSiliconTileXY - 2 * siliconGuardRing;
	tkrActiveTileZ = tkrSiliconThickness;
	numberOfTKRStrips = G4int(tkrActiveTileXY / tkrSiliconPitch);

	siliconGuardRing = tkrActiveTileXY - numberOfTKRStrips * tkrSiliconPitch;
	tkrActiveTileXY = tkrSiliconTileXY - 2 * siliconGuardRing;

	tkrXStripX = tkrYStripY = tkrSiliconPitch;
	tkrYStripX = tkrXStripY = tkrActiveTileXY;
	tkrZStrip = tkrSiliconThickness;

	calSizeXY = tkrSizeXY;
	calSizeZ = 2 * numberOfCALLayers * calBarThickness;

	calBarX = calSizeXY;
	calBarY = calSizeXY / (numberOfCALBars);
	calBarZ = calBarThickness;

	actSizeXY = tkrSizeXY + 2 * acdTKRDistance + 2 * acdThickness;
	actSizeZ = acdThickness;

	acl1SizeX = tkrSizeXY + 2 * acdTKRDistance + acdThickness;
	acl1SizeY = acdThickness;
	acl1SizeZ = tkrSizeZ + calSizeZ + acdTKRDistance + calTKRDistance;

	acl2SizeX = acdThickness;
	acl2SizeY = tkrSizeXY + 2 * acdTKRDistance + acdThickness;
	acl2SizeZ = tkrSizeZ + calSizeZ + acdTKRDistance + calTKRDistance;

	payloadSizeZ = 1.1 * (acl1SizeZ + actSizeZ);
	payloadSizeXY = (actSizeXY);

	worldSizeZ = 1.5 * payloadSizeZ;
	worldSizeXY = 1.5 * payloadSizeXY;
}
#endif
