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

//GEANT4 - Depth-of-Interaction enabled Positron emission tomography (PET) advanced example 

//Contributors

// Abdella M. Ahmed (1, 2), Andrew Chacon (1, 2), Harley Rutherford (1, 2),
// Hideaki Tashima (3), Go Akamatsu (3), Akram Mohammadi (3), Eiji Yoshida (3), Taiga Yamaya (3)
// Susanna Guatelli (2), and Mitra Safavi-Naeini (1, 2)

// (1) Australian Nuclear Science and Technology Organisation, Australia
// (2) University of Wollongong, Australia
// (3) National Institute of Radiological Sciences, Japan



#ifndef doiPETDetectorConstruction_h
#define doiPETDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4Element.hh"
#include "doiPETGlobalParameters.hh"

class G4Material;
class G4VPhysicalVolume;
class G4LogicalVolume;
//
class doiPETAnalysis;
class doiPETDetectorConstructionMessenger;

/// Detector construction class to define materials and geometry.
///
/// Crystals are positioned in Ring, with an appropriate rotation matrix. 
/// Several copies of Ring are placed in the full detector.

class doiPETDetectorConstruction : public G4VUserDetectorConstruction
{
public:
	doiPETDetectorConstruction();
	virtual ~doiPETDetectorConstruction();

public:
	virtual G4VPhysicalVolume* Construct();

	void ConstructPhantom(G4LogicalVolume*);
	void ChangePhantom(G4String choice);
	void SetPhantomPosition(G4ThreeVector);
	void SetPhantomRadius(G4double);
	void SetPhantomLength(G4double);


private:
	void DefineMaterials();
	doiPETDetectorConstructionMessenger* fDetectorMessenger;
	doiPETAnalysis* pAnalysis;

	G4LogicalVolume* phantom_logicalV;
	G4VPhysicalVolume* phantom_physicalV;
	// G4LogicalVolume* gelatin_logicalV;
	// G4VPhysicalVolume* gelatin_physicalV;

	//wolrd
	G4LogicalVolume* world_logicalV;
	G4VPhysicalVolume* world_physicalV;

	//detector block
	G4LogicalVolume* blockDetector_logicalV;
	G4VPhysicalVolume* blockDetector_physicalV;

	//air volume to fill the detector block
	G4LogicalVolume* airBox_logicalV;
	G4VPhysicalVolume* airBox_physicalV;

	//crystals
	G4LogicalVolume* crystal_logicalV;
	G4VPhysicalVolume* crystal_physicalV;


	//water
	G4LogicalVolume* water_logicalV;
	G4VPhysicalVolume* water_physicalV;

	//lung
	G4LogicalVolume* lung_logicalV;
	G4VPhysicalVolume* lung_physicalV;

	//test line phantom
	G4LogicalVolume* test_logicalV;
	G4VPhysicalVolume* test_physicalV;

	G4LogicalVolume* lung_logicalV_PMMA;
	G4VPhysicalVolume* lung_physicalVPMMA;

	//cold regions
	// G4LogicalVolume* coldRegion_logicalV;
	// G4VPhysicalVolume* coldRegion_physicalV;

	//
	//Surrounding PMMA for hot sphere
	G4LogicalVolume* hotSpherePMMA_logicalV;
	G4VPhysicalVolume* hotSpherePMMA_physicalV;

	//hot water phantom (activity is distributed)
	G4LogicalVolume* hotSphereWater_logicalV;
	G4VPhysicalVolume* hotSphereWater_physicalV;

	//surrounding PMMA cold sphere
	G4LogicalVolume* coldSpherePMMA_logicalV;
	G4VPhysicalVolume* coldSpherePMMA_physicalV;

	//cold Water phantom in the cold PMMA sphere 
	G4LogicalVolume* coldSphereWater_logicalV;
	G4VPhysicalVolume* coldSphereWater_physicalV;

	//fillable polyethylene phantom for sensitivity
	G4LogicalVolume* phantomPE_logicalV;
	G4VPhysicalVolume* phantomPE_physicalV;

	//Image quality phantom for small animal NEMA NU-4
	G4LogicalVolume* waterPhantom_logicalV;
	G4VPhysicalVolume* WaterPhantom_physicalV;

	G4LogicalVolume* rod_phantom_logicalV;
	G4VPhysicalVolume* rod_phantom_physicalV;

	G4LogicalVolume* chamberPMMA_logicalV;
	G4VPhysicalVolume* chamberPMMA_physicalV;

	//
	G4LogicalVolume* chamberWater_logicalV;
	G4VPhysicalVolume* chamberWater_physicalV;

	//
	G4LogicalVolume* chamberAir_logicalV;
	G4VPhysicalVolume* chamberAir_physicalV;


	//Dimension of the sphere
	G4double spherePositionX, spherePositionY; // spherePositionZ;
	G4double sphereDiameter;
	G4double distanceFromCenter;
	G4int numberOfSpheres;
	G4double sphereWallThickness;
	G4double zOffsetSpherePhantom;

	G4String PhantomType;

	//materials
	G4Material* air;
	G4Material* pmma;
	G4Material* water;
	G4Material* polyethylene;
	G4Material* polyethylene_NEMA;
	// G4Material* inflatedLung;
	G4Material* polystyrene;
	G4Material* Aluminum;

	//elements for GSO
	G4Element*  O;
	G4Element* Si;
	G4Element* Gd;
	G4Material* GSO;

	G4Material* crystalMaterial;
	//G4Material* phantomMaterial;

	G4bool  fCheckOverlaps;
	G4bool isotopes;


	//size of world
	G4double worldSizeX;
	G4double worldSizeY;
	G4double worldSizeZ;


	//The following is moved to doiPETGlobalParameters.hh
	//G4int numberOfCrystal_DOI;
	//G4int numberOfCrystal_tangential;
	//G4int numberOfCrystal_axial;

	////
	//G4double sizeOfCrystal_DOI; 
	//G4double sizeOfCrystal_tangential;
	//G4double sizeOfCrystal_axial;

	////
	//G4double crystalGap_DOI;
	//G4double crystalGap_tangential;
	//G4double crystalGap_axial;

	G4double sizeOfAirBox_DOI;
	G4double sizeOfAirBox_axial;
	G4double sizeOfAirBox_tangential;


	G4double sizeOfBlockDetector_DOI;
	G4double sizeOfBlockDetector_axial;
	G4double sizeOfBlockDetector_tangential;

	//G4double AluminumCoverThickness;


	//G4int numberOfPETDetector;
	//G4int numberOfRings;


	//G4double scannerRadius;
	G4double thetaDetector; //The azimuthal angle for arranging the detector in the PET ring 
	//G4double ringGap;
	G4int blockIndex;
	// G4int AlCase_Index;
	G4int crystalIndex;

	//detector position
	G4double detectorPositionX;
	G4double detectorPositionY;
	G4double detectorPositionZ;

	//crystal position
	G4double crystalPositionX; 
	G4double crystalPositionY;
	G4double crystalPositionZ;



	G4ThreeVector phantomPosition;

	//
	G4double phantomRadius;
	G4double phantomLength;

	//Phantom dimension for rectangular box (placed for therapy study)
	// G4double phantomSizeX, phantomSizeY, phantomSizeZ;

	//the following is to make the body phantom 
	G4double yOffsetBodyPhantom;
	G4double zOffsetBodyPhantom;
	G4double lengthOfBodyPhantom; //Interior length ( = 180m mm) + wallthickness (= 3mm)
	G4double radiusOfBodyPhantom;
	G4double wallThicknessOfBodyPhantom;
	G4double radiusOfLungPhantom;

	//Test phantom defnition. The phantom has the same as that that of NECR phantom except
	G4double hieghtOfTestPhantom; 
	G4double diameterOfTestPhantom;

	//To the cylindrical phantom to make the body phantom
	G4double radiusOfSmallcyl;
	G4double boxWidth;
	G4double boxHeight;


	//Image quality phantom for small animals
	G4double waterPhantomRadius;
	G4double waterPhantomLength;

	G4double rodPhantomLength;
	G4double rodDiameter;
	G4int numberOfRods;

	//Declare position for the rod phantoms
	G4double rodPositionX, rodPositionY, rodPositionZ;

	//Declare position for cold region chanmbers
	G4double chamberPositionX, chamberPositionY, chamberPositionZ;
	G4double chamberPhantomLength;
	G4double chamberDiameter;
	G4double wallThicknessOfChamber;


};

/////////////////////////////////////////////////////////////////////////////////////////

#endif
