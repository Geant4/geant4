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

//GEANT4 - Depth-of-Interaction enabled Positron emission tomography (PET) advanced example 

//Authors and contributors

// Author list to be updated, with names of co-authors and contributors from National Institute of Radiological Sciences (NIRS)

// Abdella M. Ahmed (1, 2), Andrew Chacon (1, 2), Harley Rutherford (1, 2),
// Hideaki Tashima (3), Go Akamatsu (3), Akram Mohammadi (3), Eiji Yoshida (3), Taiga Yamaya (3)
// Susanna Guatelli (2), and Mitra Safavi-Naeini (1, 2)

// (1) Australian Nuclear Science and Technology Organisation, Australia
// (2) University of Wollongong, Australia
// (3) National Institute of Radiological Sciences, Japan



//A whole-body PET scanner with depth-of-intercation (DOI) detector is 
//defined in the doiPETDetectorConstruction class and varius types phantoms 
//are also defined based on the NEMA NU-2 2012 standard protocol for
//performance evaluation PET scanners. 
//The scanner's geometrical specifications (like number of crystals, number of rings ...) are defined as a global variable in the doiPETGlobalParameters.hh
//If the phantom is other there the phantoms defined, exception will be thrown. If you want to run the simulation without phantom (for any reason),
//you can comment out ConstructPhantom(world_logicalV) method. 

#include "doiPETDetectorConstruction.hh"
#include "doiPETAnalysis.hh"
#include "doiPETDetectorConstructionMessenger.hh"
#include "doiPETGlobalParameters.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
//
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Sphere.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"


////////////////////////////////////////////////////////////////////////////////////////

doiPETDetectorConstruction::doiPETDetectorConstruction()
	: G4VUserDetectorConstruction(),
	fCheckOverlaps(true)
{
	DefineMaterials();
	fDetectorMessenger = new doiPETDetectorConstructionMessenger(this);

	//Initialize variables. The following variables can be changed via the run.mac file.
	PhantomType = "NEMA_phantom_NECR";
	phantomRadius = 100 * mm;
	phantomLength = 700 * mm;
	phantomPosition = G4ThreeVector(0,0,0);
}

//////////////////////////////////////////////////////////////////////////////////////

doiPETDetectorConstruction::~doiPETDetectorConstruction()
{
	delete fDetectorMessenger;
}

/////////////////////////////////// Define Materials /////////////////////////////////////////////////////

void doiPETDetectorConstruction::DefineMaterials()
{
	G4NistManager* nist = G4NistManager::Instance();

	//Define air
	air = nist->FindOrBuildMaterial("G4_AIR");

	//Define PMMA
	pmma  = nist->FindOrBuildMaterial("G4_PLEXIGLASS"); //Default 1.19 g/cm3

	//Define water
	water  = nist->FindOrBuildMaterial("G4_WATER");

	//Defining polyethylene from NIST and modifying the density
	polyethylene = nist->BuildMaterialWithNewDensity("polyethylene","G4_POLYETHYLENE",0.959*g/cm3);
	polyethylene->GetIonisation()->SetMeanExcitationEnergy(56*eV);

	//polyethylene_NEMA, defualt density 0.94 g/cm3, and excitation energy = 57.4 eV
	polyethylene_NEMA = nist->FindOrBuildMaterial("G4_POLYETHYLENE");


	//Define expanded polystyrene by modifiying the density to mimic lung phantom used in phantom experiment
	polystyrene = nist->BuildMaterialWithNewDensity( "polystyrene","G4_POLYSTYRENE",0.3*g/cm3); 

	isotopes = false;

	//Defile Aluminum material for the detetor cover
	Aluminum = nist->FindOrBuildMaterial("G4_Al", isotopes);

	//Define elements for the GSO  crystal (scintillator) 
	O = nist->FindOrBuildElement("O" , isotopes); 
	Si = nist->FindOrBuildElement("Si", isotopes);
	Gd = nist->FindOrBuildElement("Gd", isotopes);  


	//define GSO crystal for PET detector
	GSO = new G4Material("GSO", 6.7*g/cm3, 3);
	GSO->AddElement(Gd, 2);
	GSO->AddElement(Si, 1);
	GSO->AddElement(O,  5); 
	crystalMaterial   = nist->FindOrBuildMaterial("GSO");
}

////////////////////////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* doiPETDetectorConstruction::Construct()
{  


	//size of the world
	worldSizeX = 2 * m;//0.5 *mm
	worldSizeY = 2 * m;//0.5*mm
	worldSizeZ = 4. * m;


	// Define a solid shape to describe the world volume
	G4Box* solid_world = new G4Box("world", worldSizeX/2., worldSizeY/2., worldSizeZ/2.);

	// Define a logical volume for the world volume
	world_logicalV = new G4LogicalVolume(solid_world, air, "world_logicalV", 0,0,0);


	// Define the physical world volume
	world_physicalV = new G4PVPlacement(0,G4ThreeVector(),world_logicalV, "world_physicalV", 0, false, 0, fCheckOverlaps);
	world_logicalV->SetVisAttributes (G4VisAttributes::Invisible);

	//NOTE!!!
	//The scanner specification (like size and number of crystals in each detctor) are given in the "doiPETGlobalParameters.hh" header file.

	//Each block detector is identified with its unique number, blockIndex.
	blockIndex = 0;

	//Each crystal is identified with its unique number, crystalIndex
	crystalIndex = 0;

	//This is to place the phantom. 

	ConstructPhantom(world_logicalV);

	//Define air volume (box) to fill the detector block. Crystal elements (scintillators) is then placed.
	sizeOfAirBox_DOI = (numberOfCrystal_DOI * sizeOfCrystal_DOI) + (numberOfCrystal_DOI - 1)*crystalGap_DOI;
	sizeOfAirBox_axial = (numberOfCrystal_axial * sizeOfCrystal_axial) + (numberOfCrystal_axial - 1)*crystalGap_axial;
	sizeOfAirBox_tangential = (numberOfCrystal_tangential * sizeOfCrystal_tangential) + (numberOfCrystal_tangential - 1)*crystalGap_tangential;

	//This will pass the size of the detector block so that the position of the PMT will be calculated with respect to the axis of the detector block.
	//see doiPETAnalysis.cc
	pAnalysis = doiPETAnalysis::GetInstance();
	pAnalysis->GetSizeOfDetector(sizeOfAirBox_DOI,sizeOfAirBox_tangential, sizeOfAirBox_axial);

	G4cout<<"size of crytal element: "<<sizeOfCrystal_tangential<<" "<<sizeOfCrystal_axial<<" "<<sizeOfCrystal_DOI<<G4endl;
	G4cout<<"Size of detector block (without Al cover): "<<sizeOfAirBox_tangential<<" "<<sizeOfAirBox_axial<<" "<<sizeOfAirBox_DOI<<G4endl;



	//Define the size of the detector block. 
	sizeOfBlockDetector_DOI = sizeOfAirBox_DOI + AluminumCoverThickness;
	sizeOfBlockDetector_axial = sizeOfAirBox_axial + AluminumCoverThickness;
	sizeOfBlockDetector_tangential = sizeOfAirBox_tangential + AluminumCoverThickness;



	//Define solid shape for the detector block
	G4Box* blockDetector = new G4Box("blockDetector",sizeOfBlockDetector_DOI/2,sizeOfBlockDetector_tangential/2,sizeOfBlockDetector_axial/2);

	//Define the logical volume for the detector block
	blockDetector_logicalV = new G4LogicalVolume(blockDetector,Aluminum,"blockDetector_logicalV", 0,0,0);



	//Define air (box) inside the detector block. Crystal elements will be placed in it.
	G4Box* airBox = new G4Box("airBox", sizeOfAirBox_DOI/2, sizeOfAirBox_tangential/2,sizeOfAirBox_axial/2);

	//Define the logical volume
	airBox_logicalV = new G4LogicalVolume(airBox,air,"airBox_logicalV", 0,0,0);

	//Define its physical volume and place it inside the detector block
	airBox_physicalV = new G4PVPlacement (0,G4ThreeVector(0,0,0),airBox_logicalV,"airBox_physicalV", blockDetector_logicalV,false,0,fCheckOverlaps);


	///////////////////////////////////////// Arrange the PET ring and place the PET detectors in the ring(s) ////////////////////////////////////

	for(G4int Ring = 0; Ring< numberOfRings; Ring++)
	{
		//place the detectors in a ring along the axial direction. Note that the ring gap between two adjcent rings is measured from scintillator to scintillator.  It does not include the Aluminum thickness cover.

		detectorPositionZ = (Ring-((G4double)numberOfRings)/2 + 0.5)*(sizeOfBlockDetector_axial + ringGap - AluminumCoverThickness);

		for(G4int i = 0; i<numberOfDetector_perRing; i++)
		{
			//The azimuthal angle to arrange the detectors in a ring
			thetaDetector = (double)(i*twopi/numberOfDetector_perRing);

			//The radius of the scanner is measured from opposing crystal (scintillator) faces. It does not include the Aluminum thickness cover.
			detectorPositionX = (scannerRadius + sizeOfBlockDetector_DOI/2 - AluminumCoverThickness/2)*std::cos(thetaDetector);
			detectorPositionY = (scannerRadius + sizeOfBlockDetector_DOI/2 - AluminumCoverThickness/2)*std::sin(thetaDetector);

			//Define the rotation matrix for correct placement of detetors
			G4RotationMatrix rotm_PET = G4RotationMatrix();
			rotm_PET.rotateZ(thetaDetector);
			G4ThreeVector uz_PET = G4ThreeVector(detectorPositionX,detectorPositionY,detectorPositionZ);
			G4Transform3D transform = G4Transform3D(rotm_PET,uz_PET);

			//Define the physical volume of the detectors.
			blockDetector_physicalV = new G4PVPlacement (transform,blockDetector_logicalV,"blockDetector_physicalV", world_logicalV,false,blockIndex,fCheckOverlaps);
			blockIndex++;
			//G4cout<<Ring<<" "<<detectorPositionX- ((sizeOfBlockDetector_DOI - AluminumCoverThickness)/2)*cos(thetaDetector)<<" "<<detectorPositionY- ((sizeOfBlockDetector_DOI- AluminumCoverThickness)/2)*sin(thetaDetector)<<" "<<detectorPositionZ<<G4endl;

		}
	}

	//Define the solid crystal
	G4VSolid* CrystalSolid = new G4Box("Crystal", sizeOfCrystal_DOI/2., sizeOfCrystal_tangential/2., sizeOfCrystal_axial/2.);

	//Define the local volume of the crystal
	crystal_logicalV = new G4LogicalVolume(CrystalSolid,crystalMaterial,"Crystal_logicalV", 0,0,0);

	//Place the crystals inside the detectors and give them a unique number with crystalIndex.
	for(G4int i_DOI = 0; i_DOI<numberOfCrystal_DOI; i_DOI++){
		crystalPositionX=(i_DOI-((G4double)numberOfCrystal_DOI)/2 + 0.5)*(sizeOfCrystal_DOI + crystalGap_DOI);
		for(G4int i_axial=0; i_axial< numberOfCrystal_axial;i_axial++){
			crystalPositionZ = (i_axial-((G4double)numberOfCrystal_axial)/2 + 0.5)*(sizeOfCrystal_axial + crystalGap_axial);
			for(G4int i_tan=0; i_tan<numberOfCrystal_tangential;i_tan++){
				crystalPositionY=(i_tan-((G4double)numberOfCrystal_tangential)/2 + 0.5)*(sizeOfCrystal_tangential + crystalGap_tangential);

				//G4cout<<crystalIndex<<" "<<crystalPositionX<<" "<<crystalPositionY<<" "<<crystalPositionZ<<G4endl;
				//place the crystal inside the block detector. 
				crystal_physicalV = new G4PVPlacement (0, G4ThreeVector (crystalPositionX,crystalPositionY,crystalPositionZ), crystal_logicalV, "Crystal_physicalV", airBox_logicalV,false,crystalIndex/*,fCheckOverlaps*/);
				crystalIndex++;
			}
		}
	}

	//******************  Visualization *****************************//

	//visualization for the block detector
	G4VisAttributes* blockDetectorVisAtt;
	blockDetectorVisAtt = new G4VisAttributes(G4Colour(1,1.0,1.0));
	blockDetectorVisAtt->SetVisibility (true);
	//blockDetectorVisAtt->SetForceWireframe (true);
	blockDetector_logicalV->SetVisAttributes (blockDetectorVisAtt);
	//blockDetector_logicalV->SetVisAttributes (G4VisAttributes::Invisible);

	//visualization for the the box filled with air
	G4VisAttributes* airBoxVisAtt;
	airBoxVisAtt = new G4VisAttributes(G4Colour(1,1.0,1.0));
	airBoxVisAtt->SetVisibility (true);
	airBoxVisAtt->SetForceWireframe (true);
	airBox_logicalV->SetVisAttributes (airBoxVisAtt);
	airBox_logicalV->SetVisAttributes (G4VisAttributes::Invisible);

	//visualization for the crystal
	G4VisAttributes* crystalVisAtt;
	crystalVisAtt = new G4VisAttributes(G4Colour(0.88,0.55,1.0));
	//crystalVisAtt->SetVisibility (true);
	crystalVisAtt->SetForceWireframe (true);
	crystal_logicalV->SetVisAttributes (crystalVisAtt);
	crystal_logicalV->SetVisAttributes (G4VisAttributes::Invisible);


	//always return the physical World
	return world_physicalV;
}


/////////////////////////////////////////////// Construct Phantom //////////////////////////////////////////////

void doiPETDetectorConstruction::ConstructPhantom(G4LogicalVolume* worldLogical)
{
	G4cout<<"------------------================- " <<PhantomType<<" Phantom has been used  -==================-----------\n"<<G4endl;

	//The following phantoms are defined based on NEMA NU-2 Standards Publication (Performance Measurements of Positron Emission Tomographs)
	if(PhantomType == "NEMA_imageQualityPhantom_wholeBody"){

		//The body phantom (also called image quality phantom) can presisely be created by the G4UninSolid method from (1) one big half cylinder (with a diameter of 150 mm), 
		//(2) two small full (or quarter) cylinders (with 80mm), and (3) a box with a height width of 140 mm and height of 80mm. 
		//All the material to surround (hold) the water is PMMA.

		//offset position for the body phantom. The center of the body phantom is mover by 35 mm in the +y direction
		yOffsetBodyPhantom =  35*mm;

		//The body phantom is moved by 70 mm in the z direction so that the sphere (cold and hot spheres) are coplanar with the middle slice of the scanner
		zOffsetBodyPhantom =70*mm;

		//length (in the z-direction) of the body phantom. 
		lengthOfBodyPhantom = 180*mm; //Interior length ( = 180m mm) + wallthickness (= 3mm)

		//outer radius of the body phantom. (inner radius + wallthinkness)
		radiusOfBodyPhantom = 147 * mm;

		//wall thickness of the surrounding pmma phantom
		wallThicknessOfBodyPhantom = 3 * mm;

		//
		radiusOfSmallcyl = 77 * mm;
		boxWidth = 140*mm;// along the x-axis
		boxHeight = 77 * mm; // along the y-axis


		/************************* Surrounding PMMA ***********************************/

		//define a half cylinder
		G4Tubs* halfcyl_pmma = new G4Tubs("halfcyl_pmma", 0*mm,(radiusOfBodyPhantom + wallThicknessOfBodyPhantom),(lengthOfBodyPhantom + wallThicknessOfBodyPhantom)/2,0*deg,180*deg);

		//define two full small cylinder. (It can be also be a quarter or half). 
		G4Tubs* cyl1_pmma = new G4Tubs("cyl_pmma", 0*mm,radiusOfSmallcyl + wallThicknessOfBodyPhantom, (lengthOfBodyPhantom + wallThicknessOfBodyPhantom)/2,0*deg,360*deg);
		G4Tubs* cyl2_pmma = new G4Tubs("cyl_pmma", 0*mm,radiusOfSmallcyl + wallThicknessOfBodyPhantom, (lengthOfBodyPhantom + wallThicknessOfBodyPhantom)/2,0*deg,360*deg);

		//define a box
		G4Box*  box_pmma = new G4Box("box_pmma", boxWidth/2, (boxHeight + wallThicknessOfBodyPhantom)/2, (lengthOfBodyPhantom + wallThicknessOfBodyPhantom)/2);

		//a translation vector for the small cylinder with respect to the half cylinder at position (70 mm,0,0)
		G4ThreeVector translation_cyl1 = G4ThreeVector(boxWidth/2, 0, 0);

		//a translation vector for the small cylinder with respect to the half cylinder at position (-70 mm,0,0)
		G4ThreeVector translation_cyl2 = G4ThreeVector(-boxWidth/2, 0, 0);

		//translation for the box with respect to half cylinder
		G4ThreeVector translation_box = G4ThreeVector(0,-(boxHeight + wallThicknessOfBodyPhantom)/2, 0);

		//define union1_pmma by uniting the solids of halfcyl and cyl1
		G4UnionSolid* union1_pmma = new G4UnionSolid("union1",halfcyl_pmma,cyl1_pmma,0,translation_cyl1);

		//define union2_pmma by uniting the solids of union1_pmma and cyl2
		G4UnionSolid* union2_pmma = new G4UnionSolid("union2",union1_pmma,cyl2_pmma,0,translation_cyl2);

		//********** Now define the solid volume of the surrounding PMMA body phantom *********//
		//define phantom by uniting the solids of union2_pmma and box_pmma 
		G4UnionSolid* phantom = new G4UnionSolid("phantom",union2_pmma,box_pmma,0,translation_box);

		//Define the logical volume of the body phantom
		phantom_logicalV = new G4LogicalVolume(phantom, pmma, "phantom_logicalV",0,0,0);

		//Define the physical volume of the body phantom
		phantom_physicalV = new G4PVPlacement(0, G4ThreeVector(0,-yOffsetBodyPhantom,-(lengthOfBodyPhantom + 2*wallThicknessOfBodyPhantom)/2 + zOffsetBodyPhantom), phantom_logicalV, "pmma_phantom_physicalV", worldLogical, false, 0, fCheckOverlaps); 



		//****************************** Water inside the PMMA phantom ********************************/

		//define a half cylinder
		G4Tubs* halfcyl_water = new G4Tubs("halfcyl_water", 0*mm,radiusOfBodyPhantom,lengthOfBodyPhantom/2,0*deg,180*deg);

		//define a full small cylinder phantom. (It can be also be a quarter or half).
		G4Tubs* cyl1_water = new G4Tubs("cyl_water", 0*mm,radiusOfSmallcyl, lengthOfBodyPhantom/2,0*deg,360*deg);
		G4Tubs* cyl2_water = new G4Tubs("cyl_water", 0*mm,radiusOfSmallcyl, lengthOfBodyPhantom/2,0*deg,360*deg);

		//define a box
		G4Box*  box_water = new G4Box("box_water", boxWidth/2, boxHeight/2, lengthOfBodyPhantom/2);

		//a translation vector for the small cylinder with respect to the half cylinder at position (70 mm,0,0)
		G4ThreeVector translation_cyl1_water = G4ThreeVector(boxWidth/2, 0, 0);

		//a translation vector for the small cylinder with respect to the half cylinder at position (-70 mm,0,0)
		G4ThreeVector translation_cyl2_water = G4ThreeVector(-boxWidth/2, 0, 0);

		//translation for the box with respect to half cylinder
		G4ThreeVector translation_box_water = G4ThreeVector(0,-boxHeight/2, 0);

		//define union1_water by uniting the solids of halfcyl and cyl1
		G4UnionSolid* union1_water = new G4UnionSolid("union1",halfcyl_water,cyl1_water,0,translation_cyl1_water);

		//define union2_water by uniting the solids of union1_water and cyl2
		G4UnionSolid* union2_water = new G4UnionSolid("union2",union1_water,cyl2_water,0,translation_cyl2_water);

		//********** Now define the solid volume of the body phantom to be filled with water *********//
		//define phantom by uniting the solids of union2_water and box_water 
		G4UnionSolid* phantom_water = new G4UnionSolid("phantom_water",union2_water,box_water,0,translation_box_water);

		//Define the logical volume of the body phantom
		water_logicalV = new G4LogicalVolume(phantom_water, water, "phantom_logicalV",0,0,0);

		//Define the physical volume of the body phantom
		water_physicalV = new G4PVPlacement(0, G4ThreeVector(0,0,0), water_logicalV, "phantom_physicalV", phantom_logicalV, false, 0, fCheckOverlaps); 

		//
		/////////////////////////////// lung phantom /////////////////////////////////////////////////////

		//A lung phantom (with low density material) is inserted in the body phantom. 
		G4double wallThiknessOfLungPhantom = 4 *mm;
		radiusOfLungPhantom = 21 *mm;

		//define surrounding pmma phantom for the lung
		//Define the solid shape for the lung phantom
		G4Tubs* phantom_lungPMMA = new G4Tubs("Phantom_lung", 0*mm, radiusOfLungPhantom + wallThiknessOfLungPhantom,lengthOfBodyPhantom/2,0*deg,360*deg);

		//Define the logical volume of for the lung phantom
		lung_logicalV_PMMA = new G4LogicalVolume(phantom_lungPMMA, pmma, "coldRegion_logicalV",0,0,0);

		//Define the physical volume for the lung phantom. The center of the lung phantom is moved (in the y-axis) by the same distance as that of the body phantom.
		lung_physicalVPMMA = new G4PVPlacement(0, G4ThreeVector(0,yOffsetBodyPhantom,0), lung_logicalV_PMMA, "coldRegion_physicalV", water_logicalV, false, 0, fCheckOverlaps);


		//Define the solid shape for the lung phantom
		G4Tubs* phantom_lung = new G4Tubs("Phantom_lung", 0*mm,radiusOfLungPhantom,lengthOfBodyPhantom/2,0*deg,360*deg);

		//Define the logical volume of for the lung phantom
		lung_logicalV = new G4LogicalVolume(phantom_lung, polystyrene, "coldRegion_logicalV",0,0,0);

		//Define the physical volume for the lung phantom and place it in the phantom_lungPMMA.
		lung_physicalV = new G4PVPlacement(0, G4ThreeVector(0,0,0), lung_logicalV, "coldRegion_physicalV", lung_logicalV_PMMA, false, 0, fCheckOverlaps);
		//lung_logicalV->SetVisAttributes (G4VisAttributes::Invisible);


		////////////////////////////////// Test phantom ////////////////////////////////////////////////////////
		//This phantom has the same characteristics with that of NECR phantom except its axial position.

		hieghtOfTestPhantom = 705* mm;
		diameterOfTestPhantom = 203 * mm;
		G4double zOffset_testPhantom = 0*mm;//this is to make some gap between the phantoms if necessary

		//Define the solid shape for the test phantom
		G4Tubs* phantom_test = new G4Tubs("Phantom", 0*mm,diameterOfTestPhantom/2,hieghtOfTestPhantom/2,0*deg,360*deg);

		//Define the logical volume of the test phantom
		test_logicalV = new G4LogicalVolume(phantom_test, polyethylene_NEMA, "phantom_logicalV",0,0,0);

		//Define the physical volume of the test phantom. The test phantom is placed next to the body phantom. 
		test_physicalV = new G4PVPlacement(0,G4ThreeVector(0,0,hieghtOfTestPhantom/2+zOffsetBodyPhantom + zOffset_testPhantom), test_logicalV, "phantom_physicalV", worldLogical, false, 0,fCheckOverlaps);
		//test_logicalV->SetVisAttributes (G4VisAttributes::Invisible);

		////////////////////////////////// Six Spherical phantoms for hot and cold lesions (sources) placed inside the body phatom ///////////////////////////////////

		//Define diameter from the center of the body phantom to the center of the spheres
		distanceFromCenter = 114.4*mm;

		//Define total number of spherical phantoms (both hot and cold spheres)
		numberOfSpheres = 6;

		//Define the wall thickness of the spherical phantoms. It should be less than or equal to 1 mm according to the NEMA NU-2 protocol.
		sphereWallThickness = 1*mm;

		//The centers of the sphere phantoms are placed 68 mm from the endplate of the body phantom so that the plane through the centers of the spheres is coplanar with to the middile slice of the scanner
		zOffsetSpherePhantom = lengthOfBodyPhantom/2 + wallThicknessOfBodyPhantom- zOffsetBodyPhantom;

		//Place the spherical phantoms in the body phantom
		for(G4int i = 0; i < numberOfSpheres; i++){
			if(i == 0) sphereDiameter = 37*mm;//cold sphere
			if(i == 1) sphereDiameter = 10*mm;// hot sphere
			if(i == 2) sphereDiameter = 13*mm;// hot sphere
			if(i == 3) sphereDiameter = 17*mm;// hot sphere
			if(i == 4) sphereDiameter = 22*mm;// hot sphere
			if(i == 5) sphereDiameter = 28*mm;// cold sphere

			spherePositionX = distanceFromCenter/2*std::cos(((G4double)i*twopi/numberOfSpheres));
			spherePositionY = distanceFromCenter/2*std::sin(((G4double)i*twopi/numberOfSpheres));
			//zOffsetBodyPhantom = 0;

			//hot sphere phantoms
			if(i>0 && i <5){

				//Surrounding PMMA phantom for the hot sphere
				//Define the solid shape for the surrounding PMMA phantom for the hot sphere pahntom
				G4Sphere* hotSpherePMMA = new G4Sphere("hotSphere", 0., sphereDiameter/2 + sphereWallThickness, 0*deg, 360*deg, 0*deg, 180*deg);

				//Deifne the logical volume of the surrounding PMMA for the hot sphere phantom
				hotSpherePMMA_logicalV = new G4LogicalVolume(hotSpherePMMA, pmma , "hotSphere_logicalV",0,0,0);

				//Deifne the physical volume of the surrounding PMMA for the hot sphere phantom
				hotSpherePMMA_physicalV = new G4PVPlacement(0, G4ThreeVector(spherePositionX,spherePositionY+yOffsetBodyPhantom,zOffsetSpherePhantom), hotSpherePMMA_logicalV, "hotSphere_physicalV", water_logicalV, false, i,fCheckOverlaps);


				//Defining and placing the water in the surrounding PMMA for hot sphere
				//Define the solid shape of the water phantom for the hot sphere phantom
				G4Sphere* hotSphereWater = new G4Sphere("hotSphere", 0., sphereDiameter/2, 0*deg, 360*deg, 0*deg, 180*deg);

				//Define the logical volume of the water phatom for the hot sphere
				hotSphereWater_logicalV = new G4LogicalVolume(hotSphereWater, water , "hotSphere_logicalV",0,0,0);

				//Define the physical volume of the water phatom for the hot sphere
				hotSphereWater_physicalV = new G4PVPlacement(0, G4ThreeVector(0,0,0), hotSphereWater_logicalV, "phantom_physicalV", hotSpherePMMA_logicalV, false, i,fCheckOverlaps);
				//G4cout<<"Hot sphere "<<i<<" is placed. "<< " sphereDiameter = " <<sphereDiameter <<" "<< "Position " <<spherePositionX<<" "<< spherePositionY<<", "<<sphereDiameter<<G4endl;
				//hotSpherePMMA_logicalV->SetVisAttributes (G4VisAttributes::Invisible);
				//hotSphereWater_logicalV->SetVisAttributes (G4VisAttributes::Invisible);
			}

			//Cold sphere phantoms
			if(i==0 || i==5){
				//Surrounding PMMA phantom for the cold sphere
				//Define the solid shape for the surrounding PMMA phantom for the cold sphere pahntom
				G4Sphere* coldSpherePMMA = new G4Sphere("coldSphere", 0., sphereDiameter/2 + sphereWallThickness, 0*deg, 360*deg, 0*deg, 180*deg);

				//Deifne the logical volume of the surrounding PMMA for the cold sphere phantom
				coldSpherePMMA_logicalV = new G4LogicalVolume(coldSpherePMMA, pmma , "coldRegion_logicalV",0,0,0);

				//Deifne the physical volume of the surrounding PMMA for the cold sphere phantom
				coldSpherePMMA_physicalV = new G4PVPlacement(0, G4ThreeVector(spherePositionX,spherePositionY+yOffsetBodyPhantom,zOffsetSpherePhantom), coldSpherePMMA_logicalV, "coldRegion_physicalV", water_logicalV, false, i,fCheckOverlaps);


				//Defining and placing the water in the surrounding PMMA for cold sphere
				//Define the solid shape of the water phantom for the cold sphere phantom
				G4Sphere* coldSphereWater = new G4Sphere("coldSphere", 0., sphereDiameter/2, 0*deg, 360*deg, 0*deg, 180*deg);

				//Define the logical volume of the water phatom for the cold sphere
				coldSphereWater_logicalV = new G4LogicalVolume(coldSphereWater, water , "coldRegion_logicalV",0,0,0);

				//Define the physical volume of the water phatom for the cold sphere
				coldSphereWater_physicalV = new G4PVPlacement(0, G4ThreeVector(0,0,0), coldSphereWater_logicalV, "coldRegion_physicalV", coldSpherePMMA_logicalV, false, i,fCheckOverlaps);
				//G4cout<<"Cold sphere "<<i<<" is placed. "<< " sphereDiameter = " <<sphereDiameter <<" "<< "Position " <<spherePositionX<<" "<< spherePositionY<<" "<<sphereDiameter<<G4endl;
				//coldSpherePMMA_logicalV->SetVisAttributes (G4VisAttributes::Invisible);
				//coldSphereWater_logicalV->SetVisAttributes (G4VisAttributes::Invisible);
			}				
		}
	}

	else if(PhantomType == "NEMA_imageQualityPhantom_smallAnimal"){
		//The follwoing is for NEMA NU-2 image quality phantom for small animal
		//To see the details of the phantom, please see: http://www.qrm.de/content/pdf/QRM-MicroPET-IQ.pdf

		//Outside radius of PMMMA 
		phantomRadius = 16.75*mm;// dia=33.5*mm;

		//Outside length of PMMA
		phantomLength = 63*mm;

		//Dimension of water phantom to be filled with activity (hot region).
		waterPhantomRadius = 15*mm; 
		waterPhantomLength = 30*mm;


		//Dimension of rod phantom (hot region)
		rodPhantomLength = 20*mm;

		//There are five rod phantoms with different diameters
		numberOfRods = 5;

		//Distance from the center of the cylinder to the center of the rods
		distanceFromCenter = 7*mm;	

		//surrounding PMMA phantom.
		//Define PMMA solid
		G4Tubs* phantom = new G4Tubs("Phantom", 0*mm,phantomRadius,phantomLength/2,0*deg,360*deg);

		//Define PMMA logical volume
		phantom_logicalV = new G4LogicalVolume(phantom, pmma, "phantom_logicalV",0,0,0);

		//Define PMMA physical volume
		phantom_physicalV = new G4PVPlacement(0,G4ThreeVector(0,0,0), phantom_logicalV, "pmma_phantom_physicalV", worldLogical, false, 0,fCheckOverlaps);

		//The rods are placed at one end of the surrounding PMMA cylinderical phantom
		rodPositionZ = -waterPhantomLength/2;

		for(G4int i = 0; i < numberOfRods; i++){

			if(i == 0) rodDiameter = 1 * mm;
			if(i == 1) rodDiameter = 2 * mm;
			if(i == 2) rodDiameter = 3 * mm;
			if(i == 3) rodDiameter = 4 * mm;
			if(i == 4) rodDiameter = 5 * mm;

			rodPositionX = distanceFromCenter*std::cos(((G4double)i*twopi/numberOfRods));
			rodPositionY = distanceFromCenter*std::sin(((G4double)i*twopi/numberOfRods));

			//Define rod phantom 
			G4Tubs* rod_phantom = new G4Tubs("Phantom", 0*mm,rodDiameter/2,rodPhantomLength/2,0*deg,360*deg);

			//Define rod phantom logical volume
			rod_phantom_logicalV = new G4LogicalVolume(rod_phantom, water, "rod_phantom_logicalV",0,0,0);

			//Define rod phantom physical volume
			rod_phantom_physicalV = new G4PVPlacement(0,G4ThreeVector(rodPositionX,rodPositionY,rodPositionZ), rod_phantom_logicalV, "phantom_physicalV", phantom_logicalV, false, 0, fCheckOverlaps);
		}

		//Out dimensions of the surrounding PMMA for cold region chambers
		chamberPhantomLength = 15*mm;
		chamberDiameter = 10*mm;

		//Wall thickness of the surrounding PMMA phantom for the cold region chambers
		wallThicknessOfChamber = 1*mm;

		//The centers of the cold chambers is 
		chamberPositionX = distanceFromCenter;
		chamberPositionY = 0*mm;
		chamberPositionZ =  waterPhantomLength/2 - chamberPhantomLength/2;

		//hot region filled with water
		G4Tubs* water_phantom = new G4Tubs("Phantom", 0*mm,waterPhantomRadius,waterPhantomLength/2,0*deg,360*deg);

		waterPhantom_logicalV = new G4LogicalVolume(water_phantom, water, "waterPhantom_logicalV",0,0,0);

		//place the phantom at one end of the PMMA phantom
		WaterPhantom_physicalV = new G4PVPlacement(0,G4ThreeVector(0,0,rodPhantomLength/2), waterPhantom_logicalV, "phantom_physicalV", phantom_logicalV, false, 0,fCheckOverlaps);
		//waterPhantom_logicalV->SetVisAttributes (G4VisAttributes::Invisible);

		//define the surrounding PMMA chamber
		G4Tubs* chamberPMMA = new G4Tubs("chamber_phantom1", 0*mm,chamberDiameter/2,chamberPhantomLength/2,0*deg,360*deg);

		//define the logical volume of the surrounding PMMA chamber
		chamberPMMA_logicalV = new G4LogicalVolume(chamberPMMA, pmma, "chamberPMMA_logicalV",0,0,0);

		//define the physical volume of the surrounding PMMA chamber
		chamberPMMA_physicalV = new G4PVPlacement(0,G4ThreeVector(chamberPositionX, chamberPositionY, chamberPositionZ), chamberPMMA_logicalV, "pmma_phantom_physicalV", waterPhantom_logicalV, false, 0,fCheckOverlaps);


		//Two cold region chambers: one of them filled with (non-radioactive) water and the other with air and is placed inside a hot region
		//chamber one filled with water (cold region)
		//Define the cold region chamber phantom to be filled with water
		G4Tubs* chamberWater = new G4Tubs("chamber_phantom1", 0*mm,chamberDiameter/2 - wallThicknessOfChamber,chamberPhantomLength/2 - wallThicknessOfChamber,0*deg,360*deg);

		//Define the logical volume of the cold region phantom
		chamberWater_logicalV = new G4LogicalVolume(chamberWater, water, "chamber_phantom_logicalV",0,0,0);

		//Define the physical volume of the cold region phantom
		chamberWater_physicalV = new G4PVPlacement(0,G4ThreeVector(0,0,0), chamberWater_logicalV, "coldRegion_physicalV", chamberPMMA_logicalV, false, 0,fCheckOverlaps);

		//chamber2 filled with air 
		//Place the surrounding chamber at another location (at -chamberPositionX)
		chamberPMMA_physicalV = new G4PVPlacement(0,G4ThreeVector(-chamberPositionX, chamberPositionY, chamberPositionZ), chamberPMMA_logicalV, "pmma_phantom_physicalV", waterPhantom_logicalV, false, 0,fCheckOverlaps);

		//Define the logical volume of the cold region phantom to be filled with air
		G4Tubs* chamberAir = new G4Tubs("chamber_phantom1", 0*mm,chamberDiameter/2 - wallThicknessOfChamber,chamberPhantomLength/2 - wallThicknessOfChamber,0*deg,360*deg);

		//Define its logical volume 
		chamberAir_logicalV = new G4LogicalVolume(chamberAir, air, "chamber_phantom_logicalV",0,0,0);

		//Define the physical volume of air-filled chamber 
		chamberAir_physicalV = new G4PVPlacement(0,G4ThreeVector(0,0,0), chamberAir_logicalV, "coldRegion_physicalV", chamberPMMA_logicalV, false, 0,fCheckOverlaps);

	}

	else if(PhantomType == "NEMA_phantom_NECR"){
		//////////////////////////////   Phantom for NECR //////////////////////////////////////
		//The phantom is 203 mm in diameter and 700 mm in height and is placed at the center of the AFOV

		//Define its solid shape 
		G4Tubs* phantom = new G4Tubs("Phantom", 0*mm,phantomRadius,phantomLength/2,0*deg,360*deg);

		//Define its logical volume
		phantom_logicalV = new G4LogicalVolume(phantom, polyethylene_NEMA, "phantom_logicalV",0,0,0);

		//Define its physical volume
		phantom_physicalV = new G4PVPlacement(0,G4ThreeVector(0,0,0), phantom_logicalV, "phantom_physicalV", worldLogical, false, 0,fCheckOverlaps);
	}
	else if(PhantomType == "Phantom_sensitivity"){

		///////////////////////////////////////// Sensitiviy phantom  /////////////////////////////////////////////
		//There are six diffrent sizes of the sensitivity phantom which differe in their diameter. 
		//The size of the phantom can be changed via the macro file
		//The hight of sensitivity phantom is 700 mm

		//Define the solid shape for the sensitivity phantom which surrounds the fillable polyethylene phantom
		G4Tubs* phantom = new G4Tubs("Phantom_sensitivity", 0 *mm,phantomRadius,phantomLength/2,0*deg,360*deg);

		//Define its logical volume. The Matrial of the tubes are Aluminum
		phantom_logicalV = new G4LogicalVolume(phantom, Aluminum, "phantom_logicalV",0,0,0);

		//Define its physical volume
		phantom_physicalV = new G4PVPlacement(0, phantomPosition, phantom_logicalV, "phantom_physicalV", worldLogical, false, 0,fCheckOverlaps);


		//Define the inner most polyethylene (PE) tube with a diameter of 3 mm. This size will not be changed. This is the phantom that holds the source.
		G4Tubs* phantomPE = new G4Tubs("Phantom_sensitivity", 0 *mm,1.5*mm,phantomLength/2,0*deg,360*deg);

		//Define its logical volume. The Matrial of the tubes are polyethylene
		phantomPE_logicalV = new G4LogicalVolume(phantomPE, polyethylene_NEMA, "phantom_logicalV",0,0,0);

		//Define its physical volume
		phantomPE_physicalV = new G4PVPlacement(0, G4ThreeVector(0,0,0), phantomPE_logicalV, "phantom_physicalV", phantom_logicalV, false, 0,fCheckOverlaps);


		G4cout<<"---------------Phantom dimension: "<<"radius (mm) = "<<phantomRadius<<" "<<"Length (mm) = "<<phantomLength<<G4endl;
		G4cout<<"---------------Phantom position (mm) : "<<phantomPosition<<G4endl;
	}
	else if(PhantomType == "Phantom_spatialResolution"){
		////////////////// phantom for point source (For spatial resolution evaluation) /////////////////////////////////
		//The position of the point source phantom can be changed in the macro file. 
		//It has cylindrical shape and its radius and length can be set in the macro file

		//Define its sold shape
		G4Tubs* phantom = new G4Tubs("Phantom_point", 0., phantomRadius,phantomLength/2,0*deg,360*deg);

		//Define its logical volume
		phantom_logicalV = new G4LogicalVolume(phantom, polyethylene_NEMA , "phantom_logicalV",0,0,0);

		//place the phantom
		phantomPosition = G4ThreeVector(phantomPosition.x(),phantomPosition.y(),phantomPosition.z());

		//Define its physical volume
		phantom_physicalV = new G4PVPlacement(0, phantomPosition, phantom_logicalV, "phantom_physicalV", worldLogical, false, 0,fCheckOverlaps);
		G4cout<<"---------------Phantom dimension: "<<"radius (mm) = "<<phantomRadius<<" "<<", Length (mm) = "<<phantomLength<<G4endl;
		G4cout<<"---------------Phantom position: "<<phantomPosition<< "mm"<<G4endl;	
	}

	else if (PhantomType == "Normalization")
	{
		//The normalization phantom is a hallow cylindrical phantom to represent a (rotating) line source. The source is confined in the phantom.
		//The thickness of the phantom is 3 mm, and its diameter is 350 mm. 

		//Define the solid shape. 
		G4Tubs* phantom = new G4Tubs("Phantom", phantomRadius - 3*mm, phantomRadius,phantomLength/2,0*deg,360*deg);

		//Define its logical volume.
		phantom_logicalV = new G4LogicalVolume(phantom, polyethylene_NEMA, "phantom_logicalV",0,0,0);

		//Define its physical volume
		phantom_physicalV = new G4PVPlacement(0, phantomPosition, phantom_logicalV, "phantom_physicalV", worldLogical, false, 0,fCheckOverlaps);


		G4cout<<"---------------Phantom dimension: "<<"radius (mm) = "<<phantomRadius<<" "<<"Length = "<<phantomLength<<G4endl;
		G4cout<<"---------------Phantom position : "<<phantomPosition<< "mm"<<G4endl;

	}
	//-------------------------------------================- No Phantom Set -==================------------------
	else
	{
		G4cerr << "****************************************\nERROR: Phantom Remains: " << PhantomType <<"\n****************************************" << G4endl;
		exit(0);

	}
	///////////Visualization of phantoms///////////

	G4VisAttributes* phantomVisAtt;
	phantomVisAtt = new G4VisAttributes(G4Colour(0.88,0.55,1.0));
	phantomVisAtt->SetVisibility (true);
	phantomVisAtt->SetForceWireframe (true);
	phantom_logicalV->SetVisAttributes (phantomVisAtt);
	//phantom_logicalV->SetVisAttributes (G4VisAttributes::Invisible);	
}

//Change the type of the phantom via .mac file
void doiPETDetectorConstruction::ChangePhantom(G4String NewPhantomtype)
{
	PhantomType = NewPhantomtype;
}

//Change position of the phantom via .mac file
void doiPETDetectorConstruction::SetPhantomPosition(G4ThreeVector NewphantomPosition)
{
	phantomPosition = NewphantomPosition;
}

//Change the radius of the phantom via .mac file
void doiPETDetectorConstruction::SetPhantomRadius(G4double newPhantomRadius){
	phantomRadius = newPhantomRadius;
}

//Change the length of the phantom via .mac file
void doiPETDetectorConstruction::SetPhantomLength(G4double newPhantomLength){
	phantomLength = newPhantomLength;
}



