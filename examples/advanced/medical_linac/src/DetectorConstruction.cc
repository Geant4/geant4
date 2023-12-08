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
// Code developed by
// Silvia Pozzi (1), silvia.pozzi@iss.it
// Barbara Caccia (1), barbara.caccia@iss.it
// Carlo Mancini Terracciano (2), carlo.mancini.terracciano@roma1.infn.it
// (1) Istituto Superiore di Sanita' and INFN Roma, Italy
// (2) Univ. La Sapienza and INFN Roma, Italy

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4ProductionCuts.hh>
#include <G4NistManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VisAttributes.hh>
#include <G4Box.hh>
#include <G4VSolid.hh>
#include <G4Cons.hh>
#include <G4Tubs.hh>
#include <G4BooleanSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4SDManager.hh>
#include <G4GlobalMagFieldMessenger.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4RunManager.hh>
#include <G4GeometryManager.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4SolidStore.hh>
#include <G4Transform3D.hh>
#include <CLHEP/Vector/Rotation.h>


G4VPhysicalVolume* DetectorConstruction::Construct()
{
	detectorMessenger = new DetectorMessenger(this);

	//// Saturn 43
	G4Material* air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

	G4double worldSizeX = 3.5 * m;
	G4double worldSizeY = 3.5 * m;
	G4double worldSizeZ = 3.5 * m;

	G4VSolid* worldBox = new G4Box("world", worldSizeX/2., worldSizeY/2., worldSizeZ/2.);

	G4LogicalVolume* worldLog = new G4LogicalVolume(worldBox, air, "world");

	G4VisAttributes* visAttr = new G4VisAttributes();
	visAttr->SetVisibility(false);
	worldLog->SetVisAttributes(visAttr);

	world_phys = new G4PVPlacement(nullptr, {}, worldLog, "world", nullptr, false, 0);

	ConstructMaterials();
	ConstructAccelerator();
	ConstructPhantom();

	return world_phys;
}

void DetectorConstruction::ConstructPhantom()
{
	G4Material* water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

	G4ThreeVector phantomHalfDim = {phantomSideDim/2., phantomSideDim/2., phantomSideDim/2.};

	G4VSolid* Phantom = new G4Box("Phantom", phantomHalfDim.getX(), phantomHalfDim.getY(), phantomHalfDim.getZ());
	G4LogicalVolume* Phantom_log = new G4LogicalVolume(Phantom, water, "Phantom_log", 0, 0, 0);

	G4VSolid* spessore = new G4Box("Spessore", phantomHalfDim.getX(), phantomHalfDim.getY(), 1.25*mm);
	G4LogicalVolume* spessore_log = new G4LogicalVolume(spessore, water, "spessore_log", 0, 0, 0);
	new G4PVPlacement(nullptr, {0.,0., 1.25*mm}, "spessore_phys", spessore_log,  world_phys, false, 1);

	// axes origin on the face of the phantom
	G4ThreeVector phantomPos = {0., 0., phantomHalfDim.getZ()+voxelDepthDim/2};

	phantom_phys = new G4PVPlacement(nullptr, phantomPos,
				 "phantom_phys", Phantom_log, world_phys, false, 1);

	PhantomSegmentation(Phantom_log);

	// Region for cuts
	G4Region *regVol= new G4Region("fullWaterPhantomR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*mm);
	regVol->SetProductionCuts(cuts);

	Phantom_log -> SetRegion(regVol);
	regVol->AddRootLogicalVolume(Phantom_log);

	//--------- Visualization attributes -------------------------------
	/// Phantom
	G4VisAttributes* simpleH2OVisAtt= new G4VisAttributes(G4Colour::Cyan());
	simpleH2OVisAtt->SetVisibility(true);
	simpleH2OVisAtt->SetForceSolid(false);
	Phantom_log->SetVisAttributes(simpleH2OVisAtt);

}

G4Material* DetectorConstruction::GetMaterial(const G4String materialName)
{
	G4Material* material = nullptr;

		if (materialName == "kapton")
		{
			material = mat_Kapton;
		}
		else if (materialName == "steel")
		{
			material = mat_Ssteel;
		}
		else if (materialName == "xc10")
		{
			material = mat_XC10;
		}
		else if (materialName == "wnicu")
		{
			material = mat_WNICU;
		}
		else
		{
			G4cerr << "* DetectorConstruction::GetMaterial: material " <<
					materialName << " not found" << G4endl;
		}

		return material;
}

void DetectorConstruction::ConstructMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();

	G4Material* el_H  = nist -> FindOrBuildMaterial("G4_H");
	G4Material* el_C  = nist -> FindOrBuildMaterial("G4_C");
	G4Material* el_N  = nist -> FindOrBuildMaterial("G4_N");
	G4Material* el_O  = nist -> FindOrBuildMaterial("G4_O");
	G4Material* el_Si = nist -> FindOrBuildMaterial("G4_Si");
	G4Material* el_Cr = nist -> FindOrBuildMaterial("G4_Cr");
	G4Material* el_Mn = nist -> FindOrBuildMaterial("G4_Mn");
	G4Material* el_Fe = nist -> FindOrBuildMaterial("G4_Fe");
	G4Material* el_Ni = nist -> FindOrBuildMaterial("G4_Ni");
	G4Material* el_Cu = nist -> FindOrBuildMaterial("G4_Cu");
	G4Material* el_W  = nist -> FindOrBuildMaterial("G4_W");

	mat_Ssteel = new G4Material("StainlessSteel", 7.8 *g/cm3, 6);
	mat_Ssteel -> AddMaterial(el_Fe, 0.6898);
	mat_Ssteel -> AddMaterial(el_Cr, 0.18);
	mat_Ssteel -> AddMaterial(el_Ni, 0.10);
	mat_Ssteel -> AddMaterial(el_Mn, 0.02);
	mat_Ssteel -> AddMaterial(el_Si, 0.01);
	mat_Ssteel -> AddMaterial(el_C,  0.0002);

	mat_XC10 = new G4Material("CARBON_STEEL", 7.8 *g/cm3, 3);
	mat_XC10 -> AddMaterial(el_Fe, 0.993);
	mat_XC10 -> AddMaterial(el_Mn, 0.006);
	mat_XC10 -> AddMaterial(el_C,  0.001);

	mat_WNICU = new G4Material("Denal(WNICU)", 16.8*g/cm3, 3);
	mat_WNICU -> AddMaterial(el_W,  0.905);
	mat_WNICU -> AddMaterial(el_Ni, 0.07);
	mat_WNICU -> AddMaterial(el_Cu, 0.025);

	mat_Kapton = new G4Material("Kapton", 1.42*g/cm3, 4);
	mat_Kapton -> AddMaterial(el_C, 69.1133 *perCent);
	mat_Kapton -> AddMaterial(el_O, 20.9235 *perCent);
	mat_Kapton -> AddMaterial(el_N, 7.3270 *perCent);
	mat_Kapton -> AddMaterial(el_H, 2.6362 *perCent);

}

void DetectorConstruction::ConstructAccelerator()
{
	G4bool checkOverlap = true;

	G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

	accHalfSize = {600.*mm, 600.*mm, 600.*mm};
	G4ThreeVector initialCentre = {0.*mm, 0.*mm, -std::abs(sourceToSkinDistance)};

	G4Box* accWorldBox = new G4Box("accWorld", accHalfSize.getX(), accHalfSize.getY(), accHalfSize.getZ());
	G4LogicalVolume* accWorldLV = new G4LogicalVolume(accWorldBox, Vacuum, "accWorldL", 0, 0, 0);
	accWorld_phys = new G4PVPlacement(nullptr, initialCentre, "acceleratorBox", accWorldLV, world_phys, false, checkOverlap);

	G4VisAttributes* simpleAlSVisAtt = new G4VisAttributes(G4Colour::White());
	simpleAlSVisAtt -> SetVisibility(false);
	accWorldLV -> SetVisAttributes(simpleAlSVisAtt);

	ConstructTarget();
	ConstructPrimaryCollimator();
	ConstructVacuumWindow();
	ConstructFlatteningFilter();
	ConstructIonizationChamber();
	ConstructJawsX();
	ConstructJawsY();
}

void DetectorConstruction::ConstructTarget()
{
	G4bool checkOverlap = true;

	G4ThreeVector targetCentre, boxHalfSide, tubeCentre;
	G4double tubeRadius, tubeHalfLengthZ;
	targetCentre.set(0, 0, -3.5*mm);
	boxHalfSide.set(5.*mm, 5.*mm, 7.5*mm);

	tubeRadius 		= 3.*mm;
	tubeHalfLengthZ = 5.5*mm; // 11 mm tot

	G4Material* diskMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ti");
	G4Material* targetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

	// Physical volumes
	G4Box* box = new G4Box("targetBox", boxHalfSide.getX(), boxHalfSide.getY(), boxHalfSide.getZ());
	G4Tubs* tube = new G4Tubs("targetTube", 0.*mm, tubeRadius, tubeHalfLengthZ, 0.*deg, 360.*deg);

	G4SubtractionSolid* BoxMinusTube = new G4SubtractionSolid("TargetSolid", box, tube, nullptr, {0.,0.,-2*mm});
	G4LogicalVolume* BoxMinusTubeLV = new G4LogicalVolume(BoxMinusTube, targetMaterial, "BoxMinusTubeLV",0,0,0);
	new G4PVPlacement(0, targetCentre, "TargetPV", BoxMinusTubeLV, accWorld_phys, false, checkOverlap);

	G4Tubs* disk = new G4Tubs("targetDisk", 0.*mm, 4.*mm, 0.025*mm, 0.*deg, 360.*deg);
	G4LogicalVolume* diskLV = new G4LogicalVolume(disk, diskMaterial, "targetDiskLV",0,0,0);
	new G4PVPlacement(0, {0.,0.,-15.*mm}, "diskPV", diskLV, accWorld_phys, false, checkOverlap);

	// Visualization
	G4VisAttributes* simpleAlSVisAtt;
	simpleAlSVisAtt = new G4VisAttributes(G4Colour::Grey());
	simpleAlSVisAtt -> SetVisibility(true);
	simpleAlSVisAtt -> SetForceSolid(true);
	BoxMinusTubeLV  -> SetVisAttributes(simpleAlSVisAtt);
	diskLV -> SetVisAttributes(simpleAlSVisAtt);

	// Region for cuts
	G4Region *regVol;
	regVol = new G4Region("TargetR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(1.*cm);
	regVol -> SetProductionCuts(cuts);
	BoxMinusTubeLV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(BoxMinusTubeLV);
	diskLV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(diskLV);
}

void DetectorConstruction::ConstructPrimaryCollimator()
{
	G4bool checkOverlap = true;

	G4double bigTube1HalfLengthZ, bigTube1Radius, bigTube1Z;
	G4double bigTube2HalfLengthZ, bigTube2Radius, bigTube2Z;
	G4double bigTube3HalfLengthZ, bigTube3Radius, bigTube3Z;

	G4double tube1Radius, tube1HalfLengthZ, tube2Radius, tube2HalfLengthZ;

	G4double cone1RadiusMin, cone1RadiusMax, cone1HalfLengthZ;
	G4double cone2RadiusMin, cone2RadiusMax, cone2HalfLengthZ;
	G4double cone3RadiusMin, cone3RadiusMax, cone3HalfLengthZ;

	// upper principal tube = bigTube1
	bigTube1Radius 		= 141./2.*mm;
	bigTube1HalfLengthZ = 79./2. * mm;
	bigTube1Z      		= 5. * mm + bigTube1HalfLengthZ;
	tube1Radius 		= 34./2. * mm;
	tube1HalfLengthZ 	= 15.86/2. * mm;
	cone1RadiusMin 		= 34./2.*mm; // upper cone
	cone1RadiusMax 		= 54./2.*mm;
	cone1HalfLengthZ 	= (bigTube1HalfLengthZ - tube1HalfLengthZ)/2.*mm; //(79.-15.86)/2.*mm;

	// middle principal tube = bigTube2
	bigTube2Radius 		= 141./2. * mm;
	bigTube2HalfLengthZ = 57.5/2. * mm;
	bigTube2Z      		= bigTube1Z + bigTube1HalfLengthZ + bigTube2HalfLengthZ;
	tube2Radius 		= 53./2. * mm;
	tube2HalfLengthZ 	= 1.35/2. * mm;
	cone2RadiusMin 		= 53./2.*mm; // middle cone
	cone2RadiusMax 		= 81./2.*mm;
	cone2HalfLengthZ 	= bigTube2HalfLengthZ - tube2HalfLengthZ; // (57.50 - 1.35)/2.*mm;

	// lower principal tube = bigTube3
	bigTube3Radius 		= 241./2.*mm;
	bigTube3HalfLengthZ = 40.50/2.*mm;
	bigTube3Z 			= bigTube2Z + bigTube2HalfLengthZ + 19.5 + bigTube3HalfLengthZ;
	cone3RadiusMin 		= 88.06/2.*mm; // lower cone
	cone3RadiusMax 		= 109.76/2.*mm;
	cone3HalfLengthZ 	= 40.50/2.*mm;

	G4Material* upperTubeMaterial = GetMaterial("xc10");
	G4Material* middleTubeMaterial = GetMaterial("wnicu");
	G4Material* lowerTubeMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");

	// Physical volumes
	G4Tubs* bigTube1 = new G4Tubs("PriCollBigTube1", 0.*mm, bigTube1Radius, bigTube1HalfLengthZ, 0.*deg, 360.*deg);
	G4Tubs* tube1 = new G4Tubs("PriCollTube1", 0.*mm, tube1Radius, tube1HalfLengthZ, 0.*deg, 360.*deg);
	G4Cons* cone1 = new G4Cons("PriCollCone1", 0., cone1RadiusMin, 0., cone1RadiusMax, cone1HalfLengthZ, 0, 360.*deg);

	G4Tubs* bigTube2 = new G4Tubs("PriCollBigTube2", 0.*mm, bigTube2Radius, bigTube2HalfLengthZ, 0.*deg, 360.*deg);
	G4Tubs* tube2 = new G4Tubs("PriCollTube2", 0.*mm, tube2Radius, tube2HalfLengthZ, 0.*deg, 360.*deg);
	G4Cons* cone2 = new G4Cons("PriCollCone2", 0., cone2RadiusMin, 0., cone2RadiusMax, cone2HalfLengthZ, 0, 360.*deg);

	G4Tubs* bigTube3 = new G4Tubs("PriCollBigTube3", 0.*mm, bigTube3Radius, bigTube3HalfLengthZ, 0.*deg, 360.*deg);
	G4Cons* cone3 = new G4Cons("PriCollCone3", 0., cone3RadiusMin, 0., cone3RadiusMax, cone3HalfLengthZ, 0, 360.*deg);

	G4UnionSolid* tube1AndCone1 = new G4UnionSolid("PriCollTube1AndCone1",
				tube1, cone1, nullptr, {0., 0., (tube1HalfLengthZ + cone1HalfLengthZ)});
	G4SubtractionSolid* bigTube1NotTube1AndCone1 = new G4SubtractionSolid(
			"PriColltube1NotTube1AndCone1", bigTube1, tube1AndCone1, 0,
			{0., 0., (-bigTube1HalfLengthZ + tube1HalfLengthZ)});
	G4LogicalVolume* bigTube1NotTube1AndCone1LV = new G4LogicalVolume(
			bigTube1NotTube1AndCone1, upperTubeMaterial, "PriCollBigTube1NotTube1AndCone1LV", 0, 0, 0);

	G4UnionSolid* tube2AndCone2 = new G4UnionSolid("PriCollTube2AndCone2",
				tube2, cone2, nullptr, {0., 0., (tube2HalfLengthZ + cone2HalfLengthZ)});
	G4SubtractionSolid* bigTube2NotTube2AndCone2 = new G4SubtractionSolid(
			"PriCollBigTube2NotTube2Cone2", bigTube2, tube2AndCone2, 0,
			{0., 0., (-bigTube2HalfLengthZ + tube2HalfLengthZ)});
	G4LogicalVolume* bigTube2NotTube2AndCone2LV  = new G4LogicalVolume(
			bigTube2NotTube2AndCone2, middleTubeMaterial, "PriCollTube2NotTube2Cone2LV", 0, 0, 0);

	G4SubtractionSolid* bigTube3NotCone3 = new G4SubtractionSolid("PriCollTube4NotCone3", bigTube3, cone3);
	G4LogicalVolume* bigTube3NotCone3LV  = new G4LogicalVolume(bigTube3NotCone3,
			lowerTubeMaterial, "PriCollBigTube3NotCone3LV",0,0,0);

	new G4PVPlacement(nullptr, {0, 0, bigTube1Z}, "PriCollUpperPV",
			bigTube1NotTube1AndCone1LV, accWorld_phys, false, checkOverlap);
	new G4PVPlacement(nullptr, {0, 0, bigTube2Z}, "PriCollMiddlePV",
			bigTube2NotTube2AndCone2LV, accWorld_phys, false, checkOverlap);
	new G4PVPlacement(nullptr, {0, 0, bigTube3Z}, "PriCollLowerPV",
			bigTube3NotCone3LV, accWorld_phys, false, checkOverlap);

	// Visualization
	G4VisAttributes* simpleAlSVisAtt1 = new G4VisAttributes(G4Colour::Green());
	G4VisAttributes* simpleAlSVisAtt2 = new G4VisAttributes(G4Colour::Red());
	G4VisAttributes* simpleAlSVisAtt3 = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAtt1 -> SetVisibility(true);
	simpleAlSVisAtt1 -> SetForceSolid(false);
	simpleAlSVisAtt1 -> SetForceWireframe(true);
	bigTube1NotTube1AndCone1LV -> SetVisAttributes(simpleAlSVisAtt1);  //primary collimator
	simpleAlSVisAtt2 -> SetVisibility(true);
	simpleAlSVisAtt2 -> SetForceSolid(false);
	simpleAlSVisAtt2 -> SetForceWireframe(true);
	bigTube2NotTube2AndCone2LV -> SetVisAttributes(simpleAlSVisAtt2); //primary collimator
	simpleAlSVisAtt3 -> SetVisibility(true);
	simpleAlSVisAtt3 -> SetForceSolid(false);
	simpleAlSVisAtt3 -> SetForceWireframe(true);
	bigTube3NotCone3LV -> SetVisAttributes(simpleAlSVisAtt3); //secondary collimator

	// Region for cuts
	G4Region* regVol = new G4Region("primaryCollimatorR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(1.*cm);
	regVol->SetProductionCuts(cuts);
	bigTube1NotTube1AndCone1LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(bigTube1NotTube1AndCone1LV);
	bigTube2NotTube2AndCone2LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(bigTube2NotTube2AndCone2LV);
	bigTube3NotCone3LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(bigTube3NotCone3LV);
}

void DetectorConstruction::ConstructFlatteningFilter()
{
	G4double tube1HalfLengthZ; //tube1Radius,
	G4double cone1RadiusMin, cone1RadiusMax, cone1HalfLengthZ;
	G4double cone2RadiusMin, cone2RadiusMax, cone2HalfLengthZ;
	G4double cone3RadiusMin, cone3RadiusMax, cone3HalfLengthZ, ffZ;

	tubeFFRadius 	 = 108./2.*mm; // top principal tube
	tube1HalfLengthZ = 7.5/2.*mm;

	cone1RadiusMin 	 = 54./2.*mm; // top cone
	cone1RadiusMax   = 76./2.*mm;
	cone1HalfLengthZ = 13.7/2.*mm;

	cone2RadiusMin   = 8./2.*mm; // mid cone
	cone2RadiusMax   = 54./2.*mm;
	cone2HalfLengthZ = (44.3-13.7)/2.*mm;

	cone3RadiusMin   = 0.000001*mm; // bottom cone
	cone3RadiusMax   = 8./2.*mm;
	cone3HalfLengthZ = (46.8-44.3)/2.*mm;
	tubeFFFirstFaceZ   = 149.50*mm;

	ffZ = tubeFFFirstFaceZ + tube1HalfLengthZ; // the half point is located ad the centre of the tube1 because of the solids unions and translations

	G4Material* ffMaterial = GetMaterial("steel");

	// Physical volumes
	G4Tubs* tube1 = new G4Tubs("ffTube1", 0.*mm, tubeFFRadius, tube1HalfLengthZ, 0.*deg, 360.*deg);
	G4Cons* cone1 = new G4Cons("ffCone1", 0., cone1RadiusMax, 0., cone1RadiusMin, cone1HalfLengthZ, 0, 360.*deg);
	G4Cons* cone2 = new G4Cons("ffCone2", 0., cone2RadiusMax, 0., cone2RadiusMin, cone2HalfLengthZ, 0, 360.*deg);
	G4Cons* cone3 = new G4Cons("ffCone3", 0., cone3RadiusMax, 0., cone3RadiusMin, cone3HalfLengthZ, 0, 360.*deg);

	G4double pos = (cone1HalfLengthZ - tube1HalfLengthZ - 1.)*mm; // == (13.7/2. - 7.5/2. -1.)*mm
	G4UnionSolid *tube1AndCone1 = new G4UnionSolid("ffTube1AndCone1", tube1,
			cone1, 0, {0., 0., pos});

	pos += cone1HalfLengthZ + cone2HalfLengthZ;
	G4UnionSolid* tubeCone1AndCone2 = new G4UnionSolid("fftubeCone1AndCone2",
			tube1AndCone1, cone2, 0, {0., 0., pos});

	pos += cone2HalfLengthZ + cone3HalfLengthZ;
	G4UnionSolid* tubeCone12AndCone3 = new G4UnionSolid("fftubeCone12AndCone3",
			tubeCone1AndCone2, cone3, 0, {0., 0., pos});

	G4LogicalVolume* ffLV = new G4LogicalVolume(tubeCone12AndCone3, ffMaterial, "ffLV",0,0,0);

	new G4PVPlacement(0, {0,0,ffZ}, "ffPV", ffLV, accWorld_phys, false, 0);

	// Visualization
	G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::Magenta());
	simpleAlSVisAtt->SetVisibility(true);
	simpleAlSVisAtt->SetForceSolid(true);
	ffLV->SetVisAttributes(simpleAlSVisAtt);

	// Region for cuts
	G4Region* regVol;
	regVol= new G4Region("ffR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);
	ffLV->SetRegion(regVol);
	regVol->AddRootLogicalVolume(ffLV);
}

void DetectorConstruction::ConstructIonizationChamber()
{
	G4bool checkOverlap = true;

	G4double tubeRadius, tubeA1Z, tubeA2Z, tubeA3Z, tubeA4Z, tubeA5Z, tubeA6Z;
	G4double kaptonThickness, AlThickness1, AlThickness8;

	kaptonThickness = 0.025/2. *mm;
	AlThickness1 	= 1.6e-5/2. *cm;
	AlThickness8 	= 8e-6/2. *cm;

	tubeRadius = 110./2.*mm; // upper principal tube
	tubeA1Z = (202.5 + AlThickness8)*mm;
	tubeA2Z = (204.5 + AlThickness1)*mm;
	tubeA3Z = (206.5 + AlThickness8)*mm;
	tubeA4Z = (207.5 + AlThickness8)*mm;
	tubeA5Z = (209.5 + AlThickness1)*mm;
	tubeA6Z = (211.5 + AlThickness8)*mm;

	G4double pos1 = (AlThickness1 + kaptonThickness)*mm;
	G4double pos8 = (AlThickness8 + kaptonThickness)*mm;

	G4Material* el_Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
	G4Material* Kapton = GetMaterial("kapton");

	// Physical volumes
	G4Tubs* tubeAl8 = new G4Tubs("ICTube1A", 0.*mm, tubeRadius, AlThickness8, 0.*deg, 360.*deg);
	G4Tubs* tubeK = new G4Tubs("ICTube1B", 0.*mm, tubeRadius, kaptonThickness, 0.*deg, 360.*deg);
	G4Tubs* tubeAl1 = new G4Tubs("ICTube2A", 0.*mm, tubeRadius, AlThickness1, 0.*deg, 360.*deg);

	G4LogicalVolume* IC_Al8_LV = new G4LogicalVolume(tubeAl8, el_Al, "IC_Al8_LV", 0, 0, 0);
	G4LogicalVolume* IC_Al1_LV = new G4LogicalVolume(tubeAl1, el_Al, "IC_Al1_LV", 0, 0, 0);
	G4LogicalVolume* IC_K_LV   = new G4LogicalVolume(tubeK, Kapton, "IC_Al_LV",  0, 0, 0);

	new G4PVPlacement(0, {0,0,tubeA1Z}, "IC_AL8_PV1", IC_Al8_LV, accWorld_phys, false, checkOverlap);
	new G4PVPlacement(0, {0,0,tubeA1Z + pos8}, "IC_K_PV1", IC_K_LV, accWorld_phys, false, checkOverlap);

	new G4PVPlacement(0, {0,0,tubeA2Z}, "IC_AL1_PV2", IC_Al1_LV, accWorld_phys, false, checkOverlap);
	new G4PVPlacement(0, {0,0,tubeA2Z + pos1}, "IC_K_PV2", IC_K_LV, accWorld_phys, false, checkOverlap);

	new G4PVPlacement(0, {0,0,tubeA3Z}, "IC_AL8_PV3", IC_Al8_LV, accWorld_phys, false, checkOverlap);
	new G4PVPlacement(0, {0,0,tubeA3Z + pos8}, "IC_K_PV3", IC_K_LV, accWorld_phys, false, checkOverlap);

	new G4PVPlacement(0, {0,0,tubeA4Z}, "IC_AL8_PV4", IC_Al8_LV, accWorld_phys, false, checkOverlap);
	new G4PVPlacement(0, {0,0,tubeA4Z + pos8}, "IC_K_PV4", IC_K_LV, accWorld_phys, false, checkOverlap);

	new G4PVPlacement(0, {0,0,tubeA5Z}, "IC_AL1_PV5", IC_Al1_LV, accWorld_phys, false, checkOverlap);
	new G4PVPlacement(0, {0,0,tubeA5Z + pos1}, "IC_K_PV5", IC_K_LV, accWorld_phys, false, checkOverlap);

	new G4PVPlacement(0, {0,0,tubeA6Z}, "IC_AL8_PV6", IC_Al8_LV, accWorld_phys, false, checkOverlap);
	new G4PVPlacement(0, {0,0,tubeA6Z + pos8}, "IC_K_PV6", IC_K_LV, accWorld_phys, false, checkOverlap);

	// Visualization
	G4VisAttributes* simpleAlSVisAttAl1 = new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAttAl1 -> SetVisibility(true);
	simpleAlSVisAttAl1 -> SetForceSolid(true);
	IC_Al1_LV -> SetVisAttributes(simpleAlSVisAttAl1);

	G4VisAttributes* simpleAlSVisAttAl8 = new G4VisAttributes(G4Colour::Magenta());
	simpleAlSVisAttAl8 -> SetVisibility(true);
	simpleAlSVisAttAl8 -> SetForceSolid(true);
	IC_Al8_LV -> SetVisAttributes(simpleAlSVisAttAl8);

	G4VisAttributes* simpleAlSVisAttK = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAttK -> SetVisibility(true);
	simpleAlSVisAttK -> SetForceSolid(true);
	IC_K_LV -> SetVisAttributes(simpleAlSVisAttK);

	// Region for cuts
	G4Region * regVol = new G4Region("IonChamberR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(0.1*mm);
	regVol -> SetProductionCuts(cuts);

	IC_Al1_LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(IC_Al1_LV);
	IC_Al8_LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(IC_Al8_LV);
	IC_K_LV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(IC_K_LV);
}

void DetectorConstruction::ConstructVacuumWindow()
{
	G4double tubeRadius, tubeHalfLengthZ, tubeZ;

	tubeRadius 		= 116.53/2.*mm; // upper principal tube
	tubeHalfLengthZ = 2./2.*mm;
	tubeZ 			= 215.75*mm + tubeHalfLengthZ;

	G4Material* el_Al= G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

	// Physical volumes
	G4Tubs* tube = new G4Tubs("VacWindowTube", 0.*mm, tubeRadius, tubeHalfLengthZ, 0.*deg, 360.*deg);
	G4LogicalVolume* tubeLV = new G4LogicalVolume(tube, el_Al, "VacWindowTubeLV",0,0,0);

	new G4PVPlacement(0, {0,0,tubeZ}, "VacWindowTubePV", tubeLV, accWorld_phys, false, 0);

	// Visualization
	G4VisAttributes* simpleAlSVisAtt = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAtt -> SetVisibility(true);
	simpleAlSVisAtt -> SetForceSolid(true);
	tubeLV -> SetVisAttributes(simpleAlSVisAtt);

	// Region for cuts
	G4Region* regVol= new G4Region("VacWindowR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts->SetProductionCut(0.1*cm);
	regVol->SetProductionCuts(cuts);
	tubeLV -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(tubeLV);
}
void DetectorConstruction::ConstructJawsX()
{
	G4bool checkOverlap = false;

	G4double boxSide 		 = 101.*mm;
	G4double box1HalfLengthZ =   3.*mm;
	G4double box2HalfLengthZ =  35.*mm;
	G4double box3HalfLengthZ =  35.*mm;
	G4double box4HalfLengthZ =  27.*mm;
	G4double box5HalfLengthZ =  10.*mm;
	G4double box6HalfLengthZ =   5.*mm;

	G4Material* el_Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
	G4Material* XC10 = GetMaterial("xc10");
	G4Material* WNICU = GetMaterial("wnicu");

	G4ThreeVector boxJawXSide = {boxSide, boxSide, 212.*mm};

	G4Box* boxJaw1X = new G4Box("Jaw1Xbox", boxJawXSide.getX()/2., boxJawXSide.getY()/2., boxJawXSide.getZ()/2.);
	G4LogicalVolume* boxJaw1XLV = new G4LogicalVolume(boxJaw1X,
			G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "Jaw1XboxLV", 0, 0, 0);
	G4Box* boxJaw2X = new G4Box("Jaw2Xbox", boxJawXSide.getX()/2., boxJawXSide.getY()/2., boxJawXSide.getZ()/2.);
	G4LogicalVolume* boxJaw2XLV = new G4LogicalVolume(boxJaw2X,
				G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "Jaw2XboxLV", 0, 0, 0);

	jaw1XInitialPos.setX(boxJawXSide.getX()/2.);
	jaw1XInitialPos.setY(0.);
	jaw1XInitialPos.setZ(275.5*mm + boxJawXSide.getZ()/2.);

	G4RotationMatrix* cRotationJaw1X = new G4RotationMatrix();
	cRotationJaw1X -> rotateY(std::fabs(std::atan(jawAperture/2/(-(sourceToSkinDistance+10*cm)))));

	G4ThreeVector position = *cRotationJaw1X * jaw1XInitialPos;
	*cRotationJaw1X = cRotationJaw1X -> inverse();

	boxJaw1X_phys =
		new G4PVPlacement (cRotationJaw1X, position,
				"Jaw1XPV", boxJaw1XLV, accWorld_phys, false, 0, checkOverlap);

	jaw2XInitialPos.setX(-jaw1XInitialPos.getX());
	jaw2XInitialPos.setY(jaw1XInitialPos.getY());
	jaw2XInitialPos.setZ(jaw1XInitialPos.getZ());

	G4RotationMatrix* cRotationJaw2X = new G4RotationMatrix();
	cRotationJaw2X -> rotateY(-std::fabs(std::atan(jawAperture/2/(-(sourceToSkinDistance+10*cm)))));

	position = *cRotationJaw2X * jaw2XInitialPos;
	*cRotationJaw2X = cRotationJaw2X -> inverse();

	boxJaw2X_phys =
				new G4PVPlacement (cRotationJaw2X, position,
						"Jaw2XPV", boxJaw2XLV, accWorld_phys, false, 0, 0);

	// Physical volumes
	G4Box* box1 = new G4Box("Jaws1XBox1", boxSide/2., boxSide/2., box1HalfLengthZ/2.);
	G4Box* box2 = new G4Box("Jaws1XBox2", boxSide/2., boxSide/2., box2HalfLengthZ/2.);
	G4Box* box3 = new G4Box("Jaws1XBox3", boxSide/2., boxSide/2., box3HalfLengthZ/2.);
	G4Box* box4 = new G4Box("Jaws1XBox4", boxSide/2., boxSide/2., box4HalfLengthZ/2.);
	G4Box* box5 = new G4Box("Jaws1XBox5", boxSide/2., boxSide/2., box5HalfLengthZ/2.);
	G4Box* box6 = new G4Box("Jaws1XBox6", boxSide/2., boxSide/2., box6HalfLengthZ/2.);
	G4LogicalVolume* box1LV1 = new G4LogicalVolume(box1, XC10,  "Jaws1XLV1", 0, 0, 0);
	G4LogicalVolume* box2LV1 = new G4LogicalVolume(box2, el_Pb, "Jaws1XLV2", 0, 0, 0);
	G4LogicalVolume* box3LV1 = new G4LogicalVolume(box3, WNICU, "Jaws1XLV3", 0, 0, 0);
	G4LogicalVolume* box4LV1 = new G4LogicalVolume(box4, el_Pb, "Jaws1XLV4", 0, 0, 0);
	G4LogicalVolume* box5LV1 = new G4LogicalVolume(box5, el_Pb, "Jaws1XLV5", 0, 0, 0);
	G4LogicalVolume* box6LV1 = new G4LogicalVolume(box6, WNICU, "Jaws1XLV6", 0, 0, 0);

	G4LogicalVolume* box1LV2 = new G4LogicalVolume(box1, XC10,  "Jaws2XLV1", 0, 0, 0);
	G4LogicalVolume* box2LV2 = new G4LogicalVolume(box2, el_Pb, "Jaws2XLV2", 0, 0, 0);
	G4LogicalVolume* box3LV2 = new G4LogicalVolume(box3, WNICU, "Jaws2XLV3", 0, 0, 0);
	G4LogicalVolume* box4LV2 = new G4LogicalVolume(box4, el_Pb, "Jaws2XLV4", 0, 0, 0);
	G4LogicalVolume* box5LV2 = new G4LogicalVolume(box5, el_Pb, "Jaws2XLV5", 0, 0, 0);
	G4LogicalVolume* box6LV2 = new G4LogicalVolume(box6, WNICU, "Jaws2XLV6", 0, 0, 0);

	G4ThreeVector centre = {0., 0., 0.};
	G4double zCentreCurrentBox = -boxJawXSide.getZ()/2.+box1HalfLengthZ/2.*mm;

	centre.setZ(zCentreCurrentBox);
	new G4PVPlacement(nullptr, centre, "Jaws1XPV1", box1LV1, boxJaw1X_phys, false, 0, checkOverlap);
	new G4PVPlacement(nullptr, centre, "Jaws2XPV1", box1LV2, boxJaw2X_phys, false, 0, checkOverlap);

	zCentreCurrentBox += (box1HalfLengthZ+box2HalfLengthZ)/2.*mm;
	centre.setZ(zCentreCurrentBox);
	new G4PVPlacement(nullptr, centre, "Jaws1XPV2", box2LV1, boxJaw1X_phys, false, 0, checkOverlap);
	new G4PVPlacement(nullptr, centre, "Jaws2XPV2", box2LV2, boxJaw2X_phys, false, 0, checkOverlap);

	zCentreCurrentBox += (box2HalfLengthZ+box3HalfLengthZ)/2.*mm;
	centre.setZ(zCentreCurrentBox);
	new G4PVPlacement(nullptr, centre, "Jaws1XPV3", box3LV1, boxJaw1X_phys, false, 0, checkOverlap);
	new G4PVPlacement(nullptr, centre, "Jaws2XPV3", box3LV2, boxJaw2X_phys, false, 0, checkOverlap);

	zCentreCurrentBox += (box3HalfLengthZ+box4HalfLengthZ)/2.*mm;
	centre.setZ(zCentreCurrentBox);
	new G4PVPlacement(nullptr, centre, "Jaws1XPV4", box4LV1, boxJaw1X_phys, false, 0, checkOverlap);
	new G4PVPlacement(nullptr, centre, "Jaws2XPV4", box4LV2, boxJaw2X_phys, false, 0, checkOverlap);

	zCentreCurrentBox += (box4HalfLengthZ+box5HalfLengthZ)/2.+97.*mm;
	centre.setZ(zCentreCurrentBox);
	new G4PVPlacement(nullptr, centre, "Jaws1XPV5", box5LV1, boxJaw1X_phys, false, 0, checkOverlap);
	new G4PVPlacement(nullptr, centre, "Jaws2XPV5", box5LV2, boxJaw2X_phys, false, 0, checkOverlap);

	zCentreCurrentBox += (box5HalfLengthZ+box6HalfLengthZ)/2.;
	centre.setZ(zCentreCurrentBox);
	new G4PVPlacement(nullptr, centre, "Jaws1XPV6", box6LV1, boxJaw1X_phys, false, 0, checkOverlap);
	new G4PVPlacement(nullptr, centre, "Jaws2XPV6", box6LV2, boxJaw2X_phys, false, 0, checkOverlap);

	// Visualization
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour::Yellow());
	visAtt -> SetVisibility(false);
	visAtt -> SetForceSolid(false);

	G4VisAttributes* simpleAlSVisAttPb = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAttPb -> SetVisibility(true);
	simpleAlSVisAttPb -> SetForceSolid(true);

	G4VisAttributes* simpleAlSVisAttXC10 = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAttXC10 -> SetVisibility(true);
	simpleAlSVisAttXC10 ->SetForceSolid(true);

	G4VisAttributes* simpleAlSVisAttWNICU = new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAttWNICU ->SetVisibility(true);
	simpleAlSVisAttWNICU ->SetForceSolid(true);

	box1LV1 -> SetVisAttributes(simpleAlSVisAttXC10);
	box2LV1 -> SetVisAttributes(simpleAlSVisAttPb);
	box3LV1 -> SetVisAttributes(simpleAlSVisAttWNICU);
	box4LV1 -> SetVisAttributes(simpleAlSVisAttPb);
	box5LV1 -> SetVisAttributes(simpleAlSVisAttPb);
	box6LV1 -> SetVisAttributes(simpleAlSVisAttWNICU);

	box1LV2 -> SetVisAttributes(simpleAlSVisAttXC10);
	box2LV2 -> SetVisAttributes(simpleAlSVisAttPb);
	box3LV2 -> SetVisAttributes(simpleAlSVisAttWNICU);
	box4LV2 -> SetVisAttributes(simpleAlSVisAttPb);
	box5LV2 -> SetVisAttributes(simpleAlSVisAttPb);
	box6LV2 -> SetVisAttributes(simpleAlSVisAttWNICU);

	boxJaw1XLV -> SetVisAttributes(visAtt);
	boxJaw2XLV -> SetVisAttributes(visAtt);

	// Region for cuts
	G4Region* regVol = new G4Region("JawsXR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(2.*cm);
	regVol -> SetProductionCuts(cuts);
	box1LV1 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box1LV1);
	box2LV1 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box2LV1);
	box3LV1 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box3LV1);
	box4LV1 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box4LV1);
	box5LV1 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box5LV1);
	box6LV1 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box6LV1);

	box1LV2 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box1LV2);
	box2LV2 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box2LV2);
	box3LV2 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box3LV2);
	box4LV2 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box4LV2);
	box5LV2 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box5LV2);
	box6LV2 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box6LV2);
}

void DetectorConstruction::ConstructJawsY()
{
	G4bool checkOverlap = false;

	G4double boxSide 		 = 101.*mm;
	G4double box1HalfLengthZ =  15.*mm;
	G4double box2HalfLengthZ =  31.*mm;
	G4double box3HalfLengthZ =  35.*mm;
	G4double box4HalfLengthZ =  21.*mm;

	G4Material* el_Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
	G4Material* XC10 = GetMaterial("xc10");
	G4Material* WNICU = GetMaterial("wnicu");

	G4ThreeVector boxJawYSide = {boxSide, boxSide, 219.5*mm}; // 219, perché 219.5??

	G4Box* boxJaw1Y = new G4Box("Jaw1Ybox", boxJawYSide.getX()/2., boxJawYSide.getY()/2., boxJawYSide.getZ()/2.);
	G4LogicalVolume* boxJaw1YLV = new G4LogicalVolume(boxJaw1Y,
			G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "Jaw1YboxLV", 0, 0, 0);
	G4Box* boxJaw2Y = new G4Box("Jaw2Ybox", boxJawYSide.getX()/2., boxJawYSide.getY()/2., boxJawYSide.getZ()/2.);
		G4LogicalVolume* boxJaw2YLV = new G4LogicalVolume(boxJaw2Y,
				G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "Jaw2YboxLV", 0, 0, 0);

	jaw1YInitialPos.setX(0.);
	jaw1YInitialPos.setY(boxJawYSide.getX()/2.);
	jaw1YInitialPos.setZ(248.5*mm + boxJawYSide.getZ()/2.);

	G4RotationMatrix* cRotationJaw1Y = new G4RotationMatrix();
	cRotationJaw1Y -> rotateX(-std::fabs(std::atan(jawAperture/2/-(-(sourceToSkinDistance+10*cm)))));

	G4ThreeVector position = *cRotationJaw1Y * jaw1YInitialPos;
	*cRotationJaw1Y = cRotationJaw1Y -> inverse();

	boxJaw1Y_phys =
			new G4PVPlacement (cRotationJaw1Y, position,
					"Jaw1YPV", boxJaw1YLV, accWorld_phys, false, 0, checkOverlap);

	jaw2YInitialPos.setX(jaw1YInitialPos.getX());
	jaw2YInitialPos.setY(-jaw1YInitialPos.getY());
	jaw2YInitialPos.setZ(jaw1YInitialPos.getZ());

	G4RotationMatrix* cRotationJaw2Y = new G4RotationMatrix();
	cRotationJaw2Y -> rotateX(std::fabs(std::atan(jawAperture/2/(-(sourceToSkinDistance+10*cm)))));

	position = *cRotationJaw2Y * jaw2YInitialPos;
	*cRotationJaw2Y = cRotationJaw2Y -> inverse();

	boxJaw2Y_phys =
			new G4PVPlacement (cRotationJaw2Y, position,
					"Jaw2YPV", boxJaw2YLV, accWorld_phys, false, 0, checkOverlap);

	// Physical volumes
	G4Box *box1 = new G4Box("Jaws1YBox1", boxSide/2., boxSide/2., box1HalfLengthZ/2.);
	G4Box *box2 = new G4Box("Jaws1YBox2", boxSide/2., boxSide/2., box2HalfLengthZ/2.);
	G4Box *box3 = new G4Box("Jaws1YBox3", boxSide/2., boxSide/2., box3HalfLengthZ/2.);
	G4Box *box4 = new G4Box("Jaws1YBox4", boxSide/2., boxSide/2., box4HalfLengthZ/2.);
	G4LogicalVolume* box1LV1 = new G4LogicalVolume(box1, XC10,  "Jaws1YLV1", 0,0,0);
	G4LogicalVolume* box2LV1 = new G4LogicalVolume(box2, el_Pb, "Jaws1YLV2", 0,0,0);
	G4LogicalVolume* box3LV1 = new G4LogicalVolume(box3, WNICU, "Jaws1YLV3", 0,0,0);
	G4LogicalVolume* box4LV1 = new G4LogicalVolume(box4, el_Pb, "Jaws1YLV4", 0,0,0);

	G4LogicalVolume* box1LV2 = new G4LogicalVolume(box1, XC10,  "Jaws2YLV1", 0,0,0);
	G4LogicalVolume* box2LV2 = new G4LogicalVolume(box2, el_Pb, "Jaws2YLV2", 0,0,0);
	G4LogicalVolume* box3LV2 = new G4LogicalVolume(box3, WNICU, "Jaws2YLV3", 0,0,0);
	G4LogicalVolume* box4LV2 = new G4LogicalVolume(box4, el_Pb, "Jaws2YLV4", 0,0,0);

	G4ThreeVector centre = {0., 0., 0.};
	G4double zCentreCurrentBox = -boxJawYSide.getZ()/2.+box1HalfLengthZ/2.*mm;

	centre.setZ(zCentreCurrentBox);
	new G4PVPlacement(nullptr, centre, "Jaws1YPV1", box1LV1, boxJaw1Y_phys, false, 0, checkOverlap);
	new G4PVPlacement(nullptr, centre, "Jaws2YPV1", box1LV2, boxJaw2Y_phys, false, 0, checkOverlap);

	zCentreCurrentBox += (box1HalfLengthZ+box2HalfLengthZ)/2.*mm + 117.5*mm;
	centre.setZ(zCentreCurrentBox);
	new G4PVPlacement(nullptr, centre, "Jaws1YPV2", box2LV1, boxJaw1Y_phys, false, 0, checkOverlap);
	new G4PVPlacement(nullptr, centre, "Jaws2YPV2", box2LV2, boxJaw2Y_phys, false, 0, checkOverlap);

	zCentreCurrentBox += (box2HalfLengthZ+box3HalfLengthZ)/2.*mm;
	centre.setZ(zCentreCurrentBox);
	new G4PVPlacement(nullptr, centre, "Jaws1YPV3", box3LV1, boxJaw1Y_phys, false, 0, checkOverlap);
	new G4PVPlacement(nullptr, centre, "Jaws2YPV3", box3LV2, boxJaw2Y_phys, false, 0, checkOverlap);

	zCentreCurrentBox += (box3HalfLengthZ+box4HalfLengthZ)/2.*mm;
	centre.setZ(zCentreCurrentBox);
	new G4PVPlacement(nullptr, centre, "Jaws1YPV4", box4LV1, boxJaw1Y_phys, false, 0, checkOverlap);
	new G4PVPlacement(nullptr, centre, "Jaws2YPV4", box4LV2, boxJaw2Y_phys, false, 0, checkOverlap);

	// Visualization
	G4VisAttributes* visAtt = new G4VisAttributes(G4Colour::Yellow());
	visAtt -> SetVisibility(false);
	visAtt -> SetForceSolid(false);

	G4VisAttributes* simpleAlSVisAttPb = new G4VisAttributes(G4Colour::Blue());
	simpleAlSVisAttPb -> SetVisibility(true);
	simpleAlSVisAttPb -> SetForceSolid(true);

	G4VisAttributes* simpleAlSVisAttXC10 = new G4VisAttributes(G4Colour::Green());
	simpleAlSVisAttXC10 -> SetVisibility(true);
	simpleAlSVisAttXC10 -> SetForceSolid(true);

	G4VisAttributes* simpleAlSVisAttWNICU = new G4VisAttributes(G4Colour::Red());
	simpleAlSVisAttWNICU -> SetVisibility(true);
	simpleAlSVisAttWNICU -> SetForceSolid(true);

	box1LV1 -> SetVisAttributes(simpleAlSVisAttWNICU);
	box2LV1 -> SetVisAttributes(simpleAlSVisAttPb);
	box3LV1 -> SetVisAttributes(simpleAlSVisAttWNICU);
	box4LV1 -> SetVisAttributes(simpleAlSVisAttPb);

	box1LV2 -> SetVisAttributes(simpleAlSVisAttWNICU);
	box2LV2 -> SetVisAttributes(simpleAlSVisAttPb);
	box3LV2 -> SetVisAttributes(simpleAlSVisAttWNICU);
	box4LV2 -> SetVisAttributes(simpleAlSVisAttPb);

	boxJaw1YLV -> SetVisAttributes(visAtt);
	boxJaw2YLV -> SetVisAttributes(visAtt);

	// Region for cuts
	G4Region* regVol = new G4Region("JawsYR");
	G4ProductionCuts* cuts = new G4ProductionCuts;
	cuts -> SetProductionCut(2.*cm);
	regVol -> SetProductionCuts(cuts);

	box1LV1 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box1LV1);
	box2LV1 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box2LV1);
	box3LV1 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box3LV1);
	box4LV1 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box4LV1);

	box1LV2 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box1LV2);
	box2LV2 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box2LV2);
	box3LV2 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box3LV2);
	box4LV2 -> SetRegion(regVol);
	regVol -> AddRootLogicalVolume(box4LV2);
}

void DetectorConstruction::PhantomSegmentation(G4LogicalVolume* Phantom_log)
{

	G4ThreeVector phantomHalfDim = {phantomSideDim/2., phantomSideDim/2., phantomSideDim/2.};

	nSideCells = CheckPhantomSegmentation(phantomSideDim/voxelSideDim);
	nDepthCells = CheckPhantomSegmentation(phantomSideDim/voxelDepthDim);

	G4Material* water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

	G4String yRepName("RepY");
	G4VSolid* solYRep = new G4Box(yRepName, phantomHalfDim.getX(),
			phantomHalfDim.getY()/nSideCells, phantomHalfDim.getZ());
	G4LogicalVolume* logYRep = new G4LogicalVolume(solYRep, water, yRepName);
	new G4PVReplica(yRepName, logYRep, Phantom_log, kYAxis, nSideCells,
			2.*phantomHalfDim.getY()/nSideCells);

	G4String xRepName("RepX");
	G4VSolid* solXRep = new G4Box(xRepName, phantomHalfDim.getX()/nSideCells,
			phantomHalfDim.getY()/nSideCells, phantomHalfDim.getZ());
	G4LogicalVolume* logXRep = new G4LogicalVolume(solXRep, water, xRepName);
	new G4PVReplica(xRepName, logXRep, logYRep, kXAxis, nSideCells,
			2.*phantomHalfDim.getX()/nSideCells);

	G4String zVoxName("phantomSens");
	G4VSolid* solVoxel = new G4Box(zVoxName, phantomHalfDim.getX()/nSideCells,
			phantomHalfDim.getY()/nSideCells, phantomHalfDim.getZ()/nDepthCells);

	LVPhantomSens = new G4LogicalVolume(solVoxel, water, zVoxName); // This is the Sensitive Volume
	new G4PVReplica(zVoxName, LVPhantomSens, logXRep, kZAxis, nDepthCells,
			2.*phantomHalfDim.getZ()/nDepthCells);

}

G4int DetectorConstruction::CheckPhantomSegmentation(G4int nCells)
{
	// I want an odd number of voxels
	// check whether nCells is integer
	if ( nCells != (int)nCells )
	{
		nCells = std::ceil(nCells);
	}
	// check whether nCells is even
	if ( std::fmod(nCells,2) == 0)
	{
		nCells += 1;
	}

	return nCells;
}

G4double DetectorConstruction::GetAccOriginPosition()
{
	return sourceToSkinDistance;
}

G4double DetectorConstruction::GetFFilterRadius()
{
	return tubeFFRadius;
}

G4double DetectorConstruction::GetFFilterZ()
{
	return tubeFFFirstFaceZ;
}

G4double DetectorConstruction::GetVoxelDepthDim()
{
	return voxelDepthDim;
}

G4int DetectorConstruction::GetNumberSideCells() const
{
	return nSideCells;
}

G4int DetectorConstruction::GetNumberDepthCells() const
{
	return nDepthCells;
}

G4int DetectorConstruction::GetPhantomDepth() const
{
	return phantomSideDim;
}

void DetectorConstruction::SetJaws(G4double value)
{
	jawAperture = value;
}

void DetectorConstruction::SetTargetPosition(G4double value)
{
	sourceToSkinDistance = value;
}

void DetectorConstruction::SetPhantomSide(G4double value)
{
	phantomSideDim = value;
}

void DetectorConstruction::SetVoxelSide(G4double value)
{
	voxelSideDim = value;
}

void DetectorConstruction::SetVoxelDepth(G4double value)
{
	voxelDepthDim = value;
}

void DetectorConstruction::UpdateGeometry(G4String string, G4double value)
{

	if ( string == "fieldSide" )
	{
		//
		G4double halfFieldSide = value/2;
		G4ThreeVector position;
		G4RotationMatrix* rMatrix1X = new G4RotationMatrix();

		rMatrix1X -> rotateY(std::fabs(std::atan(halfFieldSide/-sourceToSkinDistance)));
		position = *rMatrix1X * jaw1XInitialPos;
		boxJaw1X_phys -> SetTranslation(position);
		*rMatrix1X = rMatrix1X -> inverse();
		boxJaw1X_phys -> SetRotation(rMatrix1X);

		G4RotationMatrix* rMatrix2X = new G4RotationMatrix();
		rMatrix2X -> rotateY(-std::fabs(std::atan(halfFieldSide/-sourceToSkinDistance)));

		position = *rMatrix2X * jaw2XInitialPos;
		boxJaw2X_phys -> SetTranslation(position);
		*rMatrix2X = rMatrix2X -> inverse();
		boxJaw2X_phys -> SetRotation(rMatrix2X);

		G4RotationMatrix* rMatrix1Y = new G4RotationMatrix();
		rMatrix1Y -> rotateX(-std::fabs(std::atan(halfFieldSide/-sourceToSkinDistance)));
		position = *rMatrix1Y * jaw1YInitialPos;
		boxJaw1Y_phys -> SetTranslation(position);
		*rMatrix1Y = rMatrix1Y -> inverse();
		boxJaw1Y_phys -> SetRotation(rMatrix1Y);

		G4RotationMatrix* rMatrix2Y = new G4RotationMatrix();
		rMatrix2Y -> rotateX(std::fabs(std::atan(halfFieldSide/-sourceToSkinDistance)));

		position = *rMatrix2Y * jaw2YInitialPos;
		boxJaw2Y_phys -> SetTranslation(position);
		*rMatrix2Y = rMatrix2Y -> inverse();
		boxJaw2Y_phys -> SetRotation(rMatrix2Y);
	}
	else if ( string == "sourceToSkinDistance")
	{
		accWorld_phys -> SetTranslation({0., 0., value});
	}
	else if (string == "phantomSide")
	{
		G4Box* box = (G4Box *) phantom_phys -> GetLogicalVolume() -> GetSolid();

        // change phantom side
		box -> SetXHalfLength(phantomSideDim/2.);
		box -> SetYHalfLength(phantomSideDim/2.);
		box -> SetZHalfLength(phantomSideDim/2.);
		phantom_phys -> SetTranslation({0., 0., phantomSideDim/2.});

		// clear daughters and rebuild them
		phantom_phys -> GetLogicalVolume() -> ClearDaughters();
		PhantomSegmentation(phantom_phys -> GetLogicalVolume());
	}
	else if (string == "voxelSide")
	{
		phantom_phys -> GetLogicalVolume() -> ClearDaughters();
		PhantomSegmentation(phantom_phys -> GetLogicalVolume());
	}
	else if (string == "voxelDepth")
	{
		phantom_phys -> GetLogicalVolume() -> ClearDaughters();
		PhantomSegmentation(phantom_phys -> GetLogicalVolume());
	}
	else
		G4cerr << "*** DetectorMessenger::UpdateGeometry: Command not found" << G4endl;

	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}


