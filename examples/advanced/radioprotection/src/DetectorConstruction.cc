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
// Authors: Susanna Guatelli and Francesco Romano
// susanna@uow.edu.au, francesco.romano@ct.infn.it

// Modified by Jacopo Magini: j.magini@surrey.ac.uk

// Modified by David Bolst: db001@uowmail.edu.au
// Added geometry of "bridge" microdosimeter (ConstructSiliconBridgeDetector())

#include "DetectorConstruction.hh"
#include "globals.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
//#include "G4SubtractionSolid.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "SensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4NistManager.hh"

DetectorConstruction::DetectorConstruction(AnalysisManager* analysis_manager, G4String detector)
{
	analysis = analysis_manager;
	
	detectorType = detector;
}

DetectorConstruction::~DetectorConstruction(){

}

G4VPhysicalVolume* DetectorConstruction::Construct()
{

	if( detectorType == "Diamond" ) return ConstructDiamondDetector();
	
	else if( detectorType == "MicroDiamond" ) return ConstructMicroDiamondDetector();
	
	else if( detectorType == "Silicon" ) return ConstructSiliconDetector();
	
	else if( detectorType == "SiliconBridge" ) return ConstructSiliconBridgeDetector();
	
	else
	{
		G4cout << "ERROR: " << detectorType << " is not an allowed detector type. ";
		G4cout << "Did you change some code in radioprotection.cc and/or DetectorMessenger.cc ?" << G4endl;
		
		return 0;
	}
}

G4VPhysicalVolume* DetectorConstruction::ConstructDiamondDetector()
{

//Define each individual element
//Define Nitrogen
 G4double A = 14.01 * g/mole;
 G4double Z = 7;
 G4Element* elN = new G4Element ("Nitrogen", "N", Z, A);

//Define Oxygen
 A = 16.0 * g/mole;
 Z = 8;
 G4Element* elO = new G4Element ("Oxygen", "O", Z, A);

//Define Hydrogen 
 A = 1.01 * g/mole;
 Z = 1;
 G4Element* elH = new G4Element ("Hydrogen", "H", Z, A);

//Define Boron
 A = 10.8 * g/mole;
 Z = 5;
 G4Element* elB = new G4Element ("Boron", "B", Z, A);

//Define Carbon
 A = 12.01 * g/mole;
 Z = 6;
 G4Element* elC = new G4Element ("Carbon", "C", Z, A);

//Define Air   
 G4Material* Air = new G4Material("Air", 1.29*mg/cm3, 2);
 Air -> AddElement(elN, 70*perCent);
 Air -> AddElement(elO, 30*perCent);

//Define diamond
 A = 12.01 * g/mole;
 Z = 6;
 G4Material* diamond = new G4Material("diamond", Z, A, 3.515*g/cm3);
			
//Define dopant (boron doped diamond)
 G4Material* dopant = new G4Material("dopant", 3.514*g/cm3, 2);
 dopant -> AddElement(elC, 99.9994*perCent);
 dopant -> AddElement(elB, 0.0006*perCent);

 //Define Aluminium contacts (AlContact)
 A = 26.981 * g/mole;
 Z = 13;
 G4Material* AlContact = new G4Material("AlContact", Z, A, 2.7 *g/cm3);

 //Define Gold contact (AuContact)
 A = 196.97 * g/mole;
 Z = 79;
 G4Material* AuContact = new G4Material("AuContact", Z, A, 19.3 *g/cm3);
 
 //Define PMMA (C502H8)
 // NIST reference 
 G4Material* PMMA = new G4Material("PMMA", 1.19*g/cm3, 3);
 PMMA -> AddElement(elC, 5);
 PMMA -> AddElement(elO, 2);
 PMMA -> AddElement(elH, 8);

 //define water
 G4Material* water = new G4Material("water", 1*g/cm3, 2);
 water -> AddElement(elH, 2);
 water -> AddElement(elO, 1);
	
 //Define Vacuum
 G4double vacuumDensity = 1.e-25 *g/cm3;
 G4double pressure = 3.e-18*pascal;
 G4double temperature = 2.73*kelvin;
 G4Material* vacuum = new G4Material("Galactic", Z=1., A=1.01*g/mole,
			         vacuumDensity,kStateGas,temperature,pressure);

 //Define volumes
 // World volume  has size 1cm
 G4double worldx = 0.5 * m;  //half length!!!!
 G4double worldy = 0.5 * m;
 G4double worldz = 0.5 * m;

 // World volume, containing all geometry
 G4Box* world = new G4Box("world_box", worldx, worldy, worldz);

 G4LogicalVolume* logical_world = new G4LogicalVolume(world, vacuum, "world_log", 0,0,0);

 //set the logical world volume invisible
 logical_world -> SetVisAttributes(G4VisAttributes::GetInvisible());

 G4VPhysicalVolume* physical_world = new G4PVPlacement(0,
						       G4ThreeVector(),
						       logical_world, 
							"world_phys",
							0, 
							false, 
							0);

	
 // Define the geometry of the diamond microdosimeter
 // mother volume of the detector components
 G4double DiaVol_x = 300*micrometer;
 G4double DiaVol_y = 240*micrometer;
 G4double DiaVol_z = 150*micrometer; 

 G4Box* DiaVol_box = new G4Box("DiaVol_box",DiaVol_x,DiaVol_y,DiaVol_z);

 G4LogicalVolume* logical_DiaVol = new G4LogicalVolume(DiaVol_box, diamond, "DiaVol_log", 0,0,0);

 new G4PVPlacement(0, G4ThreeVector(0,0,0), logical_DiaVol,"DiaVol_phys",
			  logical_world, 
			  false, 0, true);

 //VacBlock for contact placement
 G4double vacblock_x = 300*um;
 G4double vacblock_y = 240*um;
 G4double vacblock_z = 0.25*um; 

 G4Box* vacblock_box = new G4Box("vacblock_box",vacblock_x,vacblock_y,vacblock_z);

 G4LogicalVolume* logical_vacblock = new G4LogicalVolume(vacblock_box, vacuum, "vacblock_log", 0,0,0);

 new G4PVPlacement(0, 
	           G4ThreeVector(0,0,DiaVol_z - vacblock_z),
		   logical_vacblock,
	           "vacblock_phys",
	           logical_DiaVol, 
	           false, 
	           0, true);
//Bdl in DiaVol
 G4double Bdl_x = 300*micrometer;
 G4double Bdl_y = 240*micrometer;
 G4double Bdl_z = 0.69*micrometer; 
	
 G4Box* Bdl_box = new G4Box("Bdl_box",Bdl_x,Bdl_y,Bdl_z);

 G4LogicalVolume* logical_Bdl = new G4LogicalVolume(Bdl_box, dopant, "Bdl_log", 0,0,0);

 new G4PVPlacement(0, 
		   G4ThreeVector(0,0,DiaVol_z - Bdl_z - vacblock_z- vacblock_z),
		   logical_Bdl,
		   "Bdl_phys",
	           logical_DiaVol,   //mother volume 
                   false, 
		   0, true);

 //Diamond SV
 G4double SV_x = 75*um;
 G4double SV_y = 75*um;
 G4double SV_z = 0.69*um; 

 G4Box* SV_box = new G4Box("SV_box",SV_x,SV_y,SV_z);

 G4LogicalVolume* logical_SV = new G4LogicalVolume(SV_box, diamond, "SV_log", 0,0,0);

 new G4PVPlacement(0, G4ThreeVector(-45*um,105*um,0*um), logical_SV,"SV_phys1",
		    logical_Bdl,false, 0, true);

 new G4PVPlacement(0, G4ThreeVector(165*um,105*um,0*um), logical_SV,"SV_phys2",
		   logical_Bdl, false, 0, true);

 new G4PVPlacement(0, G4ThreeVector(-45*um,-105*um,0*um),logical_SV,"SV_phys3", 
		   logical_Bdl, false, 0, true);

 new G4PVPlacement(0, G4ThreeVector(165*um,-105*um,0*um),logical_SV,"SV_phys4",
		   logical_Bdl, false, 0, true);

//Al strips
//10 nm thickness 
 G4double AlStrip_x = 240*um;
 G4double AlStrip_y = 240*um;
 G4double AlStrip_z = vacblock_z; 

 G4Box* AlStrip = new G4Box("AlStrip",AlStrip_x,AlStrip_y,AlStrip_z);

 G4LogicalVolume* logical_AlStrip = new G4LogicalVolume(AlStrip, AlContact, "AlStrip_log", 0,0,0);

 new G4PVPlacement(0, G4ThreeVector(60*um,0,0), logical_AlStrip, "AlStrip_phys",
                   logical_vacblock, false, 0, true);

//gold cylinder in vacblock
 G4double innerRadiusOfTheTube1 = 0.*um;
 G4double outerRadiusOfTheTube1 = 45.*um;
 G4double heightOfTheTube1 = 10*nm;
 G4double startAngleOfTheTube1 = 0.*deg;
 G4double spanningAngleOfTheTube1 = 360.*deg;

 G4Tubs* GoldCylinder1 = new G4Tubs("GoldCylinder1", innerRadiusOfTheTube1, 
                 		    outerRadiusOfTheTube1,
                 		    heightOfTheTube1,
                 		    startAngleOfTheTube1, 
                 		    spanningAngleOfTheTube1);
	
 G4LogicalVolume* logical_GoldCylinder1 = new G4LogicalVolume(GoldCylinder1, AuContact, "GoldCylinder1_log", 0,0,0);

 new G4PVPlacement(0,G4ThreeVector(-245*um,0,-vacblock_z + heightOfTheTube1),
	           		   logical_GoldCylinder1,
		   		   "GoldCylinder1_phys",
		   		   logical_vacblock, false, 0, true);

//gold contacts
 G4double innerRadiusOfTheTube2 = 0.*um;
 G4double outerRadiusOfTheTube2 = 45.*um;
 G4double heightOfTheTube2 = Bdl_z;
 G4double startAngleOfTheTube2 = 0.*deg;
 G4double spanningAngleOfTheTube2 = 360.*deg;

 G4Tubs* GoldCylinder2 = new G4Tubs("GoldCylinder2",
                 		    innerRadiusOfTheTube2, 
                 		    outerRadiusOfTheTube2,
                 		    heightOfTheTube2,
                 		    startAngleOfTheTube2, 
                 		    spanningAngleOfTheTube2);
	
 G4LogicalVolume* logical_GoldCylinder2 = new G4LogicalVolume(GoldCylinder2, AuContact, "GoldCylinder2_log", 0,0,0);

 new G4PVPlacement(0, G4ThreeVector(-245*um,0,0), logical_GoldCylinder2, "GoldCylinder2_phys",
		   logical_Bdl, false, 0, true);

//gold cylinder in DiaVol
 G4double innerRadiusOfTheTube3 = 0.*um;
 G4double outerRadiusOfTheTube3 = 45.*um;
 G4double heightOfTheTube3 = 75.*um -heightOfTheTube2 - heightOfTheTube1 ;
 G4double startAngleOfTheTube3 = 0.*deg;
 G4double spanningAngleOfTheTube3 = 360.*deg;

 G4Tubs* GoldCylinder3 = new G4Tubs("GoldCylinder3",
                 		    innerRadiusOfTheTube3, 
                 		    outerRadiusOfTheTube3,
                 		    heightOfTheTube3,
                 		    startAngleOfTheTube3, 
                 		    spanningAngleOfTheTube3);
	
G4LogicalVolume* logical_GoldCylinder3 = new G4LogicalVolume(GoldCylinder3, AuContact, "GoldCylinder3_log", 0,0,0);

new G4PVPlacement(0, G4ThreeVector(-245*um,0,DiaVol_z - heightOfTheTube3 - Bdl_z - Bdl_z - vacblock_z- vacblock_z),
			   	   logical_GoldCylinder3,
			   	   "GoldCylinder3_phys",
			   	   logical_DiaVol, 
			   	   false, 
			   	   0, true);

// Visualisation attributes

        logical_DiaVol -> SetVisAttributes(G4VisAttributes(G4Colour(255,255,255))); //white
	logical_Bdl -> SetVisAttributes(G4VisAttributes(G4Colour(0,255,0)));        //green
	
	G4VisAttributes vis_SV(G4Colour(198, 226, 255));
	vis_SV.SetForceSolid(true);
	logical_SV -> SetVisAttributes(vis_SV);
        logical_vacblock -> SetVisAttributes(G4VisAttributes::GetInvisible());		
	logical_AlStrip -> SetVisAttributes(G4VisAttributes(G4Colour(0, 255, 255)));//cyan
	
	G4VisAttributes vis_GoldCylinder1(G4Colour(255, 255, 0));                    
	vis_GoldCylinder1.SetForceAuxEdgeVisible(true);
	logical_GoldCylinder1 -> SetVisAttributes(vis_GoldCylinder1);
	
	G4VisAttributes vis_GoldCylinder2(G4Colour(255, 255, 0));                    
	vis_GoldCylinder2.SetForceAuxEdgeVisible(true);
	logical_GoldCylinder2 -> SetVisAttributes(vis_GoldCylinder2); 
	
	G4VisAttributes vis_GoldCylinder3(G4Colour(255, 255, 0));                    
	vis_GoldCylinder3.SetForceAuxEdgeVisible(true);
	logical_GoldCylinder3 -> SetVisAttributes(vis_GoldCylinder3); 
        
return physical_world; 

}

G4VPhysicalVolume* DetectorConstruction::ConstructMicroDiamondDetector()
{
	//Define each individual element
	//Define Nitrogen
	G4double A = 14.01 * g/mole;
	G4double Z = 7;
	G4Element* elN = new G4Element ("Nitrogen", "N", Z, A);

	//Define Oxygen
	A = 16.0 * g/mole;
	Z = 8;
	G4Element* elO = new G4Element ("Oxygen", "O", Z, A);

	//Define Boron
	A = 10.8 * g/mole;
	Z = 5;
	G4Element* elB = new G4Element ("Boron", "B", Z, A);

	//Define Carbon
	A = 12.01 * g/mole;
	Z = 6;
	G4Element* elC = new G4Element ("Carbon", "C", Z, A);

	//Define Air   
	G4Material* Air = new G4Material("Air", 1.29*mg/cm3, 2);
	Air -> AddElement(elN, 70*perCent);
	Air -> AddElement(elO, 30*perCent);

	//Define diamond
	A = 12.01 * g/mole;
	Z = 6;
	G4Material* diamond = new G4Material("diamond", Z, A, 3.515*g/cm3);
			
	//Define p-type diamond (boron doped diamond)
	G4Material* p_diamond = new G4Material("p_diamond", 3.514*g/cm3, 2);
	// Boron concentration used is 1e20 cm-3, considering the diamond density and a Boron atomic weight of 10.811u
	p_diamond -> AddElement(elC, 99.94887*perCent);
	p_diamond -> AddElement(elB, 0.05113*perCent);

	//Define chromium contact
	G4Material* chromium = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cr");
	
	//Define Vacuum
	G4double vacuumDensity = 1.e-25 *g/cm3;
	G4double pressure = 3.e-18*pascal;
	G4double temperature = 2.73*kelvin;
	G4Material* vacuum = new G4Material("Galactic", Z=1., A=1.01*g/mole,
			         vacuumDensity,kStateGas,temperature,pressure);

	//Define volumes
	// World volume  has size 1cm
	G4double worldx = 0.5 * m;  //half length!!!!
	G4double worldy = 0.5 * m;
	G4double worldz = 0.5 * m;

	// World volume, containing all geometry
	G4Box* world = new G4Box("world_box", worldx, worldy, worldz);

	G4LogicalVolume* logical_world = new G4LogicalVolume(world, vacuum, "world_log", 0,0,0);

	//set the logical world volume invisible
	logical_world -> SetVisAttributes(G4VisAttributes::GetInvisible());

	G4VPhysicalVolume* physical_world = new G4PVPlacement(0,
								G4ThreeVector(),
								logical_world, 
								"world_phys",
								0, 
								false, 
								0);

	// sentive volume
	G4double SVside = 100.*um /2.;
	G4double SVthickness = 8.*um /2.;
	G4double SVspacing = 200.*um; //edge-edge distance
	
	G4Box* SV_box = new G4Box("SV_box", SVside, SVside, SVthickness);

	G4LogicalVolume* logical_SV = new G4LogicalVolume(SV_box, diamond, "SV_log", 0,0,0);
	
	G4VisAttributes SVcolour(G4Colour(0.5, 0.5, 0.5));
	SVcolour.SetForceSolid(true);
	logical_SV -> SetVisAttributes(SVcolour);

	// chromium front-electrode
	G4double feThickness = 50.*nm /2.; // front-electrode thickness

	G4Box* fe_box = new G4Box("frontElec_box", SVside, SVside, feThickness);

	G4LogicalVolume* logical_fe = new G4LogicalVolume(fe_box, chromium, "frontElec_log", 0,0,0);
	
	G4VisAttributes fe_colour(G4Colour::Brown());
	fe_colour.SetForceSolid(false);
	logical_fe -> SetVisAttributes(fe_colour);
	
	// p-type diamond
	G4double pDthickness = 1.*um /2.; // p-type diamond back-electrode thickness
	
	G4Box* pD_box = new G4Box("pDiam_box", SVside, SVside, pDthickness);
	
	G4LogicalVolume* logical_pD = new G4LogicalVolume(pD_box, p_diamond, "pDiam_box", 0,0,0);
	
	G4VisAttributes pDcolour(G4Colour::Blue());
	pDcolour.SetForceSolid(false);
	
	logical_pD -> SetVisAttributes(pDcolour);

	// put them in place
	G4ThreeVector SVposition = {0., 0., SVthickness};	//position of the first edge (left)
	G4ThreeVector fePosition = {0., 0., 2.*SVthickness + feThickness};
	G4ThreeVector pDposition = {0., 0., -pDthickness};
	
	G4double SVposition_x[4] = { -3.*SVside -1.5*SVspacing, -SVside -0.5*SVspacing, +SVside +0.5*SVspacing, +3.*SVside +1.5*SVspacing };

	std::ostringstream PVName;

	for( int i=0; i<4; i++)
	{	
		// sensitive volume
		SVposition[0] = SVposition_x[i];
		PVName << "SV_phys" << i;
		new G4PVPlacement(0, SVposition, logical_SV, PVName.str(),
					logical_world,
					false, 0, true);
		PVName.str("");	//reset the string
		
		// chromium front-electrode
		PVName << "frontElec_phys" << i;
		fePosition[0] = SVposition[0];
		new G4PVPlacement(0, fePosition, logical_fe, PVName.str(),
					logical_world,
					false, 0, true);
		PVName.str("");
		
		// p-type diamond back-electrode
		PVName << "pD_phys" << i;
		pDposition[0] = SVposition[0];
		new G4PVPlacement(0, pDposition, logical_pD, PVName.str(),
					logical_world,
					false, 0, true);
		PVName.str("");		
	}

	// HPHT diamond substrate (only one big substrate for simplicity)
	G4double subs_x = 2.*mm /2.;
	G4double subs_y = 0.5*mm /2.; 
	G4double sub_z = 300.*micrometer /2.; 
	
	G4Box* sub_box = new G4Box("sub_box", subs_x, subs_y, sub_z);
	
	G4LogicalVolume* logical_sub = new G4LogicalVolume(sub_box, diamond, "sub_log", 0,0,0);
	
	G4ThreeVector subPosition = {0,0, -2.*pDthickness -sub_z};
	
	new G4PVPlacement(0, subPosition, logical_sub, "sub_phys",
				logical_world,
				false, 0, true);

	G4VisAttributes subColour(G4Colour(0.5, 0.5, 0.5));
	subColour.SetForceSolid(false);
	logical_sub -> SetVisAttributes(subColour);
	
	return physical_world;
}

G4VPhysicalVolume* DetectorConstruction::ConstructSiliconDetector()
{
	//load NIST database
	G4NistManager* nist = G4NistManager::Instance();
	nist->SetVerbose(1);

	//Define each individual element
	//Define Nitrogen
	G4double A = 14.01 * g/mole;
	G4double Z = 7;
	G4Element* elN = new G4Element ("Nitrogen", "N", Z, A);

	//Define Oxygen
	A = 16.0 * g/mole;
	Z = 8;
	G4Element* elO = new G4Element ("Oxygen", "O", Z, A);

	//Define Hydrogen 
	A = 1.01 * g/mole;
	Z = 1;
	G4Element* elH = new G4Element ("Hydrogen", "H", Z, A);

	//Define Carbon
	A = 12.01 * g/mole;
	Z = 6;
	G4Element* elC = new G4Element ("Carbon", "C", Z, A);
	
	//Define Air   
	G4Material* Air = new G4Material("Air", 1.29*mg/cm3, 2);
	Air -> AddElement(elN, 70*perCent);
	Air -> AddElement(elO, 30*perCent);
 
	//Define PMMA (C502H8)
	// NIST reference 
	G4Material* PMMA = new G4Material("PMMA", 1.19*g/cm3, 3);
	PMMA -> AddElement(elC, 5);
	PMMA -> AddElement(elO, 2);
	PMMA -> AddElement(elH, 8);

	//define materials
	G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");
	G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	
	//Define Vacuum
	 G4double vacuumDensity = 1.e-25 *g/cm3;
	 G4double pressure = 3.e-18*pascal;
	 G4double temperature = 2.73*kelvin;
	 G4Material* vacuum = new G4Material("Galactic", Z=1., A=1.01*g/mole,
						 vacuumDensity,kStateGas,temperature,pressure);

	 //Define volumes
	 // World volume  has size 1cm
	 G4double worldx = 0.5 * m;  //half length!!!!
	 G4double worldy = 0.5 * m;
	 G4double worldz = 0.5 * m;

	 // World volume, containing all geometry
	 G4Box* world = new G4Box("world_box", worldx, worldy, worldz);

	 G4LogicalVolume* logical_world = new G4LogicalVolume(world, vacuum, "world_log", 0,0,0);

	 //set the logical world volume invisible
	 logical_world -> SetVisAttributes(G4VisAttributes::GetInvisible());

	 G4VPhysicalVolume* physical_world = new G4PVPlacement(0,
								   G4ThreeVector(),
								   logical_world, 
								"world_phys",
								0, 
								false, 
								0);
		
	// I need to set the size of the SV now, because some other parameters depend on it
	G4double SV_thick = 25.*um /2.; 
	
	// PMMA
	G4double PMMA_x = 200.*um /2.;
	G4double PMMA_y = 200.*um /2.;
	G4double PMMA_z = SV_thick;
	
	G4Box* PMMA_box = new G4Box("PMMA_box", PMMA_x, PMMA_y, PMMA_z);
	
	G4LogicalVolume* logical_PMMA = new G4LogicalVolume(PMMA_box, PMMA, "PMMA_log", 0,0,0);
	
	new G4PVPlacement(0, G4ThreeVector(), logical_PMMA, "PMMA_phys",
					logical_world,
					false, 0, true);
	
	logical_PMMA -> SetVisAttributes(G4VisAttributes(G4Colour(0., 1., 0.)));
	
	// sensitive volumes
	G4double SV_radius = SV_thick *2;	// full length
	
	G4Tubs* SV_cyl = new G4Tubs("SV_cyl", 0., SV_radius, SV_thick, 0.*deg, 360.*deg);
	//G4RotationMatrix* cylRot = new G4RotationMatrix;
	//cylRot->rotateY(M_PI/2.*rad);
		
	G4LogicalVolume* logical_SV = new G4LogicalVolume(SV_cyl, silicon, "SV_log", 0,0,0);
	
	G4VisAttributes SVcolour(G4Colour(0.5, 0.5, 0.5));
	SVcolour.SetForceSolid(true);
	logical_SV -> SetVisAttributes(SVcolour);
	
	G4ThreeVector SVposition;	//if(volumeName != "SV_phys1")
		
	G4double SVspacing = 10.*um;	// distance between the edges of two SV
	
	SVposition = { +SVspacing/2. +SV_radius, +SVspacing/2. +SV_radius, 0. };
	new G4PVPlacement(0, SVposition, logical_SV, "SV_phys1",
						logical_PMMA,
						false, 0, true);
	/*new G4PVPlacement(cylRot, SVposition, logical_SV, "SV_phys1",
						logical_PMMA,
						false, 0, true);*/

	SVposition = { -SVspacing/2. -SV_radius, +SVspacing/2. +SV_radius, 0. };
	new G4PVPlacement(0, SVposition, logical_SV, "SV_phys2",
						logical_PMMA,
						false, 0, true);
	
	SVposition = { -SVspacing/2. -SV_radius, -SVspacing/2. -SV_radius, 0. };
	new G4PVPlacement(0, SVposition, logical_SV, "SV_phys3",
						logical_PMMA,
						false, 0, true);
						
	SVposition = { +SVspacing/2. +SV_radius, -SVspacing/2. -SV_radius, 0. };
	new G4PVPlacement(0, SVposition, logical_SV, "SV_phys4",
						logical_PMMA,
						false, 0, true);
	
	// Si02 layer
	G4double oxyde_x = PMMA_x;
	G4double oxyde_y = PMMA_y;
	G4double oxyde_z = 1.*um /2.;
	
	G4Box* oxyde_box = new G4Box("oxyde_box", oxyde_x, oxyde_y, oxyde_z);
	
	G4LogicalVolume* logical_oxyde = new G4LogicalVolume(oxyde_box, SiO2, "oxyde_log", 0,0,0);
	
	G4ThreeVector oxyde_position = G4ThreeVector( 0, 0, PMMA_z + oxyde_z );
	new G4PVPlacement(0, oxyde_position, logical_oxyde, "oxyde_phys",
					logical_world,
					false, 0, true);
	
	logical_oxyde -> SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 0.6)));

	return physical_world; 

}

G4VPhysicalVolume* DetectorConstruction::ConstructSiliconBridgeDetector()
{
//--------------------------------------------------------
//--------------------- MATERIALS ------------------------
//--------------------------------------------------------
	G4NistManager * nistMan = G4NistManager::Instance();
	
	//Define Water   
	//G4Material* Water  = nistMan->FindOrBuildMaterial("G4_WATER");

	//Define Polyethylene
	//G4Material* poly  = nistMan->FindOrBuildMaterial("G4_POLYETHYLENE");
	
	//Define PMMA
	//G4Material* perspex  = nistMan->FindOrBuildMaterial("G4_PLEXIGLASS"); //Default 1.19g
	
	//Define Air
	G4Material* Air  = nistMan->FindOrBuildMaterial("G4_AIR");
	
	//Define Aluminium
	G4Material* Aluminium = nistMan->FindOrBuildMaterial("G4_Al");
	
	//Define Silicon
	G4Material* silicon = nistMan->FindOrBuildMaterial("G4_Si");
	
	//Define Aluminium Oxide (this is acting as the ceremic)
	//G4Material* AlOx = nistMan->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	
	//Define SiO
	G4Material* SiO2 = nistMan->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	
	//Define Pyrex Glass
	//G4Material* PyrexGlass = nistMan->FindOrBuildMaterial("G4_Pyrex_Glass");
	
//--------------------------------------------------------
//-------------------- Vis Attributes --------------------
//--------------------------------------------------------	
	G4VisAttributes* wireFrameWhiteAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	wireFrameWhiteAtt -> SetVisibility(true);
	wireFrameWhiteAtt -> SetForceWireframe(true);
	
	G4VisAttributes* wireFramePinkAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0)); 
   	wireFramePinkAtt -> SetVisibility(true);
   	wireFramePinkAtt -> SetForceWireframe(true);
	
	G4VisAttributes* solidGreyAtt = new G4VisAttributes(G4Colour(0.5, .5, .5)); 
   	solidGreyAtt -> SetVisibility(true);
   	solidGreyAtt -> SetForceSolid(true);
	
	G4VisAttributes* solidRedAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
	solidRedAtt -> SetVisibility(true);
	solidRedAtt -> SetForceSolid(true);
	
	G4VisAttributes* solidGreenAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
	solidGreenAtt -> SetVisibility(true);
	solidGreenAtt -> SetForceSolid(true);
	
	G4VisAttributes* solidYellowAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
	solidYellowAtt -> SetVisibility(true);
	solidYellowAtt -> SetForceSolid(true);
	
	G4VisAttributes* solidBlueAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
	solidBlueAtt -> SetVisibility(true);
	solidBlueAtt -> SetForceSolid(true);
	
//--------------------------------------------------------
//----------------------- Volumes ------------------------
//--------------------------------------------------------

//-------------------------World -------------------------
    G4double worldX = 10. * cm;
    G4double worldY = 10. * cm;
    G4double worldZ = 10. * cm;
    
	G4RotationMatrix* rotMatW = new G4RotationMatrix;
	rotMatW->rotateY(0*deg);
	rotMatW->rotateX(0*deg);
	rotMatW->rotateZ(0*deg);

    G4Box* world = new G4Box("world_box", worldX/2., worldY/2., worldZ/2.);
    G4LogicalVolume* logical_world = new G4LogicalVolume(world, Air, "world_log", 0,0,0);
    G4VPhysicalVolume* physical_world = new G4PVPlacement(0,G4ThreeVector(0,0,0),logical_world, "world_phys", 0, false, 0);
	
	logical_world -> SetVisAttributes(wireFrameWhiteAtt);

//------------------Detector Mother Volume------------------	
	G4double SVheight = 10.*micrometer;
	G4double insulationThickness = 1.*micrometer;
	G4double SiOoverlayerBottomThickness = 1.7*micrometer;
	G4double AlOverlayerThickness = 1.7*micrometer;
	G4double AlOverlayerRadius = 4.*micrometer;
	G4double SiOoverlayerTopThickness = 1.43*micrometer;
	G4double SiOoverlayerTopRadius = 10.5/2. * micrometer;
	G4double overLayerThickness = AlOverlayerThickness + SiOoverlayerTopThickness;

	G4double baseSiThickness = 300.*micrometer;
	G4double detectorHeight = (SVheight + baseSiThickness + insulationThickness + SiOoverlayerBottomThickness + AlOverlayerThickness + SiOoverlayerTopThickness);

	G4double SVwidth = 30.*micrometer;
	G4double pitch = 20.*micrometer; //distance between the odd and even rows

	G4double bridgingWidth = 20.*micrometer;
	G4double bridgingLength = 15.*micrometer;
	G4double bridingHeight = SVheight;

	G4double numberOfRows = 59;
	G4double numberOfColumns = (24*3);
	G4double SVareaWidth = (numberOfColumns * (SVwidth + bridgingWidth)) - 1.*bridgingWidth;
	G4double SVareaLength = numberOfRows*(SVwidth + pitch) - 1.*pitch;
	G4double bufferWidth = 100.*micrometer;
	G4double detectorWidth = SVareaWidth + bufferWidth;
	G4double detectorLength = SVareaLength + bufferWidth;
	
//------------------Detector Mother Volume------------------	
	G4Box* detectorBox = new G4Box("detectorBox", detectorWidth/2., detectorHeight/2., detectorLength/2.);
	G4LogicalVolume* logicalDetectorReg = new G4LogicalVolume(detectorBox, Air, "detectorReg_log", 0,0,0);
	new G4PVPlacement(0, G4ThreeVector(0,0,0), logicalDetectorReg, "detectorRegPhys", logical_world , false, 0, true);

	logicalDetectorReg -> SetVisAttributes(wireFrameWhiteAtt);
	
//------------------Base Silicon Volume----------------------
	G4Box* basSiBox = new G4Box("baseSiBox", detectorWidth/2., baseSiThickness/2., detectorLength/2.);
	G4LogicalVolume* logicalbaseSi = new G4LogicalVolume(basSiBox, silicon, "baseSi_log", 0,0,0);
	new G4PVPlacement(0, G4ThreeVector(0,(-detectorHeight/2. + baseSiThickness/2.),0), logicalbaseSi, "baseSiPhys", logicalDetectorReg, false, 0, true);
	
	logicalbaseSi -> SetVisAttributes(wireFramePinkAtt);

//---------------------Insulation Volume----------------------
	G4Box* SiOBox = new G4Box("SiOBox", detectorWidth/2., insulationThickness/2., detectorLength/2.);
	G4LogicalVolume* logicalSoI = new G4LogicalVolume(SiOBox, SiO2, "SoI_log", 0,0,0);
	new G4PVPlacement(0, G4ThreeVector(0,(detectorHeight/2. -overLayerThickness - SVheight - insulationThickness/2.),0), logicalSoI, "SoIPhys", logicalDetectorReg, false, 0, true);

	logicalSoI -> SetVisAttributes(solidGreyAtt);

//----------------Sensitive Volume Region---------------------
//This volume encapsulates the SVs and the "bridging" volumes
	G4Box* SVregBox = new G4Box("SVregBox", detectorWidth/2., SVheight/2., detectorLength/2.);
	G4LogicalVolume* logicalSVreg = new G4LogicalVolume(SVregBox, Air, "SVreg_log", 0,0,0);
	new G4PVPlacement(0, G4ThreeVector(0,(detectorHeight/2. - overLayerThickness - SVheight/2.),0), logicalSVreg, "SVregPhys", logicalDetectorReg, false, 0, true);

    logicalSVreg -> SetVisAttributes(wireFrameWhiteAtt);

//----------------------Bridging Volume----------------------
//This volume connects or "bridges" the main sensitive volumes 
//together but effectively act as additional sensitive volumes
	G4Box* bridgeVolBox = new G4Box("bridgeVolBox", bridgingWidth/2., bridingHeight/2., bridgingLength/2.);
	G4LogicalVolume* logicalBridgeVol = new G4LogicalVolume(bridgeVolBox, silicon, "bridgeVol_log", 0,0,0);

	logicalBridgeVol -> SetVisAttributes(solidRedAtt);
	
//---------------------SiO Bottom Overlayer------------------
	G4Box* SiObotLayerBox = new G4Box("SiObotLayerBox", SVwidth/2., SiOoverlayerBottomThickness/2., SVwidth/2.);
	G4LogicalVolume* logicalSiObotLayer = new G4LogicalVolume(SiObotLayerBox, SiO2, "logicalSiObotLayer", 0,0,0);

	logicalSiObotLayer -> SetVisAttributes(solidGreenAtt);
	
//--------------SiO Bottom Overlayer Bridging Vol-------------
	G4Box* SiObotLayerBridgeBox = new G4Box("SiObotLayerBridgeBox", bridgingWidth/2., SiOoverlayerBottomThickness/2., bridgingLength/2.);
	G4LogicalVolume* logicalSiObotLayerBridge = new G4LogicalVolume(SiObotLayerBridgeBox, SiO2, "logicalSiObotLayer", 0,0,0);

	logicalSiObotLayerBridge -> SetVisAttributes(solidGreenAtt);

//-----------------------Al contact ---------------------------
	G4Box* AlContactBox = new G4Box("AlContactBox", AlOverlayerRadius, AlOverlayerThickness/2., AlOverlayerRadius);
	G4LogicalVolume* logicalAlContact = new G4LogicalVolume(AlContactBox, Aluminium, "logicalAlContact", 0,0,0);

	logicalAlContact -> SetVisAttributes(solidYellowAtt);

//----------------------SiO Top Overlayer----------------------	
	G4Box* SiOtopLayerBox = new G4Box("SiOtopLayerBox", SiOoverlayerTopRadius, SiOoverlayerTopThickness/2., SiOoverlayerTopRadius);
	G4LogicalVolume* logicalSiOtopLayer = new G4LogicalVolume(SiOtopLayerBox, SiO2, "logicalSiOtopLayer", 0,0,0);

	logicalSiOtopLayer -> SetVisAttributes(solidBlueAtt);

//--------------------Sensitive Volume--------------------------
	G4Box* sensitiveBridgeVolume = new G4Box("water_region_sphere", SVwidth/2., SVheight/2., SVwidth/2.);
	G4LogicalVolume* SV_log = new G4LogicalVolume(sensitiveBridgeVolume, silicon, "SV_log", 0,0,0);

	SV_log -> SetVisAttributes(solidYellowAtt);
	
//---Placing SVs, 30x30 volumes and "bridging" volumes between them

	G4bool checkSVoverlap = false;
	
	//Placing "odd" SV rows/columns
	G4int globalCount = 100000;  
	
	for (G4double j = -SVareaLength/2.; j <= SVareaLength/2.; j+=2.*(SVwidth + bridgingWidth)) //creates each row
	{
		for (G4double i = -SVareaWidth/2.; i <= SVareaWidth/2.; i+=(SVwidth + bridgingWidth)) //goes through each odd row
		{
			//G4cout << "x: " << i/um << " z: " << j/um << G4endl;
			new G4PVPlacement(0, G4ThreeVector(i,0,j), SV_log , "physSensitiveBridgeVolume", logicalSVreg, false, globalCount, checkSVoverlap);
			
			new G4PVPlacement(0, G4ThreeVector(i,((detectorHeight/2. -overLayerThickness + SiOoverlayerBottomThickness/2.)),j), logicalSiObotLayer, "physSiObotLayer", logicalDetectorReg, false, globalCount, checkSVoverlap);
			
			new G4PVPlacement(0, G4ThreeVector(i,((detectorHeight/2. -overLayerThickness + SiOoverlayerBottomThickness + SiOoverlayerTopThickness/2.)),j), logicalSiOtopLayer, "physSiOtopLayer", logicalDetectorReg, false, globalCount, checkSVoverlap);
			
			globalCount++;	
		}
	}
	
	//Placing "even" SV rows/columns
	globalCount = 200000;  
	
	for (G4double j = (-SVareaLength/2. + SVwidth + bridgingWidth); j <= SVareaLength/2.; j+=2.*(SVwidth + bridgingWidth)) //creates each row
	{
		for (G4double i = -SVareaWidth/2.; i <= SVareaWidth/2.; i+=(SVwidth + bridgingWidth)) //goes through each odd row
		{
			//G4cout << "x: " << i/micrometer << " z: " << j/micrometer << G4endl;

			new G4PVPlacement(0, G4ThreeVector(i,0,j), SV_log , "physSensitiveBridgeVolume", logicalSVreg, false, globalCount, checkSVoverlap);
			
			new G4PVPlacement(0, G4ThreeVector(i,((detectorHeight/2. -overLayerThickness + SiOoverlayerBottomThickness/2.)),j), logicalSiObotLayer, "physSiObotLayer", logicalDetectorReg, false, globalCount, checkSVoverlap);
			
			new G4PVPlacement(0, G4ThreeVector(i,((detectorHeight/2. -overLayerThickness + SiOoverlayerBottomThickness + SiOoverlayerTopThickness/2.)),j), logicalSiOtopLayer, "physSiOtopLayer", logicalDetectorReg, false, globalCount, checkSVoverlap);
			
			globalCount++;
		}
	}

	//Placing "odd" bridging volumes rows/columns
	globalCount = 300000;

	for (G4double j = -SVareaLength/2.; j < SVareaLength/2.; j+=2.*(SVwidth + bridgingWidth)) //creates each row
	{
		for (G4double i = -SVareaWidth/2.; i < SVareaWidth/2.; i+=(SVwidth + bridgingWidth)) //goes through each odd row
		{
		if ( (i + (SVwidth + bridgingWidth)) < SVareaWidth/2.)
		{
			//G4cout << globalCount << G4endl;
			new G4PVPlacement(0, G4ThreeVector((i + SVwidth/2. + bridgingWidth/2.),0,(j+ bridgingLength/2.)), logicalBridgeVol , "sen_bridge", logicalSVreg, false, globalCount, checkSVoverlap);
			
			new G4PVPlacement(0, G4ThreeVector((i + SVwidth/2. + bridgingWidth/2.),((detectorHeight/2. -overLayerThickness + SiOoverlayerBottomThickness/2.)),(j+ bridgingLength/2.)), logicalSiObotLayerBridge, "physSiObotLayerBridge", logicalDetectorReg, false, globalCount, checkSVoverlap);
			
			//G4cout << "x: " << (i + SVwidth/2. + bridgingWidth/2.)/micrometer << " z: " << j/micrometer << G4endl;
			
			globalCount++;
		}	
		}
	}

	//Placing "even" bridging rows/columns
	globalCount = 400000;  

	for (G4double j = (-SVareaLength/2. + SVwidth + bridgingWidth); j < SVareaLength/2.; j+=2.*(SVwidth + bridgingWidth)) //creates each row
	{
		for (G4double i = -SVareaWidth/2.; i < SVareaWidth/2.; i+=(SVwidth + bridgingWidth)) //goes through each odd row
		{
		
		if ( (i+ (SVwidth + bridgingWidth)) < SVareaWidth/2.)
		{	
			//G4cout << globalCount << G4endl;
			new G4PVPlacement(0, G4ThreeVector((i + SVwidth/2. + bridgingWidth/2.),0,(j+ bridgingLength/2.)), logicalBridgeVol , "sen_bridge", logicalSVreg, false, globalCount, checkSVoverlap);
			
			//G4cout << "x: " << (i + SVwidth/2. + bridgingWidth/2.)/micrometer << " z: " << j/micrometer << G4endl;
			
			new G4PVPlacement(0, G4ThreeVector((i + SVwidth/2. + bridgingWidth/2.),((detectorHeight/2. -overLayerThickness + SiOoverlayerBottomThickness/2.)),(j+ bridgingLength/2.)), logicalSiObotLayerBridge, "physSiObotLayerBridge", logicalDetectorReg, false, globalCount, checkSVoverlap);
			
			globalCount++;
		}
			
		}
	}
	
//----------------------------------------
	new G4PVPlacement(0, G4ThreeVector(0,0,0), logicalAlContact, "physAlContact", logicalSiObotLayer, false, 0, 1);



	return physical_world; 
}


void DetectorConstruction::ConstructSDandField()
{
   SensitiveDetector* SD = new SensitiveDetector("SD", "DetectorHitsCollection", analysis);
   G4SDManager::GetSDMpointer()->AddNewDetector(SD);
   SetSensitiveDetector("SV_log", SD);

	if (detectorType == "SiliconBridge")
	{
		SetSensitiveDetector("bridgeVol_log", SD);
	}
}
