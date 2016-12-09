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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//

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

DetectorConstruction::DetectorConstruction(AnalysisManager* analysis_manager)
{
analysis = analysis_manager;
}

DetectorConstruction::~DetectorConstruction(){

}

G4VPhysicalVolume* DetectorConstruction::Construct()
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

void DetectorConstruction::ConstructSDandField()
{
   SensitiveDetector* SD = new SensitiveDetector("SD", "DetectorHitsCollection", analysis);
   G4SDManager::GetSDMpointer()->AddNewDetector(SD);
   SetSensitiveDetector("SV_log", SD);


}
