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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "HadrontherapyModulator.hh"
#include "G4Transform3D.hh"
#include "G4ios.hh"
#include <fstream>
#include "G4RunManager.hh"
#include "G4NistManager.hh"

HadrontherapyModulator::HadrontherapyModulator():physiMotherMod(0),
						 solidMod0(0),         logicMod0(0),          physiMod0(0),
						 solidMod1(0),         logicMod1(0),          physiMod1(0),
						 solidMod2(0),         logicMod2(0),          physiMod2(0),
						 solidMod3(0),         logicMod3(0),          physiMod3(0),
						 solidMod4(0),         logicMod4(0),          physiMod4(0),
						 solidMod5(0),         logicMod5(0),          physiMod5(0),
						 solidMod6(0),         logicMod6(0),          physiMod6(0),
						 solidMod7(0),         logicMod7(0),          physiMod7(0),
						 solidMod8(0),         logicMod8(0),          physiMod8(0),
						 solidMod9(0),         logicMod9(0),          physiMod9(0),
						 solidMod10(0),        logicMod10(0),         physiMod10(0),
						 solidMod11(0),        logicMod11(0),         physiMod11(0),
						 solidMod12(0),        logicMod12(0),         physiMod12(0),
						 solidMod13(0),        logicMod13(0),         physiMod13(0),
						 solidMod14(0),        logicMod14(0),         physiMod14(0),
						 solidMod15(0),        logicMod15(0),         physiMod15(0),
						 solidMod16(0),        logicMod16(0),         physiMod16(0),
						 solidMod17(0),        logicMod17(0),         physiMod17(0),
						 solidMod18(0),        logicMod18(0),         physiMod18(0),
						 solidMod20(0),        logicMod20(0),         physiMod20(0),
						 solidMod21(0),        logicMod21(0),         physiMod21(0),
						 solidMod22(0),        logicMod22(0),         physiMod22(0),
						 solidMod23(0),        logicMod23(0),         physiMod23(0),
						 solidMod24(0),        logicMod24(0),         physiMod24(0),
						 solidMod25(0),        logicMod25(0),         physiMod25(0),
						 solidMod26(0),        logicMod26(0),         physiMod26(0),
						 solidMod27(0),        logicMod27(0),         physiMod27(0),
						 solidMod28(0),        logicMod28(0),         physiMod28(0),
						 solidMod29(0),        logicMod29(0),         physiMod29(0),
						 solidMod30(0),        logicMod30(0),         physiMod30(0),
						 solidMod31(0),        logicMod31(0),         physiMod31(0),
						 solidMod32(0),        logicMod32(0),         physiMod32(0),
						 solidMod33(0),        logicMod33(0),         physiMod33(0),
						 solidMod34(0),        logicMod34(0),         physiMod34(0),
						 solidMod35(0),        logicMod35(0),         physiMod35(0),
						 solidMod36(0),        logicMod36(0),         physiMod36(0),
						 solidMod37(0),        logicMod37(0),         physiMod37(0),
						 solidMod38(0),        logicMod38(0),         physiMod38(0),
						 solidMod40(0),        logicMod40(0),         physiMod40(0),
						 solidMod41(0),        logicMod41(0),         physiMod41(0),
						 solidMod42(0),        logicMod42(0),         physiMod42(0),
						 solidMod43(0),        logicMod43(0),         physiMod43(0),
						 solidMod44(0),        logicMod44(0),         physiMod44(0),
						 solidMod45(0),        logicMod45(0),         physiMod45(0),
						 solidMod46(0),        logicMod46(0),         physiMod46(0),
						 solidMod47(0),        logicMod47(0),         physiMod47(0),
						 solidMod48(0),        logicMod48(0),         physiMod48(0),
						 solidMod49(0),        logicMod49(0),         physiMod49(0),
						 solidMod50(0),        logicMod50(0),         physiMod50(0),
						 solidMod51(0),        logicMod51(0),         physiMod51(0),
						 solidMod52(0),        logicMod52(0),         physiMod52(0),
						 solidMod53(0),        logicMod53(0),         physiMod53(0),
						 solidMod54(0),        logicMod54(0),         physiMod54(0),
						 solidMod55(0),        logicMod55(0),         physiMod55(0),
						 solidMod56(0),        logicMod56(0),         physiMod56(0),
						 solidMod57(0),        logicMod57(0),         physiMod57(0),
						 solidMod58(0),        logicMod58(0),         physiMod58(0),
						 solidMod60(0),        logicMod60(0),         physiMod60(0),
						 solidMod61(0),        logicMod61(0),         physiMod61(0),
						 solidMod62(0),        logicMod62(0),         physiMod62(0),
						 solidMod63(0),        logicMod63(0),         physiMod63(0),
						 solidMod64(0),        logicMod64(0),         physiMod64(0),
						 solidMod65(0),        logicMod65(0),         physiMod65(0),
						 solidMod66(0),        logicMod66(0),         physiMod66(0),
						 solidMod67(0),        logicMod67(0),         physiMod67(0),
						 solidMod68(0),        logicMod68(0),         physiMod68(0),
						 solidMod69(0),        logicMod69(0),         physiMod69(0),
						 solidMod70(0),        logicMod70(0),         physiMod70(0),
						 solidMod71(0),        logicMod71(0),         physiMod71(0),
						 solidMod72(0),        logicMod72(0),         physiMod72(0),
						 solidMod73(0),        logicMod73(0),         physiMod73(0),
						 solidMod74(0),        logicMod74(0),         physiMod74(0),
						 solidMod75(0),        logicMod75(0),         physiMod75(0),
						 solidMod76(0),        logicMod76(0),         physiMod76(0),
						 solidMod77(0),        logicMod77(0),         physiMod77(0),
						 solidMod78(0),        logicMod78(0),         physiMod78(0) 
{
  rm = new G4RotationMatrix(); 
  G4double phi = 270. *deg;     
  rm -> rotateY(phi); 
}
/////////////////////////////////////////////////////////////////////////////
HadrontherapyModulator::~HadrontherapyModulator() 
{
  delete rm;
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyModulator::BuildModulator(G4VPhysicalVolume* motherVolume)
{
  G4bool isotopes = false;
  G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);

  // You have to uncomment the following line if you want to define a PMMA material
  //  G4Material* PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);

  G4Material* Mod0Mater = airNist;
  G4Material* ModMater = airNist; // You have to set ModMater to PMMANist if you want to change modulator material (default is air)
 
  G4double innerRadiusOfTheTube = 2.5 *cm;
  G4double outerRadiusOfTheTube = 9.5 *cm;
  G4double hightOfTheTube = 0.03*cm;

  // Mother of the modulator wheel  
  G4ThreeVector positionMotherMod = G4ThreeVector(-2160.50 *mm, 30 *mm, 50 *mm);
 
  G4Box* solidMotherMod = new G4Box("MotherMod", 12 *cm, 12 *cm, 12 *cm);
 
  G4LogicalVolume * logicMotherMod = new G4LogicalVolume(solidMotherMod, Mod0Mater,"MotherMod",0,0,0);

  physiMotherMod = new G4PVPlacement(rm,positionMotherMod,  "MotherMod", 
				     logicMotherMod,    				  
				     motherVolume,      
				     false,           
				     0);              
 
  //----------------------------------------------------------
  // Mother volume of first quarter of the modulator
  //----------------------------------------------------------

  G4double hightOfTheTube0 = 5.0 *cm;
  G4double startAngleOfTheTube0 = 45 *deg;
  G4double spanningAngleOfTheTube0 = 90 *deg;
  
  G4RotationMatrix rm2;
  rm2.rotateZ(0 *deg);
 
  G4ThreeVector positionMod0 = G4ThreeVector(0*cm,0*cm,0*cm);
  
  solidMod0 = new G4Tubs("Mod0",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube0,
			 startAngleOfTheTube0, 
			 spanningAngleOfTheTube0);
  
  logicMod0 = new G4LogicalVolume(solidMod0, Mod0Mater, "Mod0",0,0,0);
  
  physiMod0 = new G4PVPlacement(G4Transform3D(rm2, positionMod0), 
				logicMod0,    
				"Mod0",       
				logicMotherMod,  
				false,         
				0);            
  
 
  //----------------------------------------------------------
  // First modulator sclice
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube1 = 54.267*deg;
  G4double spanningAngleOfTheTube1 = 71.466*deg;

 
  G4ThreeVector positionMod1 = G4ThreeVector(0*cm,0*cm,0.51*cm);
  solidMod1 = new G4Tubs("Mod1",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube,
			 startAngleOfTheTube1, 
			 spanningAngleOfTheTube1);
			       
  logicMod1 = new G4LogicalVolume(solidMod1, ModMater, "Mod1",0,0,0);
  physiMod1 = new G4PVPlacement(0,               
				positionMod1,  
				logicMod1,     
				"Mod1",        
				logicMod0,     
				false,         
				0);            
				      
				      
  //----------------------------------------------------------
  // Second modulator slice
  //----------------------------------------------------------
  
  G4double startAngleOfTheTube2 = 57.714*deg;
  G4double spanningAngleOfTheTube2 = 64.572*deg;

   
   
  G4ThreeVector positionMod2 = G4ThreeVector(0*cm,0*cm,0.45*cm);
  solidMod2 = new G4Tubs("Mod2",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube,
			 startAngleOfTheTube2, 
			 spanningAngleOfTheTube2);

  logicMod2 = new G4LogicalVolume(solidMod2, ModMater, "Mod2",0,0,0);
  physiMod2 = new G4PVPlacement(0,             
				positionMod2,  
				logicMod2,     
				"Mod2",        
				logicMod0,     
				false,         
				0);            


  //----------------------------------------------------------
  // 3th modulator slice
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube3 = 60.478*deg;
  G4double spanningAngleOfTheTube3 = 59.044*deg;

   
   
  G4ThreeVector positionMod3 = G4ThreeVector(0*cm,0*cm,0.39*cm);
  solidMod3 = new G4Tubs("Mod3",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube,
			 startAngleOfTheTube3, 
			 spanningAngleOfTheTube3);

  logicMod3 = new G4LogicalVolume(solidMod3, ModMater, "Mod3",0,0,0);
  physiMod3 = new G4PVPlacement(0,                  
				positionMod3,  
				logicMod3,     
				"Mod3",        
				logicMod0,     
				false,         
				0);            

 
  //----------------------------------------------------------
  //
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube4 = 62.668*deg;
  G4double spanningAngleOfTheTube4 = 54.664*deg;
   
   
  G4ThreeVector positionMod4 = G4ThreeVector(0*cm,0*cm,0.33*cm);
  solidMod4 = new G4Tubs("Mod4",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube,
			 startAngleOfTheTube4, 
			 spanningAngleOfTheTube4);

  logicMod4 = new G4LogicalVolume(solidMod4, ModMater, "Mod4",0,0,0);
  physiMod4 = new G4PVPlacement(0,                   // no rotation
				positionMod4,  // at (x,y,z)
				logicMod4,     // its logical volume				  
				"Mod4",        // its name
				logicMod0,      // its mother  volume
				false,           // no boolean operations
				0);              // no particular field


  //----------------------------------------------------------
  //Quinta fetta Modulatore
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube5 = 64.814*deg;
  G4double spanningAngleOfTheTube5 = 50.372*deg;
   
   
  G4ThreeVector positionMod5 = G4ThreeVector(0*cm,0*cm,0.27*cm);
  solidMod5 = new G4Tubs("Mod5",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube,
			 startAngleOfTheTube5, 
			 spanningAngleOfTheTube5);

  logicMod5 = new G4LogicalVolume(solidMod5, ModMater, "Mod5",0,0,0);
  physiMod5 = new G4PVPlacement(0,                   // no rotation
				positionMod5,  // at (x,y,z)
				logicMod5,     // its logical volume				  
				"Mod5",        // its name
				logicMod0,      // its mother  volume
				false,           // no boolean operations
				0);              // no particular field

  
  //----------------------------------------------------------
  //Sesta fetta Modulatore
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube6 = 66.706*deg;
  G4double spanningAngleOfTheTube6 = 46.588*deg;
   
   
  G4ThreeVector positionMod6 = G4ThreeVector(0*cm,0*cm,0.21*cm);
  solidMod6 = new G4Tubs("Mod6",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube,
			 startAngleOfTheTube6, 
			 spanningAngleOfTheTube6);

  logicMod6 = new G4LogicalVolume(solidMod6, ModMater, "Mod6",0,0,0);
  physiMod6 = new G4PVPlacement(0,                   // no rotation
				positionMod6,  // at (x,y,z)
				logicMod6,     // its logical volume				  
				"Mod6",        // its name
				logicMod0,      // its mother  volume
				false,           // no boolean operations
				0);              // no particular field


  //----------------------------------------------------------
  //Settima fetta Modulatore
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube7 = 68.648*deg;
  G4double spanningAngleOfTheTube7 = 42.704*deg;

   
   
  G4ThreeVector positionMod7 = G4ThreeVector(0*cm,0*cm,0.15*cm);
  solidMod7 = new G4Tubs("Mod7",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube,
			 startAngleOfTheTube7, 
			 spanningAngleOfTheTube7);

  logicMod7 = new G4LogicalVolume(solidMod7, ModMater, "Mod7",0,0,0);
  physiMod7 = new G4PVPlacement(0,                   // no rotation
				positionMod7,  // at (x,y,z)
				logicMod7,     // its logical volume				  
				"Mod7",        // its name
				logicMod0,      // its mother  volume
				false,           // no boolean operations
				0);              // no particular field



  //----------------------------------------------------------
  //Ottava fetta Modulatore
  //----------------------------------------------------------
   
  G4double startAngleOfTheTube8 = 70.472*deg;
  G4double spanningAngleOfTheTube8 = 39.056*deg;

   
  G4ThreeVector positionMod8 = G4ThreeVector(0*cm,0*cm,0.09*cm);
  solidMod8 = new G4Tubs("Mod8",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube,
			 startAngleOfTheTube8, 
			 spanningAngleOfTheTube8);

  logicMod8 = new G4LogicalVolume(solidMod8, ModMater, "Mod8",0,0,0);
  physiMod8 = new G4PVPlacement(0,                   // no rotation
				positionMod8,  // at (x,y,z)
				logicMod8,     // its logical volume				  
				"Mod8",        // its name
				logicMod0,      // its mother  volume
				false,           // no boolean operations
				0);              // no particular field




  //----------------------------------------------------------
  //Nona fetta Modulatore
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube9 = 72.288*deg;
  G4double spanningAngleOfTheTube9 = 35.424*deg;

   
  G4ThreeVector positionMod9 = G4ThreeVector(0*cm,0*cm,0.03*cm);
  solidMod9 = new G4Tubs("Mod9",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube,
			 startAngleOfTheTube9, 
			 spanningAngleOfTheTube9);

  logicMod9 = new G4LogicalVolume(solidMod9, ModMater, "Mod9",0,0,0);
  physiMod9 = new G4PVPlacement(0,                   // no rotation
				positionMod9,  // at (x,y,z)
				logicMod9,     // its logical volume				  
				"Mod9",        // its name
				logicMod0,      // its mother  volume
				false,           // no boolean operations
				0);              // no particular field


  //----------------------------------------------------------
  //Decima fetta Modulatore
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube10 = 74.061*deg;
  G4double spanningAngleOfTheTube10 = 31.878*deg;

  
  G4ThreeVector positionMod10 = G4ThreeVector(0*cm,0*cm,-0.03*cm);
  solidMod10 = new G4Tubs("Mod10",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube10, 
			  spanningAngleOfTheTube10);

  logicMod10 = new G4LogicalVolume(solidMod10, ModMater, "Mod10",0,0,0);
  physiMod10 = new G4PVPlacement(0,                   // no rotation
				 positionMod10,  // at (x,y,z)
				 logicMod10,     // its logical volume				  
				 "Mod10",        // its name
				 logicMod0,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  // Undicesima fetta Modulatore
  //----------------------------------------------------------

  G4double startAngleOfTheTube11 = 75.793*deg;
  G4double spanningAngleOfTheTube11 = 28.414*deg;
   
   
  G4ThreeVector positionMod11 = G4ThreeVector(0*cm,0*cm,-0.09*cm);
  solidMod11 = new G4Tubs("Mod11",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube11, 
			  spanningAngleOfTheTube11);

  logicMod11 = new G4LogicalVolume(solidMod11, ModMater, "Mod11",0,0,0);
  physiMod11 = new G4PVPlacement(0,                   // no rotation
				 positionMod11,  // at (x,y,z)
				 logicMod11,     // its logical volume				  
				 "Mod11",        // its name
				 logicMod0,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Dodicesima fetta Modulatore
  //----------------------------------------------------------

  G4double startAngleOfTheTube12 = 77.579*deg;
  G4double spanningAngleOfTheTube12 = 24.842*deg;

   
   
  G4ThreeVector positionMod12 = G4ThreeVector(0*cm,0*cm,-0.15*cm);
  solidMod12 = new G4Tubs("Mod12",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube12, 
			  spanningAngleOfTheTube12);

  logicMod12 = new G4LogicalVolume(solidMod12, ModMater, "Mod12",0,0,0);
  physiMod12 = new G4PVPlacement(0,                   // no rotation
				 positionMod12,  // at (x,y,z)
				 logicMod12,     // its logical volume				  
				 "Mod12",        // its name
				 logicMod0,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  

  //----------------------------------------------------------
  //Tredicesima fetta Modulatore
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube13 = 79.273*deg;
  G4double spanningAngleOfTheTube13 = 21.454*deg;
   
  
  G4ThreeVector positionMod13 = G4ThreeVector(0*cm,0*cm,-0.21*cm);
  solidMod13 = new G4Tubs("Mod13",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube13, 
			  spanningAngleOfTheTube13);

  logicMod13 = new G4LogicalVolume(solidMod13, ModMater, "Mod13",0,0,0);
  physiMod13 = new G4PVPlacement(0,                   // no rotation
				 positionMod13,  // at (x,y,z)
				 logicMod13,     // its logical volume				  
				 "Mod13",        // its name
				 logicMod0,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Quat. fetta Modulatore
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube14 = 81.014*deg;
  G4double spanningAngleOfTheTube14 = 17.972*deg;

      
  G4ThreeVector positionMod14 = G4ThreeVector(0*cm,0*cm,-0.27*cm);
  solidMod14 = new G4Tubs("Mod14",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube14, 
			  spanningAngleOfTheTube14);

  logicMod14 = new G4LogicalVolume(solidMod14, ModMater, "Mod14",0,0,0);
  physiMod14 = new G4PVPlacement(0,                   // no rotation
				 positionMod14,  // at (x,y,z)
				 logicMod14,     // its logical volume				  
				 "Mod14",        // its name
				 logicMod0,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  
  //----------------------------------------------------------
  //Quindicesima fetta Modulatore
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube15 = 82.695*deg;
  G4double spanningAngleOfTheTube15 = 14.61*deg;

      
  G4ThreeVector positionMod15 = G4ThreeVector(0*cm,0*cm,-0.33*cm);
  solidMod15 = new G4Tubs("Mod15",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube15, 
			  spanningAngleOfTheTube15);

  logicMod15 = new G4LogicalVolume(solidMod15, ModMater, "Mod15",0,0,0);
  physiMod15 = new G4PVPlacement(0,                   // no rotation
				 positionMod15,  // at (x,y,z)
				 logicMod15,     // its logical volume				  
				 "Mod15",        // its name
				 logicMod0,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Sedic. fetta Modulatore
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube16 = 84.425*deg;
  G4double spanningAngleOfTheTube16 = 11.15*deg;

   
  G4ThreeVector positionMod16 = G4ThreeVector(0*cm,0*cm,-0.39*cm);
  solidMod16 = new G4Tubs("Mod16",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube16, 
			  spanningAngleOfTheTube16);

  logicMod16 = new G4LogicalVolume(solidMod16, ModMater, "Mod16",0,0,0);
  physiMod16 = new G4PVPlacement(0,                   // no rotation
				 positionMod16,  // at (x,y,z)
				 logicMod16,     // its logical volume				  
				 "Mod16",        // its name
				 logicMod0,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Dicias. fetta Modulatore
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube17 = 86.203*deg;
  G4double spanningAngleOfTheTube17 = 7.594*deg;

   
   
  G4ThreeVector positionMod17 = G4ThreeVector(0*cm,0*cm,-0.45*cm);
  solidMod17 = new G4Tubs("Mod17",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube17, 
			  spanningAngleOfTheTube17);

  logicMod17 = new G4LogicalVolume(solidMod17, ModMater, "Mod17",0,0,0);
  physiMod17 = new G4PVPlacement(0,                   // no rotation
				 positionMod17,  // at (x,y,z)
				 logicMod17,     // its logical volume				  
				 "Mod17",        // its name
				 logicMod0,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  //----------------------------------------------------------
  //Diciot. fetta Modulatore
  //----------------------------------------------------------
 
  G4double startAngleOfTheTube18 = 87.910*deg;
  G4double spanningAngleOfTheTube18 = 4.18*deg;

   
   
  G4ThreeVector positionMod18 = G4ThreeVector(0*cm,0*cm,-0.51*cm);
  solidMod18 = new G4Tubs("Mod18",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube18, 
			  spanningAngleOfTheTube18);

  logicMod18 = new G4LogicalVolume(solidMod18, ModMater, "Mod18",0,0,0);
  physiMod18 = new G4PVPlacement(0,                   // no rotation
				 positionMod18,  // at (x,y,z)
				 logicMod18,     // its logical volume				  
				 "Mod18",        // its name
				 logicMod0,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field



  //----------------------------------------------------------
  // Mother volume of the second modulator quarter
  //----------------------------------------------------------
  
    
  G4RotationMatrix rm20;
  rm20.rotateZ(90 *deg);
  
  G4ThreeVector positionMod20 = G4ThreeVector(0*cm,0*cm,0*cm);
  
  solidMod20 = new G4Tubs("Mod20",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube0,
			  startAngleOfTheTube0, 
			  spanningAngleOfTheTube0);
  
  logicMod20 = new G4LogicalVolume(solidMod20, Mod0Mater, "Mod0",0,0,0);
  
  
  physiMod20 = new G4PVPlacement(G4Transform3D(rm20, positionMod20), 
				 logicMod20,    
				 "Mod20",       
				 logicMotherMod, 
				 false,         
				 0);            
  

    

  //----------------------------------------------------------
  // 1st modulator slice (2nd quarter)
  //----------------------------------------------------------
 
  G4ThreeVector positionMod21 = G4ThreeVector(0*cm,0*cm,0.51*cm);
  solidMod21 = new G4Tubs("Mod21",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube1, 
			  spanningAngleOfTheTube1);
  
  logicMod21 = new G4LogicalVolume(solidMod21, ModMater, "Mod21",0,0,0);
  
  physiMod21 = new G4PVPlacement(0,               // no rotation
				 positionMod21,  // at (x,y,z)
				 logicMod21,     // its logical volume				  
				 "Mod21",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  
  
  //----------------------------------------------------------
  // 2nd modulator slice (2nd quarter)
  //----------------------------------------------------------
  
  G4ThreeVector positionMod22 = G4ThreeVector(0*cm,0*cm,0.45*cm);
  
  solidMod22 = new G4Tubs("Mod22",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube2, 
			  spanningAngleOfTheTube2);
  
  logicMod22 = new G4LogicalVolume(solidMod22, ModMater, "Mod22",0,0,0);
  
  physiMod22 = new G4PVPlacement(0,               // no rotation
				 positionMod22,  // at (x,y,z)
				 logicMod22,     // its logical volume				  
				 "Mod22",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  
  
  //----------------------------------------------------------
  // 3rd modulator slice (2nd quarter)
  //----------------------------------------------------------
 
  G4ThreeVector positionMod23 = G4ThreeVector(0*cm,0*cm,0.39*cm);

  solidMod23 = new G4Tubs("Mod23",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube3, 
			  spanningAngleOfTheTube3);

  logicMod23 = new G4LogicalVolume(solidMod23, ModMater, "Mod23",0,0,0);
  physiMod23 = new G4PVPlacement(0,                   // no rotation
				 positionMod23,  // at (x,y,z)
				 logicMod23,     // its logical volume				  
				 "Mod23",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

 
  //----------------------------------------------------------
  // 4th modulator slice (2nd quarter)
  //----------------------------------------------------------
   
   
  G4ThreeVector positionMod24 = G4ThreeVector(0*cm,0*cm,0.33*cm);

  solidMod24 = new G4Tubs("Mod24",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube4, 
			  spanningAngleOfTheTube4);

  logicMod24 = new G4LogicalVolume(solidMod24, ModMater, "Mod24",0,0,0);
 
  physiMod24 = new G4PVPlacement(0,                   // no rotation
				 positionMod24,  // at (x,y,z)
				 logicMod24,     // its logical volume				  
				 "Mod24",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field


  //----------------------------------------------------------
  //Quinta fetta Modulatore 2
  //----------------------------------------------------------
 
   
   
  G4ThreeVector positionMod25 = G4ThreeVector(0*cm,0*cm,0.27*cm);

  solidMod25 = new G4Tubs("Mod25",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube5, 
			  spanningAngleOfTheTube5);
  
  logicMod25 = new G4LogicalVolume(solidMod25, ModMater, "Mod25",0,0,0);
  physiMod25 = new G4PVPlacement(0,                   // no rotation
				 positionMod25,  // at (x,y,z)
				 logicMod25,     // its logical volume				  
				 "Mod25",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  
  //----------------------------------------------------------
  //Sesta fetta Modulatore 2
  //----------------------------------------------------------
 
   
  G4ThreeVector positionMod26 = G4ThreeVector(0*cm,0*cm,0.21*cm);
  solidMod26 = new G4Tubs("Mod26",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube6, 
			  spanningAngleOfTheTube6);

  logicMod26 = new G4LogicalVolume(solidMod26, ModMater, "Mod26",0,0,0);
  physiMod26 = new G4PVPlacement(0,                   // no rotation
				 positionMod26,  // at (x,y,z)
				 logicMod26,     // its logical volume				  
				 "Mod26",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field


  //----------------------------------------------------------
  //Settima fetta Modulatore 2
  //----------------------------------------------------------
 
      
  G4ThreeVector positionMod27 = G4ThreeVector(0*cm,0*cm,0.15*cm);
  
  solidMod27 = new G4Tubs("Mod27",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube7, 
			  spanningAngleOfTheTube7);

  logicMod27 = new G4LogicalVolume(solidMod27, ModMater, "Mod27",0,0,0);
  physiMod27 = new G4PVPlacement(0,                   // no rotation
				 positionMod27,  // at (x,y,z)
				 logicMod27,     // its logical volume				  
				 "Mod27",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field



  //----------------------------------------------------------
  //Ottava fetta Modulatore 2
  //----------------------------------------------------------
   
  

   
  G4ThreeVector positionMod28 = G4ThreeVector(0*cm,0*cm,0.09*cm);
  solidMod28 = new G4Tubs("Mod28",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube8, 
			  spanningAngleOfTheTube8);

  logicMod28 = new G4LogicalVolume(solidMod28, ModMater, "Mod28",0,0,0);
  physiMod28 = new G4PVPlacement(0,                   // no rotation
				 positionMod28,  // at (x,y,z)
				 logicMod28,     // its logical volume				  
				 "Mod28",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field




  //----------------------------------------------------------
  //Nona fetta Modulatore 3
  //----------------------------------------------------------
 

     
  G4ThreeVector positionMod29 = G4ThreeVector(0*cm,0*cm,0.03*cm);
  solidMod29 = new G4Tubs("Mod29",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube9, 
			  spanningAngleOfTheTube9);

  logicMod29 = new G4LogicalVolume(solidMod29, ModMater, "Mod29",0,0,0);
  physiMod29 = new G4PVPlacement(0,                   // no rotation
				 positionMod29,  // at (x,y,z)
				 logicMod29,     // its logical volume				  
				 "Mod29",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field


  //----------------------------------------------------------
  //Decima fetta Modulatore 3
  //----------------------------------------------------------
 
   
  G4ThreeVector positionMod30 = G4ThreeVector(0*cm,0*cm,-0.03*cm);
  solidMod30 = new G4Tubs("Mod30",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube10, 
			  spanningAngleOfTheTube10);

  logicMod30 = new G4LogicalVolume(solidMod30, ModMater, "Mod30",0,0,0);
  physiMod30 = new G4PVPlacement(0,                   // no rotation
				 positionMod30,  // at (x,y,z)
				 logicMod30,     // its logical volume				  
				 "Mod30",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  // Undicesima fetta Modulatore 3
  //----------------------------------------------------------

  G4ThreeVector positionMod31 = G4ThreeVector(0*cm,0*cm,-0.09*cm);
  solidMod31 = new G4Tubs("Mod31",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube11, 
			  spanningAngleOfTheTube11);

  logicMod31 = new G4LogicalVolume(solidMod31, ModMater, "Mod31",0,0,0);
  physiMod31 = new G4PVPlacement(0,                   // no rotation
				 positionMod31,  // at (x,y,z)
				 logicMod31,     // its logical volume				  
				 "Mod31",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Dodicesima fetta Modulatore 3
  //----------------------------------------------------------

  G4ThreeVector positionMod32 = G4ThreeVector(0*cm,0*cm,-0.15*cm);
  solidMod32 = new G4Tubs("Mod32",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube12, 
			  spanningAngleOfTheTube12);

  logicMod32 = new G4LogicalVolume(solidMod32, ModMater, "Mod32",0,0,0);
  physiMod32 = new G4PVPlacement(0,                   // no rotation
				 positionMod32,  // at (x,y,z)
				 logicMod32,     // its logical volume				  
				 "Mod32",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  

  //----------------------------------------------------------
  //Tredicesima fetta Modulatore 3
  //----------------------------------------------------------
 
  G4ThreeVector positionMod33 = G4ThreeVector(0*cm,0*cm,-0.21*cm);
  solidMod33 = new G4Tubs("Mod33",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube13, 
			  spanningAngleOfTheTube13);

  logicMod33 = new G4LogicalVolume(solidMod33, ModMater, "Mod33",0,0,0);
  physiMod33 = new G4PVPlacement(0,                   // no rotation
				 positionMod33,  // at (x,y,z)
				 logicMod33,     // its logical volume				  
				 "Mod33",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Quat. fetta Modulatore 3
  //----------------------------------------------------------
 
       
  G4ThreeVector positionMod34 = G4ThreeVector(0*cm,0*cm,-0.27*cm);
  solidMod34 = new G4Tubs("Mod34",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube14, 
			  spanningAngleOfTheTube14);

  logicMod34 = new G4LogicalVolume(solidMod34, ModMater, "Mod34",0,0,0);
  physiMod34 = new G4PVPlacement(0,                   // no rotation
				 positionMod34,  // at (x,y,z)
				 logicMod34,     // its logical volume				  
				 "Mod134",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  
  //----------------------------------------------------------
  //Quindicesima fetta Modulatore 2
  //----------------------------------------------------------
 
       
  G4ThreeVector positionMod35 = G4ThreeVector(0*cm,0*cm,-0.33*cm);
  solidMod35 = new G4Tubs("Mod35",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube15, 
			  spanningAngleOfTheTube15);

  logicMod35 = new G4LogicalVolume(solidMod35, ModMater, "Mod35",0,0,0);
  physiMod35 = new G4PVPlacement(0,                   // no rotation
				 positionMod35,  // at (x,y,z)
				 logicMod35,     // its logical volume				  
				 "Mod35",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Sedic. fetta Modulatore 2
  //----------------------------------------------------------
 
   
  G4ThreeVector positionMod36 = G4ThreeVector(0*cm,0*cm,-0.39*cm);
  solidMod36 = new G4Tubs("Mod36",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube16, 
			  spanningAngleOfTheTube16);

  logicMod36 = new G4LogicalVolume(solidMod36, ModMater, "Mod36",0,0,0);
  physiMod36 = new G4PVPlacement(0,                   // no rotation
				 positionMod36,  // at (x,y,z)
				 logicMod36,     // its logical volume				  
				 "Mod36",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Dicias. fetta Modulatore 2
  //----------------------------------------------------------
    
  G4ThreeVector positionMod37 = G4ThreeVector(0*cm,0*cm,-0.45*cm);
  solidMod37 = new G4Tubs("Mod37",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube17, 
			  spanningAngleOfTheTube17);

  logicMod37 = new G4LogicalVolume(solidMod37, ModMater, "Mod37",0,0,0);
  physiMod37 = new G4PVPlacement(0,                   // no rotation
				 positionMod37,  // at (x,y,z)
				 logicMod37,     // its logical volume				  
				 "Mod37",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  //----------------------------------------------------------
  //Diciot. fetta Modulatore 2
  //----------------------------------------------------------
    
   
  G4ThreeVector positionMod38 = G4ThreeVector(0*cm,0*cm,-0.51*cm);
  solidMod38 = new G4Tubs("Mod38",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube18, 
			  spanningAngleOfTheTube18);

  logicMod38 = new G4LogicalVolume(solidMod38, ModMater, "Mod38",0,0,0);
  physiMod38 = new G4PVPlacement(0,                   // no rotation
				 positionMod38,  // at (x,y,z)
				 logicMod38,     // its logical volume				  
				 "Mod38",        // its name
				 logicMod20,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field



  //----------------------------------------------------------
  //Volume Madre  3/4 del Modulatore  Mod 40
  //----------------------------------------------------------
  
    
  G4RotationMatrix rm40;
  rm40.rotateZ(180 *deg);
  
  G4ThreeVector positionMod40 = G4ThreeVector(0*cm,0*cm,0*cm);
  
  solidMod40 = new G4Tubs("Mod40",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube0,
			  startAngleOfTheTube0, 
			  spanningAngleOfTheTube0);
  
  logicMod40 = new G4LogicalVolume(solidMod40, Mod0Mater, "Mod40",0,0,0);
  
  
  physiMod40 = new G4PVPlacement(G4Transform3D(rm40, positionMod40), 
				 logicMod40,    // its logical volume			  
				 "Mod40",        // its name
				 logicMotherMod,  // its mother  volume
				 false,         // no boolean operations
				 0);            // no particular field
  




  //----------------------------------------------------------
  //Prima fetta Modulatore 3
  //----------------------------------------------------------
 
  G4ThreeVector positionMod41 = G4ThreeVector(0*cm,0*cm,0.51*cm);
  solidMod41 = new G4Tubs("Mod41",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube1, 
			  spanningAngleOfTheTube1);
  
  logicMod41 = new G4LogicalVolume(solidMod41, ModMater, "Mod41",0,0,0);
  
  physiMod41 = new G4PVPlacement(0,               // no rotation
				 positionMod41,  // at (x,y,z)
				 logicMod41,     // its logical volume				  
				 "Mod41",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  
  
  //----------------------------------------------------------
  //Seconda fetta Modulatore 3
  //----------------------------------------------------------
  
  G4ThreeVector positionMod42 = G4ThreeVector(0*cm,0*cm,0.45*cm);
  
  solidMod42 = new G4Tubs("Mod42",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube2, 
			  spanningAngleOfTheTube2);
  
  logicMod42 = new G4LogicalVolume(solidMod42, ModMater, "Mod42",0,0,0);
  
  physiMod42 = new G4PVPlacement(0,               // no rotation
				 positionMod42,  // at (x,y,z)
				 logicMod42,     // its logical volume				  
				 "Mod42",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  
  
  //----------------------------------------------------------
  //Terza fetta Modulatore 3
  //----------------------------------------------------------
 
  G4ThreeVector positionMod43 = G4ThreeVector(0*cm,0*cm,0.39*cm);

  solidMod43 = new G4Tubs("Mod43",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube3, 
			  spanningAngleOfTheTube3);

  logicMod43 = new G4LogicalVolume(solidMod43, ModMater, "Mod43",0,0,0);
  physiMod43 = new G4PVPlacement(0,                   // no rotation
				 positionMod43,  // at (x,y,z)
				 logicMod43,     // its logical volume				  
				 "Mod43",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

 
  //----------------------------------------------------------
  //Quarta fetta Modulatore 3
  //----------------------------------------------------------
   
   
  G4ThreeVector positionMod44 = G4ThreeVector(0*cm,0*cm,0.33*cm);

  solidMod44 = new G4Tubs("Mod44",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube4, 
			  spanningAngleOfTheTube4);

  logicMod44 = new G4LogicalVolume(solidMod44, ModMater, "Mod44",0,0,0);
 
  physiMod44 = new G4PVPlacement(0,                   // no rotation
				 positionMod44,  // at (x,y,z)
				 logicMod44,     // its logical volume				  
				 "Mod44",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field


  //----------------------------------------------------------
  //Quinta fetta Modulatore 3
  //----------------------------------------------------------
 
   
   
  G4ThreeVector positionMod45 = G4ThreeVector(0*cm,0*cm,0.27*cm);

  solidMod45 = new G4Tubs("Mod45",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube5, 
			  spanningAngleOfTheTube5);
  
  logicMod45 = new G4LogicalVolume(solidMod45, ModMater, "Mod45",0,0,0);
  physiMod45 = new G4PVPlacement(0,                   // no rotation
				 positionMod45,  // at (x,y,z)
				 logicMod45,     // its logical volume				  
				 "Mod45",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  
  //----------------------------------------------------------
  //Sesta fetta Modulatore 3
  //----------------------------------------------------------
 
   
  G4ThreeVector positionMod46 = G4ThreeVector(0*cm,0*cm,0.21*cm);
  solidMod46 = new G4Tubs("Mod46",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube6, 
			  spanningAngleOfTheTube6);

  logicMod46 = new G4LogicalVolume(solidMod46, ModMater, "Mod46",0,0,0);
  physiMod46 = new G4PVPlacement(0,                   // no rotation
				 positionMod46,  // at (x,y,z)
				 logicMod46,     // its logical volume				  
				 "Mod46",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field


  //----------------------------------------------------------
  //Settima fetta Modulatore 3
  //----------------------------------------------------------
 
      
  G4ThreeVector positionMod47 = G4ThreeVector(0*cm,0*cm,0.15*cm);
  
  solidMod47 = new G4Tubs("Mod47",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube7, 
			  spanningAngleOfTheTube7);

  logicMod47 = new G4LogicalVolume(solidMod47, ModMater, "Mod47",0,0,0);
  physiMod47 = new G4PVPlacement(0,                   // no rotation
				 positionMod47,  // at (x,y,z)
				 logicMod47,     // its logical volume				  
				 "Mod47",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field



  //----------------------------------------------------------
  //Ottava fetta Modulatore 3
  //----------------------------------------------------------
   
     
  G4ThreeVector positionMod48 = G4ThreeVector(0*cm,0*cm,0.09*cm);
  solidMod48 = new G4Tubs("Mod48",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube8, 
			  spanningAngleOfTheTube8);

  logicMod48 = new G4LogicalVolume(solidMod48, ModMater, "Mod48",0,0,0);
  physiMod48 = new G4PVPlacement(0,                   // no rotation
				 positionMod48,  // at (x,y,z)
				 logicMod48,     // its logical volume				  
				 "Mod48",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field




  //----------------------------------------------------------
  //Nona fetta Modulatore 2
  //----------------------------------------------------------
 

     
  G4ThreeVector positionMod49 = G4ThreeVector(0*cm,0*cm,0.03*cm);
  solidMod49 = new G4Tubs("Mod49",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube9, 
			  spanningAngleOfTheTube9);

  logicMod49 = new G4LogicalVolume(solidMod49, ModMater, "Mod49",0,0,0);
  physiMod49 = new G4PVPlacement(0,                   // no rotation
				 positionMod49,  // at (x,y,z)
				 logicMod49,     // its logical volume				  
				 "Mod49",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field


  //----------------------------------------------------------
  //Decima fetta Modulatore 3
  //----------------------------------------------------------
 
   
  G4ThreeVector positionMod50 = G4ThreeVector(0*cm,0*cm,-0.03*cm);
  solidMod50 = new G4Tubs("Mod50",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube10, 
			  spanningAngleOfTheTube10);

  logicMod50 = new G4LogicalVolume(solidMod50, ModMater, "Mod50",0,0,0);
  physiMod50 = new G4PVPlacement(0,                   // no rotation
				 positionMod50,  // at (x,y,z)
				 logicMod50,     // its logical volume				  
				 "Mod50",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  // Undicesima fetta Modulatore 3
  //----------------------------------------------------------

  G4ThreeVector positionMod51 = G4ThreeVector(0*cm,0*cm,-0.09*cm);
  solidMod51 = new G4Tubs("Mod51",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube11, 
			  spanningAngleOfTheTube11);

  logicMod51 = new G4LogicalVolume(solidMod51, ModMater, "Mod51",0,0,0);
  physiMod51 = new G4PVPlacement(0,                   // no rotation
				 positionMod51,  // at (x,y,z)
				 logicMod51,     // its logical volume				  
				 "Mod51",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Dodicesima fetta Modulatore 3
  //----------------------------------------------------------

  G4ThreeVector positionMod52 = G4ThreeVector(0*cm,0*cm,-0.15*cm);
  solidMod52 = new G4Tubs("Mod52",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube12, 
			  spanningAngleOfTheTube12);

  logicMod52 = new G4LogicalVolume(solidMod52, ModMater, "Mod52",0,0,0);
  physiMod52 = new G4PVPlacement(0,                   // no rotation
				 positionMod52,  // at (x,y,z)
				 logicMod52,     // its logical volume				  
				 "Mod52",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  

  //----------------------------------------------------------
  //Tredicesima fetta Modulatore 3
  //----------------------------------------------------------
 
  G4ThreeVector positionMod53 = G4ThreeVector(0*cm,0*cm,-0.21*cm);
  solidMod53 = new G4Tubs("Mod53",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube13, 
			  spanningAngleOfTheTube13);

  logicMod53 = new G4LogicalVolume(solidMod53, ModMater, "Mod53",0,0,0);
  physiMod53 = new G4PVPlacement(0,                   // no rotation
				 positionMod53,  // at (x,y,z)
				 logicMod53,     // its logical volume				  
				 "Mod53",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Quat. fetta Modulatore 3
  //----------------------------------------------------------
 
       
  G4ThreeVector positionMod54 = G4ThreeVector(0*cm,0*cm,-0.27*cm);
  solidMod54 = new G4Tubs("Mod54",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube14, 
			  spanningAngleOfTheTube14);

  logicMod54 = new G4LogicalVolume(solidMod54, ModMater, "Mod54",0,0,0);
  physiMod54 = new G4PVPlacement(0,                   // no rotation
				 positionMod54,  // at (x,y,z)
				 logicMod54,     // its logical volume				  
				 "Mod154",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  
  //----------------------------------------------------------
  //Quindicesima fetta Modulatore 3
  //----------------------------------------------------------
 
       
  G4ThreeVector positionMod55 = G4ThreeVector(0*cm,0*cm,-0.33*cm);
  solidMod55 = new G4Tubs("Mod35",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube15, 
			  spanningAngleOfTheTube15);

  logicMod55 = new G4LogicalVolume(solidMod55, ModMater, "Mod55",0,0,0);
  physiMod55 = new G4PVPlacement(0,                   // no rotation
				 positionMod55,  // at (x,y,z)
				 logicMod55,     // its logical volume				  
				 "Mod55",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Sedic. fetta Modulatore 3
  //----------------------------------------------------------
 
   
  G4ThreeVector positionMod56 = G4ThreeVector(0*cm,0*cm,-0.39*cm);
  solidMod56 = new G4Tubs("Mod56",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube16, 
			  spanningAngleOfTheTube16);

  logicMod56 = new G4LogicalVolume(solidMod56, ModMater, "Mod56",0,0,0);
  physiMod56 = new G4PVPlacement(0,                   // no rotation
				 positionMod56,  // at (x,y,z)
				 logicMod56,     // its logical volume				  
				 "Mod56",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Dicias. fetta Modulatore 3
  //----------------------------------------------------------
    
  G4ThreeVector positionMod57 = G4ThreeVector(0*cm,0*cm,-0.45*cm);
  solidMod57 = new G4Tubs("Mod57",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube17, 
			  spanningAngleOfTheTube17);

  logicMod57 = new G4LogicalVolume(solidMod57, ModMater, "Mod57",0,0,0);
  physiMod57 = new G4PVPlacement(0,                   // no rotation
				 positionMod57,  // at (x,y,z)
				 logicMod57,     // its logical volume				  
				 "Mod57",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Diciot. fetta Modulatore 3
  //----------------------------------------------------------
    
   
  G4ThreeVector positionMod58 = G4ThreeVector(0*cm,0*cm,-0.51*cm);
  solidMod58 = new G4Tubs("Mod58",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube18, 
			  spanningAngleOfTheTube18);

  logicMod58 = new G4LogicalVolume(solidMod58, ModMater, "Mod58",0,0,0);
  physiMod58 = new G4PVPlacement(0,                   // no rotation
				 positionMod58,  // at (x,y,z)
				 logicMod58,     // its logical volume				  
				 "Mod58",        // its name
				 logicMod40,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field



  //----------------------------------------------------------
  //Volume Madre  4/4 del Modulatore  Mod 60
  //----------------------------------------------------------
  
    
  G4RotationMatrix rm60;
  rm60.rotateZ(270 *deg);
  
  G4ThreeVector positionMod60 = G4ThreeVector(0*cm,0*cm,0*cm);
  
  solidMod60 = new G4Tubs("Mod60",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube0,
			  startAngleOfTheTube0, 
			  spanningAngleOfTheTube0);
  
  logicMod60 = new G4LogicalVolume(solidMod60, Mod0Mater, "Mod60",0,0,0);
  
  
  physiMod60 = new G4PVPlacement(G4Transform3D(rm60, positionMod60), 
				 logicMod60,    // its logical volume			  
				 "Mod60",        // its name
				 logicMotherMod,  // its mother  volume
				 false,         // no boolean operations
				 0);            // no particular field
  


  //----------------------------------------------------------
  //Prima fetta Modulatore 4
  //----------------------------------------------------------
 
  G4ThreeVector positionMod61 = G4ThreeVector(0*cm,0*cm,0.51*cm);
  solidMod61 = new G4Tubs("Mod61",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube1, 
			  spanningAngleOfTheTube1);
  
  logicMod61 = new G4LogicalVolume(solidMod61, ModMater, "Mod61",0,0,0);
  
  physiMod61 = new G4PVPlacement(0,               // no rotation
				 positionMod61,  // at (x,y,z)
				 logicMod61,     // its logical volume				  
				 "Mod61",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  
  
  //----------------------------------------------------------
  //Seconda fetta Modulatore 4
  //----------------------------------------------------------
  
  G4ThreeVector positionMod62 = G4ThreeVector(0*cm,0*cm,0.45*cm);
  
  solidMod62 = new G4Tubs("Mod62",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube2, 
			  spanningAngleOfTheTube2);
  
  logicMod62 = new G4LogicalVolume(solidMod62, ModMater, "Mod62",0,0,0);
  
  physiMod62 = new G4PVPlacement(0,               // no rotation
				 positionMod62,  // at (x,y,z)
				 logicMod62,     // its logical volume				  
				 "Mod62",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  
  
  //----------------------------------------------------------
  //Terza fetta Modulatore 4
  //----------------------------------------------------------
 
  G4ThreeVector positionMod63 = G4ThreeVector(0*cm,0*cm,0.39*cm);

  solidMod63 = new G4Tubs("Mod63",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube3, 
			  spanningAngleOfTheTube3);

  logicMod63 = new G4LogicalVolume(solidMod63, ModMater, "Mod63",0,0,0);
  physiMod63 = new G4PVPlacement(0,                   // no rotation
				 positionMod63,  // at (x,y,z)
				 logicMod63,     // its logical volume				  
				 "Mod63",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

 
  //----------------------------------------------------------
  //Quarta fetta Modulatore 4
  //----------------------------------------------------------
   
   
  G4ThreeVector positionMod64 = G4ThreeVector(0*cm,0*cm,0.33*cm);

  solidMod64 = new G4Tubs("Mod64",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube4, 
			  spanningAngleOfTheTube4);

  logicMod64 = new G4LogicalVolume(solidMod64, ModMater, "Mod64",0,0,0);
 
  physiMod64 = new G4PVPlacement(0,                   // no rotation
				 positionMod64,  // at (x,y,z)
				 logicMod64,     // its logical volume				  
				 "Mod64",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field


  //----------------------------------------------------------
  //Quinta fetta Modulatore 3
  //----------------------------------------------------------
 
   
   
  G4ThreeVector positionMod65 = G4ThreeVector(0*cm,0*cm,0.27*cm);

  solidMod65 = new G4Tubs("Mod65",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube5, 
			  spanningAngleOfTheTube5);
  
  logicMod65 = new G4LogicalVolume(solidMod65, ModMater, "Mod65",0,0,0);
  physiMod65 = new G4PVPlacement(0,                   // no rotation
				 positionMod65,  // at (x,y,z)
				 logicMod65,     // its logical volume				  
				 "Mod65",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  
  //----------------------------------------------------------
  //Sesta fetta Modulatore 4
  //----------------------------------------------------------
 
   
  G4ThreeVector positionMod66 = G4ThreeVector(0*cm,0*cm,0.21*cm);
  solidMod66 = new G4Tubs("Mod66",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube6, 
			  spanningAngleOfTheTube6);

  logicMod66 = new G4LogicalVolume(solidMod66, ModMater, "Mod66",0,0,0);
  physiMod66 = new G4PVPlacement(0,                   // no rotation
				 positionMod66,  // at (x,y,z)
				 logicMod66,     // its logical volume				  
				 "Mod66",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field


  //----------------------------------------------------------
  //Settima fetta Modulatore 4
  //----------------------------------------------------------
 
      
  G4ThreeVector positionMod67 = G4ThreeVector(0*cm,0*cm,0.15*cm);
  
  solidMod67 = new G4Tubs("Mod67",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube7, 
			  spanningAngleOfTheTube7);

  logicMod67 = new G4LogicalVolume(solidMod67, ModMater, "Mod67",0,0,0);
  physiMod67 = new G4PVPlacement(0,                   // no rotation
				 positionMod67,  // at (x,y,z)
				 logicMod67,     // its logical volume				  
				 "Mod67",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field



  //----------------------------------------------------------
  //Ottava fetta Modulatore 4
  //----------------------------------------------------------
   
     
  G4ThreeVector positionMod68 = G4ThreeVector(0*cm,0*cm,0.09*cm);
  solidMod68 = new G4Tubs("Mod68",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube8, 
			  spanningAngleOfTheTube8);

  logicMod68 = new G4LogicalVolume(solidMod68, ModMater, "Mod68",0,0,0);
  physiMod68 = new G4PVPlacement(0,                   // no rotation
				 positionMod68,  // at (x,y,z)
				 logicMod68,     // its logical volume				  
				 "Mod68",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field




  //----------------------------------------------------------
  //Nona fetta Modulatore 4
  //----------------------------------------------------------
 

     
  G4ThreeVector positionMod69 = G4ThreeVector(0*cm,0*cm,0.03*cm);
  solidMod69 = new G4Tubs("Mod69",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube9, 
			  spanningAngleOfTheTube9);

  logicMod69 = new G4LogicalVolume(solidMod69, ModMater, "Mod69",0,0,0);
  physiMod69 = new G4PVPlacement(0,                   // no rotation
				 positionMod69,  // at (x,y,z)
				 logicMod69,     // its logical volume				  
				 "Mod69",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field


  //----------------------------------------------------------
  //Decima fetta Modulatore 4
  //----------------------------------------------------------
 
   
  G4ThreeVector positionMod70 = G4ThreeVector(0*cm,0*cm,-0.03*cm);
  solidMod70 = new G4Tubs("Mod70",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube10, 
			  spanningAngleOfTheTube10);

  logicMod70 = new G4LogicalVolume(solidMod70, ModMater, "Mod70",0,0,0);
  physiMod70 = new G4PVPlacement(0,                   // no rotation
				 positionMod70,  // at (x,y,z)
				 logicMod70,     // its logical volume				  
				 "Mod70",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  // Undicesima fetta Modulatore 4
  //----------------------------------------------------------

  G4ThreeVector positionMod71 = G4ThreeVector(0*cm,0*cm,-0.09*cm);
  solidMod71 = new G4Tubs("Mod71",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube11, 
			  spanningAngleOfTheTube11);

  logicMod71 = new G4LogicalVolume(solidMod71, ModMater, "Mod71",0,0,0);
  physiMod71 = new G4PVPlacement(0,                   // no rotation
				 positionMod71,  // at (x,y,z)
				 logicMod71,     // its logical volume				  
				 "Mod71",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Dodicesima fetta Modulatore 4
  //----------------------------------------------------------

  G4ThreeVector positionMod72 = G4ThreeVector(0*cm,0*cm,-0.15*cm);
  solidMod72 = new G4Tubs("Mod72",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube12, 
			  spanningAngleOfTheTube12);

  logicMod72 = new G4LogicalVolume(solidMod72, ModMater, "Mod72",0,0,0);
  physiMod72 = new G4PVPlacement(0,                   // no rotation
				 positionMod72,  // at (x,y,z)
				 logicMod72,     // its logical volume				  
				 "Mod72",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field
  

  //----------------------------------------------------------
  //Tredicesima fetta Modulatore 4
  //----------------------------------------------------------
 
  G4ThreeVector positionMod73 = G4ThreeVector(0*cm,0*cm,-0.21*cm);
  solidMod73 = new G4Tubs("Mod73",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube13, 
			  spanningAngleOfTheTube13);

  logicMod73 = new G4LogicalVolume(solidMod73, ModMater, "Mod73",0,0,0);
  physiMod73 = new G4PVPlacement(0,                   // no rotation
				 positionMod73,  // at (x,y,z)
				 logicMod73,     // its logical volume				  
				 "Mod73",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Quat. fetta Modulatore 4
  //----------------------------------------------------------
 
       
  G4ThreeVector positionMod74 = G4ThreeVector(0*cm,0*cm,-0.27*cm);
  solidMod74 = new G4Tubs("Mod74",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube14, 
			  spanningAngleOfTheTube14);

  logicMod74 = new G4LogicalVolume(solidMod74, ModMater, "Mod74",0,0,0);
  physiMod74 = new G4PVPlacement(0,                   // no rotation
				 positionMod74,  // at (x,y,z)
				 logicMod74,     // its logical volume				  
				 "Mod174",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  
  //----------------------------------------------------------
  //Quindicesima fetta Modulatore 4
  //----------------------------------------------------------
 
       
  G4ThreeVector positionMod75 = G4ThreeVector(0*cm,0*cm,-0.33*cm);
  solidMod75 = new G4Tubs("Mod75",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube15, 
			  spanningAngleOfTheTube15);

  logicMod75 = new G4LogicalVolume(solidMod75, ModMater, "Mod75",0,0,0);
  physiMod75 = new G4PVPlacement(0,                   // no rotation
				 positionMod75,  // at (x,y,z)
				 logicMod75,     // its logical volume				  
				 "Mod75",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Sedic. fetta Modulatore 4
  //----------------------------------------------------------
 
   
  G4ThreeVector positionMod76 = G4ThreeVector(0*cm,0*cm,-0.39*cm);
  solidMod76 = new G4Tubs("Mod76",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube16, 
			  spanningAngleOfTheTube16);

  logicMod76 = new G4LogicalVolume(solidMod76, ModMater, "Mod76",0,0,0);
  physiMod76 = new G4PVPlacement(0,                   // no rotation
				 positionMod76,  // at (x,y,z)
				 logicMod76,     // its logical volume				  
				 "Mod76",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Dicias. fetta Modulatore 4
  //----------------------------------------------------------
    
  G4ThreeVector positionMod77 = G4ThreeVector(0*cm,0*cm,-0.45*cm);
  solidMod77 = new G4Tubs("Mod57",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube17, 
			  spanningAngleOfTheTube17);

  logicMod77 = new G4LogicalVolume(solidMod77, ModMater, "Mod77",0,0,0);
  physiMod77 = new G4PVPlacement(0,                   // no rotation
				 positionMod77,  // at (x,y,z)
				 logicMod77,     // its logical volume				  
				 "Mod77",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  //----------------------------------------------------------
  //Diciot. fetta Modulatore 3
  //----------------------------------------------------------
    
   
  G4ThreeVector positionMod78 = G4ThreeVector(0*cm,0*cm,-0.51*cm);
  solidMod78 = new G4Tubs("Mod78",
			  innerRadiusOfTheTube, 
			  outerRadiusOfTheTube,
			  hightOfTheTube,
			  startAngleOfTheTube18, 
			  spanningAngleOfTheTube18);

  logicMod78 = new G4LogicalVolume(solidMod78, ModMater, "Mod78",0,0,0);
  physiMod78 = new G4PVPlacement(0,                   // no rotation
				 positionMod78,  // at (x,y,z)
				 logicMod78,     // its logical volume				  
				 "Mod78",        // its name
				 logicMod60,      // its mother  volume
				 false,           // no boolean operations
				 0);              // no particular field

  G4VisAttributes * red = new G4VisAttributes( G4Colour(1. ,0. ,0.));
  red-> SetVisibility(true);
  red-> SetForceSolid(true);
  logicMotherMod -> SetVisAttributes(G4VisAttributes::Invisible);

  logicMod0 ->SetVisAttributes(G4VisAttributes::Invisible);
  logicMod20 ->SetVisAttributes(G4VisAttributes::Invisible);
  logicMod40 ->SetVisAttributes(G4VisAttributes::Invisible);
  logicMod60 ->SetVisAttributes(G4VisAttributes::Invisible);
  logicMod1 -> SetVisAttributes(red);
  logicMod2 -> SetVisAttributes(red);
  logicMod3 -> SetVisAttributes(red);
  logicMod4 -> SetVisAttributes(red);
  logicMod5 -> SetVisAttributes(red);
  logicMod6 -> SetVisAttributes(red);
  logicMod7 -> SetVisAttributes(red);
  logicMod8 -> SetVisAttributes(red);
  logicMod9 -> SetVisAttributes(red);
  logicMod10 -> SetVisAttributes(red);
  logicMod11 -> SetVisAttributes(red);
  logicMod12 -> SetVisAttributes(red);
  logicMod13 -> SetVisAttributes(red);
  logicMod14 -> SetVisAttributes(red);
  logicMod15 -> SetVisAttributes(red);
  logicMod16 -> SetVisAttributes(red);
  logicMod17 -> SetVisAttributes(red);
  logicMod18 -> SetVisAttributes(red);
  logicMod21 -> SetVisAttributes(red);
  logicMod22 -> SetVisAttributes(red);
  logicMod23 -> SetVisAttributes(red);
  logicMod24 -> SetVisAttributes(red);
  logicMod25 -> SetVisAttributes(red);
  logicMod26 -> SetVisAttributes(red);
  logicMod27 -> SetVisAttributes(red);
  logicMod28 -> SetVisAttributes(red);
  logicMod29 -> SetVisAttributes(red);
  logicMod30 -> SetVisAttributes(red);
  logicMod31 -> SetVisAttributes(red);
  logicMod32 -> SetVisAttributes(red);
  logicMod33 -> SetVisAttributes(red);
  logicMod34 -> SetVisAttributes(red);
  logicMod35 -> SetVisAttributes(red);
  logicMod36 -> SetVisAttributes(red);
  logicMod37 -> SetVisAttributes(red);
  logicMod38 -> SetVisAttributes(red);
  logicMod41 -> SetVisAttributes(red);
  logicMod42 -> SetVisAttributes(red);
  logicMod43 -> SetVisAttributes(red);
  logicMod44 -> SetVisAttributes(red);
  logicMod45 -> SetVisAttributes(red);
  logicMod46 -> SetVisAttributes(red);
  logicMod47 -> SetVisAttributes(red);
  logicMod48 -> SetVisAttributes(red);
  logicMod49 -> SetVisAttributes(red);
  logicMod50 -> SetVisAttributes(red);
  logicMod51 -> SetVisAttributes(red);
  logicMod52 -> SetVisAttributes(red);
  logicMod53 -> SetVisAttributes(red);
  logicMod54 -> SetVisAttributes(red);
  logicMod55 -> SetVisAttributes(red);
  logicMod56 -> SetVisAttributes(red);
  logicMod57 -> SetVisAttributes(red);
  logicMod58 -> SetVisAttributes(red);
  logicMod61 -> SetVisAttributes(red);
  logicMod62 -> SetVisAttributes(red);
  logicMod63 -> SetVisAttributes(red);
  logicMod64 -> SetVisAttributes(red);
  logicMod65 -> SetVisAttributes(red);
  logicMod66 -> SetVisAttributes(red);
  logicMod67 -> SetVisAttributes(red);
  logicMod68 -> SetVisAttributes(red);
  logicMod69 -> SetVisAttributes(red);
  logicMod70 -> SetVisAttributes(red);
  logicMod71 -> SetVisAttributes(red);
  logicMod72 -> SetVisAttributes(red);
  logicMod73 -> SetVisAttributes(red);
  logicMod74 -> SetVisAttributes(red);
  logicMod75 -> SetVisAttributes(red);
  logicMod76 -> SetVisAttributes(red);
  logicMod77 -> SetVisAttributes(red);
  logicMod78 -> SetVisAttributes(red);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyModulator::SetModulatorAngle(G4double angle)
{
  G4double rotationAngle = angle;
  rm -> rotateZ(rotationAngle);
  G4cout << "!!!!!!!!!!!!! " << rotationAngle/deg << G4endl;
  physiMotherMod -> SetRotation(rm);  
  G4cout << "MODULATOR HAS BEEN ROTATED OF " << rotationAngle/deg 
	 << " deg" << G4endl;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified(); 
}


