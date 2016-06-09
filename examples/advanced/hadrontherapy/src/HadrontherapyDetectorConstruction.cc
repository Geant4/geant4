//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: HadrontherapyDetectorConstruction.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.Pablo Cirrone, Giorgio Russo and Francesco Di Rosa 
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
#include "HadrontherapyDetectorMessenger.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyCalorimeterSD.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

HadrontherapyDetectorConstruction::HadrontherapyDetectorConstruction()
  :WorldMaterial(0),
   solidTreatmentRoom(0), logicTreatmentRoom(0), physiTreatmentRoom(0),
   solidBeamLineSupport(0), logicBeamLineSupport(0), physiBeamLineSupport(0),
   solidBeamLineCover(0), logicBeamLineCover(0), physiBeamLineCover(0),
   solidBeamLineCover2(0), logicBeamLineCover2(0), physiBeamLineCover2(0),
   solidVacuumZone(0), logicVacuumZone(0), physiVacuumZone(0),
   solidFirstScatteringFoil(0), logicFirstScatteringFoil(0), physiFirstScatteringFoil(0),
   solidKaptonWindow(0), logicKaptonWindow(0), physiKaptonWindow(0),
   solidStopper(0), logicStopper(0), physiStopper(0),
   solidSecondScatteringFoil(0), logicSecondScatteringFoil(0), physiSecondScatteringFoil(0),
   solidFirstCollimator(0), logicFirstCollimator(0), physiFirstCollimator(0),
   solidHoleFirstCollimator(0), logicHoleFirstCollimator(0), physiHoleFirstCollimator(0),
   solidFirstCollimatorModulatorBox(0), logicFirstCollimatorModulatorBox(0),
   physiFirstCollimatorModulatorBox(0),
   solidHoleFirstCollimatorModulatorBox(0),
   logicHoleFirstCollimatorModulatorBox(0),
   physiHoleFirstCollimatorModulatorBox(0),
   solidMotherMod(0),logicMotherMod(0), physiMotherMod(0), Mod0Mater(0),ModMater(0),
   solidMod0(0), logicMod0(0), physiMod0(0),
   solidMod1(0), logicMod1(0), physiMod1(0),
   solidMod2(0), logicMod2(0), physiMod2(0),
   solidMod3(0), logicMod3(0), physiMod3(0),
   solidMod4(0), logicMod4(0), physiMod4(0),
   solidMod5(0), logicMod5(0), physiMod5(0),
   solidMod6(0), logicMod6(0), physiMod6(0),
   solidMod7(0), logicMod7(0), physiMod7(0),
   solidMod8(0), logicMod8(0), physiMod8(0),
   solidMod9(0), logicMod9(0), physiMod9(0),
   solidMod10(0), logicMod10(0), physiMod10(0),
   solidMod11(0), logicMod11(0), physiMod11(0),
   solidMod12(0), logicMod12(0), physiMod12(0),
   solidMod13(0), logicMod13(0), physiMod13(0),
   solidMod14(0), logicMod14(0), physiMod14(0),
   solidMod15(0), logicMod15(0), physiMod15(0),
   solidMod16(0), logicMod16(0), physiMod16(0),
   solidMod17(0), logicMod17(0), physiMod17(0),
   solidMod18(0), logicMod18(0), physiMod18(0),
   solidSecondCollimatorModulatorBox(0),
   logicSecondCollimatorModulatorBox(0),
   physiSecondCollimatorModulatorBox(0),
   solidHoleSecondCollimatorModulatorBox(0),
   logicHoleSecondCollimatorModulatorBox(0),
   physiHoleSecondCollimatorModulatorBox(0),
   solidSecondCollimator(0),logicSecondCollimator(0),physiSecondCollimator(0),
   solidHoleSecondCollimator(0),logicHoleSecondCollimator(0),physiHoleSecondCollimator(0),
   solidFirstMonitorLayer1(0),         logicFirstMonitorLayer1(0),     physiFirstMonitorLayer1(0),
   solidFirstMonitorLayer2(0),         logicFirstMonitorLayer2(0),     physiFirstMonitorLayer2(0),
   solidFirstMonitorLayer3(0),         logicFirstMonitorLayer3(0),     physiFirstMonitorLayer3(0),
   solidFirstMonitorLayer4(0),         logicFirstMonitorLayer4(0),     physiFirstMonitorLayer4(0),

   solidSecondMonitorLayer1(0),        logicSecondMonitorLayer1(0),    physiSecondMonitorLayer1(0),
   solidSecondMonitorLayer2(0),        logicSecondMonitorLayer2(0),    physiSecondMonitorLayer2(0),
   solidSecondMonitorLayer3(0),        logicSecondMonitorLayer3(0),    physiSecondMonitorLayer3(0),
   solidSecondMonitorLayer4(0),        logicSecondMonitorLayer4(0),    physiSecondMonitorLayer4(0),

   solidThirdMonitorLayer1(0),         logicThirdMonitorLayer1(0),     physiThirdMonitorLayer1(0),
   solidThirdMonitorLayer2(0),         logicThirdMonitorLayer2(0),     physiThirdMonitorLayer2(0),
   solidThirdMonitorLayer3(0),         logicThirdMonitorLayer3(0),     physiThirdMonitorLayer3(0),
   solidThirdMonitorLayer4(0),         logicThirdMonitorLayer4(0),     physiThirdMonitorLayer4(0),

   solidNozzleSupport(0),              logicNozzleSupport(0),          physiNozzleSupport(0),  
   solidHoleNozzleSupport(0),          logicHoleNozzleSupport(0),      physiHoleNozzleSupport(0),
   solidSecondHoleNozzleSupport(0),    logicSecondHoleNozzleSupport(0),physiSecondHoleNozzleSupport(0),
   solidFinalCollimator(0),            logicFinalCollimator(0),        physiFinalCollimator(0),
   solidWaterPhantom(0),               logicWaterPhantom(0),           physiWaterPhantom(0),
   calorimeterSD(0)

{
  // create commands for interactive definition of the calorimeter
  detectorMessenger = new HadrontherapyDetectorMessenger(this);
}

// -----------------------------------------------------------------------------------------------
HadrontherapyDetectorConstruction::~HadrontherapyDetectorConstruction()
{
  delete detectorMessenger;
}

//**********************************
// DEFINITION OF THE USED MATERIALS
//**********************************

G4VPhysicalVolume* HadrontherapyDetectorConstruction::Construct()
{

 
  G4String name, symbol;             //a = ATOMIC MASS
  G4double a, z, density;            //z = ATOMIC NUMBER 

  G4int ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;


  //*************************************
  // DEFINITION OF MATERIALS AS ELEMENTS
  //*************************************


  a = 1.01*g/mole;
  G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);            //HIDROGEN

  a = 12.01*g/mole;
  G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);            // CARBON

  a = 14.01*g/mole;
  G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);            // NITROGEN

  a = 16.00*g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);            // OXIGEN

  density = 16.6*g/cm3;
  a = 180.948*g/mole;
  G4Material* Ta = new G4Material(name="Tantalio", z=73,a,density);                // TANTALUM

  density = 7.14*g/cm3;
  a = 65.38*g/mole;
  G4Material* Zn = new G4Material(name="Zinc"   , z=30, a, density);               // ZINC  

  density = 8.90*g/cm3;
  a = 63.546*g/mole;
  G4Material* Cu = new G4Material(name="Copper", z=29, a, density);                // COPPER

  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);            //ALLUMINUM 


  //**************************************
  // DEFINITION OF MATERIAL AS COMPOSIT  
  //**************************************


  density = 1.000*g/cm3;
  G4Material* Water = new G4Material(name="Water", density, ncomponents=2);         // WATER
  Water->AddElement(H, natoms=2);
  Water->AddElement(O, natoms=1);

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);           // AIR
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  density = 1.18*g/cm3;
  G4Material* PMMA = new G4Material(name="PMMA", density, ncomponents=3);           // PMMA
  PMMA->AddElement(C,  natoms=5);
  PMMA->AddElement(O,  natoms=2);
  PMMA->AddElement(H,  natoms=8);

  density = 8.40*g/cm3;
  G4Material* Brass = new G4Material(name="Brass", density, ncomponents=2);         // BRASS
  Brass->AddMaterial(Zn, fractionmass=30*perCent);
  Brass->AddMaterial(Cu, fractionmass=70*perCent);

  density = 1.43*g/cm3;
  G4Material* Kapton = new G4Material(name="Kapton",density, ncomponents=4);
  Kapton->AddElement(H, natoms=10);
  Kapton->AddElement(O, natoms=5);
  Kapton->AddElement(C, natoms=22);
  Kapton->AddElement(N, natoms=2);
  
  //*************************
  // Definition of the vacuum 
  //*************************

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material* Vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole,density,
				      kStateGas,temperature,pressure);
		     
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  Mod0Mater = Air;
  ModMater = Air;
   

  //*********************************
  //  DEFINITIONS OF ALL THE VOLUMES
  //*********************************


  //******************************
  // TREATMENT ROOM (world volume)
  //    beam line along x axis
  //******************************

  G4double TreatmentRoom_x = 200.*cm;
  G4double TreatmentRoom_y = 200.*cm;
  G4double TreatmentRoom_z = 200.*cm;
  
  G4Box* TreatmentRoom
    = new G4Box("TreatmentRoom",
		TreatmentRoom_x,
		TreatmentRoom_y,
		TreatmentRoom_z);
  logicTreatmentRoom = new G4LogicalVolume(TreatmentRoom,Air,
					   "TreatmentRoom",
					   0,
					   0,
					   0);
  physiTreatmentRoom = new G4PVPlacement(0,
					 G4ThreeVector(),
					 "TreatmentRoom",
					 logicTreatmentRoom,
					 0,
					 false,
					 0);


  //********************************
  //      BEAM LINE SUPPORT
  //********************************


  G4double BeamLineSupport_x = 1.5     *m;
  G4double BeamLineSupport_y = 20.     *mm;
  G4double BeamLineSupport_z = 600.    *mm;                               

  G4double BeamLineSupportPosition_x = -700. *mm;
  G4double BeamLineSupportPosition_y = -230. *mm; 
  G4double BeamLineSupportPosition_z = 0.    *mm;

  G4Box* BeamLineSupport = new G4Box("BeamLineSupport",	
				     BeamLineSupport_x,
				     BeamLineSupport_y,
				     BeamLineSupport_z);

  logicBeamLineSupport = new G4LogicalVolume(BeamLineSupport,
					     Al,
					     "BeamLineSupport");

  physiBeamLineSupport = new G4PVPlacement(0,
					   G4ThreeVector(BeamLineSupportPosition_x,
							 BeamLineSupportPosition_y,
							 BeamLineSupportPosition_z),
					   "BeamLineSupport",
					   logicBeamLineSupport,
					   physiTreatmentRoom,
					   false,
					   0);

  //********************************
  //  BEAM LINE COVER 1 (left panel)
  //********************************

  G4double BeamLineCover_x = 1.5    *m;
  G4double BeamLineCover_y = 750.   *mm;
  G4double BeamLineCover_z = 10.    *mm;                               

  G4double BeamLineCoverPosition_x = -700. *mm;
  G4double BeamLineCoverPosition_y = -980. *mm; 
  G4double BeamLineCoverPosition_z =  600. *mm;

  G4Box* BeamLineCover = new G4Box("BeamLineCover",	
				   BeamLineCover_x,
				   BeamLineCover_y,
				   BeamLineCover_z);

  logicBeamLineCover = new G4LogicalVolume(BeamLineCover,
					   Al,
					   "BeamLineCover");

  physiBeamLineCover = new G4PVPlacement(0,
					 G4ThreeVector(BeamLineCoverPosition_x,
						       BeamLineCoverPosition_y,
						       BeamLineCoverPosition_z),
					 "BeamLineCover",
					 logicBeamLineCover,
					 physiTreatmentRoom,
					 false,
					 0);


  //********************************
  //  BEAM LINE COVER 2 (rigth panel)
  //********************************

  G4double BeamLineCover2_x = 1.5    *m;
  G4double BeamLineCover2_y = 750.   *mm;
  G4double BeamLineCover2_z = 10.    *mm;                               

  G4double BeamLineCover2Position_x = -700. *mm;
  G4double BeamLineCover2Position_y = -980. *mm; 
  G4double BeamLineCover2Position_z =  -600. *mm;

  G4Box* BeamLineCover2 = new G4Box("BeamLineCover2",	
				    BeamLineCover2_x,
				    BeamLineCover2_y,
				    BeamLineCover2_z);

  logicBeamLineCover2 = new G4LogicalVolume(BeamLineCover2,
					    Al,
					    "BeamLineCover2");

  physiBeamLineCover2 = new G4PVPlacement(0,
					  G4ThreeVector(BeamLineCover2Position_x,
							BeamLineCover2Position_y,
							BeamLineCover2Position_z),
					  "BeamLineCover2",
					  logicBeamLineCover2,
					  physiTreatmentRoom,
					  false,
					  0);

  //********************************
  // VACUUM ZONE
  // Before the beam exits
  // in air it traverses a vacuum zone 
  // in which a first scattering foil 
  // may be placed    
  //********************************

  G4double VacuumZone_x = 60.5325 *mm;
  G4double VacuumZone_y = 52.5    *mm;
  G4double VacuumZone_z = 52.5    *mm;                               

  G4double VacuumZonePosition_x = -1939.4675 *mm;
  G4double VacuumZonePosition_y = 0.         *mm; 
  G4double VacuumZonePosition_z = 0.         *mm;

  G4Box* VacuumZone= new G4Box("VacuumZone",	
			       VacuumZone_x,
			       VacuumZone_y,
			       VacuumZone_z);

  logicVacuumZone = new G4LogicalVolume(VacuumZone,
					Vacuum,
					"VacuumZone");

  physiVacuumZone = new G4PVPlacement(0,
				      G4ThreeVector(VacuumZonePosition_x,
						    VacuumZonePosition_y,
						    VacuumZonePosition_z),
				      "VacuumZone",
				      logicVacuumZone,
				      physiTreatmentRoom,
				      false,
				      0);


  //**************************
  // FIRST SCATTERING FOIL 
  // placed in the Vaccum Zone
  //**************************

  G4double FirstScatteringFoil_x = 0.0075 *mm;
  G4double FirstScatteringFoil_y = 52.5   *mm;
  G4double FirstScatteringFoil_z = 52.5   *mm;

  G4double FirstScatteringFoilPosition_x = -59.525 *mm;
  G4double FirstScatteringFoilPosition_y =  0      *mm;
  G4double FirstScatteringFoilPosition_z =  0      *mm;

  G4Box* FirstScatteringFoil = new G4Box("FirstScatteringFoil",
					 FirstScatteringFoil_x,
					 FirstScatteringFoil_y,
					 FirstScatteringFoil_z);

  logicFirstScatteringFoil = new G4LogicalVolume(FirstScatteringFoil,
						 Ta,
						 "FirstScatteringFoil");


  physiFirstScatteringFoil = new G4PVPlacement(0,
					       G4ThreeVector(FirstScatteringFoilPosition_x,
							     FirstScatteringFoilPosition_y,
							     FirstScatteringFoilPosition_z),
					       "FirstScatteringFoil",
					       logicFirstScatteringFoil,
					       physiVacuumZone,
					       false,
					       0);            
   

  //*********************
  // KAPTON WINDOW
  // Is the point in which 
  // the beam exits in air  
  //**********************

  solidKaptonWindow = new G4Box("KaptonWindow",
				0.025*mm,
				5.25*cm,
				5.25*cm);

  logicKaptonWindow = new G4LogicalVolume(solidKaptonWindow,
					  Kapton,
					  "KaptonWindow");

  physiKaptonWindow = new G4PVPlacement(0,
					G4ThreeVector(60.5075*mm,
						      0.*cm,
						      0.*cm),
					"KaptonWindow",
					logicKaptonWindow,
					physiVacuumZone,
					false,
					0); 


  // ************************
  // STOPPER
  // a small cilinder 
  // on the beam axis                     
  // ************************

  G4double phi = 90. *deg;     // Matrix definition for a 90 deg rotation. Also used for other volumes       
  G4RotationMatrix rm;               
  rm.rotateY(phi);

  G4double innerRadiusStopper = 0.*cm;
  G4double outerRadiusStopper = 2*mm;
  G4double hightStopper = 3.5*mm;
  G4double startAngleStopper = 0.*deg;
  G4double spanningAngleStopper = 360.*deg;

  G4double StopperPosition_x = -1707.435*mm;
  G4double StopperPosition_y = 0.*m;
  G4double StopperPosition_z = 0.*m;

  solidStopper = new G4Tubs("Stopper",
			    innerRadiusStopper,
			    outerRadiusStopper,
			    hightStopper,
			    startAngleStopper,
			    spanningAngleStopper);
  
  logicStopper = new G4LogicalVolume(solidStopper,
				     Brass,
				     "Stopper",
				     0,
				     0,
				     0);

  physiStopper = new G4PVPlacement(G4Transform3D(rm,
						 G4ThreeVector(StopperPosition_x,
							       StopperPosition_y,
							       StopperPosition_z)),
				   "Stopper",
				   logicStopper,
				   physiTreatmentRoom,
				   false,
				   0);

  // **********************
  // SECOND SCATTERING FOIL
  // **********************

  G4double SecondScatteringFoil_x = 0.0125 *mm;  
  G4double SecondScatteringFoil_y = 52.5   *mm;
  G4double SecondScatteringFoil_z = 52.5   *mm;

  G4double SecondScatteringFoilPosition_x = -1703.9225 *mm;
  G4double SecondScatteringFoilPosition_y =  0         *mm;
  G4double SecondScatteringFoilPosition_z =  0         *mm;                                                         
  
  solidSecondScatteringFoil = new G4Box("SecondScatteringFoil",
					SecondScatteringFoil_x,
					SecondScatteringFoil_y,
					SecondScatteringFoil_z);

  logicSecondScatteringFoil = new G4LogicalVolume(solidSecondScatteringFoil,
						  Ta,
						  "SecondScatteringFoil");

  physiSecondScatteringFoil = new G4PVPlacement(0,
						G4ThreeVector(SecondScatteringFoilPosition_x,
							      SecondScatteringFoilPosition_y,
							      SecondScatteringFoilPosition_z),
						"SeconScatteringFoil",
						logicSecondScatteringFoil,
						physiTreatmentRoom,
						false,
						0);


  // ********************
  // FIRST COLLIMATOR
  // ********************

  G4double FirstCollimator_x = 20    *mm;
  G4double FirstCollimator_y = 100   *mm;
  G4double FirstCollimator_z = 100   *mm;

  G4double FirstCollimatorPosition_x = -1683.91   *mm;
  G4double FirstCollimatorPosition_y =  0         *mm;
  G4double FirstCollimatorPosition_z =  0         *mm;

  solidFirstCollimator = new G4Box("FirstCollimator",
				   FirstCollimator_x,
				   FirstCollimator_y,
				   FirstCollimator_z);

  logicFirstCollimator = new G4LogicalVolume(solidFirstCollimator,
					     PMMA,
					     "FirstCollimator");

  physiFirstCollimator = new G4PVPlacement(0,
					   G4ThreeVector(FirstCollimatorPosition_x,
							 FirstCollimatorPosition_y,
							 FirstCollimatorPosition_z),
					   "FirstCollimator",
					   logicFirstCollimator,
					   physiTreatmentRoom,
					   false,
					   0);

  //**************************
  // HOLE OF FIRST COLLIMATOR
  //**************************

  G4double innerRadiusHoleFirstCollimator   = 0.*mm;
  G4double outerRadiusHoleFirstCollimator   = 15.*mm;
  G4double hightHoleFirstCollimator         = 20.*mm;
  G4double startAngleHoleFirstCollimator    = 0.*deg;
  G4double spanningAngleHoleFirstCollimator = 360.*deg;


  solidHoleFirstCollimator = new G4Tubs("HoleFirstCollimator",
					innerRadiusHoleFirstCollimator,
					outerRadiusHoleFirstCollimator,
					hightHoleFirstCollimator,
					startAngleHoleFirstCollimator,
					spanningAngleHoleFirstCollimator);

  logicHoleFirstCollimator = new G4LogicalVolume(solidHoleFirstCollimator,
						 Air,
						 "HoleFirstCollimator",
						 0,
						 0,
						 0);

  physiHoleFirstCollimator = new G4PVPlacement(G4Transform3D(rm,
							     G4ThreeVector()),
					       "HoleFirstCollimator",
					       logicHoleFirstCollimator,
					       physiFirstCollimator,
					       false,
					       0); 
  

  // ******************
  // SECOND COLLIMATOR
  // ******************

  G4double SecondCollimator_x = 20    *mm;
  G4double SecondCollimator_y = 100   *mm;
  G4double SecondCollimator_z = 100   *mm;

  G4double SecondCollimatorPosition_x = -779.91   *mm;
  G4double SecondCollimatorPosition_y =  0         *mm;
  G4double SecondCollimatorPosition_z =  0         *mm;


  solidSecondCollimator = new G4Box("SecondCollimator",
				    SecondCollimator_x,
				    SecondCollimator_y,
				    SecondCollimator_z);

  logicSecondCollimator = new G4LogicalVolume(solidSecondCollimator,
					      PMMA,
					      "SecondCollimator");
 
  physiSecondCollimator = new G4PVPlacement(0,
					    G4ThreeVector(SecondCollimatorPosition_x,
							  SecondCollimatorPosition_y,
							  SecondCollimatorPosition_z),
					    "SecondCollimator",
					    logicSecondCollimator,
					    physiTreatmentRoom,
					    false,
					    0);
  
  //**************************
  // HOLE OF SECOND COLLIMATOR
  //**************************

  G4double innerRadiusHoleSecondCollimator   = 0.*mm;
  G4double outerRadiusHoleSecondCollimator   = 15.*mm;
  G4double hightHoleSecondCollimator         = 20.*mm;
  G4double startAngleHoleSecondCollimator    = 0.*deg;
  G4double spanningAngleHoleSecondCollimator = 360.*deg;

  solidHoleSecondCollimator = new G4Tubs("HoleSecondCollimator",
					 innerRadiusHoleSecondCollimator,
					 outerRadiusHoleSecondCollimator,
					 hightHoleSecondCollimator,
					 startAngleHoleSecondCollimator,
					 spanningAngleHoleSecondCollimator);

  logicHoleSecondCollimator = new G4LogicalVolume(solidHoleSecondCollimator,
						  Air,
						  "HoleSecondCollimator",
						  0,
						  0,
						  0);

  physiHoleSecondCollimator = new G4PVPlacement(G4Transform3D(rm,
							      G4ThreeVector()),
						"HoleSecondCollimator",
						logicHoleSecondCollimator,
						physiSecondCollimator,
						false,
						0); 
  

  // *************************
  //            FIRST
  //  COLLIMATOR MODULATOR BOX
  // **************************

  G4double FirstCollimatorModulatorBox_x = 10    *mm;
  G4double FirstCollimatorModulatorBox_y = 200   *mm;
  G4double FirstCollimatorModulatorBox_z = 200   *mm;

  G4double FirstCollimatorModulatorBoxPosition_x = -1411.91   *mm;
  G4double FirstCollimatorModulatorBoxPosition_y =  0         *mm;
  G4double FirstCollimatorModulatorBoxPosition_z =  0         *mm;

  solidFirstCollimatorModulatorBox = new G4Box("FirstCollimatorModulatorBox",
					       FirstCollimatorModulatorBox_x,
					       FirstCollimatorModulatorBox_y,
					       FirstCollimatorModulatorBox_z);


  logicFirstCollimatorModulatorBox = new G4LogicalVolume(solidFirstCollimatorModulatorBox,
							 Al,
							 "FirstCollimatorModulatorBox");
  physiFirstCollimatorModulatorBox = new G4PVPlacement(0,
						       G4ThreeVector(FirstCollimatorModulatorBoxPosition_x,
								     FirstCollimatorModulatorBoxPosition_y,
								     FirstCollimatorModulatorBoxPosition_z),
						       "FirstCollimatorModulatorBox",
						       logicFirstCollimatorModulatorBox,
						       physiTreatmentRoom,
						       false,
						       0);

  // *************************************
  //   HOLE OF THE FIRST COLLLIMATOR 
  //                MODULATOR  BOX
  //  ************************************

  G4double innerRadiusHoleFirstCollimatorModulatorBox     = 0.   *mm;
  G4double outerRadiusHoleFirstCollimatorModulatorBox     = 31.  *mm;
  G4double hightHoleFirstCollimatorModulatorBox           = 10.  *mm;
  G4double startAngleHoleFirstCollimatorModulatorBox      = 0.   *deg;
  G4double spanningAngleHoleFirstCollimatorModulatorBox   = 360. *deg;

  solidHoleFirstCollimatorModulatorBox  = new G4Tubs("HoleFirstCollimatorModulatorBox",
						     innerRadiusHoleFirstCollimatorModulatorBox,
						     outerRadiusHoleFirstCollimatorModulatorBox,
						     hightHoleFirstCollimatorModulatorBox ,
						     startAngleHoleFirstCollimatorModulatorBox,
						     spanningAngleHoleFirstCollimatorModulatorBox);

  logicHoleFirstCollimatorModulatorBox = new G4LogicalVolume(solidHoleFirstCollimatorModulatorBox,
							     Air,
							     "HoleFirstCollimatorModulatorBox",
							     0,
							     0,
							     0);

  physiHoleFirstCollimatorModulatorBox = new G4PVPlacement(G4Transform3D(rm,
									 G4ThreeVector()),
							   "HoleFirstCollimatorModulatorBox",
							   logicHoleFirstCollimatorModulatorBox,
							   physiFirstCollimatorModulatorBox,
							   false,
							   0); 
  
  // *************************
  //            SECOND
  //  COLLIMATOR MODULATOR BOX
  // **************************

  G4double SecondCollimatorModulatorBox_x = 10    *mm;
  G4double SecondCollimatorModulatorBox_y = 200   *mm;
  G4double SecondCollimatorModulatorBox_z = 200   *mm;

  G4double SecondCollimatorModulatorBoxPosition_x = -841.91   *mm;
  G4double SecondCollimatorModulatorBoxPosition_y =  0         *mm;
  G4double SecondCollimatorModulatorBoxPosition_z =  0         *mm;

  solidSecondCollimatorModulatorBox= new G4Box("SecondCollimatorModulatorBox",
					       SecondCollimatorModulatorBox_x,
					       SecondCollimatorModulatorBox_y,
					       SecondCollimatorModulatorBox_z);

  logicSecondCollimatorModulatorBox = new G4LogicalVolume(solidSecondCollimatorModulatorBox,
							  Al,
							  "SecondCollimatorModulatorBox");

  physiSecondCollimatorModulatorBox = new G4PVPlacement(0,
							G4ThreeVector(SecondCollimatorModulatorBoxPosition_x,
								      SecondCollimatorModulatorBoxPosition_y,
								      SecondCollimatorModulatorBoxPosition_z),
							"SecondCollimatorModulatorBox",
							logicSecondCollimatorModulatorBox,
							physiTreatmentRoom,
							false,
							0);
  
  // *************************************
  //   HOLE OF THE SECOND  COLLLIMATOR 
  //                MODULATOR  BOX
  //  ************************************

  G4double innerRadiusHoleSecondCollimatorModulatorBox     = 0.   *mm;
  G4double outerRadiusHoleSecondCollimatorModulatorBox     = 31.  *mm;
  G4double hightHoleSecondCollimatorModulatorBox           = 10.  *mm;
  G4double startAngleHoleSecondCollimatorModulatorBox      = 0.   *deg;
  G4double spanningAngleHoleSecondCollimatorModulatorBox   = 360. *deg;

  solidHoleSecondCollimatorModulatorBox  = new G4Tubs("HoleSecondCollimatorModulatorBox",
						      innerRadiusHoleSecondCollimatorModulatorBox,
						      outerRadiusHoleSecondCollimatorModulatorBox,
						      hightHoleSecondCollimatorModulatorBox ,
						      startAngleHoleSecondCollimatorModulatorBox,
						      spanningAngleHoleSecondCollimatorModulatorBox);

  logicHoleSecondCollimatorModulatorBox = new G4LogicalVolume(solidHoleSecondCollimatorModulatorBox,
							      Air,
							      "HoleSecondCollimatorModulatorBox",
							      0,
							      0,
							      0);

  physiHoleSecondCollimatorModulatorBox = new G4PVPlacement(G4Transform3D(rm,
									  G4ThreeVector()),
							    "HoleSecondCollimatorModulatorBox",
							    logicHoleSecondCollimatorModulatorBox,
							    physiSecondCollimatorModulatorBox,
							    false,
							    0); 

  // ************************
  // FIRST MONITOR CHAMBER
  // ************************
  
  // Each chamber consist of 9 mm of air in a box
  // that has two layers one of kapton and one
  // of copper
  
  solidFirstMonitorLayer1 = new G4Box("FirstMonitorLayer1",
				      4.525022*mm,
				      10*cm,
				      10*cm);

  logicFirstMonitorLayer1 = new G4LogicalVolume(solidFirstMonitorLayer1,
						Kapton,
						"FirstMonitorLayer1");

  physiFirstMonitorLayer1 = new G4PVPlacement(0,
					      G4ThreeVector(-517.384978*mm,0.*cm,0.*cm),
					      "FirstMonitorLayer1",
					      logicFirstMonitorLayer1,
					      physiTreatmentRoom,
					      false,
					      0);
  
  solidFirstMonitorLayer2 = new G4Box("FirstMonitorLayer2",
				      0.000011*mm,
				      10*cm,
				      10*cm);

  logicFirstMonitorLayer2 = new G4LogicalVolume(solidFirstMonitorLayer2,
						Cu,
						"FirstMonitorLayer2");

  physiFirstMonitorLayer2 = new G4PVPlacement(0,
					      G4ThreeVector(-4.500011*mm,0.*cm,0.*cm),
					      "FirstMonitorLayer2",
					      logicFirstMonitorLayer2,
					      physiFirstMonitorLayer1,
					      false,
					      0);
  
  solidFirstMonitorLayer3 = new G4Box("FirstMonitorLayer3",
				      4.5*mm,
				      10*cm,
				      10*cm);

  logicFirstMonitorLayer3 = new G4LogicalVolume(solidFirstMonitorLayer3,
						Air,
						"FirstMonitorLayer3");

  physiFirstMonitorLayer3 = new G4PVPlacement(0,
					      G4ThreeVector(0.*mm,0.*cm,0.*cm),
					      "MonitorLayer3",
					      logicFirstMonitorLayer3,
					      physiFirstMonitorLayer1,
					      false,
					      0);
    
  solidFirstMonitorLayer4 = new G4Box("FirstMonitorLayer4",
				      0.000011*mm,
				      10*cm,
				      10*cm);

  logicFirstMonitorLayer4 = new G4LogicalVolume(solidFirstMonitorLayer4,
						Cu,
						"FirstMonitorLayer4");

  physiFirstMonitorLayer4 = new G4PVPlacement(0,
					      G4ThreeVector(4.500011*mm,0.*cm,0.*cm),
					      "FirstMonitorLayer4",
					      logicFirstMonitorLayer4,
					      physiFirstMonitorLayer1,
					      false,
					      0);

  // ************************
  // SECOND MONITOR CHAMBER
  // ************************

  solidSecondMonitorLayer1 = new G4Box("SecondMonitorLayer1",
				       4.525022*mm,
				       10*cm,
				       10*cm);

  logicSecondMonitorLayer1 = new G4LogicalVolume(solidSecondMonitorLayer1,
						 Kapton,
						 "SecondMonitorLayer1");

  physiSecondMonitorLayer1 = new G4PVPlacement(0,
					       G4ThreeVector(-386.334934 *mm,0.*cm,0.*cm),
					       "SecondMonitorLayer1",
					       logicSecondMonitorLayer1,
					       physiTreatmentRoom,
					       false,
					       0);
  
  solidSecondMonitorLayer2 = new G4Box("SecondMonitorLayer2",
				       0.000011*mm,
				       10*cm,
				       10*cm);

  logicSecondMonitorLayer2 = new G4LogicalVolume(solidSecondMonitorLayer2,
						 Cu,
						 "SecondMonitorLayer2");

  physiSecondMonitorLayer2 = new G4PVPlacement(0,
					       G4ThreeVector(-4.500011*mm,0.*cm,0.*cm),
					       "SecondMonitorLayer2",
					       logicSecondMonitorLayer2,
					       physiSecondMonitorLayer1,
					       false,
					       0);
  
  solidSecondMonitorLayer3 = new G4Box("SecondMonitorLayer3",
				       4.5*mm,
				       10*cm,
				       10*cm);

  logicSecondMonitorLayer3 = new G4LogicalVolume(solidSecondMonitorLayer3,
						 Air,
						 "SecondMonitorLayer3");

  physiSecondMonitorLayer3 = new G4PVPlacement(0,
					       G4ThreeVector(0.*mm,0.*cm,0.*cm),
					       "MonitorLayer3",
					       logicSecondMonitorLayer3,
					       physiSecondMonitorLayer1,
					       false,
					       0);
    
  solidSecondMonitorLayer4 = new G4Box("SecondMonitorLayer4",
				       0.000011*mm,
				       10*cm,
				       10*cm);

  logicSecondMonitorLayer4 = new G4LogicalVolume(solidSecondMonitorLayer4,
						 Cu,
						 "SecondMonitorLayer4");

  physiSecondMonitorLayer4 = new G4PVPlacement(0,
					       G4ThreeVector(4.500011*mm,0.*cm,0.*cm),
					       "SecondMonitorLayer4",
					       logicSecondMonitorLayer4,
					       physiSecondMonitorLayer1,
					       false,
					       0);


  //************************
  // THIRD MONITOR CHAMBER
  // ************************

  solidThirdMonitorLayer1 = new G4Box("ThirdMonitorLayer1",
				      4.525022*mm,
				      10*cm,
				      10*cm);

  logicThirdMonitorLayer1 = new G4LogicalVolume(solidThirdMonitorLayer1,
						Kapton,
						"ThirdMonitorLayer1");

  physiThirdMonitorLayer1 = new G4PVPlacement(0,
					      G4ThreeVector(-257.28489 *mm,0.*cm,0.*cm),
					      "ThirdMonitorLayer1",
					      logicThirdMonitorLayer1,
					      physiTreatmentRoom,
					      false,
					      0);
  
  solidThirdMonitorLayer2 = new G4Box("ThirdMonitorLayer2",
				      0.000011*mm,
				      10*cm,
				      10*cm);

  logicThirdMonitorLayer2 = new G4LogicalVolume(solidThirdMonitorLayer2,
						Cu,
						"ThirdMonitorLayer2");

  physiThirdMonitorLayer2 = new G4PVPlacement(0,
					      G4ThreeVector(-4.500011*mm,0.*cm,0.*cm),
					      "ThirdMonitorLayer2",
					      logicThirdMonitorLayer2,
					      physiThirdMonitorLayer1,
					      false,
					      0);
 
  solidThirdMonitorLayer3 = new G4Box("ThirdMonitorLayer3",
				      4.5*mm,
				      10*cm,
				      10*cm);

  logicThirdMonitorLayer3 = new G4LogicalVolume(solidThirdMonitorLayer3,
						Air,
						"ThirdMonitorLayer3");

  physiThirdMonitorLayer3 = new G4PVPlacement(0,
					      G4ThreeVector(0.*mm,0.*cm,0.*cm),
					      "MonitorLayer3",
					      logicThirdMonitorLayer3,
					      physiThirdMonitorLayer1,
					      false,
					      0);
  
  solidThirdMonitorLayer4 = new G4Box("ThirdMonitorLayer4",
				      0.000011*mm,
				      10*cm,
				      10*cm);

  logicThirdMonitorLayer4 = new G4LogicalVolume(solidThirdMonitorLayer4,
						Cu,
						"ThirdMonitorLayer4");

  physiThirdMonitorLayer4 = new G4PVPlacement(0,
					      G4ThreeVector(4.500011*mm,0.*cm,0.*cm),
					      "ThirdMonitorLayer4",
					      logicThirdMonitorLayer4,
					      physiSecondMonitorLayer1,
					      false,
					      0);

  // ***************
  // NOZZLE SUPPORT
  // ***************

  G4double NozzleSupport_x =  29.5 *mm;
  G4double NozzleSupport_y =  180. *mm;
  G4double NozzleSupport_z =  180. *mm;

  G4double NozzleSupportPosition_x = 647.59 *mm;
  G4double NozzleSupportPosition_y = 0. *mm;
  G4double NozzleSupportPosition_z = 0. *mm;

  solidNozzleSupport = new G4Box("NozzlSupport",
				 NozzleSupport_x,
				 NozzleSupport_y,
				 NozzleSupport_z);

  logicNozzleSupport = new G4LogicalVolume(solidNozzleSupport,
					   PMMA,
					   "NozzleSupport");

  physiNozzleSupport = new G4PVPlacement(0,
					 G4ThreeVector(NozzleSupportPosition_x,
						       NozzleSupportPosition_y,
						       NozzleSupportPosition_z),
					 "NozzleSupport",
					 logicNozzleSupport,
					 physiTreatmentRoom,
					 false,
					 0);

  // *************************
  // HOLE OF THE BOX SUPPORT 
  // It rapresents the nozzle 
  // in which will be inserted 
  // an air cilindrer
  // *************************
 
  G4double innerRadiusHoleNozzleSupport     = 18.    *mm;
  G4double outerRadiusHoleNozzleSupport     = 21.5   *mm;
  G4double hightHoleNozzleSupport           = 185.   *mm;
  G4double startAngleHoleNozzleSupport      = 0.     *deg;
  G4double spanningAngleHoleNozzleSupport   = 360.   *deg;

  G4double HoleNozzleSupportPosition_x      = 773.09 *mm;
  G4double HoleNozzleSupportPosition_y      = 0.     *mm;
  G4double HoleNozzleSupportPosition_z      = 0.     *mm;

  solidHoleNozzleSupport = new G4Tubs("HoleNozzleSupport",
				      innerRadiusHoleNozzleSupport,
				      outerRadiusHoleNozzleSupport,
				      hightHoleNozzleSupport,
				      startAngleHoleNozzleSupport,
				      spanningAngleHoleNozzleSupport);

  logicHoleNozzleSupport = new G4LogicalVolume(solidHoleNozzleSupport,
					       Brass,
					       "HoleNozzleSupport",
					       0,
					       0,
					       0);

  physiHoleNozzleSupport = new G4PVPlacement(G4Transform3D(rm,
							   G4ThreeVector(HoleNozzleSupportPosition_x,
									 HoleNozzleSupportPosition_y,
									 HoleNozzleSupportPosition_z)),
					     "HoleNozzleSupport",
					     logicHoleNozzleSupport,
					     physiTreatmentRoom,
					     false,
					     0); 
  
  // *********************************************
  // Second Hole of the BOX SUPPORT
  //  it is a block of air into the BoxNozzle
  //  
  // *********************************************

  G4double innerRadiusSecondHoleNozzleSupport   =    0.   *mm;
  G4double outerRadiusSecondHoleNozzleSupport   =    18.  *mm;
  G4double hightSecondHoleNozzleSupport         =    29.5 *mm;
  G4double startAngleSecondHoleNozzleSupport    =    0.   *deg;
  G4double spanningAngleSecondHoleNozzleSupport =    360. *deg;

  solidSecondHoleNozzleSupport = new G4Tubs("SecondHoleNozzleSupport",
					    innerRadiusSecondHoleNozzleSupport,
					    outerRadiusSecondHoleNozzleSupport,
					    hightSecondHoleNozzleSupport,
					    startAngleSecondHoleNozzleSupport,
					    spanningAngleSecondHoleNozzleSupport);

  logicSecondHoleNozzleSupport = new G4LogicalVolume(solidSecondHoleNozzleSupport,
						     Air,
						     "SecondHoleNozzleSupport",
						     0,
						     0,
						     0);

  physiSecondHoleNozzleSupport = new G4PVPlacement(G4Transform3D(rm,
								 G4ThreeVector()),
						   "SecondHoleNozzleSupport",
						   logicSecondHoleNozzleSupport,
						   physiNozzleSupport,
						   false,
						   0); 

  //******************************
  //     FINAL COLLIMATOR 25 mm
  //*****************************

  G4double innerRadiusFinalCollimator   = 12.5    *mm;
  G4double outerRadiusFinalCollimator   = 21.5    *mm;
  G4double hightFinalCollimator         = 3.5     *mm;
  G4double startAngleFinalCollimator    = 0.      *deg;
  G4double spanningAngleFinalCollimator = 360.    *deg;

  G4double FinalCollimatorPosition_x    = 961.59  *mm;  
  G4double FinalCollimatorPosition_y    = 0.      *mm; 
  G4double FinalCollimatorPosition_z    = 0.      *mm;


  solidFinalCollimator = new G4Tubs("FinalCollimator",
				    innerRadiusFinalCollimator,
				    outerRadiusFinalCollimator,
				    hightFinalCollimator,
				    startAngleFinalCollimator,
				    spanningAngleFinalCollimator);

  logicFinalCollimator = new G4LogicalVolume(solidFinalCollimator,
					     Brass,
					     "FinalCollimator",
					     0,
					     0,
					     0);

  physiFinalCollimator = new G4PVPlacement(G4Transform3D(rm,
							 G4ThreeVector(FinalCollimatorPosition_x,
								       FinalCollimatorPosition_y,
								       FinalCollimatorPosition_z)),
					   "FinalCollimator",
					   logicFinalCollimator,
					   physiTreatmentRoom,
					   false,
					   0); 




 
  // ---------------------------------------------------------------------------
  //
  // Modulator
  //   
  G4double innerRadiusOfTheTube = 2.5*cm;
  G4double outerRadiusOfTheTube = 9.5*cm;
  G4double hightOfTheTube = 0.03*cm;


  //------------------------------
  // Mother volume of modulator
  //------------------------------

  G4ThreeVector positionMotherMod = G4ThreeVector(-1011.91 *mm, 30 *mm, 50 *mm);
  solidMotherMod = new G4Box("MotherMod",5 *cm, 5 *cm,2 *cm);
  logicMotherMod = new G4LogicalVolume(solidMotherMod, Mod0Mater,"MotherMod",0,0,0);
  physiMotherMod = new G4PVPlacement(G4Transform3D(rm,positionMotherMod),      // rotation
				     logicMotherMod,     // its logical volume				  
				     "MotherMod",        // its name
				     logicTreatmentRoom,      // its mother  volume
				     false,           // no boolean operations
				     0);              // no particular field 
 
  
  //----------------------------------------------------------
  //Volume Madre  1/4 del Modulatore
  //----------------------------------------------------------

  G4double hightOfTheTube0 = 0.54*cm;
  G4double startAngleOfTheTube0 = 45*deg;
  G4double spanningAngleOfTheTube0 = 90*deg;

   
  G4ThreeVector positionMod0 = G4ThreeVector(0*cm,0*cm,0*cm);
  solidMod0 = new G4Tubs("Mod0",
			 innerRadiusOfTheTube, 
			 outerRadiusOfTheTube,
			 hightOfTheTube0,
			 startAngleOfTheTube0, 
			 spanningAngleOfTheTube0);

  logicMod0 = new G4LogicalVolume(solidMod0, Mod0Mater, "Mod0",0,0,0);
 


  G4RotationMatrix rm2;
  phi = 4 *deg;
  G4double angle = 90 *deg;
  for (G4int j = 0; j < 4 ; j ++)

    { 
      rm2.rotateZ(phi);
       
      physiMod0 = new G4PVPlacement(G4Transform3D(rm2, positionMod0), 
				    logicMod0,    // its logical volume			  
				    "Mod0",        // its name
				    logicMotherMod,  // its mother  volume
				    false,         // no boolean operations
				    j);            // no particular field
	   

      phi = phi + angle;

    }


  //----------------------------------------------------------
  //Prima fetta Modulatore
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
  physiMod1 = new G4PVPlacement(0,               // no rotation
				positionMod1,  // at (x,y,z)
				logicMod1,     // its logical volume				  
				"Mod1",        // its name
				logicMod0,      // its mother  volume
				false,           // no boolean operations
				0);              // no particular field


  //----------------------------------------------------------
  //Seconda fetta Modulatore
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
  physiMod2 = new G4PVPlacement(0,               // no rotation
				positionMod2,  // at (x,y,z)
				logicMod2,     // its logical volume				  
				"Mod2",        // its name
				logicMod0,      // its mother  volume
				false,           // no boolean operations
				0);              // no particular field


  //----------------------------------------------------------
  //Terza fetta Modulatore
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
  physiMod3 = new G4PVPlacement(0,                   // no rotation
				positionMod3,  // at (x,y,z)
				logicMod3,     // its logical volume				  
				"Mod3",        // its name
				logicMod0,      // its mother  volume
				false,           // no boolean operations
				0);              // no particular field

 
  //----------------------------------------------------------
  //Quarta fetta Modulatore
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

  //********************************
  //         WATER PHANTOM
  // is a cube full of water
  // in wich it is possible 
  // insert the detectors for 
  // depth and lateral dose profiles 
  // reconstruction
  //********************************

  G4double WaterPhantomDimension  = 200.    *mm;
  G4double WaterPhantomPosition_x = 1248.59 *mm;
  G4double WaterPhantomPosition_y = 0.0     *mm;
  G4double WaterPhantomPosition_z = 0.0     *mm;

  G4Box* WaterPhantom = new G4Box("WaterPhantom",
				  WaterPhantomDimension,
				  WaterPhantomDimension,
				  WaterPhantomDimension);

  logicWaterPhantom = new G4LogicalVolume(WaterPhantom,
					  Water,
					  "WaterPhantom");
  physiWaterPhantom = new G4PVPlacement(0,
					G4ThreeVector(WaterPhantomPosition_x,
						      WaterPhantomPosition_y,
						      WaterPhantomPosition_z),
					"WaterPhantom",
					logicWaterPhantom,
					physiTreatmentRoom,
					false,
					0);

  // *************************
  //     DOSEMETER
  // It is a detector placed
  // inside the water phantom
  // *************************

  G4int NbOfLayer = 20000;   //Number of the layer of the dosimeter

  G4double innerRadiusDosemeter         = 0.   *mm;         
  G4double outerRadiusDosemeter         = 3.  *mm;
  G4double hightDosemeter               = 0.001 *mm;
  G4double startAngleDosemeter          = 0.   *deg;
  G4double spanningAngleDosemeter       = 360. *deg;

  G4double DosemeterPosition_x          = (- WaterPhantomDimension + hightDosemeter);

  solidDosemeter = new G4Tubs("Dosemeter",
			      innerRadiusDosemeter,
			      outerRadiusDosemeter,
			      hightDosemeter,
			      startAngleDosemeter,
			      spanningAngleDosemeter);
 
  logicDosemeter = new G4LogicalVolume(solidDosemeter,
				       Water,
				       "Dosemeter",
				       0,
				       0,
				       0);
 
  if (NbOfLayer != 0)
    {
     
      for (G4int j = 0; j < NbOfLayer; j ++)
       
	{ 
	 
	  physiDosemeter = new G4PVPlacement(G4Transform3D(rm,
							   G4ThreeVector(DosemeterPosition_x,
									 0.,
									 0.)),				   
					     "Dosemeter",
					     logicDosemeter,
					     physiWaterPhantom,
					     false,
					     j);

	  DosemeterPosition_x = DosemeterPosition_x +  hightDosemeter*2.;
 
	}
    }
  else
    {
      physiDosemeter = new G4PVPlacement(G4Transform3D(rm,
						       G4ThreeVector(DosemeterPosition_x,
								     0.,
								     0.)),				   
					 "Dosemeter",
					 logicDosemeter,
					 physiWaterPhantom,
					 false,
					 0);
    }
 

  // **************
  // Cut by Region
  // **************
  G4Region* aRegion = new G4Region("Dosemeter");
  logicDosemeter->SetRegion(aRegion);
  aRegion->AddRootLogicalVolume(logicDosemeter);

  // *******************************                              
  // Sensitive Detectors: Dosimeter 
  // *******************************

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
 
  if(!calorimeterSD)
    {      							
      calorimeterSD = new HadrontherapyCalorimeterSD("CalorSD",this);
      SDman->AddNewDetector( calorimeterSD );
    }
  if (logicDosemeter)
    logicDosemeter->SetSensitiveDetector(calorimeterSD);
 
  // *************************                                        
  // Visualization attributes
  // *************************
  
  G4VisAttributes * yellow = new G4VisAttributes( G4Colour(1., 1., 0. ));
  yellow-> SetVisibility(true);
  yellow-> SetForceSolid(true);

  G4VisAttributes * gray = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5 ));
  gray-> SetVisibility(true);
  gray-> SetForceSolid(true);

  G4VisAttributes * red = new G4VisAttributes( G4Colour(1. ,0. ,0.));
  red-> SetVisibility(true);
  red-> SetForceSolid(true);

  G4VisAttributes * blue = new G4VisAttributes( G4Colour(0. ,0. ,1.));
  blue -> SetVisibility(true);
  blue -> SetForceSolid(true);

  G4VisAttributes * green = new G4VisAttributes( G4Colour(0. ,1. ,0.));
  green -> SetVisibility(true);
  green -> SetForceSolid(true);

  G4VisAttributes * white = new G4VisAttributes( G4Colour());
  white -> SetVisibility(true);
  white -> SetForceSolid(true);

  G4VisAttributes * brawn = new G4VisAttributes( G4Colour(0.8, 0.5, 0.35));
  brawn -> SetVisibility(true);
  brawn -> SetForceSolid(true);

  logicTreatmentRoom->SetVisAttributes (G4VisAttributes::Invisible);
  logicBeamLineSupport                 -> SetVisAttributes(gray);
  logicBeamLineCover                   -> SetVisAttributes(blue);
  logicBeamLineCover2                  -> SetVisAttributes(blue);
  logicFirstCollimatorModulatorBox     -> SetVisAttributes(blue);
  logicHoleFirstCollimatorModulatorBox ->SetVisAttributes (white);
  logicSecondCollimatorModulatorBox    -> SetVisAttributes(blue);
  logicFirstMonitorLayer1              -> SetVisAttributes(brawn);
  logicSecondMonitorLayer1             -> SetVisAttributes(brawn);
  logicThirdMonitorLayer1              -> SetVisAttributes(brawn);
  logicNozzleSupport                   -> SetVisAttributes(white);
  logicHoleNozzleSupport               -> SetVisAttributes(yellow); 
  logicFinalCollimator                 -> SetVisAttributes(red);
  //logicWaterPhantom                    -> SetVisAttributes(red);
  logicDosemeter                       -> SetVisAttributes(G4VisAttributes::Invisible);
  logicFirstCollimator                 -> SetVisAttributes(white);
  logicSecondCollimator                -> SetVisAttributes(white); 
  logicMotherMod -> SetVisAttributes( G4VisAttributes::Invisible );
  logicMod0 -> SetVisAttributes ( G4VisAttributes::Invisible );
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

  // *********************************
  // always return the physical World
  // *********************************
 
  return physiTreatmentRoom;
}

// --------------------------------------------------------------
void HadrontherapyDetectorConstruction::SetModulatorAngle(G4double val)
{
  // change Gap thickness and recompute the calorimeter parameters

  ModulatorAngle = val;
  delete   physiMotherMod;

  G4double alpha = 90 *deg;
  G4RotationMatrix rm3;
  rm3.rotateY(alpha);
  rm3.rotateX(ModulatorAngle);
  G4ThreeVector positionMotherMod = G4ThreeVector(-1011.91 *mm, 30 *mm, 50 *mm);
     
  physiMotherMod = new G4PVPlacement(G4Transform3D(rm3,positionMotherMod),      // rotation
				     logicMotherMod,     // its logical volume				  
				     "MotherMod",        // its name
				     logicTreatmentRoom,      // its mother  volume
				     false,           // no boolean operations
				     0);              // no particular field 
	  
  G4cout << "MODULATOR HAS BEEN ROTATED OF   " << ModulatorAngle << "   GRADIENTS" << G4endl;
  
   
 


  G4RunManager::GetRunManager()-> GeometryHasBeenModified();

}
// --------------------------------------------------------------



