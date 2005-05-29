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
// $Id: HadrontherapyDetectorConstruction.cc,v 3.0, September 2004;
// -----------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// -----------------------------------------------------------
// Code developed by:
//
// G. Candiano, G.A.P. Cirrone, F. Di Rosa, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// ----------------------------------------------------------
//
// Beam Line

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "HadrontherapyMaterial.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "HadrontherapyBeamLine.hh"
#include "G4Material.hh"

HadrontherapyBeamLine::HadrontherapyBeamLine(G4VPhysicalVolume* motherVolume):
  physiBeamLineSupport(0), 
  physiBeamLineCover(0), 
  physiBeamLineCover2(0),
  
  FirstScatteringFoil(0),
  physiFirstScatteringFoil(0),

  physiKaptonWindow(0), 
  solidStopper(0), 
  physiStopper(0), 
  SecondScatteringFoil(0),
  physiSecondScatteringFoil(0),
  
  physiFirstCollimator(0),
  solidRangeShifterBox(0), 
  logicRangeShifterBox(0),
  physiRangeShifterBox(0),
  
  physiSecondCollimator(0),
  physiFirstCollimatorModulatorBox(0),
  physiHoleFirstCollimatorModulatorBox(0),
  
  physiSecondCollimatorModulatorBox(0),
  physiHoleSecondCollimatorModulatorBox(0), physiFirstMonitorLayer1(0), physiFirstMonitorLayer2(0),
  physiFirstMonitorLayer3(0), physiFirstMonitorLayer4(0), physiSecondMonitorLayer1(0), physiSecondMonitorLayer2(0),
  physiSecondMonitorLayer3(0), physiSecondMonitorLayer4(0),physiThirdMonitorLayer1(0), physiThirdMonitorLayer2(0),
  physiThirdMonitorLayer3(0), physiThirdMonitorLayer4(0), physiNozzleSupport(0), physiHoleNozzleSupport(0),
  physiSecondHoleNozzleSupport(0),solidFinalCollimator(0),
  physiFinalCollimator(0)
{
  mother = motherVolume;
  material = new HadrontherapyMaterial();  
  
  G4double defaultRangeShifterBoxPosition_x = -2530.5 *mm;
  RangeShifterBoxPosition_x = defaultRangeShifterBoxPosition_x; 

  G4double defaultRangeShifterBox_x = 5 *mm; 
  RangeShifterBox_x = defaultRangeShifterBox_x;

  G4double defaultFirstScatteringFoil_x = 0.0075 *mm;
  FirstScatteringFoil_x = defaultFirstScatteringFoil_x; 
  
  G4double defaultSecondScatteringFoil_x = 0.0125 *mm;  
  SecondScatteringFoil_x = defaultSecondScatteringFoil_x;

  G4double defaultouterRadiusStopper = 2 *mm;
  outerRadiusStopper = defaultouterRadiusStopper; 

  G4double defaultinnerRadiusFinalCollimator = 12.5 *mm;
  innerRadiusFinalCollimator = defaultinnerRadiusFinalCollimator;

  RangeShifterBox_y = 176. *mm; 
  RangeShifterBox_z = 176. *mm;  

  RangeShifterBoxPosition_y = 0 *mm;
  RangeShifterBoxPosition_z = 0 *mm;

  G4Material* air = material -> GetMat("Air") ;
  RSMat = air;
}

HadrontherapyBeamLine::~HadrontherapyBeamLine()
{
 delete material;
}

void HadrontherapyBeamLine::HadrontherapyBeamLineSupport()
{
// ---------------------
// BEAM LINE SUPPORT
//--------------------

G4double BeamLineSupport_x = 1.5 *m;
G4double BeamLineSupport_y = 20. *mm;
G4double BeamLineSupport_z = 600.*mm;                               

G4double BeamLineSupportPosition_x = -1948.59 *mm;
G4double BeamLineSupportPosition_y = -230. *mm; 
G4double BeamLineSupportPosition_z = 0.*mm;

G4Material* Al = material -> GetMat("MatAluminum"); 

G4Box* BeamLineSupport = new G4Box("BeamLineSupport", 
				   BeamLineSupport_x, 
				   BeamLineSupport_y, 
				   BeamLineSupport_z);

G4LogicalVolume* logicBeamLineSupport = new G4LogicalVolume(BeamLineSupport, 
					   Al, 
					   "BeamLineSupport");
physiBeamLineSupport = new G4PVPlacement(0, G4ThreeVector(BeamLineSupportPosition_x, 
							  BeamLineSupportPosition_y,
							  BeamLineSupportPosition_z),
					 "BeamLineSupport", 
					 logicBeamLineSupport, 
					 mother, false, 0);

G4VisAttributes * gray = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5 ));
gray-> SetVisibility(true);
gray-> SetForceSolid(true);
logicBeamLineSupport -> SetVisAttributes(gray);

//----------------------
//  BEAM LINE COVER 1 (left panel)

G4double BeamLineCover_x = 1.5    *m;
G4double BeamLineCover_y = 750.   *mm;
G4double BeamLineCover_z = 10.    *mm;                               

G4double BeamLineCoverPosition_x = -1948.59 *mm;
G4double BeamLineCoverPosition_y = -980. *mm; 
G4double BeamLineCoverPosition_z =  600. *mm;

G4Box* BeamLineCover = new G4Box("BeamLineCover",
				 BeamLineCover_x, 
				 BeamLineCover_y, 
				 BeamLineCover_z);
G4LogicalVolume* logicBeamLineCover = new G4LogicalVolume(BeamLineCover, 
					 Al, 
					 "BeamLineCover");

physiBeamLineCover = new G4PVPlacement(0, G4ThreeVector(BeamLineCoverPosition_x,
							BeamLineCoverPosition_y,
							BeamLineCoverPosition_z),
				       "BeamLineCover", 
				       logicBeamLineCover, 
				       mother,
				       false, 
				       0);

// --------------------------
//  BEAM LINE COVER 2 (rigth panel)

G4double BeamLineCover2_x = 1.5    *m;
G4double BeamLineCover2_y = 750.   *mm;
G4double BeamLineCover2_z = 10.    *mm;                               

G4double BeamLineCover2Position_x = -1948.59 *mm;
G4double BeamLineCover2Position_y = -980. *mm; 
G4double BeamLineCover2Position_z =  -600. *mm;

G4Box* BeamLineCover2 = new G4Box("BeamLineCover2", 
				   BeamLineCover2_x, 
				   BeamLineCover2_y, 
				   BeamLineCover2_z);

G4LogicalVolume* logicBeamLineCover2 = new G4LogicalVolume(BeamLineCover2, 
					  Al, 
					  "BeamLineCover2");

physiBeamLineCover2 = new G4PVPlacement(0, G4ThreeVector(BeamLineCover2Position_x,
							 BeamLineCover2Position_y,
							 BeamLineCover2Position_z),
					"BeamLineCover2", 
					logicBeamLineCover2, 
					mother, 
					false, 
					0);

G4VisAttributes * blue = new G4VisAttributes( G4Colour(0. ,0. ,1.));
blue -> SetVisibility(true);
blue -> SetForceSolid(true);

logicBeamLineCover  -> SetVisAttributes(blue);
logicBeamLineCover2 -> SetVisAttributes(blue); 
}
void HadrontherapyBeamLine::HadrontherapyBeamScatteringFoils()
{
// ---------------------
 // VACUUM ZONE
 //-----------

G4double VacuumZone_x = 60.5325 *mm;
G4double VacuumZone_y = 52.5    *mm;
G4double VacuumZone_z = 52.5    *mm;                               

G4double VacuumZonePosition_x = -3188.05750 *mm;
G4double VacuumZonePosition_y = 0.         *mm; 
G4double VacuumZonePosition_z = 0.         *mm;

G4Material* Vacuum = material -> GetMat("Galactic");

G4Box* VacuumZone= new G4Box("VacuumZone", VacuumZone_x, VacuumZone_y, VacuumZone_z);
G4LogicalVolume* logicVacuumZone = new G4LogicalVolume(VacuumZone, Vacuum, "VacuumZone");

physiVacuumZone = new G4PVPlacement(0, G4ThreeVector(VacuumZonePosition_x, 
						     VacuumZonePosition_y, 
						     VacuumZonePosition_z),
				    "VacuumZone", 
				    logicVacuumZone, 
				    mother, 
				    false, 
				    0);

// ------------------------
// FIRST SCATTERING FOIL 

//***default first scattering foil-x***


//***************************************************

G4double FirstScatteringFoil_y = 52.5   *mm;
G4double FirstScatteringFoil_z = 52.5   *mm;

G4double FirstScatteringFoilPosition_x = -59.525 *mm;
G4double FirstScatteringFoilPosition_y =  0      *mm;
G4double FirstScatteringFoilPosition_z =  0      *mm;

G4Material* Ta = material -> GetMat("MatTantalum");

FirstScatteringFoil = new G4Box("FirstScatteringFoil", 
				       FirstScatteringFoil_x, 
				       FirstScatteringFoil_y, 
				       FirstScatteringFoil_z);

G4LogicalVolume* logicFirstScatteringFoil = new G4LogicalVolume(FirstScatteringFoil, 
					       Ta, 
					       "FirstScatteringFoil");
physiFirstScatteringFoil = new G4PVPlacement(0, G4ThreeVector(FirstScatteringFoilPosition_x,
							      FirstScatteringFoilPosition_y,
							      FirstScatteringFoilPosition_z),
					     "FirstScatteringFoil", 
					     logicFirstScatteringFoil, 
					     physiVacuumZone, 
					     false, 
					     0); 

// ---------------------
// KAPTON WINDOW
//-----------
G4Material* Kapton = material -> GetMat("Kapton") ;
G4Box* solidKaptonWindow = new G4Box("KaptonWindow", 0.025*mm, 5.25*cm, 5.25*cm);
G4LogicalVolume* logicKaptonWindow = new G4LogicalVolume(solidKaptonWindow, 
					Kapton, 
					"KaptonWindow");
physiKaptonWindow = new G4PVPlacement(0, G4ThreeVector(60.5075*mm, 0.*cm, 0.*cm), 
				      "KaptonWindow", 
				      logicKaptonWindow,
				      physiVacuumZone, 
				      false, 
				      0); 
G4VisAttributes * white = new G4VisAttributes( G4Colour());
white -> SetVisibility(true);
white -> SetForceSolid(true);

logicKaptonWindow    -> SetVisAttributes(white);

// --------------------
// STOPPER
//----------------------
//****default outer Radius stopper****

//********************************************

G4double phi = 90. *deg;     
// Matrix definition for a 90 deg rotation. Also used for other volumes       
G4RotationMatrix rm;               
rm.rotateY(phi);

G4double innerRadiusStopper = 0.*cm;
G4double hightStopper = 3.5*mm;
G4double startAngleStopper = 0.*deg;
G4double spanningAngleStopper = 360.*deg;

G4double StopperPosition_x = -2956.02 *mm;
G4double StopperPosition_y = 0.*m;
G4double StopperPosition_z = 0.*m;

solidStopper = new G4Tubs("Stopper", 
			  innerRadiusStopper, 
			  outerRadiusStopper, 
			  hightStopper, 
			  startAngleStopper, 
			  spanningAngleStopper);

G4Material* Brass = material -> GetMat("Brass") ;

G4LogicalVolume* logicStopper = new G4LogicalVolume(solidStopper, Brass, "Stopper", 0, 0, 0);

physiStopper = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(StopperPosition_x, 
								 StopperPosition_y, 
								 StopperPosition_z)),
				 "Stopper", 
				 logicStopper, 
				 mother, 
				 false, 
				 0);

G4VisAttributes * red = new G4VisAttributes( G4Colour(1. ,0. ,0.));
red-> SetVisibility(true);
red-> SetForceSolid(true);
logicStopper -> SetVisAttributes(red);

 // ------------------------
 // SECOND SCATTERING FOIL

//***default second scattering foil-x*********************


//*********************************************************
  
G4double SecondScatteringFoil_y = 52.5   *mm;
G4double SecondScatteringFoil_z = 52.5   *mm;

G4double SecondScatteringFoilPosition_x = -2952.51 *mm;
G4double SecondScatteringFoilPosition_y =  0         *mm;
G4double SecondScatteringFoilPosition_z =  0         *mm;                                                         
  
SecondScatteringFoil = new G4Box("SecondScatteringFoil", 
				      SecondScatteringFoil_x, 
				      SecondScatteringFoil_y,
				      SecondScatteringFoil_z);

G4LogicalVolume* logicSecondScatteringFoil = new G4LogicalVolume(SecondScatteringFoil, 
						Ta, 
						"SecondScatteringFoil");
physiSecondScatteringFoil = new G4PVPlacement(0, G4ThreeVector(SecondScatteringFoilPosition_x,
							       SecondScatteringFoilPosition_y,
							       SecondScatteringFoilPosition_z),
					      "SeconScatteringFoil", 
					      logicSecondScatteringFoil, 
					      mother, 
					      false, 

					      0);
logicSecondScatteringFoil   ->   SetVisAttributes(white);

}

void HadrontherapyBeamLine::HadrontherapyBeamCollimators()
{
 // ---------------------
 // FIRST COLLIMATOR

G4double FirstCollimator_x = 20    *mm;
G4double FirstCollimator_y = 100   *mm;
G4double FirstCollimator_z = 100   *mm;

G4double FirstCollimatorPosition_x = -2932.5   *mm;
G4double FirstCollimatorPosition_y =  0         *mm;
G4double FirstCollimatorPosition_z =  0         *mm;

 G4Material* PMMA = material -> GetMat("PMMA");

G4Box* solidFirstCollimator = new G4Box("FirstCollimator", 
				 FirstCollimator_x, 
				 FirstCollimator_y, 
				 FirstCollimator_z);


G4LogicalVolume* logicFirstCollimator = new G4LogicalVolume(solidFirstCollimator, 
					   PMMA, 
					   "FirstCollimator");

physiFirstCollimator = new G4PVPlacement(0, G4ThreeVector(FirstCollimatorPosition_x,
							  FirstCollimatorPosition_y,
							  FirstCollimatorPosition_z),
					 "FirstCollimator", 
					 logicFirstCollimator, 
					 mother, 
					 false, 
					 0);

// ----------------------------
// HOLE OF FIRST COLLIMATOR
//-----------------------------

G4double innerRadiusHoleFirstCollimator   = 0.*mm;
G4double outerRadiusHoleFirstCollimator   = 15.*mm;
G4double hightHoleFirstCollimator         = 20.*mm;
G4double startAngleHoleFirstCollimator    = 0.*deg;
G4double spanningAngleHoleFirstCollimator = 360.*deg;


G4Tubs* solidHoleFirstCollimator = new G4Tubs("HoleFirstCollimator", 
				      innerRadiusHoleFirstCollimator, 
				      outerRadiusHoleFirstCollimator,
				      hightHoleFirstCollimator, 
				      startAngleHoleFirstCollimator, 
				      spanningAngleHoleFirstCollimator);

G4Material* Air = material -> GetMat("Air") ;

G4LogicalVolume* logicHoleFirstCollimator = new G4LogicalVolume(solidHoleFirstCollimator, 
					       Air, 
					       "HoleFirstCollimator", 
					       0, 0, 0);
G4double phi = 90. *deg;     
// Matrix definition for a 90 deg rotation. Also used for other volumes       
G4RotationMatrix rm;               
rm.rotateY(phi);

physiHoleFirstCollimator = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()), 
					     "HoleFirstCollimator",
					     logicHoleFirstCollimator, 
					     physiFirstCollimator, 
					     false, 
					     0);   

//RANGE SHIFTER

  solidRangeShifterBox = new G4Box("RangeShifterBox",
                                            RangeShifterBox_x,
                                            RangeShifterBox_y,
                                            RangeShifterBox_z);

 logicRangeShifterBox = new G4LogicalVolume(solidRangeShifterBox,
                                             RSMat,
                                            "RangeShifterBox");
                                                                                              
 physiRangeShifterBox = new G4PVPlacement(0,
                                          G4ThreeVector(RangeShifterBoxPosition_x,                                                                        RangeShifterBoxPosition_y,  
                                                        RangeShifterBoxPosition_z),
                                          "RangeShifterBox",
					  logicRangeShifterBox,
                                          mother, 
                                          false,
					  0);                        

G4VisAttributes * yellow = new G4VisAttributes( G4Colour(1., 1., 0. ));
yellow-> SetVisibility(true);
yellow-> SetForceSolid(true);
 logicRangeShifterBox    ->    SetVisAttributes(yellow);
// ----------------------- 
// SECOND COLLIMATOR
//---------------------

G4double SecondCollimator_x = 20    *mm;
G4double SecondCollimator_y = 100   *mm;
G4double SecondCollimator_z = 100   *mm;

G4double SecondCollimatorPosition_x = -2028.5   *mm;
G4double SecondCollimatorPosition_y =  0         *mm;
G4double SecondCollimatorPosition_z =  0         *mm;

G4Box* solidSecondCollimator = new G4Box("SecondCollimator", 
				  SecondCollimator_x, 
				  SecondCollimator_y, 
				  SecondCollimator_z);

G4LogicalVolume* logicSecondCollimator = new G4LogicalVolume(solidSecondCollimator, 
					    PMMA, 
					    "SecondCollimator");
physiSecondCollimator = new G4PVPlacement(0, G4ThreeVector(SecondCollimatorPosition_x,
							   SecondCollimatorPosition_y,
							   SecondCollimatorPosition_z),
					  "SecondCollimator", 
					  logicSecondCollimator, 
					  mother, 
					  false, 
					  0);

// --------------------------------
// HOLE OF SECOND COLLIMATOR

G4double innerRadiusHoleSecondCollimator   = 0.*mm;
G4double outerRadiusHoleSecondCollimator   = 15.*mm;
G4double hightHoleSecondCollimator         = 20.*mm;
G4double startAngleHoleSecondCollimator    = 0.*deg;
G4double spanningAngleHoleSecondCollimator = 360.*deg;

G4Tubs* solidHoleSecondCollimator = new G4Tubs("HoleSecondCollimator", 
				       innerRadiusHoleSecondCollimator, 
				       outerRadiusHoleSecondCollimator,
				       hightHoleSecondCollimator, 
				       startAngleHoleSecondCollimator, 
				       spanningAngleHoleSecondCollimator);
G4LogicalVolume* logicHoleSecondCollimator = new G4LogicalVolume(solidHoleSecondCollimator,Air, 
						 "HoleSecondCollimator", 
						 0, 0, 0);

physiHoleSecondCollimator = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()), 
					      "HoleSecondCollimator",
					      logicHoleSecondCollimator, 
					      physiSecondCollimator, 
					      false, 
					      0); 
  
// ----------------------------
//            FIRST
//  COLLIMATOR MODULATOR BOX

G4double FirstCollimatorModulatorBox_x = 10.    *mm;
G4double FirstCollimatorModulatorBox_y = 200.   *mm;
G4double FirstCollimatorModulatorBox_z = 200.   *mm;

G4double FirstCollimatorModulatorBoxPosition_x = -2660.5   *mm;
G4double FirstCollimatorModulatorBoxPosition_y =  0.         *mm;
G4double FirstCollimatorModulatorBoxPosition_z =  0.         *mm;

G4Material* Al = material -> GetMat("MatAluminum"); 

G4Box* solidFirstCollimatorModulatorBox = new G4Box("FirstCollimatorModulatorBox", 
                                                     FirstCollimatorModulatorBox_x,
					             FirstCollimatorModulatorBox_y, 
                                                     FirstCollimatorModulatorBox_z);
G4LogicalVolume* logicFirstCollimatorModulatorBox = new G4LogicalVolume(solidFirstCollimatorModulatorBox, 
                                                                        Al, "FirstCollimatorModulatorBox");
physiFirstCollimatorModulatorBox = new G4PVPlacement(0, G4ThreeVector(FirstCollimatorModulatorBoxPosition_x,
								      FirstCollimatorModulatorBoxPosition_y,
								      FirstCollimatorModulatorBoxPosition_z),
						     "FirstCollimatorModulatorBox", logicFirstCollimatorModulatorBox,
						     mother, false, 0);

 // ---------------------------------------
 //   HOLE OF THE FIRST COLLLIMATOR 
 //                MODULATOR  BOX

G4double innerRadiusHoleFirstCollimatorModulatorBox     = 0.   *mm;
G4double outerRadiusHoleFirstCollimatorModulatorBox     = 31.  *mm;
G4double hightHoleFirstCollimatorModulatorBox           = 10.  *mm;
G4double startAngleHoleFirstCollimatorModulatorBox      = 0.   *deg;
G4double spanningAngleHoleFirstCollimatorModulatorBox   = 360. *deg;

G4Tubs* solidHoleFirstCollimatorModulatorBox  = new G4Tubs("HoleFirstCollimatorModulatorBox",
						   innerRadiusHoleFirstCollimatorModulatorBox,
						   outerRadiusHoleFirstCollimatorModulatorBox,
						   hightHoleFirstCollimatorModulatorBox ,
						   startAngleHoleFirstCollimatorModulatorBox,
						   spanningAngleHoleFirstCollimatorModulatorBox);
G4LogicalVolume* logicHoleFirstCollimatorModulatorBox = new G4LogicalVolume(solidHoleFirstCollimatorModulatorBox,
							   Air, "HoleFirstCollimatorModulatorBox", 0, 0, 0);
physiHoleFirstCollimatorModulatorBox = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
							 "HoleFirstCollimatorModulatorBox",
							 logicHoleFirstCollimatorModulatorBox,
							 physiFirstCollimatorModulatorBox, false, 0); 
// ------------------------------
//            SECOND
//  COLLIMATOR MODULATOR BOX

G4double SecondCollimatorModulatorBox_x = 10    *mm;
G4double SecondCollimatorModulatorBox_y = 200   *mm;
G4double SecondCollimatorModulatorBox_z = 200   *mm;

G4double SecondCollimatorModulatorBoxPosition_x = -2090.5   *mm;
G4double SecondCollimatorModulatorBoxPosition_y =  0         *mm;
G4double SecondCollimatorModulatorBoxPosition_z =  0         *mm;

G4Box* solidSecondCollimatorModulatorBox= new G4Box("SecondCollimatorModulatorBox",
					     SecondCollimatorModulatorBox_x,
					     SecondCollimatorModulatorBox_y,
					     SecondCollimatorModulatorBox_z);
G4LogicalVolume* logicSecondCollimatorModulatorBox = new G4LogicalVolume(solidSecondCollimatorModulatorBox,Al, 
                                                                         "SecondCollimatorModulatorBox");

physiSecondCollimatorModulatorBox = new G4PVPlacement(0, G4ThreeVector(SecondCollimatorModulatorBoxPosition_x,
								       SecondCollimatorModulatorBoxPosition_y,
								       SecondCollimatorModulatorBoxPosition_z),
						      "SecondCollimatorModulatorBox", logicSecondCollimatorModulatorBox,
						      mother, false, 0);
  
 // -------------------------------------------
 //   HOLE OF THE SECOND  COLLLIMATOR 
 //                MODULATOR  BOX

G4double innerRadiusHoleSecondCollimatorModulatorBox     = 0.   *mm;
G4double outerRadiusHoleSecondCollimatorModulatorBox     = 31.  *mm;
G4double hightHoleSecondCollimatorModulatorBox           = 10.  *mm;
G4double startAngleHoleSecondCollimatorModulatorBox      = 0.   *deg;
G4double spanningAngleHoleSecondCollimatorModulatorBox   = 360. *deg;

G4Tubs* solidHoleSecondCollimatorModulatorBox  = new G4Tubs("HoleSecondCollimatorModulatorBox",
						    innerRadiusHoleSecondCollimatorModulatorBox,
						    outerRadiusHoleSecondCollimatorModulatorBox,
						    hightHoleSecondCollimatorModulatorBox ,
						    startAngleHoleSecondCollimatorModulatorBox,
						    spanningAngleHoleSecondCollimatorModulatorBox);

G4LogicalVolume* logicHoleSecondCollimatorModulatorBox = new G4LogicalVolume(solidHoleSecondCollimatorModulatorBox, Air,
							    "HoleSecondCollimatorModulatorBox", 0, 0, 0);
physiHoleSecondCollimatorModulatorBox = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
							  "HoleSecondCollimatorModulatorBox",
							  logicHoleSecondCollimatorModulatorBox,
							  physiSecondCollimatorModulatorBox, false, 0); 

G4VisAttributes * blue = new G4VisAttributes( G4Colour(0. ,0. ,1.));
blue -> SetVisibility(true);
blue -> SetForceSolid(true);

logicFirstCollimator   ->    SetVisAttributes(yellow);
logicSecondCollimator    ->    SetVisAttributes(yellow);
logicFirstCollimatorModulatorBox    ->    SetVisAttributes(blue);
logicSecondCollimatorModulatorBox    ->    SetVisAttributes(blue);

}
void HadrontherapyBeamLine::HadrontherapyBeamMonitoring()
{
 // ----------------------------
 // FIRST MONITOR CHAMBER
  
// Each chamber consist of 9 mm of air in a box
// that has two layers one of kapton and one
// of copper
 G4Material* Kapton = material -> GetMat("Kapton") ;
 G4Material* Cu = material -> GetMat("MatCopper");
 G4Material* Air = material -> GetMat("Air") ;

G4Box* solidFirstMonitorLayer1 = new G4Box("FirstMonitorLayer1", 4.525022*mm, 10*cm, 10*cm);
G4LogicalVolume* logicFirstMonitorLayer1 = new G4LogicalVolume(solidFirstMonitorLayer1, Kapton, "FirstMonitorLayer1");
physiFirstMonitorLayer1 = new G4PVPlacement(0,G4ThreeVector(-1765.97498 *mm,0.*cm,0.*cm),
					    "FirstMonitorLayer1", logicFirstMonitorLayer1, mother, false, 0);

G4Box* solidFirstMonitorLayer2 = new G4Box("FirstMonitorLayer2", 0.000011*mm, 10*cm, 10*cm);
G4LogicalVolume* logicFirstMonitorLayer2 = new G4LogicalVolume(solidFirstMonitorLayer2, Cu, "FirstMonitorLayer2");
physiFirstMonitorLayer2 = new G4PVPlacement(0, G4ThreeVector(-4.500011*mm,0.*cm,0.*cm),
					    "FirstMonitorLayer2", logicFirstMonitorLayer2, physiFirstMonitorLayer1,
					    false, 0);
  
G4Box* solidFirstMonitorLayer3 = new G4Box("FirstMonitorLayer3", 4.5*mm, 10*cm, 10*cm);
G4LogicalVolume* logicFirstMonitorLayer3 = new G4LogicalVolume(solidFirstMonitorLayer3, Air, "FirstMonitorLayer3");
physiFirstMonitorLayer3 = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*cm,0.*cm), "MonitorLayer3",
					    logicFirstMonitorLayer3, physiFirstMonitorLayer1, false, 0);
    
G4Box* solidFirstMonitorLayer4 = new G4Box("FirstMonitorLayer4", 0.000011*mm, 10*cm, 10*cm);
G4LogicalVolume* logicFirstMonitorLayer4 = new G4LogicalVolume(solidFirstMonitorLayer4, Cu, "FirstMonitorLayer4");
physiFirstMonitorLayer4 = new G4PVPlacement(0, G4ThreeVector(4.500011*mm,0.*cm,0.*cm), "FirstMonitorLayer4",
					    logicFirstMonitorLayer4, physiFirstMonitorLayer1, false, 0);
 // -----------------------------
 // SECOND MONITOR CHAMBER

G4Box* solidSecondMonitorLayer1 = new G4Box("SecondMonitorLayer1", 4.525022*mm, 10*cm, 10*cm);
G4LogicalVolume* logicSecondMonitorLayer1 = new G4LogicalVolume(solidSecondMonitorLayer1, Kapton, "SecondMonitorLayer1");
physiSecondMonitorLayer1 = new G4PVPlacement(0, G4ThreeVector(-1634.92493 *mm,0.*cm,0.*cm),
					    "SecondMonitorLayer1", logicSecondMonitorLayer1,mother, false, 0);

G4Box* solidSecondMonitorLayer2 = new G4Box("SecondMonitorLayer2", 0.000011*mm, 10*cm,10*cm);
G4LogicalVolume* logicSecondMonitorLayer2 = new G4LogicalVolume(solidSecondMonitorLayer2, Cu, "SecondMonitorLayer2");
physiSecondMonitorLayer2 = new G4PVPlacement(0, G4ThreeVector(-4.500011*mm,0.*cm,0.*cm), "SecondMonitorLayer2",
					     logicSecondMonitorLayer2, physiSecondMonitorLayer1, false, 0);
  
G4Box* solidSecondMonitorLayer3 = new G4Box("SecondMonitorLayer3", 4.5*mm, 10*cm, 10*cm);
G4LogicalVolume* logicSecondMonitorLayer3 = new G4LogicalVolume(solidSecondMonitorLayer3, Air, "SecondMonitorLayer3");
physiSecondMonitorLayer3 = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*cm,0.*cm), "MonitorLayer3",
					     logicSecondMonitorLayer3, physiSecondMonitorLayer1, false, 0);
    
G4Box* solidSecondMonitorLayer4 = new G4Box("SecondMonitorLayer4", 0.000011*mm, 10*cm,10*cm);
G4LogicalVolume* logicSecondMonitorLayer4 = new G4LogicalVolume(solidSecondMonitorLayer4, Cu, "SecondMonitorLayer4");
physiSecondMonitorLayer4 = new G4PVPlacement(0, G4ThreeVector(4.500011*mm,0.*cm,0.*cm), "SecondMonitorLayer4",
					     logicSecondMonitorLayer4, physiSecondMonitorLayer1, false, 0);

// ---------------------------
// THIRD MONITOR CHAMBER

G4Box* solidThirdMonitorLayer1 = new G4Box("ThirdMonitorLayer1", 4.525022*mm, 10*cm, 10*cm);
G4LogicalVolume* logicThirdMonitorLayer1 = new G4LogicalVolume(solidThirdMonitorLayer1, Kapton, "ThirdMonitorLayer1");
physiThirdMonitorLayer1 = new G4PVPlacement(0, G4ThreeVector(-1505.87489 *mm,0.*cm,0.*cm),
					    "ThirdMonitorLayer1", logicThirdMonitorLayer1, mother, false, 0);

G4Box* solidThirdMonitorLayer2 = new G4Box("ThirdMonitorLayer2", 0.000011*mm, 10*cm, 10*cm);
G4LogicalVolume* logicThirdMonitorLayer2 = new G4LogicalVolume(solidThirdMonitorLayer2, Cu, "ThirdMonitorLayer2");
physiThirdMonitorLayer2 = new G4PVPlacement(0, G4ThreeVector(-4.500011*mm,0.*cm,0.*cm), "ThirdMonitorLayer2",
					    logicThirdMonitorLayer2, physiThirdMonitorLayer1, false, 0);
 
G4Box* solidThirdMonitorLayer3 = new G4Box("ThirdMonitorLayer3", 4.5*mm, 10*cm, 10*cm);
G4LogicalVolume* logicThirdMonitorLayer3 = new G4LogicalVolume(solidThirdMonitorLayer3, Air,"ThirdMonitorLayer3");
physiThirdMonitorLayer3 = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*cm,0.*cm), "MonitorLayer3",
					    logicThirdMonitorLayer3, physiThirdMonitorLayer1, false, 0);
  
G4Box* solidThirdMonitorLayer4 = new G4Box("ThirdMonitorLayer4", 0.000011*mm, 10*cm, 10*cm);
G4LogicalVolume* logicThirdMonitorLayer4 = new G4LogicalVolume(solidThirdMonitorLayer4, Cu, "ThirdMonitorLayer4");
physiThirdMonitorLayer4 = new G4PVPlacement(0, G4ThreeVector(4.500011*mm,0.*cm,0.*cm), "ThirdMonitorLayer4",
					    logicThirdMonitorLayer4, physiSecondMonitorLayer1, false, 0);
   }
 void HadrontherapyBeamLine::HadrontherapyBeamNozzle()
   {

// ---------------------
// NOZZLE SUPPORT
//--------------------
G4double NozzleSupport_x =  29.5 *mm;
G4double NozzleSupport_y =  180. *mm;
G4double NozzleSupport_z =  180. *mm;

G4double NozzleSupportPosition_x = -601.00 *mm;
G4double NozzleSupportPosition_y = 0. *mm;
G4double NozzleSupportPosition_z = 0. *mm;

G4Material* PMMA = material -> GetMat("PMMA");
G4Material* Brass = material -> GetMat("Brass") ;
G4Material* Air = material -> GetMat("Air") ;

G4double phi = 90. *deg;     
// Matrix definition for a 90 deg rotation. Also used for other volumes       
G4RotationMatrix rm;               
rm.rotateY(phi);

G4Box* solidNozzleSupport = new G4Box("NozzlSupport", NozzleSupport_x, NozzleSupport_y, NozzleSupport_z);
G4LogicalVolume* logicNozzleSupport = new G4LogicalVolume(solidNozzleSupport, PMMA, "NozzleSupport");
physiNozzleSupport = new G4PVPlacement(0, G4ThreeVector(NozzleSupportPosition_x, NozzleSupportPosition_y, 
                                                        NozzleSupportPosition_z),
			       "NozzleSupport", logicNozzleSupport, mother, false, 0);

// --------------------------
// HOLE OF THE BOX SUPPORT 
 
G4double innerRadiusHoleNozzleSupport     = 18.    *mm;
G4double outerRadiusHoleNozzleSupport     = 21.5   *mm;
G4double hightHoleNozzleSupport           = 185.   *mm;
G4double startAngleHoleNozzleSupport      = 0.     *deg;
G4double spanningAngleHoleNozzleSupport   = 360.   *deg;

G4double HoleNozzleSupportPosition_x      = -475.5 *mm;
G4double HoleNozzleSupportPosition_y      = 0.     *mm;
G4double HoleNozzleSupportPosition_z      = 0.     *mm;

G4Tubs* solidHoleNozzleSupport = new G4Tubs("HoleNozzleSupport", innerRadiusHoleNozzleSupport, outerRadiusHoleNozzleSupport,
				    hightHoleNozzleSupport, startAngleHoleNozzleSupport, spanningAngleHoleNozzleSupport);
G4LogicalVolume* logicHoleNozzleSupport = new G4LogicalVolume(solidHoleNozzleSupport, Brass, "HoleNozzleSupport", 0, 0, 0);
physiHoleNozzleSupport = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(HoleNozzleSupportPosition_x,
								       HoleNozzleSupportPosition_y,
								       HoleNozzleSupportPosition_z)),
					   "HoleNozzleSupport", logicHoleNozzleSupport, mother, false, 0); 

 //---------------------------------------
 // Second Hole of the BOX SUPPORT

G4double innerRadiusSecondHoleNozzleSupport   =    0.   *mm;
G4double outerRadiusSecondHoleNozzleSupport   =    18.  *mm;
G4double hightSecondHoleNozzleSupport         =    29.5 *mm;
G4double startAngleSecondHoleNozzleSupport    =    0.   *deg;
G4double spanningAngleSecondHoleNozzleSupport =    360. *deg;

G4Tubs* solidSecondHoleNozzleSupport = new G4Tubs("SecondHoleNozzleSupport", innerRadiusSecondHoleNozzleSupport,
					  outerRadiusSecondHoleNozzleSupport, hightSecondHoleNozzleSupport,
					  startAngleSecondHoleNozzleSupport, spanningAngleSecondHoleNozzleSupport);
G4LogicalVolume* logicSecondHoleNozzleSupport = new G4LogicalVolume(solidSecondHoleNozzleSupport, Air, "SecondHoleNozzleSupport", 0, 0, 0);
physiSecondHoleNozzleSupport = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()), "SecondHoleNozzleSupport",
						 logicSecondHoleNozzleSupport, physiNozzleSupport, false, 0); 


G4VisAttributes * yellow = new G4VisAttributes( G4Colour(1., 1., 0. ));
yellow-> SetVisibility(true);
yellow-> SetForceSolid(true);
logicHoleNozzleSupport    ->    SetVisAttributes(yellow); 
   }

 void HadrontherapyBeamLine::HadrontherapyBeamFinalCollimator()
   {

 // --------------------------------
 //     FINAL COLLIMATOR 25 mm
 //-------------------------
//*****default inner radius final collimator********************


//**************************************************************

G4double outerRadiusFinalCollimator   = 21.5    *mm;
G4double hightFinalCollimator         = 3.5     *mm;
G4double startAngleFinalCollimator    = 0.      *deg;
G4double spanningAngleFinalCollimator = 360.    *deg;
G4double FinalCollimatorPosition_x    = -287.0  *mm;  
G4double FinalCollimatorPosition_y    = 0.      *mm; 
G4double FinalCollimatorPosition_z    = 0.      *mm;

G4double phi = 90. *deg;     

// Matrix definition for a 90 deg rotation. Also used for other volumes       
G4RotationMatrix rm;               
rm.rotateY(phi);

G4Material* Brass = material -> GetMat("Brass");

solidFinalCollimator = new G4Tubs("FinalCollimator", innerRadiusFinalCollimator, outerRadiusFinalCollimator,
				  hightFinalCollimator, 
				  startAngleFinalCollimator, spanningAngleFinalCollimator);

G4LogicalVolume* logicFinalCollimator = new G4LogicalVolume(solidFinalCollimator, 
							    Brass, "FinalCollimator", 0, 0, 0);

physiFinalCollimator = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(FinalCollimatorPosition_x,
									 FinalCollimatorPosition_y,
									 FinalCollimatorPosition_z)),
					 "FinalCollimator", logicFinalCollimator, mother, false, 0); 

 
G4VisAttributes * yellow = new G4VisAttributes( G4Colour(1., 1., 0. ));
yellow-> SetVisibility(true);
yellow-> SetForceSolid(true);
logicFinalCollimator -> SetVisAttributes(yellow); 
   }


void HadrontherapyBeamLine::setRangeShifterXPos(G4double val)
{
  G4double value = val;
  physiRangeShifterBox -> SetTranslation(G4ThreeVector
					 (value, RangeShifterBoxPosition_y, RangeShifterBoxPosition_z)); 
  G4cout << "The Range Shifter is translated to"<< value/mm <<"mm along the X axis" <<G4endl;
}

void HadrontherapyBeamLine::setRangeShifterX(G4double val)
{
  G4double value = val;
  solidRangeShifterBox -> SetXHalfLength(value) ;
  G4cout << "RangeShifter size X (mm): "<<solidRangeShifterBox -> GetXHalfLength()/mm
         << G4endl;
}

void HadrontherapyBeamLine::SetFirstScatteringFoil(G4double val)
{ 
 G4double value = val;
 FirstScatteringFoil -> SetXHalfLength(value);
 G4cout <<"The X size of the first scattering foil is (mm):"<<  
   FirstScatteringFoil -> GetXHalfLength()/mm 
     << G4endl;
}

void HadrontherapyBeamLine::SetSecondScatteringFoil(G4double val)
{
  G4double value = val; 
  SecondScatteringFoil -> SetXHalfLength(value);
  G4cout <<"The X size of the second scattering foil is (mm):"<<  
  SecondScatteringFoil -> GetXHalfLength()/mm 
     << G4endl;
}

void HadrontherapyBeamLine::SetOuterRadiusStopper(G4double val)
{ 
  G4double value = val;
 solidStopper -> SetOuterRadius(value);
 G4cout << "OuterRadius od the Stopper is (mm):"
        << solidStopper -> GetOuterRadius()/mm
        << G4endl;
}

void HadrontherapyBeamLine::SetInnerRadiusFinalCollimator(G4double val)
{
  // G4double value = val;
 solidFinalCollimator -> SetInnerRadius(val);
 G4cout<<"Inner Radius of the final collimator is (mm):"
       <<solidFinalCollimator -> GetInnerRadius()/mm
       <<G4endl; 
}
void HadrontherapyBeamLine::SetRSMaterial(G4String materialChoice)
{
 G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
 
 if (pttoMaterial) 
    {
      RSMat  = pttoMaterial;
      logicRangeShifterBox -> SetMaterial(pttoMaterial);  
    }
}
