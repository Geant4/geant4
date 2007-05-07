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
// $Id: HadrontherapyBeamLine.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

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
#include "G4SubtractionSolid.hh"

HadrontherapyBeamLine::HadrontherapyBeamLine(G4VPhysicalVolume* motherVolume):
  physiBeamLineSupport(0),
  firstScatteringFoil(0), physiFirstScatteringFoil(0), physiKaptonWindow(0), 
  solidStopper(0), physiStopper(0), 
  secondScatteringFoil(0), physiSecondScatteringFoil(0),
  physiFirstCollimator(0), solidRangeShifterBox(0), logicRangeShifterBox(0),
  physiRangeShifterBox(0), physiSecondCollimator(0),
  physiFirstCollimatorModulatorBox(0),
  physiHoleFirstCollimatorModulatorBox(0),
  physiSecondCollimatorModulatorBox(0),
  physiHoleSecondCollimatorModulatorBox(0), 
  physiFirstMonitorLayer1(0), physiFirstMonitorLayer2(0),
  physiFirstMonitorLayer3(0), physiFirstMonitorLayer4(0),
  physiSecondMonitorLayer1(0), physiSecondMonitorLayer2(0),
  physiSecondMonitorLayer3(0), physiSecondMonitorLayer4(0),
  physiThirdMonitorLayer1(0), physiThirdMonitorLayer2(0),
  physiThirdMonitorLayer3(0), physiThirdMonitorLayer4(0),
  physiNozzleSupport(0), physiHoleNozzle(0),
  solidFinalCollimator(0),
  physiFinalCollimator(0)
{
  mother = motherVolume;
  material = new HadrontherapyMaterial();  

  G4double defaultFirstScatteringFoilXSize = 0.0075 *mm;
  firstScatteringFoilXSize = defaultFirstScatteringFoilXSize; 
  
  G4double defaultOuterRadiusStopper = 2 *mm;
  outerRadiusStopper = defaultOuterRadiusStopper; 

  G4double defaultSecondScatteringFoilXSize = 0.0125 *mm;  
  secondScatteringFoilXSize = defaultSecondScatteringFoilXSize;

  G4double defaultRangeShifterXSize = 5. *mm; 
  rangeShifterXSize = defaultRangeShifterXSize;

  G4double defaultRangeShifterXPosition = -2530.5 *mm;
  rangeShifterXPosition = defaultRangeShifterXPosition; 

  G4double defaultinnerRadiusFinalCollimator = 12.5 *mm;
  innerRadiusFinalCollimator = defaultinnerRadiusFinalCollimator;

  G4Material* air = material -> GetMat("Air") ;
  RSMat = air;
}

HadrontherapyBeamLine::~HadrontherapyBeamLine()
{
  delete material;
}

void HadrontherapyBeamLine::HadrontherapyBeamLineSupport()
{
  // ------------------//
  // Beam line support //
  //-------------------//

  const G4double beamLineSupportXSize = 1.5*m;
  const G4double beamLineSupportYSize = 20.*mm;
  const G4double beamLineSupportZSize = 600.*mm;                               

  const G4double beamLineSupportXPosition = -1948.59 *mm;
  const G4double beamLineSupportYPosition = -230. *mm; 
  const G4double beamLineSupportZPosition = 0.*mm;

  G4Material* Al = material -> GetMat("MatAluminum"); 

  G4Box* beamLineSupport = new G4Box("BeamLineSupport", 
				     beamLineSupportXSize, 
				     beamLineSupportYSize, 
				     beamLineSupportZSize);

  G4LogicalVolume* logicBeamLineSupport = new G4LogicalVolume(beamLineSupport, 
							      Al, 
							      "BeamLineSupport");
  physiBeamLineSupport = new G4PVPlacement(0, G4ThreeVector(beamLineSupportXPosition, 
							    beamLineSupportYPosition,
							    beamLineSupportZPosition),
					   "BeamLineSupport", 
					   logicBeamLineSupport, 
					   mother, false, 0);
 
  // Visualisation attributes of the beam line support
  G4VisAttributes * gray = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5 ));
  gray-> SetVisibility(true);
  gray-> SetForceSolid(true);
  logicBeamLineSupport -> SetVisAttributes(gray);

  
}

void HadrontherapyBeamLine::HadrontherapyBeamScatteringFoils()
{
  // ------------//
  // Vacuum area //
  //-------------//

  const G4double vacuumZoneXSize = 60.5325 *mm;
  const G4double vacuumZoneYSize = 52.5 *mm;
  const G4double vacuumZoneZSize = 52.5 *mm;                               

  const G4double vacuumZoneXPosition = -3188.05750 *mm;

  G4Material* vacuum = material -> GetMat("Galactic");

  G4Box* vacuumZone = new G4Box("VacuumZone", vacuumZoneXSize, vacuumZoneYSize, vacuumZoneZSize);

  G4LogicalVolume* logicVacuumZone = new G4LogicalVolume(vacuumZone, vacuum, "VacuumZone");

  G4VPhysicalVolume* physiVacuumZone = new G4PVPlacement(0, G4ThreeVector(vacuumZoneXPosition, 0., 0.),
				      "VacuumZone", 
				      logicVacuumZone, 
				      mother, 
				      false, 
				      0);
  // ----------------------//
  // First scattering foil //
  // ----------------------//

  const G4double firstScatteringFoilYSize = 52.5   *mm;
  const G4double firstScatteringFoilZSize = 52.5   *mm;

  const G4double firstScatteringFoilXPosition = -59.525 *mm;

  G4Material* Ta = material -> GetMat("MatTantalum");

  firstScatteringFoil = new G4Box("FirstScatteringFoil", 
				  firstScatteringFoilXSize, 
				  firstScatteringFoilYSize, 
				  firstScatteringFoilZSize);

  G4LogicalVolume* logicFirstScatteringFoil = new G4LogicalVolume(firstScatteringFoil, 
								  Ta, 
								  "FirstScatteringFoil");

  physiFirstScatteringFoil = new G4PVPlacement(0, G4ThreeVector(firstScatteringFoilXPosition,
								0.,
								0.),
					       "FirstScatteringFoil", 
					       logicFirstScatteringFoil, 
					       physiVacuumZone, 
					       false, 
					       0); 

  // --------------//
  // Kapton Window //
  //---------------//

  G4Material* kapton = material -> GetMat("Kapton") ;

  const G4double kaptonWindowXSize = 0.025*mm;
  const G4double kaptonWindowYSize = 5.25*cm;
  const G4double kaptonWindowZSize = 5.25*cm;
  const G4double kaptonWindowXPosition = 60.5075*mm;
  
  G4Box* solidKaptonWindow = new G4Box("KaptonWindow", 
				       kaptonWindowXSize,
				       kaptonWindowYSize,
				       kaptonWindowZSize);

  G4LogicalVolume* logicKaptonWindow = new G4LogicalVolume(solidKaptonWindow, 
							   kapton, 
							   "KaptonWindow");

  physiKaptonWindow = new G4PVPlacement(0, G4ThreeVector(kaptonWindowXPosition, 0., 0.), 
					"KaptonWindow", 
					logicKaptonWindow,
					physiVacuumZone, 
					false, 
					0); 

  G4VisAttributes * white = new G4VisAttributes( G4Colour());
  white -> SetVisibility(true);
  white -> SetForceSolid(true);

  logicKaptonWindow -> SetVisAttributes(white);

  // --------//
  // Stopper //
  //---------//
  
  G4double phi = 90. *deg;     
  // Matrix definition for a 90 deg rotation with respect to Y axis
  G4RotationMatrix rm;               
  rm.rotateY(phi);

  const G4double innerRadiusStopper = 0.*cm;
  const G4double hightStopper = 3.5*mm;
  const G4double startAngleStopper = 0.*deg;
  const G4double spanningAngleStopper = 360.*deg;

  const G4double stopperXPosition = -2956.04 *mm;
  const G4double stopperYPosition = 0.*m;
  const G4double stopperZPosition = 0.*m;

  solidStopper = new G4Tubs("Stopper", 
			    innerRadiusStopper, 
			    outerRadiusStopper, 
			    hightStopper, 
			    startAngleStopper, 
			    spanningAngleStopper);

  G4Material* brass = material -> GetMat("Brass") ;

  G4LogicalVolume* logicStopper = new G4LogicalVolume(solidStopper, brass, "Stopper", 0, 0, 0);

  physiStopper = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(stopperXPosition, 
								   stopperYPosition, 
								   stopperZPosition)),
				   "Stopper", 
				   logicStopper, 
				   mother, 
				   false, 
				   0);

  G4VisAttributes * red = new G4VisAttributes(G4Colour(1. ,0. ,0.));
  red-> SetVisibility(true);
  red-> SetForceSolid(true);
  logicStopper -> SetVisAttributes(red);

  // -----------------------//
  // Second scattering foil //
  // -----------------------//
  
  const G4double secondScatteringFoilYSize = 52.5   *mm;
  const G4double secondScatteringFoilZSize = 52.5   *mm;

  const G4double secondScatteringFoilXPosition = -2952.52 *mm;
  const G4double secondScatteringFoilYPosition =  0         *mm;
  const G4double secondScatteringFoilZPosition =  0         *mm;                                                         
  
  secondScatteringFoil = new G4Box("SecondScatteringFoil", 
				   secondScatteringFoilXSize, 
				   secondScatteringFoilYSize,
				   secondScatteringFoilZSize);

  G4LogicalVolume* logicSecondScatteringFoil = new G4LogicalVolume(secondScatteringFoil, 
								   Ta, 
								   "SecondScatteringFoil");

  physiSecondScatteringFoil = new G4PVPlacement(0, G4ThreeVector(secondScatteringFoilXPosition,
								 secondScatteringFoilYPosition,
								 secondScatteringFoilZPosition),
						"SeconScatteringFoil", 
						logicSecondScatteringFoil, 
						mother, 
						false, 
						0);

  logicSecondScatteringFoil -> SetVisAttributes(white);
}

void HadrontherapyBeamLine::HadrontherapyBeamCollimators()
{
  // -----------------//
  // First collimator //
  // -----------------//

  const G4double firstCollimatorXSize = 20.*mm;
  const G4double firstCollimatorYSize = 100.*mm;
  const G4double firstCollimatorZSize = 100.*mm;

  const G4double firstCollimatorXPosition = -2932.5*mm;
  const G4double firstCollimatorYPosition = 0.*mm;
  const G4double firstCollimatorZPosition = 0.*mm;

  G4Material* PMMA = material -> GetMat("PMMA");

  G4Box* solidFirstCollimator = new G4Box("FirstCollimator", 
					  firstCollimatorXSize, 
					  firstCollimatorYSize, 
					  firstCollimatorZSize);


  G4LogicalVolume* logicFirstCollimator = new G4LogicalVolume(solidFirstCollimator, 
							      PMMA, 
							      "FirstCollimator");

  physiFirstCollimator = new G4PVPlacement(0, G4ThreeVector(firstCollimatorXPosition,
							    firstCollimatorYPosition,
							    firstCollimatorZPosition),
					   "FirstCollimator", 
					   logicFirstCollimator, 
					   mother, 
					   false, 
					   0);

  // ----------------------------//
  // Hole of the first collimator//
  //-----------------------------//

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
  // --------------//
  // Range shifter //
  // --------------// 

  const G4double rangeShifterYSize = 176. *mm; 
  const G4double rangeShifterZSize = 176. *mm;  

  solidRangeShifterBox = new G4Box("RangeShifterBox",
				   rangeShifterXSize,
				   rangeShifterYSize,
				   rangeShifterZSize);

  logicRangeShifterBox = new G4LogicalVolume(solidRangeShifterBox,
                                             RSMat,
					     "RangeShifterBox");
                                                                                              
  physiRangeShifterBox = new G4PVPlacement(0,
					   G4ThreeVector(rangeShifterXPosition, 0., 0.),                                               
					   "RangeShifterBox",
					   logicRangeShifterBox,
					   mother, 
					   false,
					   0);                        

  G4VisAttributes * yellow = new G4VisAttributes(G4Colour(1., 1., 0. ));
  yellow-> SetVisibility(true);
  yellow-> SetForceSolid(true);
  logicRangeShifterBox -> SetVisAttributes(yellow);

  // ------------------// 
  // Second collimator //
  //-------------------.//

  const G4double secondCollimatorXPosition = -2028.5*mm;
  const G4double secondCollimatorYPosition =  0*mm;
  const G4double secondCollimatorZPosition =  0*mm;

  const G4double secondCollimatorXSize = 20.*mm;
  const G4double secondCollimatorYSize = 100.*mm;
  const G4double secondCollimatorZSize = 100.*mm;
 
  G4Box* solidSecondCollimator = new G4Box("SecondCollimator", 
					  secondCollimatorXSize, 
					  secondCollimatorYSize, 
					  secondCollimatorZSize);


  G4LogicalVolume* logicSecondCollimator = new G4LogicalVolume(solidSecondCollimator, 
							      PMMA, 
							      "SecondCollimator");

  physiSecondCollimator = new G4PVPlacement(0, G4ThreeVector(secondCollimatorXPosition,
							    secondCollimatorYPosition,
							    secondCollimatorZPosition),
					   "SecondCollimator", 
					   logicSecondCollimator, 
					   mother, 
					   false, 
					   0);



  // ------------------------------//
  // Hole of the second collimator //
  // ------------------------------//

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


  G4LogicalVolume* logicHoleSecondCollimator = new G4LogicalVolume(solidHoleSecondCollimator, 
								  Air, 
								  "HoleSecondCollimator", 
								  0, 0, 0);
  G4double phi2 = 90. *deg;     
  // Matrix definition for a 90 deg rotation. Also used for other volumes       
  G4RotationMatrix rm2;               
  rm2.rotateY(phi2);

  physiHoleSecondCollimator = new G4PVPlacement(G4Transform3D(rm2, G4ThreeVector()), 
					       "HoleSecondCollimator",
					       logicHoleSecondCollimator, 
					       physiSecondCollimator, 
					       false, 
					       0);   


  
  // ---------------------------------//
  //  First Collimator modulator box //
  // ---------------------------------//

  const G4double firstCollimatorModulatorXSize = 10.*mm;
  const G4double firstCollimatorModulatorYSize = 200.*mm;
  const G4double firstCollimatorModulatorZSize = 200.*mm;

  const G4double firstCollimatorModulatorXPosition = -2660.5*mm;
  const G4double firstCollimatorModulatorYPosition = 0.*mm;
  const G4double firstCollimatorModulatorZPosition = 0.*mm;

  G4Material* Al = material -> GetMat("MatAluminum"); 

  G4Box* solidFirstCollimatorModulatorBox = new G4Box("FirstCollimatorModulatorBox", 
						      firstCollimatorModulatorXSize,
						      firstCollimatorModulatorYSize, 
						      firstCollimatorModulatorZSize);

  G4LogicalVolume* logicFirstCollimatorModulatorBox = new G4LogicalVolume(solidFirstCollimatorModulatorBox, 
									  Al, "FirstCollimatorModulatorBox");

  physiFirstCollimatorModulatorBox = new G4PVPlacement(0, G4ThreeVector(firstCollimatorModulatorXPosition,
									firstCollimatorModulatorYPosition,
									firstCollimatorModulatorZPosition),
						       "FirstCollimatorModulatorBox",
						       logicFirstCollimatorModulatorBox,
						       mother, false, 0);

  // ----------------------------------//
  //   Hole of the first modulator box //
  // ----------------------------------//
  const G4double innerRadiusHoleFirstCollimatorModulatorBox = 0.*mm;
  const G4double outerRadiusHoleFirstCollimatorModulatorBox = 31.*mm;
  const G4double hightHoleFirstCollimatorModulatorBox = 10.*mm;
  const G4double startAngleHoleFirstCollimatorModulatorBox = 0.*deg;
  const G4double spanningAngleHoleFirstCollimatorModulatorBox = 360.*deg;

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
  // -------------------------------------------//
  //       Second collimator modulator box      //
  // -------------------------------------------//
 
  
  const G4double secondCollimatorModulatorXSize = 10.*mm;
  const G4double secondCollimatorModulatorYSize = 200.*mm;
  const G4double secondCollimatorModulatorZSize = 200.*mm;

  const G4double secondCollimatorModulatorXPosition = -2090.5 *mm;
  const G4double secondCollimatorModulatorYPosition = 0.*mm;
  const G4double secondCollimatorModulatorZPosition = 0.*mm;


  G4Box* solidSecondCollimatorModulatorBox = new G4Box("SecondCollimatorModulatorBox", 
						      secondCollimatorModulatorXSize,
						      secondCollimatorModulatorYSize, 
						      secondCollimatorModulatorZSize);

  G4LogicalVolume* logicSecondCollimatorModulatorBox = new G4LogicalVolume(solidSecondCollimatorModulatorBox, 
									  Al, "SecondCollimatorModulatorBox");

  physiSecondCollimatorModulatorBox = new G4PVPlacement(0, G4ThreeVector(secondCollimatorModulatorXPosition,
									secondCollimatorModulatorYPosition,
									secondCollimatorModulatorZPosition),
						       "SecondCollimatorModulatorBox",
						       logicSecondCollimatorModulatorBox,
						       mother, false, 0);
 

 
   // -------------------------------// 
  //   Hole of the second collimator modulator box //
  // -------------------------------//

  const G4double innerRadiusHoleSecondCollimatorModulatorBox = 0.*mm;
  const G4double outerRadiusHoleSecondCollimatorModulatorBox = 31.*mm;
  const G4double hightHoleSecondCollimatorModulatorBox = 10.*mm;
  const G4double startAngleHoleSecondCollimatorModulatorBox = 0.*deg;
  const G4double spanningAngleHoleSecondCollimatorModulatorBox = 360.*deg;

  G4Tubs* solidHoleSecondCollimatorModulatorBox  = new G4Tubs("HoleSecondCollimatorModulatorBox",
							     innerRadiusHoleSecondCollimatorModulatorBox,
							     outerRadiusHoleSecondCollimatorModulatorBox,
							     hightHoleSecondCollimatorModulatorBox ,
							     startAngleHoleSecondCollimatorModulatorBox,
							     spanningAngleHoleSecondCollimatorModulatorBox);

  G4LogicalVolume* logicHoleSecondCollimatorModulatorBox = new G4LogicalVolume(solidHoleSecondCollimatorModulatorBox,
									      Air, "HoleSecondCollimatorModulatorBox", 0, 0, 0);

  physiHoleSecondCollimatorModulatorBox = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
							   "HoleSecondCollimatorModulatorBox",
							   logicHoleSecondCollimatorModulatorBox,
							   physiSecondCollimatorModulatorBox, false, 0); 

  G4VisAttributes * blue = new G4VisAttributes( G4Colour(0. ,0. ,1.));
  blue -> SetVisibility(true);
  blue -> SetForceSolid(true);

  logicFirstCollimator -> SetVisAttributes(yellow);
  logicFirstCollimatorModulatorBox -> SetVisAttributes(blue);
  logicSecondCollimatorModulatorBox -> SetVisAttributes(blue);

}

void HadrontherapyBeamLine::HadrontherapyBeamMonitoring()
{
  // ----------------------------
  // Firt monitor chamber
  // ----------------------------  

  // Each chamber consist of 9 mm of air in a box
  // that has two layers one of kapton and one
  // of copper

  G4Material* Kapton = material -> GetMat("Kapton") ;
  G4Material* Cu = material -> GetMat("MatCopper");
  G4Material* Air = material -> GetMat("Air") ;

  const G4double monitor1XSize = 4.525022*mm;
  const G4double monitor2XSize = 0.000011*mm;
  const G4double monitor3XSize = 4.5*mm;
  const G4double monitorYSize = 10.*cm; 
  const G4double monitorZSize = 10.*cm;
    
  const G4double monitor1XPosition = -1765.97498 *mm;
  const G4double monitor2XPosition = -4.500011*mm;
  const G4double monitor4XPosition = 4.500011*mm;

  G4Box* solidFirstMonitorLayer1 = new G4Box("FirstMonitorLayer1", monitor1XSize, monitorYSize, monitorZSize);

  G4LogicalVolume* logicFirstMonitorLayer1 = new G4LogicalVolume(solidFirstMonitorLayer1, Kapton, "FirstMonitorLayer1");

  physiFirstMonitorLayer1 = new G4PVPlacement(0,G4ThreeVector(monitor1XPosition,0.*cm,0.*cm),
					      "FirstMonitorLayer1", logicFirstMonitorLayer1, mother, false, 0);



  G4Box* solidFirstMonitorLayer2 = new G4Box("FirstMonitorLayer2", monitor2XSize, monitorYSize, monitorZSize);
 
  G4LogicalVolume* logicFirstMonitorLayer2 = new G4LogicalVolume(solidFirstMonitorLayer2, Cu, "FirstMonitorLayer2");
 
  physiFirstMonitorLayer2 = new G4PVPlacement(0, G4ThreeVector(monitor2XPosition,0.*cm,0.*cm),
					      "FirstMonitorLayer2", logicFirstMonitorLayer2, physiFirstMonitorLayer1,
					      false, 0);
  


  G4Box* solidFirstMonitorLayer3 = new G4Box("FirstMonitorLayer3", monitor3XSize, monitorYSize, monitorZSize);

  G4LogicalVolume* logicFirstMonitorLayer3 = new G4LogicalVolume(solidFirstMonitorLayer3, Air, "FirstMonitorLayer3");

  physiFirstMonitorLayer3 = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*cm,0.*cm), "MonitorLayer3",
					      logicFirstMonitorLayer3, physiFirstMonitorLayer1, false, 0);


    
  G4Box* solidFirstMonitorLayer4 = new G4Box("FirstMonitorLayer4", monitor2XSize, monitorYSize, monitorZSize);

  G4LogicalVolume* logicFirstMonitorLayer4 = new G4LogicalVolume(solidFirstMonitorLayer4, Cu, "FirstMonitorLayer4");

  physiFirstMonitorLayer4 = new G4PVPlacement(0, G4ThreeVector(monitor4XPosition,0.*cm,0.*cm), 
					      "FirstMonitorLayer4", logicFirstMonitorLayer4, physiFirstMonitorLayer1, false, 0);



  // ------------------------//
  // Second monitor chamber  //
  // ------------------------//
  

 G4Box* solidSecondMonitorLayer1 = new G4Box("SecondMonitorLayer1", monitor1XSize, monitorYSize, monitorZSize);

 G4LogicalVolume* logicSecondMonitorLayer1 = new G4LogicalVolume(solidSecondMonitorLayer1, Kapton, "SecondMonitorLayer1");

 physiSecondMonitorLayer1 = new G4PVPlacement(0, G4ThreeVector(-1634.92493 *mm,0.*cm,0.*cm),
					       "SecondMonitorLayer1", logicSecondMonitorLayer1,mother, false, 0);



 G4Box* solidSecondMonitorLayer2 = new G4Box("SecondMonitorLayer2", monitor2XSize, monitorYSize, monitorZSize);
 
 G4LogicalVolume* logicSecondMonitorLayer2 = new G4LogicalVolume(solidSecondMonitorLayer2, Cu, "SecondMonitorLayer2");
 
 physiSecondMonitorLayer2 = new G4PVPlacement(0, G4ThreeVector( monitor2XPosition,0.*cm,0.*cm), "SecondMonitorLayer2",
				       logicSecondMonitorLayer2, physiSecondMonitorLayer1, false, 0);
  


 G4Box* solidSecondMonitorLayer3 = new G4Box("SecondMonitorLayer3", monitor3XSize, monitorYSize, monitorZSize);

 G4LogicalVolume* logicSecondMonitorLayer3 = new G4LogicalVolume(solidSecondMonitorLayer3, Air, "SecondMonitorLayer3");

 physiSecondMonitorLayer3 = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*cm,0.*cm), "MonitorLayer3",
					       logicSecondMonitorLayer3, physiSecondMonitorLayer1, false, 0);

  

 G4Box* solidSecondMonitorLayer4 = new G4Box("SecondMonitorLayer4", monitor2XSize, monitorYSize, monitorZSize);

 G4LogicalVolume* logicSecondMonitorLayer4 = new G4LogicalVolume(solidSecondMonitorLayer4, Cu, "SecondMonitorLayer4");

 physiSecondMonitorLayer4 = new G4PVPlacement(0, G4ThreeVector(monitor4XPosition,0.*cm,0.*cm), "SecondMonitorLayer4",
					       logicSecondMonitorLayer4, physiSecondMonitorLayer1, false, 0);

  // -----------------------// 
  // Third monitor chamber  //
  // -----------------------//


 G4Box* solidThirdMonitorLayer1 = new G4Box("ThirdMonitorLayer1", monitor1XSize, monitorYSize, monitorZSize);

 G4LogicalVolume* logicThirdMonitorLayer1 = new G4LogicalVolume(solidThirdMonitorLayer1, Kapton, "ThirdMonitorLayer1");

 physiThirdMonitorLayer1 = new G4PVPlacement(0, G4ThreeVector(-1505.87489 *mm,0.*cm,0.*cm),
					      "ThirdMonitorLayer1", logicThirdMonitorLayer1, mother, false, 0);



 G4Box* solidThirdMonitorLayer2 = new G4Box("ThirdMonitorLayer2", monitor2XSize, monitorYSize, monitorZSize);
 
 G4LogicalVolume* logicThirdMonitorLayer2 = new G4LogicalVolume(solidThirdMonitorLayer2, Cu, "ThirdMonitorLayer2");
 
 physiThirdMonitorLayer2 = new G4PVPlacement(0, G4ThreeVector(monitor2XPosition, 0.*cm,0.*cm), 
						"ThirdMonitorLayer2",
						logicThirdMonitorLayer2, 
						physiThirdMonitorLayer1, 
						false, 0);
 


 G4Box* solidThirdMonitorLayer3 = new G4Box("ThirdMonitorLayer3", monitor3XSize, monitorYSize, monitorZSize);

 G4LogicalVolume* logicThirdMonitorLayer3 = new G4LogicalVolume(solidThirdMonitorLayer3, Air, "ThirdMonitorLayer3");

 physiThirdMonitorLayer3 = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*cm,0.*cm), "MonitorLayer3",
					      logicThirdMonitorLayer3, physiThirdMonitorLayer1, false, 0);
  





 G4Box* solidThirdMonitorLayer4 = new G4Box("ThirdMonitorLayer4", monitor2XSize, monitorYSize, monitorZSize);

 G4LogicalVolume* logicThirdMonitorLayer4 = new G4LogicalVolume(solidThirdMonitorLayer4, Cu, "ThirdMonitorLayer4");
 
 physiThirdMonitorLayer4 = new G4PVPlacement(0, G4ThreeVector(monitor4XPosition,0.*cm,0.*cm), "ThirdMonitorLayer4",
					      logicThirdMonitorLayer4, physiThirdMonitorLayer1, false, 0);
}



void HadrontherapyBeamLine::HadrontherapyBeamNozzle()
{
  // ---------------//
  // Nozzle support //
  //----------------//
 
  const G4double nozzleSupportXSize = 29.50 *mm;
  const G4double nozzleSupportYSize = 180. *mm;
  const G4double nozzleSupportZSize = 180. *mm;

  const G4double nozzleSupportXPosition = -601.00 *mm;

  G4Material* PMMA = material -> GetMat("PMMA");
  G4Material* Brass = material -> GetMat("Brass") ;
  //  G4Material* Air = material -> GetMat("Air") ;

  G4double phi = 90. *deg;     

  // Matrix definition for a 90 deg rotation. Also used for other volumes       
  G4RotationMatrix rm;               
  rm.rotateY(phi);

  // G4Subtraction  
 G4Box* solidNozzleBox = new G4Box("NozzleBox", nozzleSupportXSize, nozzleSupportYSize, nozzleSupportZSize);


  const G4double innerRadiusHoleNozzle = 0.*mm;
  const G4double outerRadiusHoleNozzle = 21.5 *mm;
  const G4double hightHoleNozzle = 29.5 *mm;
  const G4double startAngleHoleNozzle = 0.*deg;
  const G4double spanningAngleHoleNozzle = 360.*deg;

  G4Tubs* solidHoleNozzle  = new G4Tubs("HoleNozzle",
					innerRadiusHoleNozzle,
					outerRadiusHoleNozzle,
					hightHoleNozzle,
					startAngleHoleNozzle,
				        spanningAngleHoleNozzle);

  G4SubtractionSolid* solidNozzleSupport = new G4SubtractionSolid("NozzleSupport",solidNozzleBox, solidHoleNozzle,
							&rm, G4ThreeVector(0 ,0,0 ));

  G4LogicalVolume* logicNozzleSupport = new G4LogicalVolume(solidNozzleSupport, PMMA, "NozzleSupport");
 
  physiNozzleSupport = new G4PVPlacement(0, G4ThreeVector(nozzleSupportXPosition,0., 0.),
					 "NozzleSupport", logicNozzleSupport, mother, false, 0);


  G4VisAttributes * blue = new G4VisAttributes( G4Colour(0. ,0. ,1.));
  blue -> SetVisibility(true);
  blue -> SetForceSolid(false);

  logicNozzleSupport -> SetVisAttributes(blue);
  
  //logicHoleNozzle -> SetVisAttributes(blue);
  
  // ---------------------------------//
  // First hole of the noozle support //
  // ---------------------------------//
  
  const G4double innerRadiusHoleNozzleSupport = 18.*mm;
  const G4double outerRadiusHoleNozzleSupport = 21.5 *mm;
  const G4double hightHoleNozzleSupport = 185.*mm;
  const G4double startAngleHoleNozzleSupport = 0.*deg;
  const G4double spanningAngleHoleNozzleSupport = 360.*deg;

  const G4double holeNozzleSupportXPosition = -475.5 *mm;
 
  G4Tubs* solidHoleNozzleSupport = new G4Tubs("HoleNozzleSupport", 
					      innerRadiusHoleNozzleSupport, 
					      outerRadiusHoleNozzleSupport,
					      hightHoleNozzleSupport, 
					      startAngleHoleNozzleSupport, 
					      spanningAngleHoleNozzleSupport);

  G4LogicalVolume* logicHoleNozzleSupport = new G4LogicalVolume(solidHoleNozzleSupport, 
								Brass, 
								"HoleNozzleSupport", 
								0, 0, 0);

  physiHoleNozzleSupport = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(holeNozzleSupportXPosition,
									      0., 0.)),
					     "HoleNozzleSupport", logicHoleNozzleSupport, mother, false, 0); 

    G4VisAttributes * yellow = new G4VisAttributes( G4Colour(1., 1., 0. ));
  yellow-> SetVisibility(true);
  yellow-> SetForceSolid(true);
  logicHoleNozzleSupport -> SetVisAttributes(yellow);

}

void HadrontherapyBeamLine::HadrontherapyBeamFinalCollimator()
{

  // -----------------------//
  //     Final collimator   //
  //------------------------//
  
  const G4double outerRadiusFinalCollimator = 21.5*mm;
  const G4double hightFinalCollimator = 3.5*mm;
  const G4double startAngleFinalCollimator = 0.*deg;
  const G4double spanningAngleFinalCollimator = 360.*deg;
  const G4double finalCollimatorXPosition = -287.0 *mm;  
  
  G4double phi = 90. *deg;     

  // Matrix definition for a 90 deg rotation. Also used for other volumes       
  G4RotationMatrix rm;               
  rm.rotateY(phi);

  G4Material* Brass = material -> GetMat("Brass");

  solidFinalCollimator = new G4Tubs("FinalCollimator", innerRadiusFinalCollimator, 
				    outerRadiusFinalCollimator,
				    hightFinalCollimator, 
				    startAngleFinalCollimator, 
				    spanningAngleFinalCollimator);

  G4LogicalVolume* logicFinalCollimator = new G4LogicalVolume(solidFinalCollimator, 
							      Brass, "FinalCollimator", 0, 0, 0);

  physiFinalCollimator = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(finalCollimatorXPosition,0.,0.)),
					   "FinalCollimator", logicFinalCollimator, mother, false, 0); 

  G4VisAttributes * yellow = new G4VisAttributes( G4Colour(1., 1., 0. ));
  yellow-> SetVisibility(true);
  yellow-> SetForceSolid(true);
  logicFinalCollimator -> SetVisAttributes(yellow); 

  
}

void HadrontherapyBeamLine::SetRangeShifterXPosition(G4double value)
{
  physiRangeShifterBox -> SetTranslation(G4ThreeVector(value, 0., 0.)); 
  G4cout << "The Range Shifter is translated to"<< value/mm <<"mm along the X axis" <<G4endl;
}

void HadrontherapyBeamLine::SetRangeShifterXSize(G4double value)
{
  solidRangeShifterBox -> SetXHalfLength(value) ;
  G4cout << "RangeShifter size X (mm): "<< ((solidRangeShifterBox -> GetXHalfLength())*2.)/mm
         << G4endl;
}

void HadrontherapyBeamLine::SetFirstScatteringFoilXSize(G4double value)
{ 
  firstScatteringFoil -> SetXHalfLength(value);
  G4cout <<"The X size of the first scattering foil is (mm):"<<  
    ((firstScatteringFoil -> GetXHalfLength())*2.)/mm 
	 << G4endl;
}

void HadrontherapyBeamLine::SetSecondScatteringFoilXSize(G4double value)
{
  secondScatteringFoil -> SetXHalfLength(value);
  G4cout <<"The X size of the second scattering foil is (mm):"<<  
  ((secondScatteringFoil -> GetXHalfLength())*2.)/mm 
	 << G4endl;
}

void HadrontherapyBeamLine::SetOuterRadiusStopper(G4double value)
{ 
  solidStopper -> SetOuterRadius(value);
  G4cout << "OuterRadius od the Stopper is (mm):"
	 << solidStopper -> GetOuterRadius()/mm
	 << G4endl;
}

void HadrontherapyBeamLine::SetInnerRadiusFinalCollimator(G4double value)
{
  solidFinalCollimator -> SetInnerRadius(value);
  G4cout<<"Inner Radius of the final collimator is (mm):"
	<< solidFinalCollimator -> GetInnerRadius()/mm
	<< G4endl; 
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
