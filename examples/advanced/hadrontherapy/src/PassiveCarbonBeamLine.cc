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
// **************************************************************************************
//
// HADRONTHERAPY:  a Geant4-based application for proton/ion-therapy studies
// _________________________________________________________________________


#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "G4NistElementBuilder.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyModulator.hh"
#include "PassiveCarbonBeamLine.hh"
#include "G4SystemOfUnits.hh"
//#include "FaradayCup.hh"

//G4bool PassiveCarbonBeamLine::doCalculation = false;
/////////////////////////////////////////////////////////////////////////////
PassiveCarbonBeamLine::PassiveCarbonBeamLine():
physicalTreatmentRoom(0), hadrontherapyDetectorConstruction(0),
physiBeamLineSupport(0), physiBeamLineCover(0), physiBeamLineCover2(0), 
physiKaptonWindow(0), 
physiFirstMonitorLayer1(0), physiFirstMonitorLayer2(0),
physiFirstMonitorLayer3(0), physiFirstMonitorLayer4(0),
physiNozzleSupport(0), physiHoleNozzleSupport(0),
physiSecondHoleNozzleSupport(0),
solidFinalCollimator(0),
physiFinalCollimator(0)
{

//***************************** PW ***************************************

  static G4String ROGeometryName = "DetectorROGeometry";
  RO = new HadrontherapyDetectorROGeometry(ROGeometryName);
  
  

  G4cout << "Going to register Parallel world...";
  RegisterParallelWorld(RO);
  G4cout << "... done" << G4endl;
//***************************** PW ***************************************
}

/////////////////////////////////////////////////////////////////////////////
PassiveCarbonBeamLine::~PassiveCarbonBeamLine()
{
	delete hadrontherapyDetectorConstruction;
}

/////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* PassiveCarbonBeamLine::Construct()
{ 
	// Sets default geometry and materials
	SetDefaultDimensions();
	
	// Construct the whole CarbonPassive Beam Line 
	ConstructPassiveCarbonBeamLine();
	
	
//***************************** PW ***************************************
  if (!hadrontherapyDetectorConstruction)

//***************************** PW ***************************************
	
	// HadrontherapyDetectorConstruction builds ONLY the phantom and the detector with its associated ROGeometry
	hadrontherapyDetectorConstruction = new HadrontherapyDetectorConstruction(physicalTreatmentRoom); 
	
//***************************** PW ***************************************

 hadrontherapyDetectorConstruction->InitializeDetectorROGeometry(RO,hadrontherapyDetectorConstruction->GetDetectorToWorldPosition());

//***************************** PW ***************************************		
	return physicalTreatmentRoom;
}

// In the following method the DEFAULTS used in the geometry of 
// passive beam line are provided
// HERE THE USER CAN CHANGE THE GEOMETRY CHARACTERISTICS OF BEAM
// LINE ELEMENTS, ALTERNATIVELY HE/SHE CAN USE THE MACRO FILE (IF A 
// MESSENGER IS PROVIDED)
//
// DEFAULT MATERIAL ARE ALSO PROVIDED	
// and COLOURS ARE ALSO DEFINED
// ----------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::SetDefaultDimensions()
{
	// Set of coulors that can be used
	white = new G4VisAttributes( G4Colour());
	white -> SetVisibility(true);
	white -> SetForceSolid(true);
	
	blue = new G4VisAttributes(G4Colour(0. ,0. ,1.));
	blue -> SetVisibility(true);
	blue -> SetForceSolid(true);
	
	gray = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5 ));
	gray-> SetVisibility(true);
	gray-> SetForceSolid(true);
	
	red = new G4VisAttributes(G4Colour(1. ,0. ,0.));
	red-> SetVisibility(true);
	red-> SetForceSolid(true);
	
	yellow = new G4VisAttributes(G4Colour(1., 1., 0. ));
	yellow-> SetVisibility(true);
	yellow-> SetForceSolid(true);
	
	green = new G4VisAttributes( G4Colour(25/255. , 255/255. ,  25/255. ));
	green -> SetVisibility(true);
	green -> SetForceSolid(true);
	
	darkGreen = new G4VisAttributes( G4Colour(0/255. , 100/255. ,  0/255. ));
	darkGreen -> SetVisibility(true);
	darkGreen -> SetForceSolid(true);
	
	darkOrange3 = new G4VisAttributes( G4Colour(205/255. , 102/255. ,  000/255. ));
	darkOrange3 -> SetVisibility(true);
	darkOrange3 -> SetForceSolid(true);
	
	skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
	skyBlue -> SetVisibility(true);
	skyBlue -> SetForceSolid(true);
	
	
	// VACUUM PIPE: first track of the beam line is inside vacuum;
	// The PIPE contains  KAPTON WINDOW
	G4double defaultVacuumZoneXSize = 80.5325 *mm;
	vacuumZoneXSize = defaultVacuumZoneXSize;
	
	G4double defaultVacuumZoneYSize = 52.5 *mm;
	vacuumZoneYSize = defaultVacuumZoneYSize;
	
	G4double defaultVacuumZoneZSize = 52.5 *mm;
	vacuumZoneZSize = defaultVacuumZoneZSize;
	
	// XXX -1775 mm (xKapton to WORLD) - 80.5075 (xKapton to vacuumZone) 
	G4double defaultVacuumZoneXPosition = -1855.5075 *mm;
	vacuumZoneXPosition = defaultVacuumZoneXPosition;
	
	
	// KAPTON WINDOW: it permits the passage of the beam from vacuum to air
	G4double defaultKaptonWindowXSize = 0.025*mm;
	kaptonWindowXSize = defaultKaptonWindowXSize;
	
	G4double defaultKaptonWindowYSize = 5.25*cm;
	kaptonWindowYSize = defaultKaptonWindowYSize;
	
	G4double defaultKaptonWindowZSize = 5.25*cm;
	kaptonWindowZSize = defaultKaptonWindowZSize;
	
	G4double defaultKaptonWindowXPosition = 80.5075*mm;
	kaptonWindowXPosition = defaultKaptonWindowXPosition;
	
	// FIRST SCATTERING FOIL: a thin foil performing a first scattering
	// of the original beam
	G4double defaultFirstScatteringFoilXSize = 0.025 *mm;
	firstScatteringFoilXSize = defaultFirstScatteringFoilXSize; 
	
	G4double defaultFirstScatteringFoilYSize = 105.0   *mm;
	firstScatteringFoilYSize = defaultFirstScatteringFoilYSize;
	
	G4double defaultFirstScatteringFoilZSize = 105   *mm;
	firstScatteringFoilZSize = defaultFirstScatteringFoilZSize;
	
	G4double defaultFirstScatteringFoilXPosition = 0.0 *mm;
	firstScatteringFoilXPosition = defaultFirstScatteringFoilXPosition;
	
	
	
	// STOPPER AND SCATTERING FOIL SIMULATED TO TEST THEIR EFFECT
	// IN THE LATERAL DOSE DISTRIBUTION
	// STOPPER: is a small cylinder able to stop the central component
	// of the beam (having a gaussian shape). It is connected to the SECON SCATTERING FOIL
	// and represent the second element of the scattering system
	G4double defaultInnerRadiusStopper = 0.*cm;
	innerRadiusStopper = defaultInnerRadiusStopper;
	
	G4double defaultHeightStopper = 7.0 *mm;
	heightStopper = defaultHeightStopper;
	
	G4double defaultStartAngleStopper = 0.*deg;
	startAngleStopper = defaultStartAngleStopper;
	
	G4double defaultSpanningAngleStopper = 360.*deg;
	spanningAngleStopper = defaultSpanningAngleStopper;
	
	G4double defaultStopperXPosition = -1675.0 *mm;
	stopperXPosition = defaultStopperXPosition;
	
	G4double defaultStopperYPosition = 0.*m;
	stopperYPosition = defaultStopperYPosition;
	
	G4double defaultStopperZPosition = 0.*m;
	stopperZPosition = defaultStopperZPosition;
	
	G4double defaultOuterRadiusStopper = 2 *mm;
	outerRadiusStopper = defaultOuterRadiusStopper; 
	
	// SECOND SCATTERING FOIL: it is another thin foil and provides the
	// final diffusion of the beam. It represents the third element of the scattering
	// system;
	G4double defaultSecondScatteringFoilXSize = 0.025 *mm;  
	secondScatteringFoilXSize = defaultSecondScatteringFoilXSize;
	
	G4double defaultSecondScatteringFoilYSize = 105.0   *mm;
	secondScatteringFoilYSize = defaultSecondScatteringFoilYSize;
	
	G4double defaultSecondScatteringFoilZSize = 105.0   *mm;
	secondScatteringFoilZSize = defaultSecondScatteringFoilZSize;
	
	G4double defaultSecondScatteringFoilXPosition = defaultStopperXPosition + defaultHeightStopper + defaultSecondScatteringFoilXSize/2;
	secondScatteringFoilXPosition = defaultSecondScatteringFoilXPosition;
	
	G4double defaultSecondScatteringFoilYPosition =  0 *mm;
	secondScatteringFoilYPosition = defaultSecondScatteringFoilYPosition;
	
	G4double defaultSecondScatteringFoilZPosition =  0 *mm;
	secondScatteringFoilZPosition = defaultSecondScatteringFoilZPosition;
	
	
	// FINAL COLLIMATOR: is the collimator giving the final transversal shape
	// of the beam 
	G4double defaultinnerRadiusFinalCollimator = 12.5 *mm;
	innerRadiusFinalCollimator = defaultinnerRadiusFinalCollimator;
	
	// DEFAULT DEFINITION OF THE MATERIALS
	// All elements and compound definition follows the NIST database
	
	// ELEMENTS
	G4bool isotopes = false;
	G4Material* aluminumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al", isotopes);
	G4Material* copperNistAsMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);
	G4Element* zincNist = G4NistManager::Instance()->FindOrBuildElement("Zn");
	G4Element* copperNist = G4NistManager::Instance()->FindOrBuildElement("Cu");
	G4Material* tantalumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ta", isotopes); 
	
	// COMPOUND
	G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
	G4Material* kaptonNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON", isotopes);
	G4Material* galacticNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic", isotopes);
	G4Material* PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);
	
	G4double d; // Density
	G4int nComponents;// Number of components 
	G4double fractionmass; // Fraction in mass of an element in a material
	
	d = 8.40*g/cm3;
	nComponents = 2;
	G4Material* brass = new G4Material("Brass", d, nComponents);  
	brass -> AddElement(zincNist, fractionmass = 30 *perCent);
	brass -> AddElement(copperNist, fractionmass = 70 *perCent);
	
//***************************** PW ***************************************

// DetectorROGeometry Material
  new G4Material("dummyMat", 1., 1.*g/mole, 1.*g/cm3);

//***************************** PW ***************************************
	
	
	// MATERIAL ASSIGNMENT
	// Support of the beam line
	beamLineSupportMaterial = aluminumNist;
	
	// Vacuum pipe
	vacuumZoneMaterial = galacticNist;
	
	// Material of the firt scattering foil
	firstScatteringFoilMaterial = tantalumNist;
	
	// Material of kapton window 
	kaptonWindowMaterial = kaptonNist;
	
	// Material of the stopper
	stopperMaterial = brass;
	
	// Material of the second scattering foil 
	secondScatteringFoilMaterial = tantalumNist; 
	
	// Materials of the monitor chamber
	layer1MonitorChamberMaterial = kaptonNist;
	layer2MonitorChamberMaterial = copperNistAsMaterial;
	layer3MonitorChamberMaterial = airNist;
	layer4MonitorChamberMaterial = copperNistAsMaterial;
	
	
	// material of the final nozzle
	nozzleSupportMaterial = PMMANist;
	holeNozzleSupportMaterial = brass;
	seconHoleNozzleSupportMaterial = airNist;
	
	// Material of the final collimator
	finalCollimatorMaterial = brass;
}

/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::ConstructPassiveCarbonBeamLine()
{ 
	// -----------------------------
	// Treatment room - World volume
	//------------------------------
	// Treatment room sizes
	const G4double worldX = 400.0 *cm;
	const G4double worldY = 400.0 *cm;
	const G4double worldZ = 400.0 *cm;
	G4bool isotopes = false;
	
	G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
	G4Box* treatmentRoom = new G4Box("TreatmentRoom",worldX,worldY,worldZ);
	G4LogicalVolume* logicTreatmentRoom = new G4LogicalVolume(treatmentRoom, 
															  airNist, 
															  "logicTreatmentRoom", 
															  0,0,0);
	physicalTreatmentRoom = new G4PVPlacement(0,
											  G4ThreeVector(),
											  "physicalTreatmentRoom", 
											  logicTreatmentRoom, 
											  0,false,0);
	
	
	// The treatment room is invisible in the Visualisation
	logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::Invisible);
	
	// Components of the Passive Carbon Beam Line
	HadrontherapyBeamLineSupport();
	VacuumToAirInterface();
	ScatteringSystem();
	HadrontherapyBeamMonitoring();
	HadrontherapyBeamNozzle();
	HadrontherapyBeamFinalCollimator();
}

/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::HadrontherapyBeamLineSupport()
{
	// ------------------//
	// BEAM LINE SUPPORT //
	//-------------------//
	const G4double beamLineSupportXSize = 1.5*m;
	const G4double beamLineSupportYSize = 20.*mm;
	const G4double beamLineSupportZSize = 600.*mm;
	
	const G4double beamLineSupportXPosition = -1745.09 *mm;
	const G4double beamLineSupportYPosition = -230. *mm; 
	const G4double beamLineSupportZPosition = 0.*mm;
	
	G4Box* beamLineSupport = new G4Box("BeamLineSupport", 
									   beamLineSupportXSize, 
									   beamLineSupportYSize, 
									   beamLineSupportZSize);
	
	G4LogicalVolume* logicBeamLineSupport = new G4LogicalVolume(beamLineSupport, 
																beamLineSupportMaterial, 
																"BeamLineSupport");
	physiBeamLineSupport = new G4PVPlacement(0, G4ThreeVector(beamLineSupportXPosition, 
															  beamLineSupportYPosition,
															  beamLineSupportZPosition),
											 "BeamLineSupport", 
											 logicBeamLineSupport, 
											 physicalTreatmentRoom, false, 0);
	
	// Visualisation attributes of the beam line support
	
	logicBeamLineSupport -> SetVisAttributes(gray);
	
	//---------------------------------//
	//  Beam line cover 1 (left panel) //
	//---------------------------------//
	const G4double beamLineCoverXSize = 1.5*m;
	const G4double beamLineCoverYSize = 750.*mm;
	const G4double beamLineCoverZSize = 10.*mm; 
	
	const G4double beamLineCoverXPosition = -1745.09 *mm;
	const G4double beamLineCoverYPosition = -1000.*mm; 
	const G4double beamLineCoverZPosition = 610.*mm;
	
	G4Box* beamLineCover = new G4Box("BeamLineCover",
									 beamLineCoverXSize, 
									 beamLineCoverYSize, 
									 beamLineCoverZSize);
	
	G4LogicalVolume* logicBeamLineCover = new G4LogicalVolume(beamLineCover,
															  beamLineSupportMaterial, 
															  "BeamLineCover");
	
	physiBeamLineCover = new G4PVPlacement(0, G4ThreeVector(beamLineCoverXPosition,
															beamLineCoverYPosition,
															beamLineCoverZPosition),
										   "BeamLineCover", 
										   logicBeamLineCover, 
										   physicalTreatmentRoom,
										   false, 
										   0);
	
	// ---------------------------------//
	//  Beam line cover 2 (rigth panel) //
	// ---------------------------------//
	// It has the same characteristic of beam line cover 1 but set in a different position
	physiBeamLineCover2 = new G4PVPlacement(0, G4ThreeVector(beamLineCoverXPosition,
															 beamLineCoverYPosition,
															 - beamLineCoverZPosition),
											"BeamLineCover2", 
											logicBeamLineCover, 
											physicalTreatmentRoom, 
											false, 
											0);
	
	
	logicBeamLineCover -> SetVisAttributes(blue);
}

/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::VacuumToAirInterface()
{
	// ------------//
	// VACUUM PIPE //
	//-------------//
	//
	// First track of the beam line is inside vacuum;
	// The PIPE contains the FIRST SCATTERING FOIL and the KAPTON WINDOW
	G4Box* vacuumZone = new G4Box("VacuumZone", 
								  vacuumZoneXSize, 
								  vacuumZoneYSize, 
								  vacuumZoneZSize);
	
	G4LogicalVolume* logicVacuumZone = new G4LogicalVolume(vacuumZone, 
														   vacuumZoneMaterial, 
														   "VacuumZone");
	
	G4VPhysicalVolume* physiVacuumZone = new G4PVPlacement(0, 
														   G4ThreeVector(vacuumZoneXPosition, 0., 0.),
														   "VacuumZone",
														   logicVacuumZone, 
														   physicalTreatmentRoom, 
														   false, 
														   0);
	
	
	
	
	
	// --------------------------//
	// THE FIRST SCATTERING FOIL //
	// --------------------------//
	// A thin foil performing a first scattering
	// of the original beam
	
	firstScatteringFoil = new G4Box("FirstScatteringFoil", 
									firstScatteringFoilXSize/2, 
									firstScatteringFoilYSize/2, 
									firstScatteringFoilZSize/2);
	
	G4LogicalVolume* logicFirstScatteringFoil = new G4LogicalVolume(firstScatteringFoil,
																	firstScatteringFoilMaterial,
																	"FirstScatteringFoil");
	
	physiFirstScatteringFoil = new G4PVPlacement(0, 
												 G4ThreeVector(firstScatteringFoilXPosition, 0.,0.),
												 "FirstScatteringFoil", 
												 logicFirstScatteringFoil, 
												 physiVacuumZone, 
												 false, 0); 
	
	logicFirstScatteringFoil -> SetVisAttributes(skyBlue);
	
	
	
	// -------------------//
	// THE KAPTON WINDOWS //
	//--------------------//
	//It permits the passage of the beam from vacuum to air
	
	G4Box* solidKaptonWindow = new G4Box("KaptonWindow", 
										 kaptonWindowXSize,
										 kaptonWindowYSize,
										 kaptonWindowZSize);
	
	G4LogicalVolume* logicKaptonWindow = new G4LogicalVolume(solidKaptonWindow, 
															 kaptonWindowMaterial,
															 "KaptonWindow");
	
	physiKaptonWindow = new G4PVPlacement(0, G4ThreeVector(kaptonWindowXPosition, 0., 0.), 
										  "KaptonWindow", logicKaptonWindow,
										  physiVacuumZone, false,	0); 
	
	logicKaptonWindow -> SetVisAttributes(darkOrange3);
}

/////////////////////////////////////////////////////////////////////////////  
void PassiveCarbonBeamLine::ScatteringSystem()
{
	// ------------//
	// THE STOPPER //
	//-------------//
	// Is a small cylinder able to stop the central component
	// of the beam (having a gaussian shape). It is connected to the SECON SCATTERING FOIL
	// and represent the second element of the scattering system
	
	G4double phi = 90. *deg;
		// Matrix definition for a 90 deg rotation with respect to Y axis
	G4RotationMatrix rm;
	rm.rotateY(phi);
	
	solidStopper = new G4Tubs("Stopper", 
							  innerRadiusStopper, 
							  outerRadiusStopper,
							  heightStopper/2, 
							  startAngleStopper,
							  spanningAngleStopper);
	
	logicStopper = new G4LogicalVolume(solidStopper, 
									   stopperMaterial, 
									   "Stopper", 
									   0, 0, 0);
	
	physiStopper = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(stopperXPosition, 
																	 stopperYPosition, 
																	 stopperZPosition)),
									 "Stopper", 
									 logicStopper, 
									 physicalTreatmentRoom, 
									 false, 
									 0);
	
	logicStopper -> SetVisAttributes(red);
	
	// ---------------------------//
	// THE SECOND SCATTERING FOIL //
	// ---------------------------//
	// It is another thin foil and provides the
	// final diffusion of the beam. It represents the third element of the scattering
	// system;
	
	secondScatteringFoil = new G4Box("SecondScatteringFoil", 
									 secondScatteringFoilXSize/2, 
									 secondScatteringFoilYSize/2,
									 secondScatteringFoilZSize/2);
	
	G4LogicalVolume* logicSecondScatteringFoil = new G4LogicalVolume(secondScatteringFoil, 
																	 secondScatteringFoilMaterial,
																	 "SecondScatteringFoil");
	
	physiSecondScatteringFoil = new G4PVPlacement(0, G4ThreeVector(secondScatteringFoilXPosition,
																   secondScatteringFoilYPosition,
																   secondScatteringFoilZPosition),
												  "SeconScatteringFoil", 
												  logicSecondScatteringFoil, 
												  physicalTreatmentRoom, 
												  false, 
												  0);
	
	logicSecondScatteringFoil -> SetVisAttributes(skyBlue);
	
}

/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::HadrontherapyBeamMonitoring()
{
	// ----------------------------
	//   THE FIRST MONITOR CHAMBER
	// ----------------------------  
	// A monitor chamber is a free-air  ionisation chamber
	// able to measure do carbon fluence during the treatment.
	// Here its responce is not simulated in terms of produced 
	// charge but only the energy losses are taked into account.
	// Each chamber consist of 9 mm of air in a box
	// that has two layers one of kapton and one
	// of copper
	const G4double monitor1XSize = 4.525022*mm;
	const G4double monitor2XSize = 0.000011*mm;
	const G4double monitor3XSize = 4.5*mm;
	const G4double monitorYSize = 10.*cm; 
	const G4double monitorZSize = 10.*cm;
		// XXX (Camera monitor size = 9.050088 mm)
	const G4double monitor1XPosition = -1450.474956 *mm;
	const G4double monitor2XPosition = -4.500011*mm;
	const G4double monitor4XPosition = 4.500011*mm;
	
	G4Box* solidFirstMonitorLayer1 = new G4Box("FirstMonitorLayer1", 
											   monitor1XSize, 
											   monitorYSize, 
											   monitorZSize);
	
	G4LogicalVolume* logicFirstMonitorLayer1 = new G4LogicalVolume(solidFirstMonitorLayer1, 
																   layer1MonitorChamberMaterial,
																   "FirstMonitorLayer1");
	
	physiFirstMonitorLayer1 = new G4PVPlacement(0,
												G4ThreeVector(monitor1XPosition,0.*cm,0.*cm),
												"FirstMonitorLayer1", 
												logicFirstMonitorLayer1, 
												physicalTreatmentRoom, 
												false, 
												0);
	
	G4Box* solidFirstMonitorLayer2 = new G4Box("FirstMonitorLayer2", 
											   monitor2XSize, 
											   monitorYSize, 
											   monitorZSize);
	
	G4LogicalVolume* logicFirstMonitorLayer2 = new G4LogicalVolume(solidFirstMonitorLayer2,
																   layer2MonitorChamberMaterial,
																   "FirstMonitorLayer2");
	
	physiFirstMonitorLayer2 = new G4PVPlacement(0, G4ThreeVector(monitor2XPosition,0.*cm,0.*cm),
												"FirstMonitorLayer2", 
												logicFirstMonitorLayer2, 
												physiFirstMonitorLayer1,
												false, 
												0);
	
	G4Box* solidFirstMonitorLayer3 = new G4Box("FirstMonitorLayer3", 
											   monitor3XSize, 
											   monitorYSize, 
											   monitorZSize);
	
	G4LogicalVolume* logicFirstMonitorLayer3 = new G4LogicalVolume(solidFirstMonitorLayer3, 
																   layer3MonitorChamberMaterial, 
																   "FirstMonitorLayer3");
	
	physiFirstMonitorLayer3 = new G4PVPlacement(0, 
												G4ThreeVector(0.*mm,0.*cm,0.*cm), 
												"MonitorLayer3",
												logicFirstMonitorLayer3, 
												physiFirstMonitorLayer1, 
												false, 
												0);
    
	G4Box* solidFirstMonitorLayer4 = new G4Box("FirstMonitorLayer4", 
											   monitor2XSize, 
											   monitorYSize, 
											   monitorZSize);
	
	G4LogicalVolume* logicFirstMonitorLayer4 = new G4LogicalVolume(solidFirstMonitorLayer4, 
																   layer4MonitorChamberMaterial, 
																   "FirstMonitorLayer4");
	
	physiFirstMonitorLayer4 = new G4PVPlacement(0, G4ThreeVector(monitor4XPosition,0.*cm,0.*cm), 
												"FirstMonitorLayer4",
												logicFirstMonitorLayer4,
												physiFirstMonitorLayer1, false, 0);
	
	logicFirstMonitorLayer3 -> SetVisAttributes(white);
	
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::HadrontherapyBeamNozzle()
{
	// ------------------------------//
	// THE FINAL TUBE AND COLLIMATOR //
	//-------------------------------//
	// The last part of the transport beam line consists of
	// a 59 mm thick PMMA slab (to stop all the diffused radiation), a 285 mm brass tube
	// (to well collimate the carbon beam) and a final collimator with 25 mm diameter
	// aperture (that provide the final trasversal shape of the beam)
	
	// -------------------//
	//     PMMA SUPPORT   //
	// -------------------//

	const G4double nozzleSupportXSize = 29.5 *mm;
	const G4double nozzleSupportYSize = 180. *mm;
	const G4double nozzleSupportZSize = 180. *mm;
	//XXX Placed at 
	const G4double nozzleSupportXPosition = -558. *mm;
	
	G4double phi = 90. *deg;     
	// Matrix definition for a 90 deg rotation. Also used for other volumes       
	G4RotationMatrix rm;               
	rm.rotateY(phi);
	
	G4Box* solidNozzleSupport = new G4Box("NozzleSupport", 
										  nozzleSupportXSize, 
										  nozzleSupportYSize, 
										  nozzleSupportZSize);
	
	G4LogicalVolume* logicNozzleSupport = new G4LogicalVolume(solidNozzleSupport, 
															  nozzleSupportMaterial, 
															  "NozzleSupport");
	
	physiNozzleSupport = new G4PVPlacement(0, G4ThreeVector(nozzleSupportXPosition,0., 0.),
										   "NozzleSupport", 
										   logicNozzleSupport, 
										   physicalTreatmentRoom, 
										   false, 
										   0);
	
	logicNozzleSupport -> SetVisAttributes(yellow);
	// -------------------//
	//     BRASS TUBE     //
	// -------------------//
	const G4double innerRadiusHoleNozzleSupport = 18.*mm;
	const G4double outerRadiusHoleNozzleSupport = 21.5 *mm;
	//XXX h/2 = 142.5 mm
	const G4double hightHoleNozzleSupportFirst = nozzleSupportXSize;
	const G4double hightHoleNozzleSupport = 113.0*mm;
	const G4double startAngleHoleNozzleSupport = 0.*deg;
	const G4double spanningAngleHoleNozzleSupport = 360.*deg;
	const G4double holeNozzleSupportXPosition = -415.5 *mm;
	G4Tubs* solidNozzleSupportHole = new G4Tubs("NozzleSupportHole1", 	innerRadiusHoleNozzleSupport, 
												outerRadiusHoleNozzleSupport,
												hightHoleNozzleSupportFirst, 
												startAngleHoleNozzleSupport, 
												spanningAngleHoleNozzleSupport);
	
	G4LogicalVolume* logicNozzleSupportHole = new G4LogicalVolume(solidNozzleSupportHole, 
															  holeNozzleSupportMaterial, 
															  "NozzleSupportHole1");
	
	physiNozzleSupportHole = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(0, 0., 0.)),
											   "HoleNozzleSupportHole1", 
											   logicNozzleSupportHole, 
											   physiNozzleSupport, false, 0);
		
	G4Tubs* solidHoleNozzleSupport = new G4Tubs("HoleNozzleSupport", 
												innerRadiusHoleNozzleSupport, 
												outerRadiusHoleNozzleSupport,
												hightHoleNozzleSupport, 
												startAngleHoleNozzleSupport, 
												spanningAngleHoleNozzleSupport);
	
	G4LogicalVolume* logicHoleNozzleSupport = new G4LogicalVolume(solidHoleNozzleSupport, 
																  holeNozzleSupportMaterial, 
																  "HoleNozzleSupport", 
																  0, 0, 0);
	
	physiHoleNozzleSupport = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(holeNozzleSupportXPosition, 0., 0.)),
											   "HoleNozzleSupport", 
											   logicHoleNozzleSupport, 
											   physicalTreatmentRoom, false, 0); 
	logicNozzleSupportHole -> SetVisAttributes(darkOrange3);		
	logicHoleNozzleSupport -> SetVisAttributes(darkOrange3);
	
	//--------------------------------------------------------------//
	// HOLE OF THE BRASS TUBE (otherwise we'll have PMMA)           //
	//--------------------------------------------------------------//
	const G4double innerRadiusSecondHoleNozzleSupport = 0.*mm;
	const G4double outerRadiusSecondHoleNozzleSupport = 18.*mm;
	const G4double hightSecondHoleNozzleSupport = 29.5 *mm;
	const G4double startAngleSecondHoleNozzleSupport = 0.*deg;
	const G4double spanningAngleSecondHoleNozzleSupport = 360.*deg;
	
	G4Tubs* solidSecondHoleNozzleSupport = new G4Tubs("SecondHoleNozzleSupport",
													  innerRadiusSecondHoleNozzleSupport,
													  outerRadiusSecondHoleNozzleSupport, 
													  hightSecondHoleNozzleSupport,
													  startAngleSecondHoleNozzleSupport, 
													  spanningAngleSecondHoleNozzleSupport);
	
	G4LogicalVolume* logicSecondHoleNozzleSupport = new G4LogicalVolume(solidSecondHoleNozzleSupport, 
																		seconHoleNozzleSupportMaterial,
																		"SecondHoleNozzleSupport",
																		0, 
																		0,
																		0);
	
	physiSecondHoleNozzleSupport = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()), 
													 "SecondHoleNozzleSupport",
													 logicSecondHoleNozzleSupport, 
													 physiNozzleSupport, 
													 false, 0); 
	
	
	logicHoleNozzleSupport -> SetVisAttributes(darkOrange3); 
}

/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::HadrontherapyBeamFinalCollimator()
{
	// -----------------------//
	//     FINAL COLLIMATOR   //
	//------------------------//
	const G4double outerRadiusFinalCollimator = 21.5*mm;
	const G4double hightFinalCollimator = 3.5*mm;
	const G4double startAngleFinalCollimator = 0.*deg;
	const G4double spanningAngleFinalCollimator = 360.*deg;
	//XXX
	const G4double finalCollimatorXPosition = -299.0 *mm;  
	
	G4double phi = 90. *deg;     
	
	// Matrix definition for a 90 deg rotation. Also used for other volumes       
	G4RotationMatrix rm;
	rm.rotateY(phi);
	
	solidFinalCollimator = new G4Tubs("FinalCollimator", 
									  innerRadiusFinalCollimator, 
									  outerRadiusFinalCollimator,
									  hightFinalCollimator, 
									  startAngleFinalCollimator, 
									  spanningAngleFinalCollimator);
	
	G4LogicalVolume* logicFinalCollimator = new G4LogicalVolume(solidFinalCollimator, 
																finalCollimatorMaterial, 
																"FinalCollimator", 
																0,
																0, 
																0);
	
	physiFinalCollimator = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(finalCollimatorXPosition,0.,0.)),
											 "FinalCollimator", logicFinalCollimator, physicalTreatmentRoom, false, 0); 
	
	logicFinalCollimator -> SetVisAttributes(yellow); 
}

