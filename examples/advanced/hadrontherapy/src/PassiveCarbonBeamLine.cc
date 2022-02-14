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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

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
#include "G4Trd.hh"
#include "PassiveCarbonBeamLineMessenger.hh"
#include "G4UnionSolid.hh"

using namespace std;

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
<<<<<<< HEAD
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
=======
    // Set of coulors that can be used
    white = new G4VisAttributes( G4Colour());
    white -> SetVisibility(true);
    white -> SetForceSolid(true);
    
    black = new G4VisAttributes( G4Colour(1., 1., 1.));
    black -> SetVisibility(true);
    black -> SetForceSolid(true);
    
    
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
    
    
    // FINAL COLLIMATOR: is the collimator giving the final transversal shape
    // of the beam
    
    G4double defaultinnerRadiusFinalCollimator = 12.5 *mm;
    innerRadiusFinalCollimator = defaultinnerRadiusFinalCollimator;
    
    // DEFAULT DEFINITION OF THE MATERIALS
    // All elements and compound definition follows the NIST database
    
    // ELEMENTS
    G4bool isotopes = false;
    aluminumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al", isotopes);
    G4Element* zincNist = G4NistManager::Instance()->FindOrBuildElement("Zn");
    G4Element* copperNist = G4NistManager::Instance()->FindOrBuildElement("Cu");
    
    // MATERIAL (Including compounds)
    copperNistMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);
    airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
    kaptonNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON", isotopes);
    galacticNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic", isotopes);
    PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);
    tantalumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ta", isotopes);
    
    G4double d; // Density
    G4int nComponents;// Number of components
    G4double fractionmass; // Fraction in mass of an element in a material
    
    d = 8.40*g/cm3;
    nComponents = 2;
    brass = new G4Material("Brass", d, nComponents);
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
    firstScatteringFoilMaterial = tantalumNist;
    
    // Material of kapton window
    kaptonWindowMaterial = kaptonNist;
    
    // Material of ripple filter
    rippleFilterMaterial = PMMANist;
    rippleFilterBoxMaterial = airNist;
    
    // Materials of the monitor chamber
    layer1MonitorChamberMaterial = kaptonNist;
    layer2MonitorChamberMaterial = copperNistMaterial;
    layer3MonitorChamberMaterial = airNist;
    layer4MonitorChamberMaterial = copperNistMaterial;
    
    // Material of the final nozzle
    nozzleSupportMaterial = PMMANist;
    holeNozzleSupportMaterial = airNist;
    seconHoleNozzleSupportMaterial = airNist;
    
    // Material of the final collimator
    brassTubeMaterial = brassTube2Material = brassTube3Material = brass;
    
    // Material of the final collimator
    finalCollimatorMaterial = brass;
    
    // Material of the PMMA collimator
    PMMACollimatorMaterial = airNist;
    
    
    G4Element* hydrogenNist = G4NistManager::Instance()->FindOrBuildElement("H");
    G4Element* oxygenNist = G4NistManager::Instance()->FindOrBuildElement("O");
    G4Element* carbonNist = G4NistManager::Instance()->FindOrBuildElement("C");
    
    G4Material* plastic = new G4Material("Plastic", 1.18*g/cm3, nComponents=3);
    plastic -> AddElement(carbonNist, 21);
    plastic -> AddElement(oxygenNist, 4);
    plastic -> AddElement(hydrogenNist, 24);
    PMMANist = plastic;
    
    
    
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
    
    airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
    treatmentRoom = new G4Box("TreatmentRoom",worldX,worldY,worldZ);
    logicTreatmentRoom = new G4LogicalVolume(treatmentRoom,
                                             airNist,
                                             "logicTreatmentRoom",
                                             0,0,0);
    physicalTreatmentRoom = new G4PVPlacement(0,
                                              G4ThreeVector(),
                                              "physicalTreatmentRoom",
                                              logicTreatmentRoom,
                                              0,false,0);
    
    
    // The treatment room is invisible in the Visualisation
    logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::GetInvisible());
    
    // Components of the Passive Carbon Beam Line
    HadrontherapyBeamLineSupport();
    VacuumToAirInterface();
    HadrontherapyBeamMonitoring();
    HadrontherapyBeamNozzle();
    HadrontherapyBeamFinalCollimator();
    HadrontherapyPMMACollimator();
   // HadrontherapyRippleFilter();
    
    
   // StopperCostruction();
    
    
   HadrontherapyRidgeFilter();
}

/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::HadrontherapyBeamLineSupport()
{
    // ------------------//
    // BEAM LINE SUPPORT //
    //-------------------//
    
    beamLineSupportXSize = 1.5*m;
    beamLineSupportYSize = 20.*mm;
    beamLineSupportZSize = 600.*mm;
    
    beamLineSupportXPosition = -1745.09 *mm;
    beamLineSupportYPosition = -230. *mm;
    beamLineSupportZPosition = 0.*mm;
    
    beamLineSupport = new G4Box("BeamLineSupport",
                                beamLineSupportXSize,
                                beamLineSupportYSize,
                                beamLineSupportZSize);
    
    logicBeamLineSupport = new G4LogicalVolume(beamLineSupport,
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
    beamLineCoverXSize = 1.5*m;
    beamLineCoverYSize = 750.*mm;
    beamLineCoverZSize = 10.*mm;
    
    beamLineCoverXPosition = -1745.09 *mm;
    beamLineCoverYPosition = -1000.*mm;
    beamLineCoverZPosition = 610.*mm;
    
    beamLineCover = new G4Box("BeamLineCover",
                              beamLineCoverXSize,
                              beamLineCoverYSize,
                              beamLineCoverZSize);
    
    logicBeamLineCover = new G4LogicalVolume(beamLineCover,
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
<<<<<<< HEAD
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
=======
    // ------------//
    // VACUUM PIPE //
    //-------------//
    //
    // First track of the beam line is inside vacuum;
    
    vacuumZoneXSize = 100 *mm;
    vacuumZoneYSize = 52.5 *mm;
    vacuumZoneZSize = 52.5 *mm;
    vacuumPipeXPosition = -1708.0 *mm;
    
    
    vacuumZone = new G4Box("VacuumZone",
                           vacuumZoneXSize/2,
                           vacuumZoneYSize/2,
                           vacuumZoneZSize/2);
    
    logicVacuumZone = new G4LogicalVolume(vacuumZone, vacuumZoneMaterial,"VacuumZone");
    
    physiVacuumZone = new G4PVPlacement(0, G4ThreeVector(vacuumPipeXPosition, 0., 0.),
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
    
    firstScatteringFoilXSize = 0.015 *mm;
    firstScatteringFoilYSize = 52.5 *mm;
    firstScatteringFoilZSize = 52.5 *mm;
    firstScatteringFoilXPosition = 0.0 *mm;
    
    firstScatteringFoil = new G4Box("FirstScatteringFoil",
                                    firstScatteringFoilXSize/2,
                                    firstScatteringFoilYSize/2,
                                    firstScatteringFoilZSize/2);
    
    logicFirstScatteringFoil = new G4LogicalVolume(firstScatteringFoil,
                                                   firstScatteringFoilMaterial,
                                                   "FirstScatteringFoil");
    
    physiFirstScatteringFoil = new G4PVPlacement(0, G4ThreeVector(firstScatteringFoilXPosition, 0.,0.),
                                                 "FirstScatteringFoil", logicFirstScatteringFoil, physiVacuumZone,
                                                 false, 0);
    
    logicFirstScatteringFoil -> SetVisAttributes(skyBlue);
    
    // -------------------//
    // THE KAPTON WINDOWS //
    //--------------------//
    //It permits the passage of the beam from vacuum to air
    
    // KAPTON WINDOW: it permits the passage of the beam from vacuum to air
    kaptonWindowXSize = 0.050 *mm;
    kaptonWindowYSize = 52.5 *mm;
    kaptonWindowZSize = 52.5 *mm;
    kaptonWindowXPosition = vacuumZoneXSize/2 - kaptonWindowXSize/2;
    
    solidKaptonWindow = new G4Box("KaptonWindow",
                                  kaptonWindowXSize/2,
                                  kaptonWindowYSize/2,
                                  kaptonWindowZSize/2);
    
    logicKaptonWindow = new G4LogicalVolume(solidKaptonWindow,
                                            kaptonWindowMaterial,
                                            "KaptonWindow");
    
    physiKaptonWindow = new G4PVPlacement(0, G4ThreeVector(kaptonWindowXPosition, 0., 0.),
                                          "KaptonWindow", logicKaptonWindow,
                                          physiVacuumZone, false,    0);
    
    logicKaptonWindow -> SetVisAttributes(darkOrange3);
}
/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::HadrontherapyRippleFilter()
{
    
    G4double defaultRippleFilterXPosition = -1638.0*mm;
    G4double ripple_position=(defaultRippleFilterXPosition);
    G4double RF_x = 200.0 * mm;
    G4double RF_y = 200.0 * mm;
    G4double RF_z = 1.4 * mm;
    G4double RFbase_z = 0.2 * mm;
    G4double RFtrd_z = RF_z - RFbase_z;
    G4double RFtrd_top = 1e-4 * mm;
    G4double RFtrd_bottom = 1.5 * mm;
    G4double distanceBetweenTrd = 0.1*mm;
    
    
    
    
    G4double theta = -90. *deg;
    // Matrix definition for a "theta" deg rotation with respect to Y axis
    G4RotationMatrix rot;
    rot.rotateY(theta);
    
    
    SolidRippleFilter= new G4Box("RippleFilter",
                                 RF_x/2 + 1*mm,
                                 RF_y/2 + 1*mm,
                                 RF_z/2 + 1*mm);
    
    LogicRippleFilter = new G4LogicalVolume(SolidRippleFilter,
                                            rippleFilterBoxMaterial,
                                            "LogicRippleFilter",
                                            0,0,0);
    
    PhysiRippleFilter = new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(ripple_position,0,0)),
                                          "PhysiRippleFilter",
                                          LogicRippleFilter,
                                          physicalTreatmentRoom,
                                          false,
                                          1,
                                          true);
    
    PhysiRippleFilter = new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(ripple_position + 10*cm,0,0)),
                                          "PhysiRippleFilter",
                                          LogicRippleFilter,
                                          physicalTreatmentRoom,
                                          false,
                                          2,
                                          true);
    
    LogicRippleFilter -> SetVisAttributes(G4VisAttributes::GetInvisible());
    
    SolidRippleFilterBase = new G4Box("RippleFilterBase",
                                      RF_x/2,
                                      RF_y/2,
                                      RFbase_z/2);
    
    LogicRippleFilterBase = new G4LogicalVolume(SolidRippleFilterBase,
                                                rippleFilterMaterial,
                                                "LogicRippleFilterBase",
                                                0,0,0);
    
    LogicRippleFilterBase -> SetVisAttributes(green);
    
    PhysiRippleFilterBase = new G4PVPlacement(0,
                                              G4ThreeVector(0, 0, -RF_z/2 + RFbase_z/2),
                                              "PhysiRippleFilter",
                                              LogicRippleFilterBase,
                                              PhysiRippleFilter,
                                              false,
                                              0,
                                              false);
    
    SolidRippleFilterTrd = new G4Trd("SolidRippleFilterTrd",
                                     RF_x/2,
                                     RF_x/2,
                                     RFtrd_bottom/2,
                                     RFtrd_top/2,
                                     RFtrd_z/2);
    
    LogicRippleFilterTrd = new G4LogicalVolume(SolidRippleFilterTrd,
                                               rippleFilterMaterial,
                                               "LogicRippleFilterTrd",
                                               0,0,0);
    
    LogicRippleFilterTrd -> SetVisAttributes(green);
    
    G4int numberOfTrd = static_cast<int>(std::floor( RF_y / (RFtrd_bottom+distanceBetweenTrd) ));
    
    G4int N = static_cast<int>( std::floor(numberOfTrd-1)/2 );
    
    G4int copyNumber = 0;
    
    for( int i = -N; i <= N; i++ )
    {
        PhysiRippleFilterTrd = new G4PVPlacement(0,
                                                 G4ThreeVector(0,
                                                               i*(RFtrd_bottom+distanceBetweenTrd),
                                                               -RF_z/2+RFbase_z+RFtrd_z/2),
                                                 "PhysiRippleFilterTrd",
                                                 LogicRippleFilterTrd,
                                                 PhysiRippleFilter,
                                                 false,
                                                 copyNumber,
                                                 false);
        
        copyNumber++;
    }
    
}
/////////////////////////////////////////////////////////////////////////////

void PassiveCarbonBeamLine::StopperCostruction(){
    
    supportFoil= new G4Box("supportFoil",25/2*um,10*cm,10*cm);
    LogicSupportFoil = new G4LogicalVolume(supportFoil,tantalumNist,"logicSupportFoil");
    PhysiSupportFoil = new G4PVPlacement(0,G4ThreeVector(-1398.*mm+25/2*um,0.*mm,0.*mm),"logicSupportFoil",LogicSupportFoil,physicalTreatmentRoom,false,0 );
    LogicSupportFoil -> SetVisAttributes(skyBlue);
    
    
    
    stopper = new G4Tubs("Stopper",0*mm,3*mm,3.5*mm,0.*deg,360.*deg);
    LogicStopper= new G4LogicalVolume(stopper,brass,"logicalStopper");
    
    G4double ti = -270. *deg;
    G4RotationMatrix rt;
    rt.rotateY(ti);
    
    PhysicStopper= new G4PVPlacement(G4Transform3D(rt,G4ThreeVector(-1398.*mm-3.5*mm,0*mm,0*mm)),"logicalStopper",LogicStopper,physicalTreatmentRoom,false,0);
    
    LogicStopper -> SetVisAttributes(red);
    
    
}

void PassiveCarbonBeamLine::HadrontherapyRidgeFilter(){
    
   
    G4double defaultRidgeXPosition= -1270.0*mm; 
    const G4double XBase = 0.5 *mm;
    const G4double YBase =6.03*cm;
    const G4double ZBase =6.03*cm;
    
    
    
    SolidRidgeBase = new G4Box("BaseRidgeFilter",
                               XBase,
                               YBase/2,
                               ZBase/2);
    
    LogicRidgeBase = new G4LogicalVolume(SolidRidgeBase,
                                         PMMANist,
                                         "BaseRidgeFilter");
    
    PhysiRidgeFilterBase = new G4PVPlacement(0,
                                             G4ThreeVector(
                                                           -1270*mm, 
                                                           0.,
                                                           0.),
                                             "BaseRidgeFilter",
                                             LogicRidgeBase,
                                             physicalTreatmentRoom,
                                             false,
                                             0);
    
    LogicRidgeBase->SetVisAttributes(red);
    
    
    SolidRidgeMother = new G4Box("MotherRidgeSOL",1.7*mm/2,1.7*mm/2,4.72*mm/2);
    LogicRidgeMother = new G4LogicalVolume(SolidRidgeMother,airNist,"MotherRidge");
    
    G4Trd* trapp1=new G4Trd("Trapp1SOL",1.7*mm/2,1.68*mm/2,1.7*mm/2,1.68*mm/2,0.22*mm/2);
    G4LogicalVolume* LogicTrapp1=new G4LogicalVolume(trapp1,PMMANist,"Trapp1");
 PhysiTrapp1=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22/2)*mm),LogicTrapp1,"Trapp1",LogicRidgeMother,false,0);
    
    G4Trd* trapp2=new G4Trd("Trapp2SOL",1.68*mm/2,1.64*mm/2,1.68*mm/2,1.64*mm/2,0.08*mm/2);
    G4LogicalVolume* LogicTrapp2=new G4LogicalVolume(trapp2,PMMANist,"Trapp2");
  PhysiTrapp2=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08/2)*mm),LogicTrapp2,"Trapp2",LogicRidgeMother,false,0);
    
    G4Trd* trapp3=new G4Trd("Trapp3SOL",1.64*mm/2,1.58*mm/2,1.64*mm/2,1.58*mm/2,0.14*mm/2);
    G4LogicalVolume* LogicTrapp3=new G4LogicalVolume(trapp3,PMMANist,"Trapp3");
 PhysiTrapp3=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14/2)*mm),LogicTrapp3,"Trapp3",LogicRidgeMother,false,0);
    
    G4Trd* trapp4=new G4Trd("Trapp4SOL",1.58*mm/2,1.50*mm/2,1.58*mm/2,1.50*mm/2,0.22*mm/2);
    G4LogicalVolume* LogicTrapp4=new G4LogicalVolume(trapp4,PMMANist,"Trapp4");
    PhysiTrapp4=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22/2)*mm),LogicTrapp4,"Trapp4",LogicRidgeMother,false,0);
    
    G4Trd* trapp5=new G4Trd("Trapp5SOL",1.50*mm/2,1.40*mm/2,1.50*mm/2,1.40*mm/2,0.31*mm/2);
    G4LogicalVolume* LogicTrapp5=new G4LogicalVolume(trapp5,PMMANist,"Trapp5");
   PhysiTrapp5=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31/2)*mm),LogicTrapp5,"Trapp5",LogicRidgeMother,false,0);
    
    G4Trd* trapp6=new G4Trd("Trapp6SOL",1.40*mm/2,1.26*mm/2,1.40*mm/2,1.26*mm/2,0.48*mm/2);
    G4LogicalVolume* LogicTrapp6=new G4LogicalVolume(trapp6,PMMANist,"Trapp6");
   PhysiTrapp6=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31+0.48/2)*mm),LogicTrapp6,"Trapp6",LogicRidgeMother,false,0);
    
    G4Trd* trapp7=new G4Trd("Trapp7SOL",1.26*mm/2,0.94*mm/2,1.26*mm/2,0.94*mm/2,1.20*mm/2);
    G4LogicalVolume* LogicTrapp7=new G4LogicalVolume(trapp7,PMMANist,"Trapp7");
PhysiTrapp7=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31+0.48+1.20/2)*mm),LogicTrapp7,"Trapp7",LogicRidgeMother,false,0);
    
    G4Trd* trapp8=new G4Trd("Trapp8SOL",0.94*mm/2,0.78*mm/2,0.94*mm/2,0.78*mm/2,0.57*mm/2);
    G4LogicalVolume* LogicTrapp8=new G4LogicalVolume(trapp8,PMMANist,"Trapp8");
    PhysiTrapp8=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31+0.48+1.20+0.57/2)*mm),LogicTrapp8,"Trapp8",LogicRidgeMother,false,0);
    
    G4Trd* trapp9=new G4Trd("Trapp9SOL",0.78*mm/2,0.66*mm/2,0.78*mm/2,0.66*mm/2,0.39*mm/2);
    G4LogicalVolume* LogicTrapp9=new G4LogicalVolume(trapp9,PMMANist,"Trapp9");
  PhysiTrapp9=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31+0.48+1.20+0.57+0.39/2)*mm),LogicTrapp9,"Trapp9",LogicRidgeMother,false,0);
    
    G4Trd* trapp10=new G4Trd("Trapp10SOL",0.66*mm/2,0.56*mm/2,0.66*mm/2,0.56*mm/2,0.29*mm/2);
    G4LogicalVolume* LogicTrapp10=new G4LogicalVolume(trapp10,PMMANist,"Trapp10");
 PhysiTrapp10=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31+0.48+1.20+0.57+0.39+0.29/2)*mm),LogicTrapp10,"Trapp10",LogicRidgeMother,false,0);
    
    G4Trd* trapp11=new G4Trd("Trapp11SOL",0.56*mm/2,0.46*mm/2,0.56*mm/2,0.46*mm/2,0.25*mm/2);
    G4LogicalVolume* LogicTrapp11=new G4LogicalVolume(trapp11,PMMANist,"Trapp11");
PhysiTrapp11=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31+0.48+1.20+0.57+0.39+0.29+0.25/2)*mm),LogicTrapp11,"Trapp11",LogicRidgeMother,false,0);
    
    G4Trd* trapp12=new G4Trd("Trapp12SOL",0.46*mm/2,0.38*mm/2,0.46*mm/2,0.38*mm/2,0.17*mm/2);
    G4LogicalVolume* LogicTrapp12=new G4LogicalVolume(trapp12,PMMANist,"Trapp12");
   PhysiTrapp12=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31+0.48+1.20+0.57+0.39+0.29+0.25+0.17/2)*mm),LogicTrapp12,"Trapp12",LogicRidgeMother,false,0);
    
    G4Trd* trapp13=new G4Trd("Trapp13SOL",0.38*mm/2,0.30*mm/2,0.38*mm/2,0.30*mm/2,0.14*mm/2);
    G4LogicalVolume* LogicTrapp13=new G4LogicalVolume(trapp13,PMMANist,"Trapp13");
   PhysiTrapp13=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31+0.48+1.20+0.57+0.39+0.29+0.25+0.17+0.14/2)*mm),LogicTrapp13,"Trapp13",LogicRidgeMother,false,0);
    
    G4Trd* trapp14=new G4Trd("Trapp14SOL",0.30*mm/2,0.24*mm/2,0.30*mm/2,0.24*mm/2,0.09*mm/2);
    G4LogicalVolume* LogicTrapp14=new G4LogicalVolume(trapp14,PMMANist,"Trapp13");
PhysiTrapp14=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31+0.48+1.20+0.57+0.39+0.29+0.25+0.17+0.14+0.09/2)*mm),LogicTrapp14,"Trapp14",LogicRidgeMother,false,0);
    
    G4Trd* trapp15=new G4Trd("Trapp14SOL",0.24*mm/2,0.18*mm/2,0.24*mm/2,0.18*mm/2,0.07*mm/2);
    G4LogicalVolume* LogicTrapp15=new G4LogicalVolume(trapp15,PMMANist,"Trapp15");
    PhysiTrapp15=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31+0.48+1.20+0.57+0.39+0.29+0.25+0.17+0.14+0.09+0.07/2)*mm),LogicTrapp15,"Trapp15",LogicRidgeMother,false,0);
    
    G4Trd* trapp16=new G4Trd("Trapp16SOL",0.18*mm/2,0.12*mm/2,0.18*mm/2,0.12*mm/2,0.10*mm/2);
    G4LogicalVolume* LogicTrapp16=new G4LogicalVolume(trapp16,PMMANist,"Trapp16");
  PhysiTrapp16=new G4PVPlacement(0,G4ThreeVector(0.,0.,(-4.72/2+0.22+0.08+0.14+0.22+0.31+0.48+1.20+0.57+0.39+0.29+0.25+0.17+0.14+0.09+0.07+0.10/2)*mm),LogicTrapp16,"Trapp16",LogicRidgeMother,false,0);
    
    LogicTrapp1->SetVisAttributes(green);
    LogicTrapp2->SetVisAttributes(green);
    LogicTrapp3->SetVisAttributes(green);
    LogicTrapp4->SetVisAttributes(green);
    LogicTrapp5->SetVisAttributes(green);
    LogicTrapp6->SetVisAttributes(green);
    LogicTrapp7->SetVisAttributes(green);
    LogicTrapp8->SetVisAttributes(green);
    LogicTrapp9->SetVisAttributes(green);
    LogicTrapp10->SetVisAttributes(green);
    LogicTrapp11->SetVisAttributes(green);
    LogicTrapp12->SetVisAttributes(green);
    LogicTrapp13->SetVisAttributes(green);
    LogicTrapp14->SetVisAttributes(green);
    LogicTrapp15->SetVisAttributes(green);
    LogicTrapp16->SetVisAttributes(green);
    G4VisAttributes* visAttr = new G4VisAttributes();
    visAttr->SetVisibility(false);
    LogicRidgeMother->SetVisAttributes(visAttr);
    
    
    
    
    
    
    G4int numberOfLayers = 30;
    G4double minZ = 30.15*mm-0.3*mm-1.70/2*mm;
    G4double minY = 30.15*mm-0.3*mm-1.70/2*mm;
    G4double sum_space = 0.3*mm+1.70*mm;
    
    
    std::vector<G4ThreeVector> singleTrapPositions;
    
    for (int i = 0; i < numberOfLayers; i++)
    {
        for (int j = 0; j < numberOfLayers; j++)
        {
            singleTrapPositions.push_back({defaultRidgeXPosition-4.72/2*mm-0.5*mm,
                minY - i*sum_space,
                minZ - j*sum_space,
                
            });
        }
    }
    
    G4double ti = -  90. *deg;
    G4RotationMatrix rt;
    rt.rotateY(ti);
    
    
    G4int peaks = numberOfLayers*numberOfLayers;
    for (int i = 0; i < peaks; i++)
    {
        
        std::ostringstream tName;tName << "singleTrap" << i;
        new G4PVPlacement(G4Transform3D(rt,
                                        singleTrapPositions[i]),
                          LogicRidgeMother,
                          "tName.str()",
                          logicTreatmentRoom,
                          0,
                          i);
        
    }
    
}





>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

/////////////////////////////////////////////////////////////////////////////  
void PassiveCarbonBeamLine::ScatteringSystem()
{
    
    // ----------------------//
    //   PMMA COLLIMATOR     //
    // ----------------------//
    PMMACollimatorSupportXSize = 25.0 *mm;
    PMMACollimatorSupportYSize = 200. *mm;
    PMMACollimatorSupportZSize = 200. *mm;
    
    
    PMMACollimatorXPosition = -1257.0 *mm; 
    
    G4double phi = 90. *deg;
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidPMMACollimatorSupport = new G4Box("PMMACollimatorSupport",
                                           PMMACollimatorSupportXSize/2,
                                           PMMACollimatorSupportYSize/2,
                                           PMMACollimatorSupportZSize/2);
    
    logicPMMACollimatorSupport = new G4LogicalVolume(solidPMMACollimatorSupport,
                                                     PMMACollimatorMaterial,
                                                     "PMMACollimatorSupport");
    
    physiPMMACollimatorSupport = new G4PVPlacement(0, G4ThreeVector(PMMACollimatorXPosition,0., 0.),
                                                   "PMMACollimatorSupport",
                                                   logicPMMACollimatorSupport,
                                                   physicalTreatmentRoom,
                                                   false,
                                                   0);
    
    
    yellow = new G4VisAttributes(G4Colour(1., 1., 0. ));
    yellow-> SetVisibility(true);
    yellow-> SetForceWireframe(true);
    
    logicPMMACollimatorSupport -> SetVisAttributes(yellow);
    
    // ----------------------//
    //     PMMA COLLIMATOR   //
    //-----------------------//
    innerRadiusPMMACollimator= 4.5 *cm;
    outerRadiusPMMACollimator = 4.6 *cm;
    hightPMMACollimator = PMMACollimatorSupportXSize;
    startAnglePMMACollimator = 0.*deg;
    spanningAnglePMMACollimator = 360.*deg;
    
    
    solidPMMACollimator = new G4Tubs("PMMACollimator",
                                     innerRadiusPMMACollimator,
                                     outerRadiusPMMACollimator,
                                     hightPMMACollimator/2,
                                     startAnglePMMACollimator,
                                     spanningAnglePMMACollimator);
    
    logicPMMACollimator = new G4LogicalVolume(solidPMMACollimator,
                                              PMMACollimatorMaterial,
                                              "PMMACollimator",
                                              0,
                                              0,
                                              0);
    
    physiPMMACollimator = new G4PVPlacement(G4Transform3D(rm,
                                                          G4ThreeVector(0,0.,0.)),
                                            "PMMACollimator",
                                            logicPMMACollimator,
                                            physiPMMACollimatorSupport,
                                            false,
                                            0);
    
    logicPMMACollimator -> SetVisAttributes(yellow);
    
    
    
    
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}

/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::HadrontherapyBeamMonitoring()
{
    // ----------------------------
    //       MONITOR CHAMBER
    // ----------------------------
    // A monitor chamber is a free-air  ionisation chamber
    // able to measure do carbon fluence during the treatment.
    // Here its responce is not simulated in terms of produced
    // charge but only the energy losses are taked into account.
    // Each chamber consist of 9 mm of air in a box
    // that has two layers one of kapton and one
    // of copper
    
    monitor1XSize = 4.525022*mm;
    monitor2XSize = 0.000011*mm;
    monitor3XSize = 4.5*mm;
    monitorYSize = 10.*cm;
    monitorZSize = 10.*cm;
    monitor1XPosition = -1239.974978 *mm;
    monitor2XPosition = -4.500011*mm;
    monitor4XPosition = 4.500011*mm;
    
    solidFirstMonitorLayer1 = new G4Box("FirstMonitorLayer1",
                                        monitor1XSize,
                                        monitorYSize,
                                        monitorZSize);
    
    logicFirstMonitorLayer1 = new G4LogicalVolume(solidFirstMonitorLayer1,
                                                  layer1MonitorChamberMaterial,
                                                  "FirstMonitorLayer1");
    
    physiFirstMonitorLayer1 = new G4PVPlacement(0,
                                                G4ThreeVector(monitor1XPosition,0.*cm,0.*cm),
                                                "FirstMonitorLayer1",
                                                logicFirstMonitorLayer1,
                                                physicalTreatmentRoom,
                                                false,
                                                0);
    
    solidFirstMonitorLayer2 = new G4Box("FirstMonitorLayer2",
                                        monitor2XSize,
                                        monitorYSize,
                                        monitorZSize);
    
    logicFirstMonitorLayer2 = new G4LogicalVolume(solidFirstMonitorLayer2,
                                                  layer2MonitorChamberMaterial,
                                                  "FirstMonitorLayer2");
    
    physiFirstMonitorLayer2 = new G4PVPlacement(0, G4ThreeVector(monitor2XPosition,0.*cm,0.*cm),
                                                "FirstMonitorLayer2",
                                                logicFirstMonitorLayer2,
                                                physiFirstMonitorLayer1,
                                                false,
                                                0);
    
    solidFirstMonitorLayer3 = new G4Box("FirstMonitorLayer3",
                                        monitor3XSize,
                                        monitorYSize,
                                        monitorZSize);
    
    logicFirstMonitorLayer3 = new G4LogicalVolume(solidFirstMonitorLayer3,
                                                  layer3MonitorChamberMaterial,
                                                  "FirstMonitorLayer3");
    
    physiFirstMonitorLayer3 = new G4PVPlacement(0,
                                                G4ThreeVector(0.*mm,0.*cm,0.*cm),
                                                "MonitorLayer3",
                                                logicFirstMonitorLayer3,
                                                physiFirstMonitorLayer1,
                                                false,
                                                0);
    
    solidFirstMonitorLayer4 = new G4Box("FirstMonitorLayer4",
                                        monitor2XSize,
                                        monitorYSize,
                                        monitorZSize);
    
    logicFirstMonitorLayer4 = new G4LogicalVolume(solidFirstMonitorLayer4,
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
    // a 59 mm thick PMMA slab (to stop all the diffused radiation), a 370 mm brass tube
    // (to well collimate the proton beam) and a final collimator with 25 mm diameter
    // aperture (that provide the final trasversal shape of the beam)
    
    // -------------------//
    //     PMMA SUPPORT   //
    // -------------------//
    
    nozzleSupportXSize = 50 *mm;
    nozzleSupportYSize = 360. *mm;
    nozzleSupportZSize = 360. *mm;
    nozzleSupportXPosition = -423.0 *mm;
    
    G4double phi = 90. *deg;
    // Matrix definition for a 90 deg rotation. Also used for other volumes
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidNozzleSupport = new G4Box("NozzleSupport",
                                   nozzleSupportXSize/2,
                                   nozzleSupportYSize/2,
                                   nozzleSupportZSize/2);
    
    logicNozzleSupport = new G4LogicalVolume(solidNozzleSupport,
                                             nozzleSupportMaterial,
                                             "NozzleSupport");
    
    physiNozzleSupport = new G4PVPlacement(0, G4ThreeVector(nozzleSupportXPosition,0., 0.),
                                           "NozzleSupport",
                                           logicNozzleSupport,
                                           physicalTreatmentRoom,
                                           false,
                                           0);
    
    yellow = new G4VisAttributes(G4Colour(1., 1., 0. ));
    yellow-> SetVisibility(true);
    yellow-> SetForceWireframe(true);
    
    //------------------------------------//
    // HOLE IN THE SUPPORT                //
    //------------------------------------//
    innerRadiusHoleNozzleSupport = 0.*mm;
    outerRadiusHoleNozzleSupport = 22.0*mm;
    hightHoleNozzleSupport = nozzleSupportXSize;
    startAngleHoleNozzleSupport = 0.*deg;
    spanningAngleHoleNozzleSupport = 360.*deg;
    
    solidHoleNozzleSupport = new G4Tubs("HoleNozzleSupport",
                                        innerRadiusHoleNozzleSupport,
                                        outerRadiusHoleNozzleSupport,
                                        hightHoleNozzleSupport/2,
                                        startAngleHoleNozzleSupport,
                                        spanningAngleHoleNozzleSupport);
    
    logicHoleNozzleSupport = new G4LogicalVolume(solidHoleNozzleSupport,
                                                 holeNozzleSupportMaterial,
                                                 "HoleNozzleSupport",
                                                 0,
                                                 0,
                                                 0);
    
    physiHoleNozzleSupport = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                               "HoleNozzleSupport",
                                               logicHoleNozzleSupport,
                                               physiNozzleSupport,
                                               false, 0);
    
    
    // --------------------------------------//
    //     BRASS TUBE 3 (beam line side)    //
    // -------------------------------------//
    innerRadiusBrassTube3 = 13.5 *mm;
    outerRadiusBrassTube3 = 21.5 *mm;
    hightBrassTube3 = 20.0 *mm;
    startAngleBrassTube3 = 0.*deg;
    spanningAngleBrassTube3 = 360.*deg;
    brassTube3XPosition = -458.0 *mm;
    
    solidBrassTube3 = new G4Tubs("BrassTube3",
                                 innerRadiusBrassTube3,
                                 outerRadiusBrassTube3,
                                 hightBrassTube3/2,
                                 startAngleBrassTube3,
                                 spanningAngleBrassTube3);
    
    logicBrassTube3 = new G4LogicalVolume(solidBrassTube3,
                                          brassTube3Material,
                                          "BrassTube3",
                                          0, 0, 0);
    
    physiBrassTube3 = new G4PVPlacement(G4Transform3D(rm,
                                                      G4ThreeVector(brassTube3XPosition,
                                                                    0.,
                                                                    0.)),
                                        "BrassTube3",
                                        logicBrassTube3,
                                        physicalTreatmentRoom,
                                        false,
                                        0);
    
    logicBrassTube3 -> SetVisAttributes(darkOrange3);
    
    // ----------------------------------------------//
    //     BRASS TUBE 2 (inside the PMMA support)    //
    // ----------------------------------------------//
    
    innerRadiusBrassTube2 = 13.5 *mm;
    outerRadiusBrassTube2 = 21.5 *mm;
    hightBrassTube2 = nozzleSupportXSize;
    startAngleBrassTube2 = 0.*deg;
    spanningAngleBrassTube2 = 360.*deg;
    
    
    solidBrassTube2 = new G4Tubs("BrassTube2",
                                 innerRadiusBrassTube2,
                                 outerRadiusBrassTube2,
                                 hightBrassTube2/2,
                                 startAngleBrassTube2,
                                 spanningAngleBrassTube2);
    
    logicBrassTube2 = new G4LogicalVolume(solidBrassTube2,
                                          brassTube2Material,
                                          "BrassTube2",
                                          0, 0, 0);
    physiBrassTube2 = new G4PVPlacement(0,
                                        G4ThreeVector(0,0.,0.),
                                        logicBrassTube2,
                                        "BrassTube2",
                                        logicHoleNozzleSupport,
                                        false,
                                        0);
    
    logicBrassTube2 -> SetVisAttributes(darkOrange3);
    
    // ---------------------------------//
    //     BRASS TUBE 1 (phantom side)    //
    // ---------------------------------//
    innerRadiusBrassTube= 18.*mm;
    outerRadiusBrassTube = 21.5 *mm;
    hightBrassTube = 208.0 *mm;
    startAngleBrassTube = 0.*deg;
    spanningAngleBrassTube = 360.*deg;
    brassTubeXPosition = -294 *mm;
    solidBrassTube = new G4Tubs("BrassTube",
                                innerRadiusBrassTube,
                                outerRadiusBrassTube,
                                hightBrassTube/2,
                                startAngleBrassTube,
                                spanningAngleBrassTube);
    
    logicBrassTube = new G4LogicalVolume(solidBrassTube,
                                         brassTubeMaterial,
                                         "BrassTube",
                                         0, 0, 0);
    
    physiBrassTube = new G4PVPlacement(G4Transform3D(rm,
                                                     G4ThreeVector(brassTubeXPosition,
                                                                   0.,
                                                                   0.)),
                                       "BrassTube",
                                       logicBrassTube,
                                       physicalTreatmentRoom,
                                       false,
                                       0);
    
    logicBrassTube -> SetVisAttributes(darkOrange3);
}

/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::HadrontherapyBeamFinalCollimator()
{
    // -----------------------//
    //     FINAL COLLIMATOR   //
    //------------------------//
    outerRadiusFinalCollimator = 21.5*mm;
    hightFinalCollimator = 7.0 *mm;
    startAngleFinalCollimator = 0.*deg;
    spanningAngleFinalCollimator = 360.*deg;
    finalCollimatorXPosition = -186.5 *mm;
    
    G4double phi = 90. *deg;
    
    // Matrix definition for a 90 deg rotation. Also used for other volumes
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidFinalCollimator = new G4Tubs("FinalCollimator",
                                      innerRadiusFinalCollimator,
                                      outerRadiusFinalCollimator,
                                      hightFinalCollimator/2,
                                      startAngleFinalCollimator,
                                      spanningAngleFinalCollimator);
    
    logicFinalCollimator = new G4LogicalVolume(solidFinalCollimator,
                                               finalCollimatorMaterial,
                                               "FinalCollimator",
                                               0,
                                               0,
                                               0);
    
    physiFinalCollimator = new G4PVPlacement(G4Transform3D(rm,
                                                           G4ThreeVector(finalCollimatorXPosition,0.,0.)),
                                             "FinalCollimator",
                                             logicFinalCollimator,
                                             physicalTreatmentRoom,
                                             false,
                                             0);
    
    logicFinalCollimator -> SetVisAttributes(yellow);
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////// MESSENGER ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void PassiveCarbonBeamLine::SetRippleFilterXPosition(G4double value)
{
    PhysiRippleFilter -> SetTranslation(G4ThreeVector(value, 0., 0.));
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout << "The Ripple Filter is translated to"<< value/mm <<"mm along the X axis" <<G4endl;
}


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

