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

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "G4NistElementBuilder.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyModulator.hh"
#include "PassiveProtonBeamLine.hh"
#include "PassiveProtonBeamLineMessenger.hh"


//G4bool PassiveProtonBeamLine::doCalculation = false;
/////////////////////////////////////////////////////////////////////////////
PassiveProtonBeamLine::PassiveProtonBeamLine():
modulator(0), physicalTreatmentRoom(0),hadrontherapyDetectorConstruction(0),
physiBeamLineSupport(0), physiBeamLineCover(0), physiBeamLineCover2(0),
firstScatteringFoil(0), physiFirstScatteringFoil(0), physiKaptonWindow(0),
solidStopper(0), physiStopper(0), secondScatteringFoil(0), physiSecondScatteringFoil(0),
physiFirstCollimator(0), solidRangeShifterBox(0), logicRangeShifterBox(0),
physiRangeShifterBox(0), physiSecondCollimator(0), physiFirstCollimatorModulatorBox(0),
physiHoleFirstCollimatorModulatorBox(0), physiSecondCollimatorModulatorBox(0),
physiHoleSecondCollimatorModulatorBox(0), physiMOPIMotherVolume(0),
physiFirstMonitorLayer1(0), physiFirstMonitorLayer2(0), physiFirstMonitorLayer3(0),
physiFirstMonitorLayer4(0), physiSecondMonitorLayer1(0), physiSecondMonitorLayer2(0),
physiSecondMonitorLayer3(0), physiSecondMonitorLayer4(0), physiNozzleSupport(0), physiBrassTube(0), solidFinalCollimator(0), physiFinalCollimator(0)
{
    // Messenger to change parameters of the passiveProtonBeamLine geometry
    passiveMessenger = new PassiveProtonBeamLineMessenger(this);
    
    //***************************** PW ***************************************
    static G4String ROGeometryName = "DetectorROGeometry";
    RO = new HadrontherapyDetectorROGeometry(ROGeometryName);
    
    G4cout << "Going to register Parallel world...";
    RegisterParallelWorld(RO);
    G4cout << "... done" << G4endl;

}
/////////////////////////////////////////////////////////////////////////////
PassiveProtonBeamLine::~PassiveProtonBeamLine()
{
    delete passiveMessenger;
    delete hadrontherapyDetectorConstruction;

}

/////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* PassiveProtonBeamLine::Construct()
{
    // Sets default geometry and materials
    SetDefaultDimensions();
    
    // Construct the whole Passive Beam Line
    ConstructPassiveProtonBeamLine();
    
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
void PassiveProtonBeamLine::SetDefaultDimensions()
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
    // The PIPE contains the FIRST SCATTERING FOIL and the KAPTON WINDOW
    G4double defaultVacuumZoneXSize = 100.0 *mm;
    vacuumZoneXSize = defaultVacuumZoneXSize;
    
    G4double defaultVacuumZoneYSize = 52.5 *mm;
    vacuumZoneYSize = defaultVacuumZoneYSize;
    
    G4double defaultVacuumZoneZSize = 52.5 *mm;
    vacuumZoneZSize = defaultVacuumZoneZSize;
    
    G4double defaultVacuumZoneXPosition = -3010.0 *mm;
    vacuumZoneXPosition = defaultVacuumZoneXPosition;
    
    // FIRST SCATTERING FOIL: a thin foil performing a first scattering
    // of the original beam
    G4double defaultFirstScatteringFoilXSize = 0.0075 *mm;
    firstScatteringFoilXSize = defaultFirstScatteringFoilXSize;
    
    G4double defaultFirstScatteringFoilYSize = 52.5   *mm;
    firstScatteringFoilYSize = defaultFirstScatteringFoilYSize;
    
    G4double defaultFirstScatteringFoilZSize = 52.5   *mm;
    firstScatteringFoilZSize = defaultFirstScatteringFoilZSize;
    
    G4double defaultFirstScatteringFoilXPosition = 0.0 *mm;
    firstScatteringFoilXPosition = defaultFirstScatteringFoilXPosition;
    
    // KAPTON WINDOW: it prmits the passage of the beam from vacuum to air
    G4double defaultKaptonWindowXSize = 0.010*mm;
    kaptonWindowXSize = defaultKaptonWindowXSize;
    
    G4double defaultKaptonWindowYSize = 5.25*cm;
    kaptonWindowYSize = defaultKaptonWindowYSize;
    
    G4double defaultKaptonWindowZSize = 5.25*cm;
    kaptonWindowZSize = defaultKaptonWindowZSize;
    
    G4double defaultKaptonWindowXPosition = 100.0*mm - defaultKaptonWindowXSize;
    kaptonWindowXPosition = defaultKaptonWindowXPosition;
    
    // STOPPER: is a small cylinder able to stop the central component
    // of the beam (having a gaussian shape). It is connected to the SECON SCATTERING FOIL
    // and represent the second element of the scattering system
    G4double defaultInnerRadiusStopper = 0.*cm;
    innerRadiusStopper = defaultInnerRadiusStopper;
    
    G4double defaultHeightStopper = 3.5*mm;
    heightStopper = defaultHeightStopper;
    
    G4double defaultStartAngleStopper = 0.*deg;
    startAngleStopper = defaultStartAngleStopper;
    
    G4double defaultSpanningAngleStopper = 360.*deg;
    spanningAngleStopper = defaultSpanningAngleStopper;
    
    G4double defaultStopperXPosition = -2705.0 *mm;
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
    G4double defaultSecondScatteringFoilXSize = 0.0125 *mm;
    secondScatteringFoilXSize = defaultSecondScatteringFoilXSize;
    
    G4double defaultSecondScatteringFoilYSize = 52.5   *mm;
    secondScatteringFoilYSize = defaultSecondScatteringFoilYSize;
    
    G4double defaultSecondScatteringFoilZSize = 52.5   *mm;
    secondScatteringFoilZSize = defaultSecondScatteringFoilZSize;
    
    G4double defaultSecondScatteringFoilXPosition = defaultStopperXPosition + defaultHeightStopper + defaultSecondScatteringFoilXSize;
    secondScatteringFoilXPosition = defaultSecondScatteringFoilXPosition;
    
    G4double defaultSecondScatteringFoilYPosition =  0 *mm;
    secondScatteringFoilYPosition = defaultSecondScatteringFoilYPosition;
    
    G4double defaultSecondScatteringFoilZPosition =  0 *mm;
    secondScatteringFoilZPosition = defaultSecondScatteringFoilZPosition;
    
    // RANGE SHIFTER: is a slab of PMMA acting as energy degreader of
    // primary beam
    
    //Default material of the range shifter
    
    G4double defaultRangeShifterXSize = 5. *mm;
    rangeShifterXSize = defaultRangeShifterXSize;
    
    G4double defaultRangeShifterYSize = 176. *mm;
    rangeShifterYSize = defaultRangeShifterYSize;
    
    G4double defaultRangeShifterZSize = 176. *mm;
    rangeShifterZSize = defaultRangeShifterZSize;
    
    G4double defaultRangeShifterXPosition = -2393.0 *mm;
    rangeShifterXPosition = defaultRangeShifterXPosition;
    
    G4double defaultRangeShifterYPosition = 0. *mm;
    rangeShifterYPosition = defaultRangeShifterYPosition;
    
    G4double defaultRangeShifterZPosition = 0. *mm;
    rangeShifterZPosition = defaultRangeShifterZPosition;
    
    // MOPI DETECTOR: two orthogonal microstrip gas detectors developed
    // by the INFN Section of Turin in collaboration with some
    // of the author of this example. It permits the
    // on-line check of the beam simmetry via the signal
    // integration of the collected charge for each strip.
    
    // Mother volume of MOPI
    
    G4double defaultMOPIMotherVolumeXSize = 12127.0 *um;
    MOPIMotherVolumeXSize = defaultMOPIMotherVolumeXSize;
    
    G4double defaultMOPIMotherVolumeYSize = 40.0 *cm;
    MOPIMotherVolumeYSize = defaultMOPIMotherVolumeYSize;
    
    G4double defaultMOPIMotherVolumeZSize = 40.0 *cm;
    MOPIMotherVolumeZSize = defaultMOPIMotherVolumeZSize;
    
    G4double defaultMOPIMotherVolumeXPosition = -1000.0 *mm;
    MOPIMotherVolumeXPosition = defaultMOPIMotherVolumeXPosition;
    
    G4double defaultMOPIMotherVolumeYPosition = 0.0 *mm;
    MOPIMotherVolumeYPosition = defaultMOPIMotherVolumeYPosition;
    
    G4double defaultMOPIMotherVolumeZPosition = 0.0 *mm;
    MOPIMotherVolumeZPosition = defaultMOPIMotherVolumeZPosition;
    
    // First Kapton Layer of MOPI
    G4double defaultMOPIFirstKaptonLayerXSize = 35 *um;
    MOPIFirstKaptonLayerXSize = defaultMOPIFirstKaptonLayerXSize;
    
    G4double defaultMOPIFirstKaptonLayerYSize = 30 *cm;
    MOPIFirstKaptonLayerYSize = defaultMOPIFirstKaptonLayerYSize;
    
    G4double defaultMOPIFirstKaptonLayerZSize = 30 *cm;
    MOPIFirstKaptonLayerZSize = defaultMOPIFirstKaptonLayerZSize;
    
    G4double defaultMOPIFirstKaptonLayerXPosition = -(MOPIMotherVolumeXSize/2 - (MOPIFirstKaptonLayerXSize/2));
    MOPIFirstKaptonLayerXPosition = defaultMOPIFirstKaptonLayerXPosition;
    
    G4double defaultMOPIFirstKaptonLayerYPosition = 0.0 *mm;
    MOPIFirstKaptonLayerYPosition = defaultMOPIFirstKaptonLayerYPosition;
    
    G4double defaultMOPIFirstKaptonLayerZPosition = 0.0 *mm;
    MOPIFirstKaptonLayerZPosition = defaultMOPIFirstKaptonLayerZPosition;
    
    //First Aluminum  Layer of MOPI
    G4double defaultMOPIFirstAluminumLayerXSize = 15 *um;
    MOPIFirstAluminumLayerXSize = defaultMOPIFirstAluminumLayerXSize;
    
    G4double defaultMOPIFirstAluminumLayerYSize = 30 *cm;
    MOPIFirstAluminumLayerYSize = defaultMOPIFirstAluminumLayerYSize;
    
    G4double defaultMOPIFirstAluminumLayerZSize = 30 *cm;
    MOPIFirstAluminumLayerZSize = defaultMOPIFirstAluminumLayerZSize;
    
    G4double defaultMOPIFirstAluminumLayerXPosition =
    MOPIFirstKaptonLayerXPosition + MOPIFirstKaptonLayerXSize/2 + MOPIFirstAluminumLayerXSize/2;
    MOPIFirstAluminumLayerXPosition = defaultMOPIFirstAluminumLayerXPosition;
    
    G4double defaultMOPIFirstAluminumLayerYPosition = 0.0 *mm;
    MOPIFirstAluminumLayerYPosition = defaultMOPIFirstAluminumLayerYPosition;
    
    G4double defaultMOPIFirstAluminumLayerZPosition = 0.0 *mm;
    MOPIFirstAluminumLayerZPosition = defaultMOPIFirstAluminumLayerZPosition;
    
    // First Air gap of MOPI
    G4double defaultMOPIFirstAirGapXSize = 6000 *um;
    MOPIFirstAirGapXSize = defaultMOPIFirstAirGapXSize;
    
    G4double defaultMOPIFirstAirGapYSize = 30 *cm;
    MOPIFirstAirGapYSize = defaultMOPIFirstAirGapYSize;
    
    G4double defaultMOPIFirstAirGapZSize = 30 *cm;
    MOPIFirstAirGapZSize = defaultMOPIFirstAirGapZSize;
    
    G4double defaultMOPIFirstAirGapXPosition =
    MOPIFirstAluminumLayerXPosition + MOPIFirstAluminumLayerXSize/2 + MOPIFirstAirGapXSize/2;
    MOPIFirstAirGapXPosition = defaultMOPIFirstAirGapXPosition;
    
    G4double defaultMOPIFirstAirGapYPosition = 0.0 *mm;
    MOPIFirstAirGapYPosition = defaultMOPIFirstAirGapYPosition;
    
    G4double defaultMOPIFirstAirGapZPosition = 0.0 *mm;
    MOPIFirstAirGapZPosition = defaultMOPIFirstAirGapZPosition;
    
    // Cathode of MOPI
    G4double defaultMOPICathodeXSize = 25.0 *um;
    MOPICathodeXSize = defaultMOPICathodeXSize;
    
    G4double defaultMOPICathodeYSize = 30.0 *cm;
    MOPICathodeYSize = defaultMOPICathodeYSize;
    
    G4double defaultMOPICathodeZSize = 30.0 *cm;
    MOPICathodeZSize = defaultMOPICathodeZSize;
    
    G4double defaultMOPICathodeXPosition =
    MOPIFirstAirGapXPosition + MOPIFirstAirGapXSize/2 + MOPICathodeXSize/2;
    MOPICathodeXPosition = defaultMOPICathodeXPosition;
    
    G4double defaultMOPICathodeYPosition = 0.0 *mm;
    MOPICathodeYPosition = defaultMOPICathodeYPosition;
    
    G4double defaultMOPICathodeZPosition = 0.0 *mm;
    MOPICathodeZPosition = defaultMOPICathodeZPosition;
    
    // Second Air gap of MOPI
    G4double defaultMOPISecondAirGapXSize = 6000 *um;
    MOPISecondAirGapXSize = defaultMOPISecondAirGapXSize;
    
    G4double defaultMOPISecondAirGapYSize = 30 *cm;
    MOPISecondAirGapYSize = defaultMOPISecondAirGapYSize;
    
    G4double defaultMOPISecondAirGapZSize = 30 *cm;
    MOPISecondAirGapZSize = defaultMOPISecondAirGapZSize;
    
    G4double defaultMOPISecondAirGapXPosition =
    MOPICathodeXPosition + MOPICathodeXSize/2 + MOPISecondAirGapXSize/2;
    MOPISecondAirGapXPosition = defaultMOPISecondAirGapXPosition;
    
    G4double defaultMOPISecondAirGapYPosition = 0.0 *mm;
    MOPISecondAirGapYPosition = defaultMOPISecondAirGapYPosition;
    
    G4double defaultMOPISecondAirGapZPosition = 0.0 *mm;
    MOPISecondAirGapZPosition = defaultMOPISecondAirGapZPosition;
    
    //Second Aluminum  Layer of MOPI
    G4double defaultMOPISecondAluminumLayerXSize = 15 *um;
    MOPISecondAluminumLayerXSize = defaultMOPISecondAluminumLayerXSize;
    
    G4double defaultMOPISecondAluminumLayerYSize = 30 *cm;
    MOPISecondAluminumLayerYSize = defaultMOPISecondAluminumLayerYSize;
    
    G4double defaultMOPISecondAluminumLayerZSize = 30 *cm;
    MOPISecondAluminumLayerZSize = defaultMOPISecondAluminumLayerZSize;
    
    G4double defaultMOPISecondAluminumLayerXPosition =
    MOPISecondAirGapXPosition + MOPISecondAirGapXSize/2 + MOPISecondAluminumLayerXSize/2;
    MOPISecondAluminumLayerXPosition = defaultMOPISecondAluminumLayerXPosition;
    
    G4double defaultMOPISecondAluminumLayerYPosition = 0.0 *mm;
    MOPISecondAluminumLayerYPosition = defaultMOPISecondAluminumLayerYPosition;
    
    G4double defaultMOPISecondAluminumLayerZPosition = 0.0 *mm;
    MOPISecondAluminumLayerZPosition = defaultMOPISecondAluminumLayerZPosition;
    
    // Second Kapton Layer of MOPI
    G4double defaultMOPISecondKaptonLayerXSize = 35 *um;
    MOPISecondKaptonLayerXSize = defaultMOPISecondKaptonLayerXSize;
    
    G4double defaultMOPISecondKaptonLayerYSize = 30 *cm;
    MOPISecondKaptonLayerYSize = defaultMOPISecondKaptonLayerYSize;
    
    G4double defaultMOPISecondKaptonLayerZSize = 30 *cm;
    MOPISecondKaptonLayerZSize = defaultMOPISecondKaptonLayerZSize;
    
    G4double defaultMOPISecondKaptonLayerXPosition =
    MOPISecondAluminumLayerXPosition + MOPISecondAluminumLayerXSize/2 + MOPISecondKaptonLayerXSize/2;
    MOPISecondKaptonLayerXPosition = defaultMOPISecondKaptonLayerXPosition;
    
    G4double defaultMOPISecondKaptonLayerYPosition = 0.0 *mm;
    MOPISecondKaptonLayerYPosition = defaultMOPISecondKaptonLayerYPosition;
    
    G4double defaultMOPISecondKaptonLayerZPosition = 0.0 *mm;
    MOPISecondKaptonLayerZPosition = defaultMOPISecondKaptonLayerZPosition;
    
    
    // FINAL COLLIMATOR: is the collimator giving the final transversal shape
    // of the beam
    G4double defaultinnerRadiusFinalCollimator = 7.5 *mm;
    innerRadiusFinalCollimator = defaultinnerRadiusFinalCollimator;
    
    // DEFAULT DEFINITION OF THE MATERIALS
    // All elements and compound definition follows the NIST database
    
    // ELEMENTS
    G4bool isotopes = false;
    G4Material* aluminumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al", isotopes);
    G4Material* tantalumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ta", isotopes);
    G4Material* copperNistAsMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);
    G4Element* zincNist = G4NistManager::Instance()->FindOrBuildElement("Zn");
    G4Element* copperNist = G4NistManager::Instance()->FindOrBuildElement("Cu");
    
    // COMPOUND
    G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
    G4Material* kaptonNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON", isotopes);
    G4Material* galacticNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic", isotopes);
    G4Material* PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);
    G4Material* mylarNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR", isotopes);
    
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
    // Range shifter
    rangeShifterMaterial = airNist;
    
    // Support of the beam line
    beamLineSupportMaterial = aluminumNist;
    
    // Vacuum pipe
    vacuumZoneMaterial = galacticNist;
    
    // Material of the fisrt scattering foil
    firstScatteringFoilMaterial = tantalumNist;
    
    // Material of kapton window
    kaptonWindowMaterial = kaptonNist;
    
    // Material of the stopper
    stopperMaterial = brass;
    
    // Material of the second scattering foil
    secondScatteringFoilMaterial = tantalumNist;
    
    // Materials of the collimators
    firstCollimatorMaterial = PMMANist;
    holeFirstCollimatorMaterial = airNist;
    
    // Box containing the modulator wheel
    modulatorBoxMaterial = aluminumNist;
    holeModulatorBoxMaterial = airNist;
    
    // Materials of the monitor chamber
    layer1MonitorChamberMaterial = kaptonNist;
    layer2MonitorChamberMaterial = copperNistAsMaterial;
    layer3MonitorChamberMaterial = airNist;
    layer4MonitorChamberMaterial = copperNistAsMaterial;
    
    // Mother volume of the MOPI detector
    MOPIMotherVolumeMaterial = airNist;
    MOPIFirstKaptonLayerMaterial = kaptonNist;
    MOPIFirstAluminumLayerMaterial = aluminumNist;
    MOPIFirstAirGapMaterial = airNist;
    MOPICathodeMaterial = mylarNist;
    MOPISecondAirGapMaterial = airNist;
    MOPISecondAluminumLayerMaterial = aluminumNist;
    MOPISecondKaptonLayerMaterial = kaptonNist;
    
    // material of the final nozzle
    nozzleSupportMaterial = PMMANist;
    brassTubeMaterial = brassTube2Material = brassTube3Material = brass;
    holeNozzleSupportMaterial = airNist;
    
    // Material of the final collimator
    finalCollimatorMaterial = brass;
}

/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::ConstructPassiveProtonBeamLine()
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
    
    // Components of the Passive Proton Beam Line
    HadrontherapyBeamLineSupport();
    HadrontherapyBeamScatteringFoils();
    HadrontherapyRangeShifter();
    HadrontherapyBeamCollimators();
    HadrontherapyBeamMonitoring();
    HadrontherapyMOPIDetector();
    HadrontherapyBeamNozzle();
    HadrontherapyBeamFinalCollimator();
    
    // The following lines construc a typical modulator wheel inside the Passive Beam line.
    // Please remember to set the nodulator material (default is air, i.e. no modulator!)
    // in the HadrontherapyModulator.cc file
    modulator = new HadrontherapyModulator();
    modulator -> BuildModulator(physicalTreatmentRoom);
}

/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::HadrontherapyBeamLineSupport()
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
    const G4double beamLineCoverZPosition = 600.*mm;
    
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
void PassiveProtonBeamLine::HadrontherapyBeamScatteringFoils()
{
   // ------------//
    // VACUUM PIPE //
    //-------------//
    //
    // First track of the beam line is inside vacuum;
    // The PIPE contains the FIRST SCATTERING FOIL and the KAPTON WINDOW
    G4Box* vacuumZone = new G4Box("VacuumZone", vacuumZoneXSize, vacuumZoneYSize, vacuumZoneZSize);
    G4LogicalVolume* logicVacuumZone = new G4LogicalVolume(vacuumZone, vacuumZoneMaterial, "VacuumZone");
    G4VPhysicalVolume* physiVacuumZone = new G4PVPlacement(0, G4ThreeVector(vacuumZoneXPosition, 0., 0.),
                                                           "VacuumZone", logicVacuumZone, physicalTreatmentRoom, false, 0);
    // --------------------------//
    // THE FIRST SCATTERING FOIL //
    // --------------------------//
    // A thin foil performing a first scattering
    // of the original beam
    firstScatteringFoil = new G4Box("FirstScatteringFoil",
                                    firstScatteringFoilXSize,
                                    firstScatteringFoilYSize,
                                    firstScatteringFoilZSize);
    
    G4LogicalVolume* logicFirstScatteringFoil = new G4LogicalVolume(firstScatteringFoil,
                                                                    firstScatteringFoilMaterial,
                                                                    "FirstScatteringFoil");
    
    physiFirstScatteringFoil = new G4PVPlacement(0, G4ThreeVector(firstScatteringFoilXPosition, 0.,0.),
                                                 "FirstScatteringFoil", logicFirstScatteringFoil, physiVacuumZone,
                                                 false, 0);
    
    logicFirstScatteringFoil -> SetVisAttributes(skyBlue);
    // -------------------//
    // THE KAPTON WINDOWS //
    //--------------------//
    //It prmits the passage of the beam from vacuum to air
 
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
                              heightStopper,
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
                                     secondScatteringFoilXSize,
                                     secondScatteringFoilYSize,
                                     secondScatteringFoilZSize);
    
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
void PassiveProtonBeamLine::HadrontherapyRangeShifter()
{
    // ---------------------------- //
    //         THE RANGE SHIFTER    //
    // -----------------------------//
    // It is a slab of PMMA acting as energy degreader of
    // primary beam
 
    
    solidRangeShifterBox = new G4Box("RangeShifterBox",
                                     rangeShifterXSize,
                                     rangeShifterYSize,
                                     rangeShifterZSize);
    
    logicRangeShifterBox = new G4LogicalVolume(solidRangeShifterBox,
                                               rangeShifterMaterial,
                                               "RangeShifterBox");
    physiRangeShifterBox = new G4PVPlacement(0,
                                             G4ThreeVector(rangeShifterXPosition, 0., 0.),
                                             "RangeShifterBox",
                                             logicRangeShifterBox,
                                             physicalTreatmentRoom,
                                             false,
                                             0);
    
    
    logicRangeShifterBox -> SetVisAttributes(yellow);
     
    
}
/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::HadrontherapyBeamCollimators()
{
    
    
    // -----------------//
    // FIRST COLLIMATOR //
    // -----------------//
    // It is a slab of PMMA with an hole in its center
    const G4double firstCollimatorXSize = 20.*mm;
    const G4double firstCollimatorYSize = 100.*mm;
    const G4double firstCollimatorZSize = 100.*mm;
    
    const G4double firstCollimatorXPosition = -2673.00*mm;
    const G4double firstCollimatorYPosition = 0.*mm;
    const G4double firstCollimatorZPosition = 0.*mm;
    
    
    G4Box* solidFirstCollimator = new G4Box("FirstCollimator",
                                            firstCollimatorXSize,
                                            firstCollimatorYSize,
                                            firstCollimatorZSize);
    
    G4LogicalVolume* logicFirstCollimator = new G4LogicalVolume(solidFirstCollimator,
                                                                firstCollimatorMaterial,
                                                                "FirstCollimator");
    
    physiFirstCollimator = new G4PVPlacement(0, G4ThreeVector(firstCollimatorXPosition,
                                                              firstCollimatorYPosition,
                                                              firstCollimatorZPosition),
                                             "FirstCollimator",
                                             logicFirstCollimator,
                                             physicalTreatmentRoom,
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
    
    G4LogicalVolume* logicHoleFirstCollimator = new G4LogicalVolume(solidHoleFirstCollimator,
                                                                    holeFirstCollimatorMaterial,
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
    // ------------------//
    // SECOND COLLIMATOR //
    //-------------------//
    // It is a slab of PMMA with an hole in its center
    const G4double secondCollimatorXPosition = -1900.00*mm;
    const G4double secondCollimatorYPosition =  0*mm;
    const G4double secondCollimatorZPosition =  0*mm;
    
    physiSecondCollimator = new G4PVPlacement(0, G4ThreeVector(secondCollimatorXPosition,
                                                               secondCollimatorYPosition,
                                                               secondCollimatorZPosition),
                                              "SecondCollimator",
                                              logicFirstCollimator,
                                              physicalTreatmentRoom,
                                              false,
                                              0);
    
    // ------------------------------//
    // Hole of the second collimator //
    // ------------------------------//
    physiHoleSecondCollimator = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                                  "HoleSecondCollimator",
                                                  logicHoleFirstCollimator,
                                                  physiSecondCollimator,
                                                  false,
                                                  0);
    
    // --------------------------------------//
    // FIRST SIDE OF THE MODULATOR BOX      //
    // --------------------------------------//
    // The modulator box is an aluminum box in which
    // the range shifter and the energy modulator are located
    // In this example only the entrance and exit
    // faces of the box are simulated.
    // Each face is an aluminum slab with an hole in its center
    
    const G4double firstCollimatorModulatorXSize = 10.*mm;
    const G4double firstCollimatorModulatorYSize = 200.*mm;
    const G4double firstCollimatorModulatorZSize = 200.*mm;
    
    const G4double firstCollimatorModulatorXPosition = -2523.00*mm;
    const G4double firstCollimatorModulatorYPosition = 0.*mm;
    const G4double firstCollimatorModulatorZPosition = 0.*mm;
    
   G4Box* solidFirstCollimatorModulatorBox = new G4Box("FirstCollimatorModulatorBox",
                                                        firstCollimatorModulatorXSize,
                                                        firstCollimatorModulatorYSize,
                                                        firstCollimatorModulatorZSize);
    
    G4LogicalVolume* logicFirstCollimatorModulatorBox = new G4LogicalVolume(solidFirstCollimatorModulatorBox,
                                                                            modulatorBoxMaterial,
                                                                            "FirstCollimatorModulatorBox");
    
    physiFirstCollimatorModulatorBox = new G4PVPlacement(0, G4ThreeVector(firstCollimatorModulatorXPosition,
                                                                          firstCollimatorModulatorYPosition,
                                                                          firstCollimatorModulatorZPosition),
                                                         "FirstCollimatorModulatorBox",
                                                         logicFirstCollimatorModulatorBox,
                                                         physicalTreatmentRoom, false, 0);
    
    // ----------------------------------------------------//
    //   Hole of the first collimator of the modulator box //
    // ----------------------------------------------------//
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
                                                                                holeModulatorBoxMaterial,
                                                                                "HoleFirstCollimatorModulatorBox",
                                                                                0, 0, 0);
    
    physiHoleFirstCollimatorModulatorBox = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                                             "HoleFirstCollimatorModulatorBox",
                                                             logicHoleFirstCollimatorModulatorBox,
                                                             physiFirstCollimatorModulatorBox, false, 0);
    
    // --------------------------------------------------//
    //       SECOND SIDE OF THE MODULATOR BOX            //
    // --------------------------------------------------//
    const G4double secondCollimatorModulatorXSize = 10.*mm;
    const G4double secondCollimatorModulatorYSize = 200.*mm;
    const G4double secondCollimatorModulatorZSize = 200.*mm;
    
    const G4double secondCollimatorModulatorXPosition = -1953.00 *mm;
    
    const G4double secondCollimatorModulatorYPosition = 0.*mm;
    const G4double secondCollimatorModulatorZPosition = 0.*mm;
    
    G4Box* solidSecondCollimatorModulatorBox = new G4Box("SecondCollimatorModulatorBox",
                                                         secondCollimatorModulatorXSize,
                                                         secondCollimatorModulatorYSize,
                                                         secondCollimatorModulatorZSize);
    
    G4LogicalVolume* logicSecondCollimatorModulatorBox = new G4LogicalVolume(solidSecondCollimatorModulatorBox,
                                                                             modulatorBoxMaterial,
                                                                             "SecondCollimatorModulatorBox");
    
    physiSecondCollimatorModulatorBox = new G4PVPlacement(0, G4ThreeVector(secondCollimatorModulatorXPosition,
                                                                           secondCollimatorModulatorYPosition,
                                                                           secondCollimatorModulatorZPosition),
                                                          "SecondCollimatorModulatorBox",
                                                          logicSecondCollimatorModulatorBox,
                                                          physicalTreatmentRoom, false, 0);
    
    // ----------------------------------------------//
    //   Hole of the second collimator modulator box //
    // ----------------------------------------------//
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
                                                                                 holeModulatorBoxMaterial,
                                                                                 "HoleSecondCollimatorModulatorBox",
                                                                                 0, 0, 0);
    
    physiHoleSecondCollimatorModulatorBox = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                                              "HoleSecondCollimatorModulatorBox",
                                                              logicHoleSecondCollimatorModulatorBox,
                                                              physiSecondCollimatorModulatorBox, false, 0);
    
    logicFirstCollimator -> SetVisAttributes(yellow);
    logicFirstCollimatorModulatorBox -> SetVisAttributes(blue);
    logicSecondCollimatorModulatorBox -> SetVisAttributes(blue);


  
  }

/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::HadrontherapyBeamMonitoring()
{
    // ----------------------------
    //   THE FIRST MONITOR CHAMBER
    // ----------------------------
    // A monitor chamber is a free-air  ionisation chamber
    // able to measure do proton fluence during the treatment.
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
    const G4double monitor1XPosition = -1262.47498 *mm;
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
    // ----------------------------//
    // THE SECOND MONITOR CHAMBER  //
    // ----------------------------//
    physiSecondMonitorLayer1 = new G4PVPlacement(0, G4ThreeVector(-1131.42493 *mm,0.*cm,0.*cm),
                                                 "SecondMonitorLayer1", logicFirstMonitorLayer1,physicalTreatmentRoom, false, 0);
    
    physiSecondMonitorLayer2 = new G4PVPlacement(0, G4ThreeVector( monitor2XPosition,0.*cm,0.*cm), "SecondMonitorLayer2",
                                                 logicFirstMonitorLayer2, physiSecondMonitorLayer1, false, 0);
    
    physiSecondMonitorLayer3 = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*cm,0.*cm), "MonitorLayer3",
                                                 logicFirstMonitorLayer3, physiSecondMonitorLayer1, false, 0);
    
    physiSecondMonitorLayer4 = new G4PVPlacement(0, G4ThreeVector(monitor4XPosition,0.*cm,0.*cm), "SecondMonitorLayer4",
                                                 logicFirstMonitorLayer4, physiSecondMonitorLayer1, false, 0);
    
    logicFirstMonitorLayer3 -> SetVisAttributes(white);
    

    
     }
/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::HadrontherapyMOPIDetector()
{

    // --------------------------------//
    //        THE MOPI DETECTOR        //
    // --------------------------------//
    // MOPI DETECTOR: two orthogonal microstrip gas detectors developed
    // by the INFN Section of Turin in collaboration with some
    // of the author of this example. It permits the
    // on-line check of the beam simmetry via the signal
    // integration of the collected charge for each strip.
    //
    // In this example it is simulated as:
    // 1. First anode: 35 mu of kapton + 15 mu of aluminum,
    // 2. First air gap: 6 mm of air,
    // 3. The cathode: 1 mu Al + 25 mu mylar + 1 mu Al
    //    (in common with the two air gap),
    // 4. Second air gap: 6 mm of air,
    // 5  Second anode: 15 mu Al + 35 mu kapton
    // Color used in the graphical output
    
    
    // Mother volume
    solidMOPIMotherVolume = new G4Box("MOPIMotherVolume",
                                      MOPIMotherVolumeXSize/2,
                                      MOPIMotherVolumeYSize/2,
                                      MOPIMotherVolumeYSize/2);
    
    logicMOPIMotherVolume = new G4LogicalVolume(solidMOPIMotherVolume,
                                                MOPIMotherVolumeMaterial,
                                                "MOPIMotherVolume");
    physiMOPIMotherVolume = new G4PVPlacement(0,
                                              G4ThreeVector(MOPIMotherVolumeXPosition,
                                                            MOPIMotherVolumeYPosition,
                                                            MOPIMotherVolumeZPosition),
                                              "MOPIMotherVolume",
                                              logicMOPIMotherVolume,
                                              physicalTreatmentRoom,
                                              false,
                                              0);
    
    // First Kapton layer
    solidMOPIFirstKaptonLayer = new G4Box("MOPIFirstKaptonLayer",
                                          MOPIFirstKaptonLayerXSize/2,
                                          MOPIFirstKaptonLayerYSize/2 ,
                                          MOPIFirstKaptonLayerZSize/2);
    
    logicMOPIFirstKaptonLayer = new G4LogicalVolume(solidMOPIFirstKaptonLayer,
                                                    MOPIFirstKaptonLayerMaterial,
                                                    "MOPIFirstKaptonLayer");
    
    physiMOPIFirstKaptonLayer = new G4PVPlacement(0,
                                                  G4ThreeVector(MOPIFirstKaptonLayerXPosition,
                                                                MOPIFirstKaptonLayerYPosition ,
                                                                MOPIFirstKaptonLayerZPosition),
                                                  "MOPIFirstKaptonLayer",
                                                  logicMOPIFirstKaptonLayer,
                                                  physiMOPIMotherVolume,
                                                  false,
                                                  0);
    
    // First Aluminum layer
    solidMOPIFirstAluminumLayer = new G4Box("MOPIFirstAluminumLayer",
                                            MOPIFirstAluminumLayerXSize/2,
                                            MOPIFirstAluminumLayerYSize/2 ,
                                            MOPIFirstAluminumLayerZSize/2);
    
    logicMOPIFirstAluminumLayer = new G4LogicalVolume(solidMOPIFirstAluminumLayer,
                                                      MOPIFirstAluminumLayerMaterial,
                                                      "MOPIFirstAluminumLayer");
    
    physiMOPIFirstAluminumLayer = new G4PVPlacement(0,
                                                    G4ThreeVector(MOPIFirstAluminumLayerXPosition,
                                                                  MOPIFirstAluminumLayerYPosition ,
                                                                  MOPIFirstAluminumLayerZPosition),
                                                    "MOPIFirstAluminumLayer",
                                                    logicMOPIFirstAluminumLayer, physiMOPIMotherVolume, false, 0);
    
    // First Air GAP
    solidMOPIFirstAirGap = new G4Box("MOPIFirstAirGap",
                                     MOPIFirstAirGapXSize/2,
                                     MOPIFirstAirGapYSize/2,
                                     MOPIFirstAirGapZSize/2);
    
    logicMOPIFirstAirGap = new G4LogicalVolume(solidMOPIFirstAirGap,
                                               MOPIFirstAirGapMaterial,
                                               "MOPIFirstAirgap");
    
    physiMOPIFirstAirGap = new G4PVPlacement(0,
                                             G4ThreeVector(MOPIFirstAirGapXPosition,
                                                           MOPIFirstAirGapYPosition ,
                                                           MOPIFirstAirGapZPosition),
                                             "MOPIFirstAirGap",
                                             logicMOPIFirstAirGap, physiMOPIMotherVolume, false, 0);
    
    
    // The Cathode
    solidMOPICathode = new G4Box("MOPICathode",
                                 MOPICathodeXSize/2,
                                 MOPICathodeYSize/2,
                                 MOPICathodeZSize/2);
    
    logicMOPICathode = new G4LogicalVolume(solidMOPICathode,
                                           MOPICathodeMaterial,
                                           "MOPICathode");
    
    physiMOPICathode = new G4PVPlacement(0,
                                         G4ThreeVector(MOPICathodeXPosition,
                                                       MOPICathodeYPosition ,
                                                       MOPICathodeZPosition),
                                         "MOPICathode",
                                         logicMOPICathode,
                                         physiMOPIMotherVolume, false, 0);
    
    // Second Air GAP
    solidMOPISecondAirGap = new G4Box("MOPISecondAirGap",
                                      MOPISecondAirGapXSize/2,
                                      MOPISecondAirGapYSize/2,
                                      MOPISecondAirGapZSize/2);
    
    logicMOPISecondAirGap = new G4LogicalVolume(solidMOPISecondAirGap,
                                                MOPISecondAirGapMaterial,
                                                "MOPISecondAirgap");
    
    physiMOPISecondAirGap = new G4PVPlacement(0,
                                              G4ThreeVector(MOPISecondAirGapXPosition,
                                                            MOPISecondAirGapYPosition ,
                                                            MOPISecondAirGapZPosition),
                                              "MOPISecondAirGap",
                                              logicMOPISecondAirGap, physiMOPIMotherVolume, false, 0);
    
    // Second Aluminum layer
    solidMOPISecondAluminumLayer = new G4Box("MOPISecondAluminumLayer",
                                             MOPISecondAluminumLayerXSize/2,
                                             MOPISecondAluminumLayerYSize/2 ,
                                             MOPISecondAluminumLayerZSize/2);
    
    logicMOPISecondAluminumLayer = new G4LogicalVolume(solidMOPISecondAluminumLayer,
                                                       MOPISecondAluminumLayerMaterial,
                                                       "MOPISecondAluminumLayer");
    
    physiMOPISecondAluminumLayer = new G4PVPlacement(0,
                                                     G4ThreeVector(MOPISecondAluminumLayerXPosition,
                                                                   MOPISecondAluminumLayerYPosition ,
                                                                   MOPISecondAluminumLayerZPosition),
                                                     "MOPISecondAluminumLayer",
                                                     logicMOPISecondAluminumLayer,
                                                     physiMOPIMotherVolume,
                                                     false,
                                                     0);
    
    // Second Kapton layer
    solidMOPISecondKaptonLayer = new G4Box("MOPISecondKaptonLayer",
                                           MOPISecondKaptonLayerXSize/2,
                                           MOPISecondKaptonLayerYSize/2 ,
                                           MOPISecondKaptonLayerZSize/2);
    
    logicMOPISecondKaptonLayer = new G4LogicalVolume(solidMOPISecondKaptonLayer,
                                                     MOPIFirstKaptonLayerMaterial,
                                                     "MOPISecondKaptonLayer");
    
    physiMOPISecondKaptonLayer = new G4PVPlacement(0,
                                                   G4ThreeVector(MOPISecondKaptonLayerXPosition,
                                                                 MOPISecondKaptonLayerYPosition ,
                                                                 MOPISecondKaptonLayerZPosition),
                                                   "MOPISecondKaptonLayer",
                                                   logicMOPISecondKaptonLayer,
                                                   physiMOPIMotherVolume,
                                                   false,
                                                   0);
    
    logicMOPIFirstAirGap -> SetVisAttributes(darkGreen);
    logicMOPISecondAirGap -> SetVisAttributes(darkGreen);
    
    
    }
/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::HadrontherapyBeamNozzle()
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
    const G4double nozzleSupportXSize = 29.5 *mm;
    const G4double nozzleSupportYSize = 180. *mm;
    const G4double nozzleSupportZSize = 180. *mm;
    
    const G4double nozzleSupportXPosition = -397.50 *mm;
    
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
    
    
    
    //------------------------------------//
    // HOLE IN THE SUPPORT                //
    //------------------------------------//
    const G4double innerRadiusHoleNozzleSupport = 0.*mm;
    const G4double outerRadiusHoleNozzleSupport = 21.5*mm;
    const G4double hightHoleNozzleSupport = 29.5 *mm;
    const G4double startAngleHoleNozzleSupport = 0.*deg;
    const G4double spanningAngleHoleNozzleSupport = 360.*deg;
    
    G4Tubs* solidHoleNozzleSupport = new G4Tubs("HoleNozzleSupport",
                                                innerRadiusHoleNozzleSupport,
                                                outerRadiusHoleNozzleSupport,
                                                hightHoleNozzleSupport,
                                                startAngleHoleNozzleSupport,
                                                spanningAngleHoleNozzleSupport);
    
    G4LogicalVolume* logicHoleNozzleSupport = new G4LogicalVolume(solidHoleNozzleSupport,
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
    
    logicHoleNozzleSupport -> SetVisAttributes(darkOrange3);
    
    // ---------------------------------//
    //     BRASS TUBE 1 (phantom side)    //
    // ---------------------------------//
    const G4double innerRadiusBrassTube= 18.*mm;
    const G4double outerRadiusBrassTube = 21.5 *mm;
    const G4double hightBrassTube = 140.5*mm;
    const G4double startAngleBrassTube = 0.*deg;
    const G4double spanningAngleBrassTube = 360.*deg;
    
    const G4double brassTubeXPosition = -227.5 *mm;
    
    G4Tubs* solidBrassTube = new G4Tubs("BrassTube",
                                        innerRadiusBrassTube,
                                        outerRadiusBrassTube,
                                        hightBrassTube,
                                        startAngleBrassTube,
                                        spanningAngleBrassTube);
    
    G4LogicalVolume* logicBrassTube = new G4LogicalVolume(solidBrassTube,
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
    
    // ----------------------------------------------//
    //     BRASS TUBE 2 (inside the PMMA support)    //
    // ----------------------------------------------//
    const G4double innerRadiusBrassTube2= 18.*mm;
    const G4double outerRadiusBrassTube2 = 21.5 *mm;
    const G4double hightBrassTube2 = 29.5*mm;
    const G4double startAngleBrassTube2 = 0.*deg;
    const G4double spanningAngleBrassTube2 = 360.*deg;
    
    
    G4Tubs* solidBrassTube2 = new G4Tubs("BrassTube2",
                                         innerRadiusBrassTube2,
                                         outerRadiusBrassTube2,
                                         hightBrassTube2,
                                         startAngleBrassTube2,
                                         spanningAngleBrassTube2);
    
    G4LogicalVolume* logicBrassTube2 = new G4LogicalVolume(solidBrassTube2,
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
    
    
    // --------------------------------------//
    //     BRASS TUBE 3 (beam line side)    //
    // -------------------------------------//
    const G4double innerRadiusBrassTube3= 18.*mm;
    const G4double outerRadiusBrassTube3 = 21.5 *mm;
    const G4double hightBrassTube3 = 10.0 *mm;
    const G4double startAngleBrassTube3 = 0.*deg;
    const G4double spanningAngleBrassTube3 = 360.*deg;
    
    const G4double brassTube3XPosition = -437 *mm;
    
    G4Tubs* solidBrassTube3 = new G4Tubs("BrassTube3",
                                         innerRadiusBrassTube3,
                                         outerRadiusBrassTube3,
                                         hightBrassTube3,
                                         startAngleBrassTube3,
                                         spanningAngleBrassTube3);
    
    G4LogicalVolume* logicBrassTube3 = new G4LogicalVolume(solidBrassTube3,
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
}

/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::HadrontherapyBeamFinalCollimator()
{
    // -----------------------//
    //     FINAL COLLIMATOR   //
    //------------------------//
    const G4double outerRadiusFinalCollimator = 21.5*mm;
    const G4double hightFinalCollimator = 3.5*mm;
    const G4double startAngleFinalCollimator = 0.*deg;
    const G4double spanningAngleFinalCollimator = 360.*deg;
    const G4double finalCollimatorXPosition = -83.5 *mm;
    
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
/////////////////////////// MESSENGER ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::SetRangeShifterXPosition(G4double value)
{
    physiRangeShifterBox -> SetTranslation(G4ThreeVector(value, 0., 0.));
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout << "The Range Shifter is translated to"<< value/mm <<"mm along the X axis" <<G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::SetRangeShifterXSize(G4double value)
{
    solidRangeShifterBox -> SetXHalfLength(value) ;
    G4cout << "RangeShifter size X (mm): "<< ((solidRangeShifterBox -> GetXHalfLength())*2.)/mm
    << G4endl;
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::SetFirstScatteringFoilXSize(G4double value)
{
    firstScatteringFoil -> SetXHalfLength(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout <<"The X size of the first scattering foil is (mm):"<<
    ((firstScatteringFoil -> GetXHalfLength())*2.)/mm
    << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::SetSecondScatteringFoilXSize(G4double value)
{
    secondScatteringFoil -> SetXHalfLength(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout <<"The X size of the second scattering foil is (mm):"<<
    ((secondScatteringFoil -> GetXHalfLength())*2.)/mm
    << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::SetOuterRadiusStopper(G4double value)
{
    solidStopper -> SetOuterRadius(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout << "OuterRadius od the Stopper is (mm):"
    << solidStopper -> GetOuterRadius()/mm
    << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::SetInnerRadiusFinalCollimator(G4double value)
{
    solidFinalCollimator -> SetInnerRadius(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout<<"Inner Radius of the final collimator is (mm):"
	<< solidFinalCollimator -> GetInnerRadius()/mm
	<< G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::SetRSMaterial(G4String materialChoice)
{
    if (G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice, false) )
    {
        if (pttoMaterial)
        {
            rangeShifterMaterial  = pttoMaterial;
            logicRangeShifterBox -> SetMaterial(pttoMaterial);
            G4cout << "The material of the Range Shifter has been changed to " << materialChoice << G4endl;
        }
    }
    else
    {
        G4cout << "WARNING: material \"" << materialChoice << "\" doesn't exist in NIST elements/materials"
	    " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
        G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl;
    }
}

/////////////////////////////////////////////////////////////////////////////
void PassiveProtonBeamLine::SetModulatorAngle(G4double value)
{
    modulator -> SetModulatorAngle(value);
    //G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
/////////////////////////////////////////////////////////////////////////////
