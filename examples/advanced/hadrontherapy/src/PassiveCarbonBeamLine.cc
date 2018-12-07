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
// Simulation of the "Zero degree" experimental beamline of INFN-LNS (Catania, Italy).

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

//G4bool PassiveCarbonBeamLine::doCalculation = false;
/////////////////////////////////////////////////////////////////////////////
PassiveCarbonBeamLine::PassiveCarbonBeamLine():
physicalTreatmentRoom(0),hadrontherapyDetectorConstruction(0),
physiBeamLineSupport(0), physiBeamLineCover(0), physiBeamLineCover2(0),
physiKaptonWindow(0),PhysiRippleFilter(0),PhysiRippleFilterBase(0),PhysiRippleFilterTrd(0),
physiFirstMonitorLayer1(0), physiFirstMonitorLayer2(0),
physiFirstMonitorLayer3(0), physiFirstMonitorLayer4(0),
physiNozzleSupport(0), physiHoleNozzleSupport(0)
{
    
    // Messenger to change parameters of the passiveCarbonBeamLine geometry
    PassiveCarbonMessenger = new PassiveCarbonBeamLineMessenger(this);
    
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
    delete PassiveCarbonMessenger;
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
    rippleFilterMaterial = airNist;
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
    HadrontherapyRippleFilter();
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
                                          physiVacuumZone, false,	0);
    
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
void PassiveCarbonBeamLine::HadrontherapyPMMACollimator()
{

    // ----------------------//
    //   PMMA COLLIMATOR     //
    // ----------------------//
    PMMACollimatorSupportXSize = 25.0 *mm;
    PMMACollimatorSupportYSize = 200. *mm;
    PMMACollimatorSupportZSize = 200. *mm;
    
    PMMACollimatorXPosition = -1082.00 *mm;
    
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
    monitor1XPosition = -1059.0 *mm;
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


/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::SetInnerRadiusFinalCollimator(G4double value)
{
    solidFinalCollimator -> SetInnerRadius(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout<<"Inner Radius of the final collimator is (mm):"
    << solidFinalCollimator -> GetInnerRadius()/mm
    << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void PassiveCarbonBeamLine::SetRippleFilterMaterial(G4String materialChoice)
{
    if (G4Material* RFMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice, false) )
    {
        if (RFMaterial)
        {
            rippleFilterMaterial  = RFMaterial;
            LogicRippleFilter -> SetMaterial(RFMaterial);
            LogicRippleFilterBase -> SetMaterial(RFMaterial);
            LogicRippleFilterTrd -> SetMaterial(RFMaterial);
            G4cout << "The material of the Ripple Filter has been changed to " << materialChoice << G4endl;
        }
    }
    else
    {
        G4cout << "WARNING: material \"" << materialChoice << "\" doesn't exist in NIST elements/materials"
        " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
        G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl;
    }
}



