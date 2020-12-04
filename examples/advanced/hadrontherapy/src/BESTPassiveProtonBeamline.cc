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
#include "BESTPassiveProtonBeamLine.hh"
#include "BESTPassiveProtonBeamLineMessenger.hh"


//G4bool PassiveProtonBeamLine::doCalculation = false;
/////////////////////////////////////////////////////////////////////////////
BESTPassiveProtonBeamLine::BESTPassiveProtonBeamLine():
modulator(0), physicalTreatmentRoom(0),hadrontherapyDetectorConstruction(0),
physiBeamLineSupport(0), physiBeamLineCover(0), physiBeamLineCover2(0),
BESTfirstScatteringFoil(0), physiBESTFirstScatteringFoil(0), physiBESTKaptonWindow(0),
solidBESTStopper(0), physiBESTStopper(0), BESTsecondScatteringFoil(0), physiBESTSecondScatteringFoil(0),
physiBESTFirstCollimator(0), solidBESTRangeShifterBox(0), logicBESTRangeShifterBox(0),
physiBESTRangeShifterBox(0), physiBESTSecondCollimator(0), physiBESTFirstCollimatorModulatorBox(0),
physiBESTHoleFirstCollimatorModulatorBox(0), physiBESTSecondCollimatorModulatorBox(0),
physiBESTHoleSecondCollimatorModulatorBox(0),
chamberPhys(0),innerchamberPhys(0),enterWindowPhys(0),enterElectrodePhys(0),kaptonLayerPhys1(0),copperLayerPhys1(0),nickelLayerPhys1(0),fFirstCavityPhys(0),centralElectrode1Phys(0), centralWindowPhys(0), centralElectrode2Phys(0),fSecondCavityPhys(0),exitElectrodePhys(0),kaptonLayerPhys2(0),copperLayerPhys2(0),nickelLayerPhys2(0),exitWindowPhys(0),physiNozzleSupport(0), physiBrassTube(0), solidFinalCollimator(0), physiFinalCollimator(0)
{
    // Messenger to change parameters of the passiveProtonBeamLine geometry
    passiveMessenger = new BESTPassiveProtonBeamLineMessenger(this);
    
    //***************************** PW ***************************************
    static G4String ROGeometryName = "DetectorROGeometry";
    RO = new HadrontherapyDetectorROGeometry(ROGeometryName);
    
    G4cout << "Going to register Parallel world...";
    RegisterParallelWorld(RO);
    G4cout << "... done" << G4endl;

}
/////////////////////////////////////////////////////////////////////////////
BESTPassiveProtonBeamLine::~BESTPassiveProtonBeamLine()
{
   // delete passiveMessenger;
    delete hadrontherapyDetectorConstruction;

}

/////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* BESTPassiveProtonBeamLine::Construct()
{
    // Sets default geometry and materials
    SetDefaultDimensions();
    
    // Construct the whole BEST Passive Beam Line
    ConstructBESTPassiveProtonBeamLine();
    
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
void BESTPassiveProtonBeamLine::SetDefaultDimensions()
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
    G4double defaultBESTVacuumZoneXSize = 100.0 *mm;
    BESTvacuumZoneXSize = defaultBESTVacuumZoneXSize;
    
    G4double defaultBESTVacuumZoneYSize = 52.5 *mm;
    BESTvacuumZoneYSize = defaultBESTVacuumZoneYSize;
    
    G4double defaultBESTVacuumZoneZSize = 52.5 *mm;
    BESTvacuumZoneZSize = defaultBESTVacuumZoneZSize;
    
    G4double defaultBESTVacuumZoneXPosition = -3010.0 *mm;
    BESTvacuumZoneXPosition = defaultBESTVacuumZoneXPosition;
    
    // FIRST SCATTERING FOIL: a thin foil performing a first scattering
    // of the original beam
    G4double defaultBESTFirstScatteringFoilXSize = 0.0075 *mm;
    BESTfirstScatteringFoilXSize = defaultBESTFirstScatteringFoilXSize;
    
    G4double defaultBESTFirstScatteringFoilYSize = 52.5   *mm;
   BESTfirstScatteringFoilYSize = defaultBESTFirstScatteringFoilYSize;
    
    G4double defaultBESTFirstScatteringFoilZSize = 52.5   *mm;
    BESTfirstScatteringFoilZSize = defaultBESTFirstScatteringFoilZSize;
    
    G4double defaultBESTFirstScatteringFoilXPosition = 0.0 *mm;
    BESTfirstScatteringFoilXPosition = defaultBESTFirstScatteringFoilXPosition;
    
    // KAPTON WINDOW: it prmits the passage of the beam from vacuum to air
    G4double defaultBESTKaptonWindowXSize = 0.010*mm;
    BESTkaptonWindowXSize = defaultBESTKaptonWindowXSize;
    
    G4double defaultBESTKaptonWindowYSize = 5.25*cm;
    BESTkaptonWindowYSize = defaultBESTKaptonWindowYSize;
    
    G4double defaultBESTKaptonWindowZSize = 5.25*cm;
    BESTkaptonWindowZSize = defaultBESTKaptonWindowZSize;
    
    G4double defaultBESTKaptonWindowXPosition = 100.0*mm - defaultBESTKaptonWindowXSize;
    BESTkaptonWindowXPosition = defaultBESTKaptonWindowXPosition;
    
    // STOPPER: is a small cylinder able to stop the central component
    // of the beam (having a gaussian shape). It is connected to the SECON SCATTERING FOIL
    // and represent the second element of the scattering system
    G4double defaultBESTInnerRadiusStopper = 0.*cm;
    BESTinnerRadiusStopper = defaultBESTInnerRadiusStopper;
    
    G4double defaultBESTHeightStopper = 4.5*mm;
    BESTheightStopper = defaultBESTHeightStopper;
    
    G4double defaultBESTStartAngleStopper = 0.*deg;
    BESTstartAngleStopper = defaultBESTStartAngleStopper;
    
    G4double defaultBESTSpanningAngleStopper = 360.*deg;
    BESTspanningAngleStopper = defaultBESTSpanningAngleStopper;
    
    G4double defaultBESTStopperXPosition = -2705.0 *mm;
    BESTstopperXPosition = defaultBESTStopperXPosition;
    
    G4double defaultBESTStopperYPosition = 0.*m;
    BESTstopperYPosition = defaultBESTStopperYPosition;
    
    G4double defaultBESTStopperZPosition = 0.*m;
    BESTstopperZPosition = defaultBESTStopperZPosition;
    
    G4double defaultBESTOuterRadiusStopper = 3 *mm;
    BESTouterRadiusStopper = defaultBESTOuterRadiusStopper;
    
    // SECOND SCATTERING FOIL: it is another thin foil and provides the
    // final diffusion of the beam. It represents the third element of the scattering
    // system;
    G4double defaultBESTSecondScatteringFoilXSize = 0.0125 *mm;
   BESTsecondScatteringFoilXSize = defaultBESTSecondScatteringFoilXSize;
    
    G4double defaultBESTSecondScatteringFoilYSize = 52.5   *mm;
    BESTsecondScatteringFoilYSize = defaultBESTSecondScatteringFoilYSize;
    
    G4double defaultBESTSecondScatteringFoilZSize = 52.5   *mm;
    BESTsecondScatteringFoilZSize = defaultBESTSecondScatteringFoilZSize;
    
    G4double defaultBESTSecondScatteringFoilXPosition = defaultBESTStopperXPosition + defaultBESTHeightStopper + defaultBESTSecondScatteringFoilXSize;
    BESTsecondScatteringFoilXPosition = defaultBESTSecondScatteringFoilXPosition;
    
    G4double defaultBESTSecondScatteringFoilYPosition =  0 *mm;
    BESTsecondScatteringFoilYPosition = defaultBESTSecondScatteringFoilYPosition;
    
    G4double defaultBESTSecondScatteringFoilZPosition =  0 *mm;
    BESTsecondScatteringFoilZPosition = defaultBESTSecondScatteringFoilZPosition;
    
    // RANGE SHIFTER: is a slab of PMMA acting as energy degreader of
    // primary beam
    
    //Default material of the range shifter
    
    G4double defaultBESTRangeShifterXSize = 5. *mm;
    BESTrangeShifterXSize = defaultBESTRangeShifterXSize;
    
    G4double defaultBESTRangeShifterYSize = 176. *mm;
    BESTrangeShifterYSize = defaultBESTRangeShifterYSize;
    
    G4double defaultBESTRangeShifterZSize = 176. *mm;
    BESTrangeShifterZSize = defaultBESTRangeShifterZSize;
    
    G4double defaultBESTRangeShifterXPosition = -2393.0 *mm;
    BESTrangeShifterXPosition = defaultBESTRangeShifterXPosition;
    
    G4double defaultBESTRangeShifterYPosition = 0. *mm;
    BESTrangeShifterYPosition = defaultBESTRangeShifterYPosition;
    
    G4double defaultBESTRangeShifterZPosition = 0. *mm;
    BESTrangeShifterZPosition = defaultBESTRangeShifterZPosition;
    
   
    
    // FINAL COLLIMATOR: is the collimator giving the final transversal shape
    // of the beam
    G4double defaultBESTinnerRadiusFinalCollimator = 7.5 *mm;
    BESTinnerRadiusFinalCollimator = defaultBESTinnerRadiusFinalCollimator;
    
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
    
G4Material* nichelNistAsMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ni", isotopes);

    


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
  CopperLayerMaterial=copperNistAsMaterial;
    NichelLayerMaterial=nichelNistAsMaterial;
    KaptonLayerMaterial=kaptonNist;
    WindowMaterial=mylarNist;
    CentralWindowMaterial=mylarNist;
    wallMaterial=aluminumNist;
    ElectrodeMaterial=aluminumNist;
    CavityMaterial=airNist;
      
    
    // material of the final nozzle
    nozzleSupportMaterial = PMMANist;
    brassTubeMaterial = brassTube2Material = brassTube3Material = brass;
    holeNozzleSupportMaterial = airNist;
    
    // Material of the final collimator
    finalCollimatorMaterial = brass;
}

/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::ConstructBESTPassiveProtonBeamLine()
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
    //logicTreatmentRoom -> SetVisAttributes(G4VisAttributes::Invisible);
    
    // Components of the BEST Passive Proton Beam Line
   BESTBeamLineSupport();
    BESTBeamScatteringFoils();
    BESTRangeShifter();
    BESTBeamCollimators();
    BESTBeamMonitoring();
    BESTBeamNozzle();
    BESTBeamFinalCollimator();
    
    // The following lines construc a typical modulator wheel inside the Passive Beam line.
    // Please remember to set the nodulator material (default is air, i.e. no modulator!)
    // in the HadrontherapyModulator.cc file
    modulator = new HadrontherapyModulator();
    modulator -> BuildModulator(physicalTreatmentRoom);
}

/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::BESTBeamLineSupport()
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
void BESTPassiveProtonBeamLine::BESTBeamScatteringFoils()
{
   // ------------//
    // VACUUM PIPE //
    //-------------//
    //
    // First track of the beam line is inside vacuum;
    // The PIPE contains the FIRST SCATTERING FOIL and the KAPTON WINDOW
    G4Box* BESTvacuumZone = new G4Box("VacuumZone", BESTvacuumZoneXSize, BESTvacuumZoneYSize, BESTvacuumZoneZSize);
    G4LogicalVolume* logicBESTVacuumZone = new G4LogicalVolume(BESTvacuumZone, vacuumZoneMaterial, "VacuumZone");
    G4VPhysicalVolume* physiBESTVacuumZone = new G4PVPlacement(0, G4ThreeVector(BESTvacuumZoneXPosition, 0., 0.),
                                                           "VacuumZone",logicBESTVacuumZone, physicalTreatmentRoom, false, 0);
    // --------------------------//
    // THE FIRST SCATTERING FOIL //
    // --------------------------//
    // A thin foil performing a first scattering
    // of the original beam
    BESTfirstScatteringFoil = new G4Box("FirstScatteringFoil",
                                    BESTfirstScatteringFoilXSize,
                                    BESTfirstScatteringFoilYSize,
                                    BESTfirstScatteringFoilZSize);
    
    G4LogicalVolume* logicBESTFirstScatteringFoil = new G4LogicalVolume(BESTfirstScatteringFoil,
                                                                    firstScatteringFoilMaterial,
                                                                    "FirstScatteringFoil");
    
    physiBESTFirstScatteringFoil = new G4PVPlacement(0, G4ThreeVector(BESTfirstScatteringFoilXPosition, 0.,0.),
                                                 "FirstScatteringFoil", logicBESTFirstScatteringFoil, physiBESTVacuumZone,
                                                 false, 0);
    
    logicBESTFirstScatteringFoil -> SetVisAttributes(skyBlue);
    // -------------------//
    // THE KAPTON WINDOWS //
    //--------------------//
    //It prmits the passage of the beam from vacuum to air
 
    G4Box* solidBESTKaptonWindow = new G4Box("BESTKaptonWindow",
                                         BESTkaptonWindowXSize,
                                         BESTkaptonWindowYSize,
                                         BESTkaptonWindowZSize);
    
    G4LogicalVolume* logicBESTKaptonWindow = new G4LogicalVolume(solidBESTKaptonWindow,
                                                             kaptonWindowMaterial,
                                                             "BESTKaptonWindow");
    
    physiBESTKaptonWindow = new G4PVPlacement(0, G4ThreeVector(BESTkaptonWindowXPosition, 0., 0.),
                                          "BESTKaptonWindow", logicBESTKaptonWindow,
                                          physiBESTVacuumZone, false,	0);
    
    logicBESTKaptonWindow -> SetVisAttributes(darkOrange3);
    
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
    
    solidBESTStopper = new G4Tubs("Stopper",
                              BESTinnerRadiusStopper,
                              BESTouterRadiusStopper,
                              BESTheightStopper,
                              BESTstartAngleStopper,
                              BESTspanningAngleStopper);
    
    logicBESTStopper = new G4LogicalVolume(solidBESTStopper,
                                       stopperMaterial,
                                       "Stopper",
                                       0, 0, 0);
    
    physiBESTStopper = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(BESTstopperXPosition,
                                                                     BESTstopperYPosition,
                                                                     BESTstopperZPosition)),
                                     "BESTStopper",
                                     logicBESTStopper,
                                     physicalTreatmentRoom,
                                     false,
                                     0);
    
    logicBESTStopper -> SetVisAttributes(red);
    
    // ---------------------------//
    // THE SECOND SCATTERING FOIL //
    // ---------------------------//
    // It is another thin foil and provides the
    // final diffusion of the beam. It represents the third element of the scattering
    // system;
    
    BESTsecondScatteringFoil = new G4Box("SecondScatteringFoil",
                                     BESTsecondScatteringFoilXSize,
                                     BESTsecondScatteringFoilYSize,
                                     BESTsecondScatteringFoilZSize);
    
    G4LogicalVolume* logicBESTSecondScatteringFoil = new G4LogicalVolume(BESTsecondScatteringFoil,
                                                                     secondScatteringFoilMaterial,
                                                                     "SecondScatteringFoil");
    
    physiBESTSecondScatteringFoil = new G4PVPlacement(0, G4ThreeVector(BESTsecondScatteringFoilXPosition,
                                                                   BESTsecondScatteringFoilYPosition,
                                                                   BESTsecondScatteringFoilZPosition),
                                                  "SeconScatteringFoil",
                                                  logicBESTSecondScatteringFoil,
                                                  physicalTreatmentRoom,
                                                  false,
                                                  0);
    
    logicBESTSecondScatteringFoil -> SetVisAttributes(skyBlue);
 
 
}
/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::BESTRangeShifter()
{
    // ---------------------------- //
    //         THE RANGE SHIFTER    //
    // -----------------------------//
    // It is a slab of PMMA acting as energy degreader of
    // primary beam
 
    
    solidBESTRangeShifterBox = new G4Box("RangeShifterBox",
                                     BESTrangeShifterXSize,
                                     BESTrangeShifterYSize,
                                     BESTrangeShifterZSize);
    
    logicBESTRangeShifterBox = new G4LogicalVolume(solidBESTRangeShifterBox,
                                               rangeShifterMaterial,
                                               "RangeShifterBox");
    physiBESTRangeShifterBox = new G4PVPlacement(0,
                                             G4ThreeVector(BESTrangeShifterXPosition, 0., 0.),
                                             "RangeShifterBox",
                                             logicBESTRangeShifterBox,
                                             physicalTreatmentRoom,
                                             false,
                                             0);
    
    
    logicBESTRangeShifterBox -> SetVisAttributes(yellow);
     
    
}
/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::BESTBeamCollimators()
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
    
    
    G4Box* solidBESTFirstCollimator = new G4Box("FirstCollimator",
                                            firstCollimatorXSize,
                                            firstCollimatorYSize,
                                            firstCollimatorZSize);
    
    G4LogicalVolume* logicBESTFirstCollimator = new G4LogicalVolume(solidBESTFirstCollimator,
                                                                firstCollimatorMaterial,
                                                                "FirstCollimator");
    
    physiBESTFirstCollimator = new G4PVPlacement(0, G4ThreeVector(firstCollimatorXPosition,
                                                              firstCollimatorYPosition,
                                                              firstCollimatorZPosition),
                                             "FirstCollimator",
                                             logicBESTFirstCollimator,
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
    
    G4Tubs* solidBESTHoleFirstCollimator = new G4Tubs("HoleFirstCollimator",
                                                  innerRadiusHoleFirstCollimator,
                                                  outerRadiusHoleFirstCollimator,
                                                  hightHoleFirstCollimator,
                                                  startAngleHoleFirstCollimator,
                                                  spanningAngleHoleFirstCollimator);
    
    G4LogicalVolume* logicBESTHoleFirstCollimator = new G4LogicalVolume(solidBESTHoleFirstCollimator,
                                                                    holeFirstCollimatorMaterial,
                                                                    "HoleFirstCollimator",
                                                                    0, 0, 0);
    G4double phi = 90. *deg;
    // Matrix definition for a 90 deg rotation. Also used for other volumes
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    physiBESTHoleFirstCollimator = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                                 "HoleFirstCollimator",
                                                 logicBESTHoleFirstCollimator,
                                                 physiBESTFirstCollimator,
                                                 false,
                                                 0);
    // ------------------//
    // SECOND COLLIMATOR //
    //-------------------//
    // It is a slab of PMMA with an hole in its center
    const G4double secondCollimatorXPosition = -1900.00*mm;
    const G4double secondCollimatorYPosition =  0*mm;
    const G4double secondCollimatorZPosition =  0*mm;
    
    physiBESTSecondCollimator = new G4PVPlacement(0, G4ThreeVector(secondCollimatorXPosition,
                                                               secondCollimatorYPosition,
                                                               secondCollimatorZPosition),
                                              "SecondCollimator",
                                              logicBESTFirstCollimator,
                                              physicalTreatmentRoom,
                                              false,
                                              0);
    
    // ------------------------------//
    // Hole of the second collimator //
    // ------------------------------//
    physiBESTHoleSecondCollimator = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                                  "HoleSecondCollimator",
                                                  logicBESTHoleFirstCollimator,
                                                  physiBESTSecondCollimator,
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
    
   G4Box* solidBESTFirstCollimatorModulatorBox = new G4Box("FirstCollimatorModulatorBox",
                                                        firstCollimatorModulatorXSize,
                                                        firstCollimatorModulatorYSize,
                                                        firstCollimatorModulatorZSize);
    
    G4LogicalVolume* logicBESTFirstCollimatorModulatorBox = new G4LogicalVolume(solidBESTFirstCollimatorModulatorBox,
                                                                            modulatorBoxMaterial,
                                                                            "FirstCollimatorModulatorBox");
    
    physiBESTFirstCollimatorModulatorBox = new G4PVPlacement(0, G4ThreeVector(firstCollimatorModulatorXPosition,
                                                                          firstCollimatorModulatorYPosition,
                                                                          firstCollimatorModulatorZPosition),
                                                         "FirstCollimatorModulatorBox",
                                                         logicBESTFirstCollimatorModulatorBox,
                                                         physicalTreatmentRoom, false, 0);
    
    // ----------------------------------------------------//
    //   Hole of the first collimator of the modulator box //
    // ----------------------------------------------------//
    const G4double innerRadiusHoleFirstCollimatorModulatorBox = 0.*mm;
    const G4double outerRadiusHoleFirstCollimatorModulatorBox = 31.*mm;
    const G4double hightHoleFirstCollimatorModulatorBox = 10.*mm;
    const G4double startAngleHoleFirstCollimatorModulatorBox = 0.*deg;
    const G4double spanningAngleHoleFirstCollimatorModulatorBox = 360.*deg;
    
    G4Tubs* solidBESTHoleFirstCollimatorModulatorBox  = new G4Tubs("HoleFirstCollimatorModulatorBox",
                                                               innerRadiusHoleFirstCollimatorModulatorBox,
                                                               outerRadiusHoleFirstCollimatorModulatorBox,
                                                               hightHoleFirstCollimatorModulatorBox ,
                                                               startAngleHoleFirstCollimatorModulatorBox,
                                                               spanningAngleHoleFirstCollimatorModulatorBox);
    
    G4LogicalVolume* logicBESTHoleFirstCollimatorModulatorBox = new G4LogicalVolume(solidBESTHoleFirstCollimatorModulatorBox,
                                                                                holeModulatorBoxMaterial,
                                                                                "HoleFirstCollimatorModulatorBox",
                                                                                0, 0, 0);
    
    physiBESTHoleFirstCollimatorModulatorBox = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                                             "HoleFirstCollimatorModulatorBox",
                                                             logicBESTHoleFirstCollimatorModulatorBox,
                                                             physiBESTFirstCollimatorModulatorBox, false, 0);
    
    // --------------------------------------------------//
    //       SECOND SIDE OF THE MODULATOR BOX            //
    // --------------------------------------------------//
    const G4double secondCollimatorModulatorXSize = 10.*mm;
    const G4double secondCollimatorModulatorYSize = 200.*mm;
    const G4double secondCollimatorModulatorZSize = 200.*mm;
    
    const G4double secondCollimatorModulatorXPosition = -1953.00 *mm;
    
    const G4double secondCollimatorModulatorYPosition = 0.*mm;
    const G4double secondCollimatorModulatorZPosition = 0.*mm;
    
    G4Box* solidBESTSecondCollimatorModulatorBox = new G4Box("SecondCollimatorModulatorBox",
                                                         secondCollimatorModulatorXSize,
                                                         secondCollimatorModulatorYSize,
                                                         secondCollimatorModulatorZSize);
    
    G4LogicalVolume* logicBESTSecondCollimatorModulatorBox = new G4LogicalVolume(solidBESTSecondCollimatorModulatorBox,
                                                                             modulatorBoxMaterial,
                                                                             "SecondCollimatorModulatorBox");
    
    physiBESTSecondCollimatorModulatorBox = new G4PVPlacement(0, G4ThreeVector(secondCollimatorModulatorXPosition,
                                                                           secondCollimatorModulatorYPosition,
                                                                           secondCollimatorModulatorZPosition),
                                                          "SecondCollimatorModulatorBox",
                                                          logicBESTSecondCollimatorModulatorBox,
                                                          physicalTreatmentRoom, false, 0);
    
    // ----------------------------------------------//
    //   Hole of the second collimator modulator box //
    // ----------------------------------------------//
    const G4double innerRadiusHoleSecondCollimatorModulatorBox = 0.*mm;
    const G4double outerRadiusHoleSecondCollimatorModulatorBox = 31.*mm;
    const G4double hightHoleSecondCollimatorModulatorBox = 10.*mm;
    const G4double startAngleHoleSecondCollimatorModulatorBox = 0.*deg;
    const G4double spanningAngleHoleSecondCollimatorModulatorBox = 360.*deg;
    
    G4Tubs* solidBESTHoleSecondCollimatorModulatorBox  = new G4Tubs("HoleSecondCollimatorModulatorBox",
                                                                innerRadiusHoleSecondCollimatorModulatorBox,
                                                                outerRadiusHoleSecondCollimatorModulatorBox,
                                                                hightHoleSecondCollimatorModulatorBox ,
                                                                startAngleHoleSecondCollimatorModulatorBox,
                                                                spanningAngleHoleSecondCollimatorModulatorBox);
    
    G4LogicalVolume* logicHBESToleSecondCollimatorModulatorBox = new G4LogicalVolume(solidBESTHoleSecondCollimatorModulatorBox,
                                                                                 holeModulatorBoxMaterial,
                                                                                 "HoleSecondCollimatorModulatorBox",
                                                                                 0, 0, 0);
    
    physiBESTHoleSecondCollimatorModulatorBox = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector()),
                                                              "HoleSecondCollimatorModulatorBox",
                                                             logicHBESToleSecondCollimatorModulatorBox,
                                                              physiBESTSecondCollimatorModulatorBox, false, 0);
    
    logicBESTFirstCollimator -> SetVisAttributes(yellow);
    logicBESTFirstCollimatorModulatorBox -> SetVisAttributes(blue);
    logicBESTSecondCollimatorModulatorBox -> SetVisAttributes(blue);


  
  }

/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::BESTBeamMonitoring()
{

//////// Double-GAP Monitor Chamber ////////////
  
    G4double fBoxOuterThickness { 100 * mm };
    G4double fBoxOuterHeight { 42 * cm };
    G4double fBoxOuterWidth { 42 * cm };
    G4double fEnterWindowThickness { 25 * um };
    G4double fCentralWindowThickness { 20 * um };
    G4double fExitWindowThickness { 25 * um };
    G4double fBoxInnerThickness { fBoxOuterThickness - fEnterWindowThickness - fExitWindowThickness };
    G4double fBoxInnerHeight { 41.95 * cm };
    G4double fBoxInnerWidth { 41.95 * cm };
    G4double fInnerHeight { 16 * cm };
    G4double fInnerWidth { 16 * cm };
    G4double fCentralElectrodeThickness { 2 * um };
    G4double fCopperLayerThickness { 5 * um };
    G4double fKaptonLayerThickness { 25 * um };
    G4double fNickelLayerThickness { 2 * um };
    G4double fFirstCavityThickness { 5 * mm };
    G4double fSecondCavityThickness { 10 * mm };

    G4double monitor1XPosition = -1059.0 *mm;
    G4String name;

G4double phi = 90. *deg;
    // Matrix definition for a 90 deg rotation. Also used for other volumes
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
  
    name = "monitorChamber-externalBox";
    G4Box* chamberSolid = new G4Box(name, fBoxOuterWidth / 2, fBoxOuterHeight / 2, fBoxOuterThickness / 2);
    G4LogicalVolume* chamberLog = new G4LogicalVolume(chamberSolid, wallMaterial, name);
    chamberPhys=new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(monitor1XPosition-10*cm,0.*cm,0.*cm)),name,chamberLog,physicalTreatmentRoom,
                                                false,
                                                0);
   // chamberLog->SetVisAttributes(green);

    name = "monitorChamber-innerBox";
   
    G4Box* innerchamberSolid = new G4Box(name, fBoxInnerWidth / 2, fBoxInnerHeight / 2, fBoxInnerThickness / 2);
    G4LogicalVolume* innerchamberLog = new G4LogicalVolume(innerchamberSolid, CavityMaterial, name);
    innerchamberPhys=new G4PVPlacement(0, G4ThreeVector(0,0,0), name,innerchamberLog, chamberPhys, false, 0);
   // innerchamberLog->SetVisAttributes(green);

    name = "monitorChamber-enterWindow";
 

    G4Box* enterWindowSolid = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fEnterWindowThickness / 2);
    G4LogicalVolume* enterWindowLog = new G4LogicalVolume(enterWindowSolid, WindowMaterial, name);
    enterWindowPhys= new G4PVPlacement(0,  G4ThreeVector((-fBoxOuterThickness + fEnterWindowThickness) / 2,0, 0), name,enterWindowLog, chamberPhys, false, 0);
    enterWindowLog->SetVisAttributes(green);

    name = "monitorChamber-enterElectrode";
    G4double enterElectrodeThickness = fCopperLayerThickness + fKaptonLayerThickness + fNickelLayerThickness;
    G4double exitElectrodeThickness = enterElectrodeThickness;
   
    G4Box* enterElectrodeSolid = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, enterElectrodeThickness / 2);
    G4LogicalVolume* enterElectrodeLog = new G4LogicalVolume(enterElectrodeSolid, CavityMaterial, name);

    enterElectrodePhys=new G4PVPlacement(0, G4ThreeVector(-(fCentralWindowThickness + enterElectrodeThickness) / 2 - (fCentralElectrodeThickness + fFirstCavityThickness),0,0),  name,enterElectrodeLog, innerchamberPhys, false, 0);
    enterElectrodeLog->SetVisAttributes(green);

    name = "monitorChamber-kaptonLayer1";
    G4double position_kaptonLayer1 =-enterElectrodeThickness+fKaptonLayerThickness / 2;
    G4Box* kaptonLayerSolid1 = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fKaptonLayerThickness / 2);
    G4LogicalVolume* kaptonLayerLog1 = new G4LogicalVolume(kaptonLayerSolid1,  KaptonLayerMaterial, name);
    kaptonLayerPhys1= new G4PVPlacement(0, G4ThreeVector(position_kaptonLayer1,0,0),  name,kaptonLayerLog1, enterElectrodePhys, false, 0);
    kaptonLayerLog1->SetVisAttributes(green);

    name = "monitorChamber-copperLayer1";
 
G4double position_copperLayer1=position_kaptonLayer1+(fKaptonLayerThickness+fCopperLayerThickness)/2;
    G4Box* copperLayerSolid1 = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fCopperLayerThickness / 2);
    G4LogicalVolume* copperLayerLog1 = new G4LogicalVolume(copperLayerSolid1, CopperLayerMaterial, name);
    copperLayerPhys1=new G4PVPlacement(0, G4ThreeVector(position_copperLayer1,0,0), name,copperLayerLog1, enterElectrodePhys, false, 0);

    name = "monitorChamber-nickelLayer1";
 
G4double position_nichelLayer1=position_copperLayer1+(fCopperLayerThickness + fNickelLayerThickness)/2;

    G4Box* nickelLayerSolid = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fNickelLayerThickness / 2);
    G4LogicalVolume* nickelLayerLog1 = new G4LogicalVolume(nickelLayerSolid,  NichelLayerMaterial, name);
    nickelLayerPhys1=new G4PVPlacement(0,  G4ThreeVector(position_nichelLayer1,0,0), name,nickelLayerLog1, enterElectrodePhys, false, 0);

    name = "monitorChamber-firstCavity";
  
  
G4double position_firstCavity=-fCentralWindowThickness/2-fCentralElectrodeThickness-(fFirstCavityThickness) / 2;
    G4Box* firstCavitySolid = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fFirstCavityThickness / 2);
    G4LogicalVolume*fFirstCavityLog = new G4LogicalVolume(firstCavitySolid, CavityMaterial, name);
    fFirstCavityPhys=new G4PVPlacement(0,  G4ThreeVector( position_firstCavity,0,0), name, fFirstCavityLog, innerchamberPhys, false, 0);

    name = "monitorChamber-centralElectrode1";
 
G4double position_centralElectrode1=(-fCentralWindowThickness-fCentralElectrodeThickness) / 2;
    G4Box* centralElectrode1Solid = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fCentralElectrodeThickness / 2);
    G4LogicalVolume* centralElectrode1Log = new G4LogicalVolume(centralElectrode1Solid, ElectrodeMaterial, name);
    centralElectrode1Phys=new G4PVPlacement(0, G4ThreeVector(position_centralElectrode1,0,0), name,centralElectrode1Log, innerchamberPhys, false, 0);

    name = "monitorChamber-centralWindow";
  //  position={0,0,0};
    G4Box* centralWindowSolid = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fCentralWindowThickness / 2);
    G4LogicalVolume* centralWindowLog = new G4LogicalVolume(centralWindowSolid, CentralWindowMaterial, name);
    centralWindowPhys =new G4PVPlacement(0, G4ThreeVector(0,0,0),  name,centralWindowLog, innerchamberPhys, false, 0);
    centralWindowLog->SetVisAttributes(green);

    name = "monitorChamber-centralElectrode2";
 
G4double position_centralElectrode2=(fCentralWindowThickness+fCentralElectrodeThickness) / 2;
    G4Box* centralElectrode2Solid = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fCentralElectrodeThickness / 2);
    G4LogicalVolume* centralElectrode2Log = new G4LogicalVolume(centralElectrode2Solid, ElectrodeMaterial, name);
   centralElectrode2Phys= new G4PVPlacement(0, G4ThreeVector(position_centralElectrode2,0,0),name, centralElectrode2Log,innerchamberPhys, false, 0);

    name = "monitorChamber-secondCavity";
 
G4double position_secondCavity=(fCentralWindowThickness)/2+fCentralElectrodeThickness+(fSecondCavityThickness) / 2;
    G4Box* secondCavitySolid = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fSecondCavityThickness / 2);
    G4LogicalVolume*fSecondCavityLog = new G4LogicalVolume(secondCavitySolid, CavityMaterial, name);
    fSecondCavityPhys=new G4PVPlacement(0, G4ThreeVector(position_secondCavity,0,0), name, fSecondCavityLog, innerchamberPhys, false, 0);

    name="monitorChamber-exitElectrode";
  
G4double position_exitElectrode=(fCentralWindowThickness)/2+fCentralElectrodeThickness+fSecondCavityThickness+(exitElectrodeThickness) / 2;
    G4Box* exitElectrodeSolid = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, exitElectrodeThickness / 2);
    G4LogicalVolume* exitElectrodeLog = new G4LogicalVolume(exitElectrodeSolid, ElectrodeMaterial, name);
    exitElectrodePhys=new G4PVPlacement(0, G4ThreeVector(position_exitElectrode,0,0), name,exitElectrodeLog,innerchamberPhys, false, 0);
    exitElectrodeLog->SetVisAttributes(green);

    name = "monitorChamber-kaptonLayer2";
   
G4double position_kaptonLayer2=(exitElectrodeThickness-fKaptonLayerThickness) / 2;
    G4Box* kaptonLayerSolid2 = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fKaptonLayerThickness / 2);
    G4LogicalVolume* kaptonLayerLog2 = new G4LogicalVolume(kaptonLayerSolid2,  KaptonLayerMaterial, name);
    kaptonLayerPhys2=new G4PVPlacement(0, G4ThreeVector(position_kaptonLayer2,0,0), name,kaptonLayerLog2, exitElectrodePhys, false, 0);

    name = "monitorChamber-copperLayer2";
  
G4double position_copperLayer2=(exitElectrodeThickness)/2-fKaptonLayerThickness-(fCopperLayerThickness)/2;
    G4Box* copperLayerSolid2 = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fCopperLayerThickness / 2);
    G4LogicalVolume* copperLayerLog2 = new G4LogicalVolume(copperLayerSolid2,  CopperLayerMaterial, name);
    copperLayerPhys2=new G4PVPlacement(0, G4ThreeVector(position_copperLayer2,0,0), name, copperLayerLog2, exitElectrodePhys, false, 0);

    name = "monitorChamber-nickelLayer2";
  
G4double position_nichelLayer2=(exitElectrodeThickness)/2-fKaptonLayerThickness-fCopperLayerThickness-fNickelLayerThickness/2;

    G4Box* nickelLayerSolid2 = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fNickelLayerThickness / 2);
    G4LogicalVolume* nickelLayerLog2 = new G4LogicalVolume(nickelLayerSolid2,  NichelLayerMaterial, name);
     nickelLayerPhys2=new G4PVPlacement(0, G4ThreeVector(position_nichelLayer2,0,0), name,nickelLayerLog2, exitElectrodePhys, false, 0);

    name = "monitorChamber-exitWindow";
 
G4double position_exitWindow=(fBoxOuterThickness - fEnterWindowThickness) / 2;
    G4Box* exitWindowSolid = new G4Box(name, fInnerWidth / 2, fInnerHeight / 2, fExitWindowThickness / 2);
    G4LogicalVolume* exitWindowLog = new G4LogicalVolume(exitWindowSolid, WindowMaterial, name);
   exitWindowPhys= new G4PVPlacement(0, G4ThreeVector(position_exitWindow,0,0), name,exitWindowLog, chamberPhys, false, 0);
    exitWindowLog->SetVisAttributes(green);

    
     }
/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::BESTBeamNozzle()
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
void BESTPassiveProtonBeamLine::BESTBeamFinalCollimator()
{
    // -----------------------//
    //     FINAL COLLIMATOR   //
    //------------------------//
    const G4double outerRadiusFinalCollimator = 21.5*mm;
    const G4double hightFinalCollimator = 4.5*mm;
    const G4double startAngleFinalCollimator = 0.*deg;
    const G4double spanningAngleFinalCollimator = 360.*deg;
    const G4double finalCollimatorXPosition = -82.5 *mm;
    
    G4double phi = 90. *deg;
    
    // Matrix definition for a 90 deg rotation. Also used for other volumes
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidFinalCollimator = new G4Tubs("FinalCollimator",
                                      BESTinnerRadiusFinalCollimator,
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


/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::SetRangeShifterXSize(G4double value)
{
    solidBESTRangeShifterBox -> SetXHalfLength(value) ;
    G4cout << "RangeShifter size X (mm): "<< ((solidBESTRangeShifterBox -> GetXHalfLength())*2.)/mm
    << G4endl;
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::SetFirstScatteringFoilXSize(G4double value)
{
    BESTfirstScatteringFoil -> SetXHalfLength(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout <<"The X size of the first scattering foil is (mm):"<<
    ((BESTfirstScatteringFoil -> GetXHalfLength())*2.)/mm
    << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::SetSecondScatteringFoilXSize(G4double value)
{
    BESTsecondScatteringFoil -> SetXHalfLength(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout <<"The X size of the second scattering foil is (mm):"<<
    ((BESTsecondScatteringFoil -> GetXHalfLength())*2.)/mm
    << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::SetOuterRadiusStopper(G4double value)
{
    solidBESTStopper -> SetOuterRadius(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout << "OuterRadius od the Stopper is (mm):"
    << solidBESTStopper -> GetOuterRadius()/mm
    << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::SetInnerRadiusFinalCollimator(G4double value)
{
    solidFinalCollimator -> SetInnerRadius(value);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout<<"Inner Radius of the final collimator is (mm):"
	<< solidFinalCollimator -> GetInnerRadius()/mm
	<< G4endl;
}

/////////////////////////////////////////////////////////////////////////////
void BESTPassiveProtonBeamLine::SetRSMaterial(G4String materialChoice)
{
    if (G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice, false) )
    {
        if (pttoMaterial)
        {
            rangeShifterMaterial  = pttoMaterial;
            logicBESTRangeShifterBox -> SetMaterial(pttoMaterial);
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
void BESTPassiveProtonBeamLine::SetModulatorAngle(G4double value)
{
    modulator -> SetModulatorAngle(value);
    //G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
/////////////////////////////////////////////////////////////////////////////
