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
#include "G4UnionSolid.hh"
#include "G4Trd.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "G4NistElementBuilder.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "TrentoPassiveProtonBeamLine.hh"
#include "TrentoPassiveProtonBeamLineMessenger.hh"

/////////////////////////////////////////////////////////////////////////////
TrentoPassiveProtonBeamLine::TrentoPassiveProtonBeamLine():
 logicTreatmentRoom(0), physicalTreatmentRoom(0),hadrontherapyDetectorConstruction(0),
physiBeamLineSupport(0), physiBeamLineCover(0), physiBeamLineCover2(0),
physiMonitorLayer1(0), physiMonitorLayer2(0),
  ScatteringFoil(0), logicScatteringFoil(0), physiScatteringFoil(0), ridgeFilterPhys(0), preCollimator(0), physiPreCollimator(0), Collimator(0), physiCollimator(0), solidAirTube(0), solidAirPreTube(0), physiAirTube(0), physiAirPreTube(0)


{
    // Messenger to change parameters of the passiveTrentoBeamLine geometry
    TrentoPassiveMessenger = new TrentoPassiveProtonBeamLineMessenger(this);
    
    //***************************** PW ***********************************
    static G4String ROGeometryName = "DetectorROGeometry";
    RO = new HadrontherapyDetectorROGeometry(ROGeometryName);
    
    G4cout << "Going to register Parallel world...";
    RegisterParallelWorld(RO);
    G4cout << "... done" << G4endl;
    //********************************************************************
    
}
/////////////////////////////////////////////////////////////////////////////
TrentoPassiveProtonBeamLine::~TrentoPassiveProtonBeamLine()
{
    delete TrentoPassiveMessenger;
    delete hadrontherapyDetectorConstruction;
}

using namespace std;

/////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* TrentoPassiveProtonBeamLine::Construct()
{
    // Sets default geometry and materials
    SetDefaultDimensions();
    
    // Construct the whole Passive Beam Line
    ConstructTrentoPassiveProtonBeamLine();
    
    //***************************** PW ***************************************
    if (!hadrontherapyDetectorConstruction)
        
        // HadrontherapyDetectorConstruction builds ONLY the phantom and the detector with its associated ROGeometry
        hadrontherapyDetectorConstruction = new HadrontherapyDetectorConstruction(physicalTreatmentRoom);
    
    
    //********************************************************************
    
    hadrontherapyDetectorConstruction->InitializeDetectorROGeometry(RO,hadrontherapyDetectorConstruction->GetDetectorToWorldPosition());
    
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
void TrentoPassiveProtonBeamLine::SetDefaultDimensions()
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
    darkOrange3 -> SetForceSolid(false);
	
    skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
    skyBlue -> SetVisibility(true);
    skyBlue -> SetForceSolid(true);
    
    
    
    // SCATTERING FOIL: it is a thin foil that provides the
    // final diffusion of the beam.
    G4double defaultScatteringFoilXSize = 0.75*cm;
    ScatteringFoilXSize = defaultScatteringFoilXSize;
    
    G4double defaultScatteringFoilYSize = 100*mm;
    ScatteringFoilYSize = defaultScatteringFoilYSize;
    
    G4double defaultScatteringFoilZSize = 100 *mm;
    ScatteringFoilZSize = defaultScatteringFoilZSize;
    
    G4double defaultScatteringFoilXPosition = -273.*cm;
    ScatteringFoilXPosition = defaultScatteringFoilXPosition;
    
    G4double defaultScatteringFoilYPosition =  0 *mm;
    ScatteringFoilYPosition = defaultScatteringFoilYPosition;
    
    G4double defaultScatteringFoilZPosition =  0 *mm;
    ScatteringFoilZPosition = defaultScatteringFoilZPosition;

    //PRE COLLIMATOR: it is a PMMA collimator
    G4double defaultPreCollimatorXHalfSide = 12.5 * cm;
    preCollimatorXHalfSide = defaultPreCollimatorXHalfSide;

    G4double defaultPreCollimatorXPosition = -50.*cm;
    preCollimatorXPosition = defaultPreCollimatorXPosition;

    //DEFAULT DEFINITION OF THE AIR TUBE
    G4double defaultYHalfSideAirTube= 5.*cm;
    YHalfSideAirTube = defaultYHalfSideAirTube;

    G4double defaultZHalfSideAirTube = 5.*cm;
    ZHalfSideAirTube = defaultZHalfSideAirTube;

    
    // DEFAULT DEFINITION OF THE MATERIALS
    // All elements and compound definition follows the NIST database
    
    // ELEMENTS
    G4bool isotopes = false;
    G4Material* aluminumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al", isotopes);
    G4Material* copperNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);
    G4Element* zincNist = G4NistManager::Instance()->FindOrBuildElement("Zn");
    G4Element* copper = G4NistManager::Instance()->FindOrBuildElement("Cu");
    G4Element* silicNist = G4NistManager::Instance()->FindOrBuildElement("Si");
    G4Element* hydrogenNist = G4NistManager::Instance()->FindOrBuildElement("H");
    G4Element* oxygenNist = G4NistManager::Instance()->FindOrBuildElement("O");
    G4Element* carbonNist = G4NistManager::Instance()->FindOrBuildElement("C");
    
    // COMPOUND
    G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
    G4Material* PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);
    G4Material* MylarNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");

    G4int nComponents; // Number of components
    G4int nAtoms; //Number of atoms
    G4double fractionmass; // Fraction in mass of an element in a material

    //Quartz
    G4Material* SiO2 = new G4Material("quartz", 2.200*g/cm3, nComponents=2);
    SiO2->AddElement(silicNist, nAtoms=1);
    SiO2->AddElement(oxygenNist , nAtoms=2);
    
    //Epoxy (for FR4 )
    G4Material* Epoxy = new G4Material("Epoxy" , 1.2*g/cm3, nComponents=2);
    Epoxy->AddElement(hydrogenNist, nAtoms=2);
    Epoxy->AddElement(carbonNist, nAtoms=2);
    
    //FR4 (Glass + Epoxy)
    G4Material* FR4 = new G4Material("FR4"  , 1.86*g/cm3, nComponents=2);
    FR4->AddMaterial(SiO2, fractionmass=0.528);
    FR4->AddMaterial(Epoxy, fractionmass=0.472);
    
    //Brass
    G4Material* brass = new G4Material("Brass", 8.40*g/cm3, nComponents=2);
    brass -> AddElement(zincNist, fractionmass = 30 *perCent);
    brass -> AddElement(copper, fractionmass = 70 *perCent);

    //Plastic
    G4Material* plastic = new G4Material("Plastic", 1.40*g/cm3, nComponents=3);
    plastic -> AddElement(carbonNist, nAtoms = 21);
    plastic -> AddElement(oxygenNist, nAtoms = 4);
    plastic -> AddElement(hydrogenNist, nAtoms = 24);
     
    
    //***************************** PW ***************************************
    
    // DetectorROGeometry Material
    new G4Material("dummyMat", 1., 1.*g/mole, 1.*g/cm3);
    
    
    // MATERIAL ASSIGNMENT
    // Support of the beam line
    beamLineSupportMaterial = aluminumNist;

    // Matreial of the monitor chamber
    layerMonitorChamberMaterial = MylarNist;
    layerDefaultMaterial = airNist;
    internalStructureMaterial =  FR4;
    FoilMaterial = copperNist;
    airgapMaterial = airNist;
    
    // Material of the scattering foil
    ScatteringFoilMaterial = copperNist;

    //Material of the ridge filter
    singleTrapMaterial = plastic;
    
    // Materials of the collimators
    preCollimatorMaterial = PMMANist;
    CollimatorMaterial = brass;

    // Material of the final brass tube
    airTubeMaterial = airTube2Material = airTube3Material = airNist;
    
}

/////////////////////////////////////////////////////////////////////////////
void TrentoPassiveProtonBeamLine::ConstructTrentoPassiveProtonBeamLine()
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
    
    G4Box* treatmentRoom = new G4Box("TreatmentRoom",
                                     worldX,
                                     worldY,
                                     worldZ);
    
    logicTreatmentRoom = new G4LogicalVolume(treatmentRoom,
                                             airNist,
                                             "logicTreatmentRoom",
                                             0,0,0);
    
    physicalTreatmentRoom = new G4PVPlacement(0,
                                              G4ThreeVector(),
                                              "physicalTreatmentRoom",
                                              logicTreatmentRoom,
                                              0,
                                              false,
                                              0);
    
    
    // The treatment room is invisible in the Visualisation
    logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::Invisible);
    
    // Components of the Passive Proton Beam Line
    HadrontherapyBeamLineSupport();
    HadrontherapyBeamMonitoring();
    HadrontherapyBeamScatteringFoils();
    HadrontherapyRidgeFilter();
    HadrontherapyBeamCollimators();
    
}

/////////////////////////////////////////////////////////////////////////////
void TrentoPassiveProtonBeamLine::HadrontherapyBeamLineSupport()
{
    // ------------------//
    // BEAM LINE SUPPORT //
    //-------------------//
    const G4double beamLineSupportXSize = 1.2*m;
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
                                             physicalTreatmentRoom,
                                             false,
                                             0);
    
    // Visualisation attributes of the beam line support
    logicBeamLineSupport -> SetVisAttributes(gray);
    
    //---------------------------------//
    //  Beam line cover 1 (left panel) //
    //---------------------------------//
    const G4double beamLineCoverXSize = 1.2*m;
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
void TrentoPassiveProtonBeamLine::HadrontherapyBeamMonitoring()
{
  // ----------------------------
  //   MONITOR CHAMBER
  // ----------------------------
  // The  monitor chamber is a ionisation chamber
  // Each chamber consist in a sequence of 2 mm of air and
  // 8 microm copper layers which contain 800 microm of Fr4 in a box
  // that has the two layer of Mylar

    
////////////////////////////////////////////
///External Structure of the Monitor Chamber
////////////////////////////////////////////
    
  const G4double monitorXSize = 12.*um;
  const G4double monitorYSize = 12.7*cm;
  const G4double monitorZSize = 12.7*cm;
  
  const G4double shift = 17.*cm;
  
  const G4double monitorlayer1XPosition = -269.*cm +shift;
  const G4double monitorlayer2XPosition = -270.4104*cm +shift;
  
  
// First Mylar layer of the monitor chamber
  G4Box* solidMonitorLayer1 = new G4Box("MonitorLayer1",
					monitorXSize/2,
					monitorYSize/2,
					monitorZSize/2);
  
  G4LogicalVolume* logicMonitorLayer1 = new G4LogicalVolume(solidMonitorLayer1,
							    layerMonitorChamberMaterial,
							    "MonitorLayer1");
    
  physiMonitorLayer1 = new G4PVPlacement(0,
                                         G4ThreeVector(monitorlayer1XPosition, 0.*cm, 0.*cm),
                                         "MonitorLayer1",
                                         logicMonitorLayer1,
                                         physicalTreatmentRoom,
                                         false,
                                         0);

// Second Mylar layer of the monitor chamber
  G4Box* solidMonitorLayer2 = new G4Box("MonitorLayer2",
                                        monitorXSize/2,
                                        monitorYSize/2,
                                        monitorZSize/2);
  
  G4LogicalVolume* logicMonitorLayer2 = new G4LogicalVolume(solidMonitorLayer2,
                                                            layerMonitorChamberMaterial,
                                                            "MonitorLayer2");
  
  physiMonitorLayer2 = new G4PVPlacement(0,
					 G4ThreeVector(monitorlayer2XPosition, 0.*cm, 0.*cm),
					 "MonitorLayer2",
					 logicMonitorLayer2,
					 physicalTreatmentRoom,
					 false,
					 0);    
  
  logicMonitorLayer1 -> SetVisAttributes(gray);
  logicMonitorLayer2 -> SetVisAttributes(gray);
  
    
///////////////////////////////////////////
// Internal Structure of the Monitor Chamber
////////////////////////////////////////////
    
  const G4double layerThickness = 816*um;
  const G4double layerXposition = -270.2684*cm;
  const G4int nofLayers = 5;
  const G4double internalStructureThickness = 800.*um;
  const G4double airGapThickness = 2.*mm;
  const G4double foilThickness = 8.*um;
  const G4double FirstCoppperLayerXPosition = -404.*um;
  const G4double SecondCoppperLayerXPosition = 404.*um;
  const G4double InternalStructureXPosition = 0.;
  
// Air Layer
  G4VSolid* SolidAirLayer= new G4Box("Layer",
                              layerThickness,
                              monitorYSize/2,
                              monitorZSize/2);
    
  new G4LogicalVolume(SolidAirLayer,
                      layerDefaultMaterial,
                      "Layer");



// First Copper Layer
  G4VSolid* SolidFirstCoppperLayer = new G4Box("SolidFirstCoppperLayer",
                                  foilThickness/2,
                                  monitorYSize/2,
                                  monitorZSize/2);
  
  G4LogicalVolume* LogicFirstCoppperLayer = new G4LogicalVolume(
                                                     SolidFirstCoppperLayer,
                                                     FoilMaterial,
                                                     "SolidFirstCoppperLayer");
  
  
// Fr4 Internal Layer
  G4VSolid* SolidInternalStructure = new G4Box("SolidInternalStructure",
                                          internalStructureThickness/2,
                                          monitorYSize/2,
                                          monitorZSize/2);
  
  G4LogicalVolume* LogicInternalStructure = new G4LogicalVolume(
                                                             SolidInternalStructure,
                                                             internalStructureMaterial,
                                                             "SolidInternalStructure");

    
    
// Second Copper Layer
  G4VSolid* SolidSecondCoppperLayer = new G4Box("SolidSecondCoppperLayer",
                                   foilThickness/2,
                                   monitorYSize/2,
                                   monitorZSize/2); // its size
  
  G4LogicalVolume* LogicSecondCoppperLayer= new G4LogicalVolume(
                                                     SolidSecondCoppperLayer,
                                                     FoilMaterial,
                                                     "SolidSecondCoppperLayer");

  
  G4RotationMatrix Ra, Rm;
  G4ThreeVector Ta;
  G4Transform3D Tr;
  G4AssemblyVolume* assemblyMonitor = new G4AssemblyVolume();

 
    Ta.setX(FirstCoppperLayerXPosition); Ta.setY(0.); Ta.setZ(0.);
    Tr = G4Transform3D(Ra,Ta);
    assemblyMonitor->AddPlacedVolume(LogicFirstCoppperLayer, Tr);
    
    Ta.setX(InternalStructureXPosition); Ta.setY(0.); Ta.setZ(0.);
    Tr = G4Transform3D(Ra,Ta);
    assemblyMonitor->AddPlacedVolume(LogicInternalStructure, Tr );
    
    Ta.setX(SecondCoppperLayerXPosition); Ta.setY(0.); Ta.setZ(0.);
    Tr = G4Transform3D(Ra,Ta);
    
    assemblyMonitor->AddPlacedVolume(LogicSecondCoppperLayer, Tr );
   
  
    for( unsigned int i = 0; i<nofLayers; i++ )
   {
     G4ThreeVector Tm(shift+layerXposition + i*(airGapThickness), 0., 0.);
     Tr = G4Transform3D(Rm,Tm);
     assemblyMonitor -> MakeImprint(logicTreatmentRoom, Tr);
   }
    
}


void TrentoPassiveProtonBeamLine::HadrontherapyBeamScatteringFoils()
{
  
  // ---------------------------//
  // SCATTERING FOIL            //
  // ---------------------------//
  // It is a metal foil and provides the
  // initial diffusion of the beam. 
  
  
  
  ScatteringFoil = new G4Box("ScatteringFoil",
			     ScatteringFoilXSize,
			     ScatteringFoilYSize,
			     ScatteringFoilZSize);
  
  
  logicScatteringFoil = new G4LogicalVolume(ScatteringFoil,
					    ScatteringFoilMaterial,
					    "ScatteringFoil");
  
  physiScatteringFoil = new G4PVPlacement(0, G4ThreeVector(ScatteringFoilXPosition,
							   ScatteringFoilYPosition,
							   ScatteringFoilZPosition),
					  "SeconScatteringFoil",
                                              logicScatteringFoil,
					  physicalTreatmentRoom,
					  false,
					  0);

  logicScatteringFoil -> SetVisAttributes(skyBlue);
}

void TrentoPassiveProtonBeamLine::HadrontherapyRidgeFilter()
{
  // ---------------------------- //
  //         THE Ridge Filter    //
  // -----------------------------//
  // Energy degreader of
  // primary beam. Made of a series of pin-substructures.
  
    
    G4double ridgeXPosition = -200*cm;
    
    const G4double XBase = 0.1 *mm;
    const G4double YBase = 38.75*mm;
    
    G4VSolid* RidgeBase = new G4Box("Base_RidgeFilter",
                                    XBase,
                                    YBase,
                                    YBase);
    
    G4LogicalVolume* RidgeBaseLog = new G4LogicalVolume(RidgeBase,
                                                        singleTrapMaterial,
                                                        "Base_RidgeFilter");
    
    new G4PVPlacement(0,
                      G4ThreeVector(
                                    -199.7*cm,
                                    0.,
                                    0.),
                      "Base_RidgeFilter",
                      RidgeBaseLog,
                      physicalTreatmentRoom,
                      false,
                      0);
    
    RidgeBaseLog->SetVisAttributes(yellow);
    
    

  G4double x0 = 1.215*mm;
  G4double x1 = 1*mm;
  G4double x2 = 0.95*mm;
  G4double x3 = 0.87*mm;
  G4double x4 = 0.76*mm;
  G4double x5 = 0.71*mm;
  G4double x6 = 0.65*mm;
  G4double x7 = 0.56*mm;
  G4double x8 = 0.529*mm;
  G4double x9 = 0.258*mm;
  G4double x10 = 0.225*mm;
  G4double x11 = 0.144*mm;
  G4double x12 = 0.055*mm;
    
  G4double heigth = 0.13*mm;
  G4double heigth0 = 0.11*mm;
  G4double heigth1 = 2.3*mm;
  G4double heigth2 = 0.65*mm;
  G4double heigth3 = 1.9*mm;
  G4double heigth4 = 0.615*mm;
  G4double heigth5 = 2.23*mm;
  G4double heigth6 = 0.3*mm;
  G4double heigth7 = 4.1*mm;
  G4double heigth8 = 0.1*mm;
  G4double heigth9 = 0.105*mm;
  G4double heigth10 = 0.065*mm;
  
  G4VSolid* Trap00 = new G4Trd("singleTrap00", x0, x1, x0, x1, heigth);
  G4VSolid* Trap0 = new G4Trd("singleTrap0", x1, x2, x1, x2, heigth0);
  G4VSolid* Trap1 = new G4Trd("singleTrap1", x2, x3, x2, x3, heigth1);
  G4VSolid* Trap2 = new G4Trd("singleTrap2", x3, x4, x3, x4, heigth2);
  G4VSolid* Trap3 = new G4Trd("singleTrap3", x4, x5, x4, x5, heigth3);
  G4VSolid* Trap4 = new G4Trd("singleTrap4", x5, x6, x5, x6, heigth4);
  G4VSolid* Trap5 = new G4Trd("singleTrap5", x6, x7, x6, x7, heigth5);
  G4VSolid* Trap6 = new G4Trd("singleTrap6", x7, x8, x7, x8, heigth6);
  G4VSolid* Trap7 = new G4Trd("singleTrap7", x8, x9, x8, x9, heigth7);
  G4VSolid* Trap8 = new G4Trd("singleTrap8", x9, x10, x9, x10, heigth8);
  G4VSolid* Trap9 = new G4Trd("singleTrap9", x10, x11, x10, x11, heigth9);
  G4VSolid* Trap10 = new G4Trd("singleTrap10", x11, x12, x11, x12, heigth10);
  
 
  G4ThreeVector tr0(0., 0., 0.24);
  G4ThreeVector tr1(0., 0., 2.54);
  G4ThreeVector tr2(0., 0., 2.55);
  G4ThreeVector tr3(0., 0., 2.845);
  G4ThreeVector tr4(0., 0., 4.4);
  
  G4RotationMatrix* yRot = new G4RotationMatrix; 
  yRot->rotateY(0); 
  G4UnionSolid* unionTrap00 = new G4UnionSolid("Trap01", Trap00, Trap0, yRot, tr0);
  G4UnionSolid* unionTrap01 = new G4UnionSolid("Trap01", unionTrap00, Trap1, yRot, tr1);
  G4UnionSolid* unionTrap23 = new G4UnionSolid("Trap23", Trap2, Trap3, yRot, tr2);
  G4UnionSolid* unionTrap45 = new G4UnionSolid("Trap45", Trap4, Trap5, yRot, tr3);
  G4UnionSolid* unionTrap67 = new G4UnionSolid("Trap67", Trap6, Trap7, yRot, tr4);
  
  G4ThreeVector tr03(0., 0., 5.09);
  G4UnionSolid* unionTrap03 = new G4UnionSolid("unionTrap03", unionTrap01, unionTrap23, yRot, tr03);
  
  G4ThreeVector tr05(0., 0., 7.935);
  G4UnionSolid* unionTrap05 = new G4UnionSolid("unionTrap05", unionTrap03, unionTrap45, yRot, tr05);
    
  G4ThreeVector tr_blocco1(0., 0., 12.335);
  G4UnionSolid* unionTrap_blocco1 = new G4UnionSolid("unionTrap_blocco1", unionTrap05, unionTrap67, yRot, tr_blocco1);
  
  G4ThreeVector tr_blocco2( 0., 0., 12.435);
  G4UnionSolid* unionTrap_blocco2 = new G4UnionSolid("unionTrap_blocco2", unionTrap_blocco1, Trap8, yRot, tr_blocco2);
  
  G4ThreeVector tr_blocco3(0., 0., 12.54);
  G4UnionSolid* unionTrap_blocco3 = new G4UnionSolid("unionTrap_blocco3", unionTrap_blocco2, Trap9, yRot, tr_blocco3);
  
  G4ThreeVector tr_tot( 0., 0., 12.605);
  G4UnionSolid* unionTrap = new G4UnionSolid("unionTrap", unionTrap_blocco3, Trap10, yRot, tr_tot);

  
  G4LogicalVolume* singleTrapLog = new G4LogicalVolume(unionTrap, singleTrapMaterial, "singleTrap");
  

  singleTrapLog->SetVisAttributes(yellow);
  
  
    G4int numberOfLayers = 31;
    G4double minZ = -37.5*mm;
    G4double minY = -37.5*mm;
    G4double sum_space = 1.25*mm;
    
    vector<G4ThreeVector> singleTrapPositions;
    
    for (int i = 0; i < numberOfLayers; i++)
      {
    	for (int j = 0; j < numberOfLayers; j++)
    	  {
	    singleTrapPositions.push_back({-0.01*cm+ridgeXPosition,
    		  minY + 2*i*sum_space,
    		  minZ + 2*j*sum_space,});
    	  }
      }
    
    G4double ti = -  90. *deg;
    G4RotationMatrix rt;
    rt.rotateY(ti);
    
    G4int peaks = numberOfLayers*numberOfLayers;
    for (int i = 0; i < peaks; i++)
      {
	
    	ostringstream tName;tName << "singleTrap" << i;
        new G4PVPlacement(G4Transform3D(rt,
                                        singleTrapPositions[i]),
                                        singleTrapLog,
                                        "tName.str()",
                                        logicTreatmentRoom,
                                        0,
                                        i);
	
      }
  
    
}

void TrentoPassiveProtonBeamLine::HadrontherapyBeamCollimators()
{
  
  // ------------------------//
  //     PRE - COLLIMATOR    //
  // -----------------------//
  // It is a plastic collimator to limit neutron dose and reduce activation
  
  const G4double preCollimatorYHalfSide = 10.*cm;
  const G4double preCollimatorZHalfSide = 10.*cm;
  
  preCollimator = new G4Box("PreCollimator",
  			     preCollimatorXHalfSide,
			    preCollimatorYHalfSide,
			    preCollimatorZHalfSide);
  
  
  G4LogicalVolume* logicPreCollimator = new G4LogicalVolume(preCollimator,
							    preCollimatorMaterial,
							    "PreCollimator");
  
  
  physiPreCollimator = new G4PVPlacement(0,
                                         G4ThreeVector(
                                                       -50.*cm,
                                                           0.,
							  0.),
					 "PreCollimator",
  					 logicPreCollimator,
                                         physicalTreatmentRoom,
                                         false,
                                         0);
  
  
  
  logicPreCollimator -> SetVisAttributes(white);
  
  
  // -----------------//
  //    COLLIMATOR     //
  // -----------------//
  // It is a brass collimator
  
  
G4double collimatorYHalfSide = 10.*cm;
G4double collimatorZHalfSide = 10.*cm;
G4double collimatorXHalfSide  = 3.25*cm;
  
G4double CollimatorXPosition = -34.25*cm;
  
  Collimator = new G4Box("Collimator",
			 collimatorXHalfSide,
			 collimatorYHalfSide,
			 collimatorZHalfSide);
  
  
  G4LogicalVolume* logicCollimator = new G4LogicalVolume(Collimator,
							 CollimatorMaterial,
							 "Collimator");
  
  
  
  physiCollimator = new G4PVPlacement(0, G4ThreeVector(CollimatorXPosition,
						       0.,
						       0.),
				      "Collimator",
				      logicCollimator,
                      physicalTreatmentRoom,
				      false,
				      0);
  
  
  
  logicCollimator -> SetVisAttributes(darkGreen);
  

    
    // ---------------------------------//
    //      AIR BOX Collimator          //
    // ---------------------------------//
    

    solidAirTube = new G4Box("AirTube",
                             collimatorXHalfSide,
                             YHalfSideAirTube,
                             ZHalfSideAirTube);
    
    G4LogicalVolume* logicAirTube = new G4LogicalVolume(solidAirTube,
                                                        airTubeMaterial,
                                                        "AirTube",
                                                        0, 0, 0);
    
    
    physiAirTube = new G4PVPlacement(0, G4ThreeVector(0,
                                                      0.,
                                                      0.),
                                     "AirTube",
                                     logicAirTube,
                                     physiCollimator,
                                     false,
                                     0);
    
    logicAirTube -> SetVisAttributes(darkOrange3);

    
    
    // // ---------------------------------//
    // //      AIR BOX PreCollimator       //
    // // ---------------------------------//

    
    solidAirPreTube = new G4Box("AirPreTube",
                                preCollimatorXHalfSide,
                                YHalfSideAirTube,
                                ZHalfSideAirTube);
    
    G4LogicalVolume* logicAirPreTube = new G4LogicalVolume(solidAirPreTube,
                                                           airTubeMaterial,
                                                           "AirPreTube",
                                                           0, 0, 0);
    
    
    physiAirPreTube = new G4PVPlacement(0,
                                        G4ThreeVector(0,
                                                      0.,
                                                      0.),
                                        "AirPreTube",
                                        logicAirPreTube,
                                        physiPreCollimator,
                                        false,
                                        0);
    
    logicAirPreTube -> SetVisAttributes(darkOrange3);



}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////// MESSENGER ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


void TrentoPassiveProtonBeamLine::SetScatteringFoilXSize(G4double value)
{
  ScatteringFoil -> SetXHalfLength(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4cout <<"The X size of the second scattering foil is (mm):"<<
      ((ScatteringFoil -> GetXHalfLength())*2.)/mm
	   << G4endl;
}


void TrentoPassiveProtonBeamLine::SetPreCollimatorXSize(G4double value)
{
  preCollimator -> SetXHalfLength(value);
  G4cout<<"The size of the pre collimator is (mm)"
	<< ((preCollimator -> GetXHalfLength())*2)/mm << G4endl;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
   
}

void TrentoPassiveProtonBeamLine::SetPreCollimatorXPosition(G4double value)
{
  physiPreCollimator -> SetTranslation(G4ThreeVector(value, 0., 0.));
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout  <<" The pre collimator is translated to" << value/mm <<"mm along the X axis" << G4endl;
   
}

void TrentoPassiveProtonBeamLine::SetAirTubeYSize(G4double value)
{
  solidAirTube -> SetYHalfLength(value);
  G4cout<<"The y side of brass tube is (mm)"
	<< ((preCollimator -> GetYHalfLength())*2) << G4endl;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();   
}


void TrentoPassiveProtonBeamLine::SetAirTubeZSize(G4double value)
{
  solidAirTube -> SetZHalfLength(value);
  G4cout<<"The z side of the brass tube is (mm)"
	<< ((Collimator -> GetZHalfLength())*2) << G4endl;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  
}


/////////////////////////////////////////////////////////////////////////////
void TrentoPassiveProtonBeamLine::SetScattererMaterial(G4String materialChoice)
{
  if (G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice, false) )
    {
      if (pttoMaterial)
        {
	  ScatteringFoilMaterial  = pttoMaterial;
	  logicScatteringFoil -> SetMaterial(pttoMaterial);
	  G4cout << "The material of the Scatterer has been changed to " << materialChoice << G4endl;
        }
    }
  else
    {
      G4cout << "WARNING: material \"" << materialChoice << "\" doesn't exist in NIST elements/materials"
	    " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
      G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl;
    }
}
