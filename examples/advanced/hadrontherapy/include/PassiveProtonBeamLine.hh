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

#ifndef PassiveProtonBeamLine_H
#define PassiveProtonBeamLine_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

class G4VPhysicalVolume;
class HadrontherapyDetectorConstruction;
class HadrontherapyModulator;
class PassiveProtonBeamLineMessenger;
class HadrontherapyDetectorROGeometry;

class PassiveProtonBeamLine : public G4VUserDetectorConstruction
{
public:

    PassiveProtonBeamLine();
    ~PassiveProtonBeamLine();
  // static G4bool doCalculation;
    
    G4VPhysicalVolume* Construct();
    //***************************** PW **************NON SERVE*************************
    
    static PassiveProtonBeamLine* GetInstance();
    
    //***************************** PW **************NON SERVE*************************
    
    void HadrontherapyBeamLineSupport();
    // Definition of the beam line support
    
    void HadrontherapyBeamScatteringFoils();
    // Definition of the first scattering foil,
    // of the Kapton window, of the stopper
    
    void HadrontherapyRangeShifter();
    // This defines the "range shifter". Is is a slab
    // (usually of PMMA" acting as energy degrader
    // of primary beam
    
    void HadrontherapyBeamCollimators();
    // Definition of the first collimator, of the range shifter,
    // of the second collimator, of the first and second
    // collimator modulators
    
    void HadrontherapyBeamMonitoring();
    // Definition of three monitor chambers
    
    void HadrontherapyMOPIDetector();
    // Construct the MOPI on-line detector
    
    void HadrontherapyBeamNozzle();
    // Definition of the beam noozle
    
    void HadrontherapyBeamFinalCollimator();
    // Definition of the final collimator
    
    // The following methods allow to change parameters
    // of some beam line components
    
    void SetRangeShifterXPosition(G4double value);
    // This method allows to move the Range Shifter along
    // the X axis
    
    void SetRangeShifterXSize(G4double halfSize);
    // This method allows to change the size of the range shifter along
    // the X axis
    
    void SetFirstScatteringFoilXSize(G4double);
    // This method allows to change the size of the first scattering foil
    // along the X axis
    
    void SetSecondScatteringFoilXSize(G4double);
    // This method allows to change the size of the second scattering foil
    // along the X axis
    
    void SetOuterRadiusStopper(G4double);
    // This method allows to change the size of the outer radius of the stopper
    
    void SetInnerRadiusFinalCollimator(G4double);
    // This method allows to change the size of the inner radius of the
    // final collimator
    
    void SetRSMaterial(G4String);
    // This method allows to change the material
    // of the range shifter
    
    void SetModulatorAngle(G4double angle);
    // This method allows moving the modulator through UI commands
    
    
private:
    static PassiveProtonBeamLine* instance;
    //passive proton line dimensions
    void SetDefaultDimensions();
    void ConstructPassiveProtonBeamLine();
    
    HadrontherapyModulator* modulator; // Pointer to the modulator
    // geometry component
    PassiveProtonBeamLineMessenger* passiveMessenger;
    G4VPhysicalVolume* physicalTreatmentRoom;
    HadrontherapyDetectorConstruction* hadrontherapyDetectorConstruction;
    
    
    G4Material* kapton;
	
    G4double vacuumZoneXSize;
    G4double vacuumZoneYSize;
    G4double vacuumZoneZSize;
    G4double vacuumZoneXPosition;
    
    G4double firstScatteringFoilXSize;
    G4double firstScatteringFoilYSize;
    G4double firstScatteringFoilZSize;
    G4double firstScatteringFoilXPosition;
    
    G4double kaptonWindowXSize;
    G4double kaptonWindowYSize;
    G4double kaptonWindowZSize;
    G4double kaptonWindowXPosition;
    
    G4double innerRadiusStopper;
    G4double heightStopper;
    G4double startAngleStopper;
    G4double spanningAngleStopper;
    G4double stopperXPosition;
    G4double stopperYPosition;
    G4double stopperZPosition;
    G4double outerRadiusStopper;
    
    G4double secondScatteringFoilXSize;
    G4double secondScatteringFoilYSize;
    G4double secondScatteringFoilZSize;
    G4double secondScatteringFoilXPosition;
    G4double secondScatteringFoilYPosition;
    G4double secondScatteringFoilZPosition;
    
    G4double rangeShifterXSize;
    G4double rangeShifterYSize;
    G4double rangeShifterZSize;
    G4double rangeShifterXPosition;
    G4double rangeShifterYPosition;
    G4double rangeShifterZPosition;
    
    
    G4VPhysicalVolume* physiBeamLineSupport;
    G4VPhysicalVolume* physiBeamLineCover;
    G4VPhysicalVolume* physiBeamLineCover2;
    G4Box* firstScatteringFoil;
    G4VPhysicalVolume* physiFirstScatteringFoil;
    G4VPhysicalVolume* physiKaptonWindow;
    
    G4Tubs* solidStopper;
    G4VPhysicalVolume* physiStopper;
    G4LogicalVolume* logicStopper;
    
    G4Box* secondScatteringFoil;
    G4VPhysicalVolume* physiSecondScatteringFoil;
    G4VPhysicalVolume* physiFirstCollimator;
    G4VPhysicalVolume* physiHoleFirstCollimator;
    G4Box* solidRangeShifterBox;
    G4LogicalVolume* logicRangeShifterBox;
    G4VPhysicalVolume* physiRangeShifterBox;
    G4VPhysicalVolume* physiSecondCollimator;
    G4VPhysicalVolume* physiHoleSecondCollimator;
    
    G4VPhysicalVolume* physiFirstCollimatorModulatorBox;
    G4VPhysicalVolume* physiHoleFirstCollimatorModulatorBox;
    
    G4VPhysicalVolume* physiSecondCollimatorModulatorBox;
    G4VPhysicalVolume* physiHoleSecondCollimatorModulatorBox;
    
    // MOPI Detector
    // Mother volume
    G4VPhysicalVolume* physiMOPIMotherVolume;
    G4LogicalVolume* logicMOPIMotherVolume;
    G4Box* solidMOPIMotherVolume;
    
    G4double MOPIMotherVolumeXSize;
    G4double MOPIMotherVolumeYSize;
    G4double MOPIMotherVolumeZSize;
    G4double MOPIMotherVolumeXPosition;
    G4double MOPIMotherVolumeYPosition;
    G4double MOPIMotherVolumeZPosition;
    
    // First Kapton layer
    G4double MOPIFirstKaptonLayerXSize;
    G4double MOPIFirstKaptonLayerYSize;
    G4double MOPIFirstKaptonLayerZSize;
    G4double MOPIFirstKaptonLayerXPosition;
    G4double MOPIFirstKaptonLayerYPosition;
    G4double MOPIFirstKaptonLayerZPosition;
    G4Box* solidMOPIFirstKaptonLayer;
    G4LogicalVolume* logicMOPIFirstKaptonLayer;
    G4VPhysicalVolume* physiMOPIFirstKaptonLayer;
    
    // First Aluminum layer
    G4double MOPIFirstAluminumLayerXSize;
    G4double MOPIFirstAluminumLayerYSize;
    G4double MOPIFirstAluminumLayerZSize;
    G4double MOPIFirstAluminumLayerXPosition;
    G4double MOPIFirstAluminumLayerYPosition;
    G4double MOPIFirstAluminumLayerZPosition;
    G4Box* solidMOPIFirstAluminumLayer;
    G4LogicalVolume* logicMOPIFirstAluminumLayer;
    G4VPhysicalVolume* physiMOPIFirstAluminumLayer;
    
    // First Air Gap
    G4double MOPIFirstAirGapXSize;
    G4double MOPIFirstAirGapYSize;
    G4double MOPIFirstAirGapZSize;
    G4double MOPIFirstAirGapXPosition;
    G4double MOPIFirstAirGapYPosition;
    G4double MOPIFirstAirGapZPosition;
    G4Box* solidMOPIFirstAirGap;
    G4LogicalVolume* logicMOPIFirstAirGap;
    G4VPhysicalVolume* physiMOPIFirstAirGap;
    
    // Cathode
    G4double MOPICathodeXSize;
    G4double MOPICathodeYSize;
    G4double MOPICathodeZSize;
    G4double MOPICathodeXPosition;
    G4double MOPICathodeYPosition;
    G4double MOPICathodeZPosition;
    G4Box* solidMOPICathode;
    G4LogicalVolume* logicMOPICathode;
    G4VPhysicalVolume* physiMOPICathode;
    
    G4VisAttributes* redWire;
    
    // First Air Gap
    G4double MOPISecondAirGapXSize;
    G4double MOPISecondAirGapYSize;
    G4double MOPISecondAirGapZSize;
    G4double MOPISecondAirGapXPosition;
    G4double MOPISecondAirGapYPosition;
    G4double MOPISecondAirGapZPosition;
    G4Box* solidMOPISecondAirGap;
    G4LogicalVolume* logicMOPISecondAirGap;
    G4VPhysicalVolume* physiMOPISecondAirGap;
    
    // First Aluminum layer
    G4double MOPISecondAluminumLayerXSize;
    G4double MOPISecondAluminumLayerYSize;
    G4double MOPISecondAluminumLayerZSize;
    G4double MOPISecondAluminumLayerXPosition;
    G4double MOPISecondAluminumLayerYPosition;
    G4double MOPISecondAluminumLayerZPosition;
    G4Box* solidMOPISecondAluminumLayer;
    G4LogicalVolume* logicMOPISecondAluminumLayer;
    G4VPhysicalVolume* physiMOPISecondAluminumLayer;
    
    // Second Kapton layer
    G4double MOPISecondKaptonLayerXSize;
    G4double MOPISecondKaptonLayerYSize;
    G4double MOPISecondKaptonLayerZSize;
    G4double MOPISecondKaptonLayerXPosition;
    G4double MOPISecondKaptonLayerYPosition;
    G4double MOPISecondKaptonLayerZPosition;
    G4Box* solidMOPISecondKaptonLayer;
    G4LogicalVolume* logicMOPISecondKaptonLayer;
    G4VPhysicalVolume* physiMOPISecondKaptonLayer;
    
    G4double innerRadiusFinalCollimator;
    G4VPhysicalVolume* mother;
    
    G4VPhysicalVolume* physiFirstMonitorLayer1;
    G4VPhysicalVolume* physiFirstMonitorLayer2;
    G4VPhysicalVolume* physiFirstMonitorLayer3;
    G4VPhysicalVolume* physiFirstMonitorLayer4;
    G4VPhysicalVolume* physiSecondMonitorLayer1;
    G4VPhysicalVolume* physiSecondMonitorLayer2;
    G4VPhysicalVolume* physiSecondMonitorLayer3;
    G4VPhysicalVolume* physiSecondMonitorLayer4;
    G4VPhysicalVolume* physiNozzleSupport;
    G4VPhysicalVolume* physiHoleNozzleSupport;
    G4VPhysicalVolume* physiBrassTube;
    G4VPhysicalVolume* physiBrassTube2;
    G4VPhysicalVolume* physiBrassTube3;
    G4Tubs* solidFinalCollimator;
    G4VPhysicalVolume* physiFinalCollimator;
    
    G4VisAttributes* blue;
    G4VisAttributes* gray;
    G4VisAttributes* white;
    G4VisAttributes* red;
    G4VisAttributes* yellow;
    G4VisAttributes* green;
    G4VisAttributes* darkGreen;
    G4VisAttributes* darkOrange3;
    G4VisAttributes* skyBlue;
    
    G4Material* rangeShifterMaterial;
    G4Material* beamLineSupportMaterial;
    G4Material* vacuumZoneMaterial;
    G4Material* firstScatteringFoilMaterial;
    G4Material* kaptonWindowMaterial;
    G4Material* stopperMaterial;
    G4Material* secondScatteringFoilMaterial;
    G4Material* firstCollimatorMaterial;
    G4Material* holeFirstCollimatorMaterial;
    G4Material* modulatorBoxMaterial;
    G4Material* holeModulatorBoxMaterial;
    G4Material* layer1MonitorChamberMaterial;
    G4Material* layer2MonitorChamberMaterial;
    G4Material* layer3MonitorChamberMaterial;
    G4Material* layer4MonitorChamberMaterial;
    G4Material* MOPIMotherVolumeMaterial;
    G4Material* MOPIFirstKaptonLayerMaterial;
    G4Material* MOPIFirstAluminumLayerMaterial;
    G4Material* MOPIFirstAirGapMaterial;
    G4Material* MOPICathodeMaterial;
    G4Material* MOPISecondAirGapMaterial;
    G4Material* MOPISecondAluminumLayerMaterial;
    G4Material* MOPISecondKaptonLayerMaterial;
    G4Material* nozzleSupportMaterial;
    G4Material* holeNozzleSupportMaterial;
    
    G4Material* brassTubeMaterial;
    G4Material* brassTube2Material;
    G4Material* brassTube3Material;
    G4Material* finalCollimatorMaterial;
    
    
    HadrontherapyDetectorROGeometry* RO;
    
    
};
#endif
