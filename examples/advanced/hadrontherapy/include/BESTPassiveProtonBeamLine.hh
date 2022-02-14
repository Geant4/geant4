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

#ifndef BESTPassiveProtonBeamLine_H
#define BESTPassiveProtonBeamLine_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

class G4VPhysicalVolume;
class HadrontherapyDetectorConstruction;
class HadrontherapyModulator;
class BESTPassiveProtonBeamLineMessenger;
class HadrontherapyDetectorROGeometry;

class BESTPassiveProtonBeamLine : public G4VUserDetectorConstruction
{
public:

    BESTPassiveProtonBeamLine();
    ~BESTPassiveProtonBeamLine();
  // static G4bool doCalculation;
    
    G4VPhysicalVolume* Construct();
    //***************************** PW **************NON SERVE*************************
    
    static BESTPassiveProtonBeamLine* GetInstance();
    
    //***************************** PW **************NON SERVE*************************
    
    void BESTBeamLineSupport();
    // Definition of the beam line support
    
    void BESTBeamScatteringFoils();
    // Definition of the first scattering foil,
    // of the Kapton window, of the stopper
    
    void BESTRangeShifter();
    // This defines the "range shifter". Is is a slab
    // (usually of PMMA" acting as energy degrader
    // of primary beam
    
    void BESTBeamCollimators();
    // Definition of the first collimator, of the range shifter,
    // of the second collimator, of the first and second
    // collimator modulators
    
    void BESTBeamMonitoring();
    // Definition of three monitor chambers
    
    
    void BESTBeamNozzle();
    // Definition of the beam noozle
    
    void BESTBeamFinalCollimator();
    // Definition of the final collimator
    
    // The following methods allow to change parameters
    // of some beam line components
    
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
    static BESTPassiveProtonBeamLine* instance;
    //passive proton line dimensions
    void SetDefaultDimensions();
    void ConstructBESTPassiveProtonBeamLine();
    
    HadrontherapyModulator* modulator; // Pointer to the modulator
    // geometry component
    BESTPassiveProtonBeamLineMessenger* passiveMessenger;
    G4VPhysicalVolume* physicalTreatmentRoom;
    HadrontherapyDetectorConstruction* hadrontherapyDetectorConstruction;
    
    
    G4Material* kapton;
	
    G4double BESTvacuumZoneXSize;
    G4double BESTvacuumZoneYSize;
    G4double BESTvacuumZoneZSize;
    G4double BESTvacuumZoneXPosition;
    
    G4double BESTfirstScatteringFoilXSize;
    G4double BESTfirstScatteringFoilYSize;
    G4double BESTfirstScatteringFoilZSize;
    G4double BESTfirstScatteringFoilXPosition;
    
    G4double BESTkaptonWindowXSize;
    G4double BESTkaptonWindowYSize;
    G4double BESTkaptonWindowZSize;
    G4double BESTkaptonWindowXPosition;
    
    G4double BESTinnerRadiusStopper;
    G4double BESTheightStopper;
    G4double BESTstartAngleStopper;
    G4double BESTspanningAngleStopper;
    G4double BESTstopperXPosition;
    G4double BESTstopperYPosition;
    G4double BESTstopperZPosition;
    G4double BESTouterRadiusStopper;
    
    G4double BESTsecondScatteringFoilXSize;
    G4double BESTsecondScatteringFoilYSize;
    G4double BESTsecondScatteringFoilZSize;
    G4double BESTsecondScatteringFoilXPosition;
    G4double BESTsecondScatteringFoilYPosition;
    G4double BESTsecondScatteringFoilZPosition;
    
    G4double BESTrangeShifterXSize;
    G4double BESTrangeShifterYSize;
    G4double BESTrangeShifterZSize;
    G4double BESTrangeShifterXPosition;
    G4double BESTrangeShifterYPosition;
    G4double BESTrangeShifterZPosition;
    
    
    G4VPhysicalVolume* physiBeamLineSupport;
    G4VPhysicalVolume* physiBeamLineCover;
    G4VPhysicalVolume* physiBeamLineCover2;
    G4Box* BESTfirstScatteringFoil;
    G4VPhysicalVolume* physiBESTFirstScatteringFoil;
    G4VPhysicalVolume* physiBESTKaptonWindow;
    
    G4Tubs* solidBESTStopper;
    G4VPhysicalVolume* physiBESTStopper;
    G4LogicalVolume* logicBESTStopper;
    
    G4Box* BESTsecondScatteringFoil;
    G4VPhysicalVolume* physiBESTSecondScatteringFoil;
    G4VPhysicalVolume* physiBESTFirstCollimator;
    G4VPhysicalVolume* physiBESTHoleFirstCollimator;
    G4Box* solidBESTRangeShifterBox;
    G4LogicalVolume* logicBESTRangeShifterBox;
    G4VPhysicalVolume* physiBESTRangeShifterBox;
    G4VPhysicalVolume* physiBESTSecondCollimator;
    G4VPhysicalVolume* physiBESTHoleSecondCollimator;
    
    G4VPhysicalVolume* physiBESTFirstCollimatorModulatorBox;
    G4VPhysicalVolume* physiBESTHoleFirstCollimatorModulatorBox;
    
    G4VPhysicalVolume* physiBESTSecondCollimatorModulatorBox;
    G4VPhysicalVolume* physiBESTHoleSecondCollimatorModulatorBox;
      
    G4double BESTinnerRadiusFinalCollimator;
    G4VPhysicalVolume* mother;


G4VPhysicalVolume* chamberPhys;
  G4VPhysicalVolume*innerchamberPhys;
  G4VPhysicalVolume*enterWindowPhys;
   G4VPhysicalVolume*enterElectrodePhys;
 G4VPhysicalVolume* kaptonLayerPhys1;
  G4VPhysicalVolume*copperLayerPhys1;
  G4VPhysicalVolume*nickelLayerPhys1;
  G4VPhysicalVolume*fFirstCavityPhys;
  G4VPhysicalVolume*centralElectrode1Phys;
  G4VPhysicalVolume*centralWindowPhys;
  G4VPhysicalVolume*centralElectrode2Phys;
G4VPhysicalVolume*fSecondCavityPhys;
  G4VPhysicalVolume*exitElectrodePhys;
  G4VPhysicalVolume* kaptonLayerPhys2;
  G4VPhysicalVolume*copperLayerPhys2;
  G4VPhysicalVolume*nickelLayerPhys2;

G4VPhysicalVolume* exitWindowPhys;

  G4Material* CopperLayerMaterial;
  G4Material* NichelLayerMaterial;
  G4Material* KaptonLayerMaterial;
  G4Material* WindowMaterial;
  G4Material* CentralWindowMaterial;
  G4Material* wallMaterial;
  G4Material* ElectrodeMaterial;
  G4Material* CavityMaterial;


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
 
    
    G4Material* nozzleSupportMaterial;
    G4Material* holeNozzleSupportMaterial;
    
    G4Material* brassTubeMaterial;
    G4Material* brassTube2Material;
    G4Material* brassTube3Material;
    G4Material* finalCollimatorMaterial;
    
    
    HadrontherapyDetectorROGeometry* RO;
    
    
};
#endif
