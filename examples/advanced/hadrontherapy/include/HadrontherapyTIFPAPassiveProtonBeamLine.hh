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

#ifndef TrentoPassiveProtonBeamLine_H
#define TrentoPassiveProtonBeamLine_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

class G4VPhysicalVolume;
class HadrontherapyDetectorConstruction;
class HadrontherapyModulator;
class TrentoPassiveProtonBeamLineMessenger;
class HadrontherapyDetectorROGeometry;

class TrentoPassiveProtonBeamLine : public G4VUserDetectorConstruction
{
public:
    
    //TrentoPassiveProtonBeamLine(G4VPhysicalVolume*);
    TrentoPassiveProtonBeamLine();
    ~TrentoPassiveProtonBeamLine();
  // static G4bool doCalculation;
    
    G4VPhysicalVolume* Construct();
    //
    static TrentoPassiveProtonBeamLine* GetInstance();
    
    //
    void HadrontherapyBeamLineSupport();
    // Definition of the beam line support
    
    void HadrontherapyBeamMonitoring();
    //Definition of the monitor chamber
    
    void HadrontherapyBeamScatteringFoils();
    // Definition of the first scattering foil made of tantalum
    
    void HadrontherapyBeamCollimators();
    // Definition of the pre collimator made of plastic
    // and the collimator in front of the detector
    
    void HadrontherapyRidgeFilter();
    // The following methods allow to change parameters
    // of some beam line components (through Messenger files)

    void SetScatteringFoilXSize(G4double);
    // This method allows to change the size of the scattering foil
    // along the X axis
    
    void SetPreCollimatorXSize(G4double);
    // This method allows to change the size of the pre collimator along
    // the Z axis
    
    void SetPreCollimatorXPosition(G4double);
    // This method allows to change the position of the pre collimator along
    // the X axis
    
    void SetAirTubeYSize(G4double);
    // This method allows to change the y side of
    // air tube
    
    void SetAirTubeZSize(G4double);
    // This method allows to change the z side of
    // the air tube
    
    void SetScattererMaterial(G4String);
    // This method permit to change
    //the material of scatterer
    
private:
    static TrentoPassiveProtonBeamLine* instance;
    //passive proton line dimensions
    void SetDefaultDimensions();
    void ConstructTrentoPassiveProtonBeamLine();
      
    // geometry component
    TrentoPassiveProtonBeamLineMessenger* TrentoPassiveMessenger;
    G4LogicalVolume* logicTreatmentRoom;
    G4VPhysicalVolume* physicalTreatmentRoom;
    HadrontherapyDetectorConstruction* hadrontherapyDetectorConstruction;
    

    G4Material* kapton;
    
    G4double ScatteringFoilXSize;
    G4double ScatteringFoilYSize;
    G4double ScatteringFoilZSize;
    G4double ScatteringFoilXPosition;
    G4double ScatteringFoilYPosition;
    G4double ScatteringFoilZPosition;
    
    G4double preCollimatorXHalfSide;
    G4double preCollimatorXPosition;
    G4double preCollimatorYPosition;
    G4double preCollimatorZPosition;
   
    G4double YHalfSideAirTube;
    G4double ZHalfSideAirTube;
    
    G4VPhysicalVolume* physiBeamLineSupport;
    G4VPhysicalVolume* physiBeamLineCover;
    G4VPhysicalVolume* physiBeamLineCover2;
    
    G4VPhysicalVolume* physiMonitorLayer1;
    G4VPhysicalVolume* physiMonitorLayer2;
    
    G4VPhysicalVolume* internalMonitorStructurePV;
    G4VPhysicalVolume* firstFoilPV;
    G4VPhysicalVolume* secondFoilPV;
    G4VPhysicalVolume* airgapPV;
    
    G4Box* ScatteringFoil;
    G4LogicalVolume* logicScatteringFoil;
    G4VPhysicalVolume* physiScatteringFoil;
    
    G4VPhysicalVolume* ridgeFilterPhys;
    
    G4Box* preCollimator;
    G4VPhysicalVolume* physiPreCollimator;
    
    G4Box* Collimator;
    G4VPhysicalVolume* physiCollimator;
    
    G4Box* solidAirTube;
    G4Box* solidAirPreTube;
    
    G4VPhysicalVolume* mother;
    G4VPhysicalVolume* physiAirTube;
    G4VPhysicalVolume* physiAirPreTube;
    
    G4VisAttributes* blue;
    G4VisAttributes* gray;
    G4VisAttributes* white;
    G4VisAttributes* red;
    G4VisAttributes* yellow;
    G4VisAttributes* green;
    G4VisAttributes* darkGreen;
    G4VisAttributes* darkOrange3;
    G4VisAttributes* skyBlue;
    
    G4Material* beamLineSupportMaterial;
    G4Material* vacuumZoneMaterial;
    G4Material* layerMonitorChamberMaterial;
    G4Material* layerDefaultMaterial;
    G4Material* FoilMaterial;
    G4Material* internalStructureMaterial;
    G4Material* airgapMaterial;
    G4Material* ScatteringFoilMaterial;
    G4Material* singleTrapMaterial;
    G4Material* CollimatorMaterial;
    G4Material* preCollimatorMaterial;
    G4Material* airTubeMaterial;
    G4Material* airTube2Material;
    G4Material* airTube3Material;
    G4Material* ridgeFilterMaterial;
    
    HadrontherapyDetectorROGeometry* RO;
    
    
};
#endif
