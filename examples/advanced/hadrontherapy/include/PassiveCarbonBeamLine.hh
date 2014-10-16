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
// PassiveCarbonBeamLine.cc; 
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy

#ifndef PassiveCarbonBeamLine_H
#define PassiveCarbonBeamLine_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

class G4VPhysicalVolume;
class HadrontherapyDetectorConstruction;
class HadrontherapyModulator;

class PassiveCarbonBeamLine : public G4VUserDetectorConstruction
{
public:

	PassiveCarbonBeamLine();
	~PassiveCarbonBeamLine();

	G4VPhysicalVolume* Construct();  
	
	void HadrontherapyBeamLineSupport();
	// Definition of the beam line support
	
	// Simulation of the scattering system for the 
	// passive spread of the beam 
	void ScatteringSystem();
	
	void VacuumToAirInterface();
	// Definition of the first scattering foil, 
	// of the Kapton window, of the stopper 
	
	void HadrontherapyBeamMonitoring();
	// Definition of three monitor chambers
	
	void HadrontherapyBeamNozzle();
	// Definition of the beam noozle
	
	void HadrontherapyBeamFinalCollimator();
	// Definition of the final collimator
	
	// The following methods allow to change parameters
	// of some beam line components

        G4Material* kapton;
        G4VisAttributes* redWire;
        G4VPhysicalVolume* mother;
	G4double firstScatteringFoilXPosition;
	G4double firstScatteringFoilYPosition;
	G4double firstScatteringFoilZPosition;
	
private:
	//passive proton line dimensions
	void SetDefaultDimensions(); 
	void ConstructPassiveCarbonBeamLine();
	
	//PassiveCarbonBeamLineMessenger* passiveMessenger;  
	G4VPhysicalVolume* physicalTreatmentRoom;
	HadrontherapyDetectorConstruction* hadrontherapyDetectorConstruction; 
  //FaradayCup *FC;  

	G4double vacuumZoneXSize;
	G4double vacuumZoneYSize;
	G4double vacuumZoneZSize;
	G4double vacuumZoneXPosition;
	
	G4double kaptonWindowXSize;
	G4double kaptonWindowYSize;
	G4double kaptonWindowZSize;
	G4double kaptonWindowXPosition;
	
	G4VPhysicalVolume* physiBeamLineSupport; 
	G4VPhysicalVolume* physiBeamLineCover; 
	G4VPhysicalVolume* physiBeamLineCover2;
	G4VPhysicalVolume* physiKaptonWindow;
	
	G4Tubs* solidStopper;
	G4VPhysicalVolume* physiStopper; 
	G4LogicalVolume* logicStopper;
	G4double innerRadiusStopper;
	G4double heightStopper;
	G4double startAngleStopper;
	G4double spanningAngleStopper;
	G4double stopperXPosition;
	G4double stopperYPosition;
	G4double stopperZPosition;
	G4double outerRadiusStopper;
	
	// First scattering foil coupled with the stopper
	G4Box* firstScatteringFoil;  
	G4VPhysicalVolume* physiFirstScatteringFoil;  
	G4double firstScatteringFoilXSize;
	G4double firstScatteringFoilYSize;
	G4double firstScatteringFoilZSize;
	
	// Scattering foil coupled with the stopper
	G4Box* secondScatteringFoil;  
	G4VPhysicalVolume* physiSecondScatteringFoil;  
	G4double secondScatteringFoilXSize;
	G4double secondScatteringFoilYSize;
	G4double secondScatteringFoilZSize;
	G4double secondScatteringFoilXPosition;
	G4double secondScatteringFoilYPosition;
	G4double secondScatteringFoilZPosition;	
	
	
	G4double innerRadiusFinalCollimator;
	
	G4VPhysicalVolume* physiFirstMonitorLayer1;
	G4VPhysicalVolume* physiFirstMonitorLayer2;
	G4VPhysicalVolume* physiFirstMonitorLayer3;
	G4VPhysicalVolume* physiFirstMonitorLayer4;
	G4VPhysicalVolume* physiNozzleSupport;
	G4VPhysicalVolume* physiHoleNozzleSupport; 
	G4VPhysicalVolume* physiNozzleSupportHole;
	G4VPhysicalVolume* physiSecondHoleNozzleSupport;
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
	
	G4Material* beamLineSupportMaterial;
	G4Material* vacuumZoneMaterial;
	G4Material* firstScatteringFoilMaterial;	
	G4Material* kaptonWindowMaterial;
	G4Material* stopperMaterial;
	G4Material* secondScatteringFoilMaterial;
	
	G4Material* layer1MonitorChamberMaterial;
	G4Material* layer2MonitorChamberMaterial;
	G4Material* layer3MonitorChamberMaterial;
	G4Material* layer4MonitorChamberMaterial;
	G4Material* nozzleSupportMaterial;
	G4Material* holeNozzleSupportMaterial;
	G4Material* seconHoleNozzleSupportMaterial;
	G4Material* finalCollimatorMaterial;
	
	HadrontherapyDetectorROGeometry* RO;
};
#endif

