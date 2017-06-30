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

#ifndef PassiveCarbonBeamLine_H
#define PassiveCarbonBeamLine_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4Trd.hh"
#include "PassiveCarbonBeamLineMessenger.hh"

class G4VPhysicalVolume;
class HadrontherapyDetectorConstruction;
class HadrontherapyDetectorROGeometry;
class PassiveCarbonBeamLineMessenger;

class PassiveCarbonBeamLine : public G4VUserDetectorConstruction
{
public:

	PassiveCarbonBeamLine();
	~PassiveCarbonBeamLine();

	G4VPhysicalVolume* Construct();  
	
	void HadrontherapyBeamLineSupport();
	// Definition of the beam line support
    
	void VacuumToAirInterface();
	
	void HadrontherapyBeamMonitoring();
	// Definition of three monitor chambers
	
	void HadrontherapyBeamNozzle();
	// Definition of the beam noozle
	
	void HadrontherapyBeamFinalCollimator();
	// Definition of the final collimator
    
    void HadrontherapyPMMACollimator();
    // Definition of the PMMA collimator
    
    void HadrontherapyRippleFilter();
    // Definition of the ripple filter
	
    void SetInnerRadiusFinalCollimator(G4double);
    // This method allows to change the size of the inner radius of the
    // final collimator
    
    void SetRippleFilterMaterial(G4String);
    // This method allows to change the material
    // of the ripple filter
    
    void SetRippleFilterXPosition(G4double);
    
	// The following methods allow to change parameters
	// of some beam line components

        G4Material* kapton;
        G4VisAttributes* redWire;
        G4VPhysicalVolume* mother;
private:
    static PassiveCarbonBeamLine* instance;
	//passive carbon line dimensions
	void SetDefaultDimensions(); 
	void ConstructPassiveCarbonBeamLine();
	
    PassiveCarbonBeamLineMessenger* PassiveCarbonMessenger;
    G4VPhysicalVolume* physicalTreatmentRoom;
	HadrontherapyDetectorConstruction* hadrontherapyDetectorConstruction;

    ///////////////////////////////////////////////////////////////////////////
    // Definitions of all volumes
    // World (experimental hall)
    G4Box* treatmentRoom;
    G4LogicalVolume* logicTreatmentRoom;
  

    // Beamline support
    G4Box* beamLineSupport;
    G4LogicalVolume* logicBeamLineSupport;
    G4VPhysicalVolume* physiBeamLineSupport;
    
    // Beamline covers
    G4Box* beamLineCover;
    G4LogicalVolume* logicBeamLineCover;
    G4VPhysicalVolume* physiBeamLineCover;
    G4VPhysicalVolume* physiBeamLineCover2;
    
    // Vacuum pipe
    G4Box* vacuumZone;
    G4LogicalVolume* logicVacuumZone;
    G4VPhysicalVolume* physiVacuumZone;
    
    // Scattering foil
    G4Box* firstScatteringFoil;
    G4LogicalVolume* logicFirstScatteringFoil;
    G4VPhysicalVolume* physiFirstScatteringFoil;
    
    // Kapton window
    G4Box* solidKaptonWindow;
    G4LogicalVolume* logicKaptonWindow;
    G4VPhysicalVolume* physiKaptonWindow;

    //Ripple filter
    G4Box* SolidRippleFilter;
    G4LogicalVolume* LogicRippleFilter;
    G4VPhysicalVolume* PhysiRippleFilter;
    
    G4Box* SolidRippleFilterBase;
    G4LogicalVolume* LogicRippleFilterBase;
    G4VPhysicalVolume* PhysiRippleFilterBase;
    
    G4Trd* SolidRippleFilterTrd;
    G4LogicalVolume* LogicRippleFilterTrd;
    G4VPhysicalVolume* PhysiRippleFilterTrd;
    
    // PMMA Collimator
    G4Box* solidPMMACollimatorSupport;
    G4LogicalVolume* logicPMMACollimatorSupport;
    G4VPhysicalVolume* physiPMMACollimatorSupport;
    
    G4Tubs* solidPMMACollimator;
    G4LogicalVolume* logicPMMACollimator;
    G4VPhysicalVolume* physiPMMACollimator;

    // Monitor chamber
    G4Box* solidFirstMonitorLayer1;
    G4LogicalVolume* logicFirstMonitorLayer1;
    G4VPhysicalVolume* physiFirstMonitorLayer1;

    G4Box* solidFirstMonitorLayer2;
    G4LogicalVolume* logicFirstMonitorLayer2;
    G4VPhysicalVolume* physiFirstMonitorLayer2;
    
    G4Box* solidFirstMonitorLayer3;
    G4LogicalVolume* logicFirstMonitorLayer3;
    G4VPhysicalVolume* physiFirstMonitorLayer3;
    
    G4Box* solidFirstMonitorLayer4;
    G4LogicalVolume* logicFirstMonitorLayer4;
    G4VPhysicalVolume* physiFirstMonitorLayer4;
    
    // Final nozzle and collimator
    G4Box* solidNozzleSupport;
    G4LogicalVolume* logicNozzleSupport;
    G4VPhysicalVolume* physiNozzleSupport;
    
    G4Tubs* solidHoleNozzleSupport;
    G4LogicalVolume* logicHoleNozzleSupport;
    G4VPhysicalVolume* physiHoleNozzleSupport;

    G4Tubs* solidBrassTube3;
    G4LogicalVolume* logicBrassTube3;
    G4VPhysicalVolume* physiBrassTube3;
    
    G4Tubs* solidBrassTube2;
    G4LogicalVolume* logicBrassTube2;
    G4VPhysicalVolume* physiBrassTube2;

    G4Tubs* solidBrassTube;
    G4LogicalVolume* logicBrassTube;
    G4VPhysicalVolume* physiBrassTube;
    
    // Final collimator
    G4Tubs* solidFinalCollimator;
    G4VPhysicalVolume* physiFinalCollimator;
    G4LogicalVolume* logicFinalCollimator;
    
    ///////////////////////////////////////////////////////////////////////////
    // Dimensions and positions of all volumes
    // Beamline support
    G4double beamLineSupportXSize;
    G4double beamLineSupportYSize;
    G4double beamLineSupportZSize;
    G4double beamLineSupportXPosition;
    G4double beamLineSupportYPosition;
    G4double beamLineSupportZPosition;
    
    // Beamline covers
    G4double beamLineCoverXSize;
    G4double beamLineCoverYSize;
    G4double beamLineCoverZSize;
    G4double beamLineCoverXPosition;
    G4double beamLineCoverYPosition;
    G4double beamLineCoverZPosition;
    
    // vacuum pipe
	G4double vacuumZoneXSize;
	G4double vacuumZoneYSize;
	G4double vacuumZoneZSize;
	G4double vacuumPipeXPosition;
    
    // Scattering foil
    G4double firstScatteringFoilXSize;
    G4double firstScatteringFoilYSize;
    G4double firstScatteringFoilZSize;
    G4double firstScatteringFoilXPosition;

    // kapton Windows
    G4double kaptonWindowXSize;
    G4double kaptonWindowYSize;
    G4double kaptonWindowZSize;
    G4double kaptonWindowXPosition;
    
    // PMMA Collimator
    G4double PMMACollimatorSupportXSize;
    G4double PMMACollimatorSupportYSize;
    G4double PMMACollimatorSupportZSize;
    G4double PMMACollimatorXPosition;
    G4double innerRadiusPMMACollimator;
    G4double outerRadiusPMMACollimator;
    G4double hightPMMACollimator;
    G4double startAnglePMMACollimator;
    G4double spanningAnglePMMACollimator;
    
    // Monitor chamber
    G4double monitor1XSize;
    G4double monitor2XSize;
    G4double monitor3XSize;
    G4double monitorYSize;
    G4double monitorZSize;
    G4double monitor1XPosition;
    G4double monitor2XPosition;
    G4double monitor4XPosition;
    
    // Final nozzle and collimator
    G4double nozzleSupportXSize;
    G4double nozzleSupportYSize;
    G4double nozzleSupportZSize;
    G4double nozzleSupportXPosition;

    G4double innerRadiusHoleNozzleSupport;
    G4double outerRadiusHoleNozzleSupport;
    G4double hightHoleNozzleSupport;
    G4double startAngleHoleNozzleSupport;
    G4double spanningAngleHoleNozzleSupport;
    
    G4double innerRadiusBrassTube3;
    G4double outerRadiusBrassTube3;
    G4double hightBrassTube3;
    G4double startAngleBrassTube3;
    G4double spanningAngleBrassTube3;
    G4double brassTube3XPosition;
    
    G4double innerRadiusBrassTube2;
    G4double outerRadiusBrassTube2;
    G4double hightBrassTube2;
    G4double startAngleBrassTube2;
    G4double spanningAngleBrassTube2;
    
    G4double innerRadiusBrassTube;
    G4double outerRadiusBrassTube;
    G4double hightBrassTube;
    G4double startAngleBrassTube;
    G4double spanningAngleBrassTube;
    G4double brassTubeXPosition;
    
    // Final collimator
    G4double outerRadiusFinalCollimator;
    G4double innerRadiusFinalCollimator;
    G4double hightFinalCollimator;
    G4double startAngleFinalCollimator;
    G4double spanningAngleFinalCollimator;
    G4double finalCollimatorXPosition;
	
	// Colors
	G4VisAttributes* blue;
	G4VisAttributes* gray;
	G4VisAttributes* white;
	G4VisAttributes* red;
	G4VisAttributes* yellow;
	G4VisAttributes* green;
	G4VisAttributes* darkGreen;
	G4VisAttributes* darkOrange3;
	G4VisAttributes* skyBlue;
    G4VisAttributes* black;
	
    // Elements, compounds and materials
    G4Material *aluminumNist;
    G4Material* copperNistMaterial;
    G4Material* airNist;
    G4Material* kaptonNist;
    G4Material* galacticNist;
    G4Material* PMMANist;
    G4Material* tantalumNist;
    G4Material* brass;
    
    G4Material* beamLineSupportMaterial;
	G4Material* vacuumZoneMaterial;
	G4Material* kaptonWindowMaterial;
    G4Material* firstScatteringFoilMaterial;

	G4Material* layer1MonitorChamberMaterial;
	G4Material* layer2MonitorChamberMaterial;
	G4Material* layer3MonitorChamberMaterial;
	G4Material* layer4MonitorChamberMaterial;
	G4Material* nozzleSupportMaterial;
	G4Material* holeNozzleSupportMaterial;
	G4Material* seconHoleNozzleSupportMaterial;
	G4Material* brassTubeMaterial;
    G4Material* brassTube2Material;
    G4Material* brassTube3Material;
    G4Material* finalCollimatorMaterial;
    G4Material* PMMACollimatorMaterial;
    G4Material* rippleFilterMaterial;
    G4Material* rippleFilterBoxMaterial;
    
	HadrontherapyDetectorROGeometry* RO;
};
#endif

