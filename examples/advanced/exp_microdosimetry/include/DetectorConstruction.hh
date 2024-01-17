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
// Authors: Susanna Guatelli and Francesco Romano
// susanna@uow.edu.au, francesco.romano@ct.infn.it

// Modified by Jacopo Magini: j.magini@surrey.ac.uk

#ifndef DetectorConstruction_H 
#define DetectorConstruction_H 1
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "AnalysisManager.hh"
#include "DetectorMessenger.hh"
#include "G4NistManager.hh"

class G4VPhysicalVolume;
class DetectorMessenger;
class G4LogicalVolume;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction(AnalysisManager* analysis, DetectorMessenger* messenger);
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();

    void ConstructSDandField();

private:
	AnalysisManager* analysis;
	DetectorMessenger* messenger;
	
	G4String detectorType;
	G4double detectorSizeWidth;
	G4double detectorSizeThickness;
	G4double secondStageSizeDim;
	G4double secondStageSizeThickness;
	G4bool usingWaterPhantom;
	G4double detectorPositionDepth;
	
	// Methods called by Construct() depending on the chosen setup
	void ConstructWorldWithWaterPhantom();
	void ConstructVacuumWorld();
	void ConstructDiamondDetector();
	void ConstructMicroDiamondDetector();
	void ConstructSiliconDetector();
	void ConstructSiliconBridgeDetector();
	void ConstructDiamondTelescope();
	
	// store these variable across various Construct*() above
	G4VPhysicalVolume* physical_world;
	G4LogicalVolume* logical_motherVolumeForDetector;
	G4Material* materialOfMotherVolume;
	
	G4NistManager* nistMan;
};
#endif
