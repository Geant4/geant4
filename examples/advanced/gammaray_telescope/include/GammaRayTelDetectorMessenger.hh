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
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDetectorMessenger  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#ifndef GammaRayTelDetectorMessenger_h
#define GammaRayTelDetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>

class GammaRayTelDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class GammaRayTelDetectorMessenger: public G4UImessenger {
public:
	explicit GammaRayTelDetectorMessenger(GammaRayTelDetectorConstruction*);

	~GammaRayTelDetectorMessenger() override;

	void SetNewValue(G4UIcommand *command, G4String newValue) override;

private:
	GammaRayTelDetectorConstruction *detector;

	G4UIdirectory *directory;

	// Converter

	G4UIcmdWithAString *converterMaterialCmd;
	G4UIcmdWithADoubleAndUnit *converterThicknessCmd;

	// Silicon Tile

	G4UIcmdWithADoubleAndUnit *siliconThicknessCmd;
	G4UIcmdWithADoubleAndUnit *siliconTileXYCmd;
	G4UIcmdWithAnInteger *numberOfSiTilesCmd;
	G4UIcmdWithADoubleAndUnit *siliconPitchCmd;

	// Tracker (TKR)

	G4UIcmdWithAnInteger *numberOfTKRLayersCmd;
	G4UIcmdWithADoubleAndUnit *layerDistanceCmd;
	G4UIcmdWithADoubleAndUnit *viewsDistanceCmd;

	// Calorimeter (CAL)

	G4UIcmdWithADoubleAndUnit *calThicknessCmd;
	G4UIcmdWithAnInteger *numberOfCALBarsCmd;
	G4UIcmdWithAnInteger *numberOfCALLayersCmd;

	// Anticoincidence (ACD)

	G4UIcmdWithADoubleAndUnit *acdThicknessCmd;

	// Total

	G4UIcmdWithADoubleAndUnit *magneticFieldCmd;
	G4UIcmdWithoutParameter *updateCmd;
};
#endif
