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
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDetectorMessenger ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelDetectorMessenger.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDetectorMessenger::GammaRayTelDetectorMessenger(GammaRayTelDetectorConstruction *GammaRayTelDet) : detector(GammaRayTelDet) {
	directory = new G4UIdirectory("/payload/");
	directory->SetGuidance("GammaRayTel payload control.");

	// converter material command

	converterMaterialCmd = new G4UIcmdWithAString("/payload/setConvMat", this);
	converterMaterialCmd->SetGuidance("Select the material of the converter.");
	converterMaterialCmd->SetParameterName("choice", false);
	converterMaterialCmd->AvailableForStates(G4State_Idle);

	// converter thickness command

	converterThicknessCmd = new G4UIcmdWithADoubleAndUnit("/payload/setConvThick", this);
	converterThicknessCmd->SetGuidance("Set the thickness of the converter.");
	converterThicknessCmd->SetParameterName("Size", false);
	converterThicknessCmd->SetRange("Size>=0.");
	converterThicknessCmd->SetUnitCategory("Length");
	converterThicknessCmd->AvailableForStates(G4State_Idle);

	// tracker silicon thickness command

	siliconThicknessCmd = new G4UIcmdWithADoubleAndUnit("/payload/setSiThick", this);
	siliconThicknessCmd->SetGuidance("Set the thickness of the silicon.");
	siliconThicknessCmd->SetParameterName("Size", false);
	siliconThicknessCmd->SetRange("Size>=0.");
	siliconThicknessCmd->SetUnitCategory("Length");
	siliconThicknessCmd->AvailableForStates(G4State_Idle);

	// tracker silicon pitch command

	siliconPitchCmd = new G4UIcmdWithADoubleAndUnit("/payload/setSiPitch", this);
	siliconPitchCmd->SetGuidance("Set the pitch of silicon strips.");
	siliconPitchCmd->SetParameterName("Size", false);
	siliconPitchCmd->SetRange("Size>=0.");
	siliconPitchCmd->SetUnitCategory("Length");
	siliconPitchCmd->AvailableForStates(G4State_Idle);

	// tracker silicon tile size command

	siliconTileXYCmd = new G4UIcmdWithADoubleAndUnit("/payload/setSiTileXY", this);
	siliconTileXYCmd->SetGuidance("Set XY dimensions of a silicon tile.");
	siliconTileXYCmd->SetParameterName("Size", false);
	siliconTileXYCmd->SetRange("Size>=0.");
	siliconTileXYCmd->SetUnitCategory("Length");
	siliconTileXYCmd->AvailableForStates(G4State_Idle);

	// tracker number of silicon tiles

	numberOfSiTilesCmd = new G4UIcmdWithAnInteger("/payload/setNbOfSiTiles", this);
	numberOfSiTilesCmd->SetGuidance("Set the number of silicon tiles.");
	numberOfSiTilesCmd->SetParameterName("NbSiTiles", false);
	numberOfSiTilesCmd->SetRange("NbSiTiles>0 && NbSiTiles<100");
	numberOfSiTilesCmd->AvailableForStates(G4State_Idle);

	// tracker number of silicon layers

	numberOfTKRLayersCmd = new G4UIcmdWithAnInteger("/payload/setNbOfTKRLayers", this);
	numberOfTKRLayersCmd->SetGuidance("Set the number of TKR layers.");
	numberOfTKRLayersCmd->SetParameterName("NbTKRLayers", false);
	numberOfTKRLayersCmd->SetRange("NbTKRLayers>0 && NbTKRLayers<30");
	numberOfTKRLayersCmd->AvailableForStates(G4State_Idle);

	// tracker layer distance

	layerDistanceCmd = new G4UIcmdWithADoubleAndUnit("/payload/setLayerDistance", this);
	layerDistanceCmd->SetGuidance("Set the distance between two layers.");
	layerDistanceCmd->SetParameterName("Size", false);
	layerDistanceCmd->SetRange("Size>=0.");
	layerDistanceCmd->SetUnitCategory("Length");
	layerDistanceCmd->AvailableForStates(G4State_Idle);

	// tracker views distance

	viewsDistanceCmd = new G4UIcmdWithADoubleAndUnit("/payload/setViewsDistance", this);
	viewsDistanceCmd->SetGuidance("Set the distance between X and Y views.");
	viewsDistanceCmd->SetParameterName("Size", false);
	viewsDistanceCmd->SetRange("Size>=0.");
	viewsDistanceCmd->SetUnitCategory("Length");
	viewsDistanceCmd->AvailableForStates(G4State_Idle);

	// calorimeter detector thickness

	calThicknessCmd = new G4UIcmdWithADoubleAndUnit("/payload/setCALThick", this);
	calThicknessCmd->SetGuidance("Set the thickness of CAL detectors.");
	calThicknessCmd->SetParameterName("Size", false);
	calThicknessCmd->SetRange("Size>=0.");
	calThicknessCmd->SetUnitCategory("Length");
	calThicknessCmd->AvailableForStates(G4State_Idle);

	// calorimeter, number of detectors

	numberOfCALBarsCmd = new G4UIcmdWithAnInteger("/payload/setNbOfCALBars", this);
	numberOfCALBarsCmd->SetGuidance("Set the number of CsI bars.");
	numberOfCALBarsCmd->SetParameterName("NbSiTiles", false);
	numberOfCALBarsCmd->SetRange("NbSiTiles>0 && NbSiTiles<100");
	numberOfCALBarsCmd->AvailableForStates(G4State_Idle);

	// calorimeter, number of layers

	numberOfCALLayersCmd = new G4UIcmdWithAnInteger("/payload/setNbOfCALLayers", this);
	numberOfCALLayersCmd->SetGuidance("Set the number of CAL layers.");
	numberOfCALLayersCmd->SetParameterName("NbCALLayers", false);
	numberOfCALLayersCmd->SetRange("NbCALLayers>0 && NbCALLayers<16");
	numberOfCALLayersCmd->AvailableForStates(G4State_Idle);

	// anticoincidence detector thickness

	acdThicknessCmd = new G4UIcmdWithADoubleAndUnit("/payload/setACDThick", this);
	acdThicknessCmd->SetGuidance("Set the thickness of ACD detectors");
	acdThicknessCmd->SetParameterName("Size", false);
	acdThicknessCmd->SetRange("Size>=0.");
	acdThicknessCmd->SetUnitCategory("Length");
	acdThicknessCmd->AvailableForStates(G4State_Idle);

	// update payload geometry

	updateCmd = new G4UIcmdWithoutParameter("/payload/update", this);
	updateCmd->SetGuidance("Update the geometry of the payload.");
	updateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
	updateCmd->SetGuidance("if you changed geometrical value(s).");
	updateCmd->AvailableForStates(G4State_Idle);

	// magnetic field

	magneticFieldCmd = new G4UIcmdWithADoubleAndUnit("/payload/setField", this);
	magneticFieldCmd->SetGuidance("Define the magnetic field.");
	magneticFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
	magneticFieldCmd->SetParameterName("Bz", false);
	magneticFieldCmd->SetUnitCategory("Magnetic flux density");
	magneticFieldCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDetectorMessenger::~GammaRayTelDetectorMessenger() {
	delete converterMaterialCmd;
	delete converterThicknessCmd;
	delete numberOfSiTilesCmd;
	delete numberOfTKRLayersCmd;
	delete siliconTileXYCmd;
	delete siliconPitchCmd;
	delete siliconThicknessCmd;
	delete layerDistanceCmd;
	delete viewsDistanceCmd;
	delete acdThicknessCmd;
	delete numberOfCALLayersCmd;
	delete numberOfCALBarsCmd;
	delete calThicknessCmd;
	delete updateCmd;
	delete magneticFieldCmd;
	delete directory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue) {
	// converter

	if (command == converterMaterialCmd) {
		detector->SetConverterMaterial(newValue);
	} else if (command == converterThicknessCmd) {
		detector->SetConverterThickness(converterThicknessCmd->GetNewDoubleValue(newValue));
    } else

    // tracker (TKR)

    if (command == siliconTileXYCmd) {
		detector->SetTKRTileSizeXY(siliconTileXYCmd->GetNewDoubleValue(newValue));
	} else if (command == siliconPitchCmd) {
		detector->SetTKRSiliconPitch(siliconPitchCmd->GetNewDoubleValue(newValue));
	} else if (command == siliconThicknessCmd) {
		detector->SetTKRSiliconThickness(siliconThicknessCmd->GetNewDoubleValue(newValue));
	} else if (command == numberOfSiTilesCmd) {
		detector->SetNbOfTKRTiles(numberOfSiTilesCmd->GetNewIntValue(newValue));
	} else if (command == numberOfTKRLayersCmd) {
		detector->SetNbOfTKRLayers(numberOfTKRLayersCmd->GetNewIntValue(newValue));
	} else if (command == layerDistanceCmd) {
		detector->SetTKRLayerDistance(layerDistanceCmd->GetNewDoubleValue(newValue));
	} else if (command == viewsDistanceCmd) {
		detector->SetTKRViewsDistance(viewsDistanceCmd->GetNewDoubleValue(newValue));
	} else

	// calorimeter (CAL)

	if (command == numberOfCALLayersCmd) {
		detector->SetNbOfCALLayers(numberOfCALLayersCmd->GetNewIntValue(newValue));
	} else if (command == numberOfCALBarsCmd) {
		detector->SetNbOfCALBars(numberOfCALBarsCmd->GetNewIntValue(newValue));
	} else if (command == calThicknessCmd) {
		detector->SetCALBarThickness(calThicknessCmd->GetNewDoubleValue(newValue));
	} else

	// anticoincidence (ACD)

	if (command == acdThicknessCmd) {
		detector->SetACDThickness(acdThicknessCmd->GetNewDoubleValue(newValue));
	} else if (command == updateCmd) {
		detector->UpdateGeometry();
	} else if (command == magneticFieldCmd) {
		detector->SetMagField(magneticFieldCmd->GetNewDoubleValue(newValue));
	}
}
