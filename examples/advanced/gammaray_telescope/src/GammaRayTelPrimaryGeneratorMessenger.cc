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
//      ------------ GammaRayTelPrimaryGeneratorMessenger  ------
//           by G.Santin, F.Longo & R.Giannitrapani (13 nov 2000)
//
// ************************************************************

#include "GammaRayTelPrimaryGeneratorMessenger.hh"
#include "GammaRayTelPrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPrimaryGeneratorMessenger::GammaRayTelPrimaryGeneratorMessenger(GammaRayTelPrimaryGeneratorAction *GammaRayTelGun) : GammaRayTelAction(GammaRayTelGun) {
	rndmCmd = new G4UIcmdWithAString("/gun/random", this);
	rndmCmd->SetGuidance("Shoot randomly the incident particle.");
	rndmCmd->SetGuidance("Choice : on(default), off");
	rndmCmd->SetParameterName("choice", true);
	rndmCmd->SetDefaultValue("on");
	rndmCmd->SetCandidates("on off");
	rndmCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	sourceTypeCmd = new G4UIcmdWithAnInteger("/gun/sourceType", this);
	sourceTypeCmd->SetGuidance("Select the type of incident flux.");
	sourceTypeCmd->SetGuidance("Choice : 0(default), 1(isotropic), 2(wide parallel beam)");
	sourceTypeCmd->SetParameterName("choice", true);
	sourceTypeCmd->SetDefaultValue((G4int) 0);
	sourceTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	vertexRadiusCmd = new G4UIcmdWithADoubleAndUnit("/gun/vertexRadius", this);
	vertexRadiusCmd->SetGuidance("Radius (and unit) of sphere for vertices of incident flux.");
	vertexRadiusCmd->SetParameterName("choice", true);
	vertexRadiusCmd->SetDefaultValue((G4double) 1. * cm);
	vertexRadiusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	spectrumTypeCmd = new G4UIcmdWithAnInteger("/gun/spectrumType", this);
	spectrumTypeCmd->SetGuidance("Select the type of incident spectrum.");
	spectrumTypeCmd->SetGuidance("Choice : 0(default), 1(), 2(E^{-gamma}), 3()");
	spectrumTypeCmd->SetParameterName("choice", true);
	spectrumTypeCmd->SetDefaultValue((G4int) 0);
	spectrumTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	sourceGenCmd = new G4UIcmdWithABool("/gun/sourceGen", this);
	sourceGenCmd->SetGuidance("Select the native Generation");
	sourceGenCmd->SetGuidance("  Choice : true(native), false(GPS)");
	sourceGenCmd->SetParameterName("choice", true);
	sourceGenCmd->SetDefaultValue((G4bool) true);
	sourceGenCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPrimaryGeneratorMessenger::~GammaRayTelPrimaryGeneratorMessenger() {
	delete rndmCmd;
	delete sourceTypeCmd;
	delete vertexRadiusCmd;
	delete spectrumTypeCmd;
	delete sourceGenCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPrimaryGeneratorMessenger::SetNewValue(G4UIcommand *command, G4String newValue) {
	if (command == rndmCmd) {
		GammaRayTelAction->SetRndmFlag(newValue);
	} else if (command == sourceTypeCmd) {
		GammaRayTelAction->SetSourceType(sourceTypeCmd->GetNewIntValue(newValue));
	} else if (command == vertexRadiusCmd) {
		GammaRayTelAction->SetVertexRadius(vertexRadiusCmd->GetNewDoubleValue(newValue));
	} else if (command == spectrumTypeCmd) {
		GammaRayTelAction->SetSpectrumType(spectrumTypeCmd->GetNewIntValue(newValue));
	} else if (command == sourceGenCmd) {
		GammaRayTelAction->SetSourceGen(sourceGenCmd->GetNewBoolValue(newValue));
	}
}
