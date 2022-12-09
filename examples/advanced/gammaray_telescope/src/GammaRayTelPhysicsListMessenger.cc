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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//
// Code based on the hadrontherapy && radioprotection advanced example

#include "GammaRayTelPhysicsListMessenger.hh"
#include "GammaRayTelPhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPhysicsListMessenger::GammaRayTelPhysicsListMessenger(GammaRayTelPhysicsList *pPhys) : pPhysicsList(pPhys) {
	physDir = new G4UIdirectory("/physics/");
	physDir->SetGuidance("Commands to activate physics models and set cuts");

	gammaCutCmd = new G4UIcmdWithADoubleAndUnit("/physics/setGCut", this);
	gammaCutCmd->SetGuidance("Set gamma cut.");
	gammaCutCmd->SetParameterName("Gcut", false);
	gammaCutCmd->SetUnitCategory("Length");
	gammaCutCmd->SetRange("Gcut>0.0");
	gammaCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	electronCutCmd = new G4UIcmdWithADoubleAndUnit("/physics/setECut", this);
	electronCutCmd->SetGuidance("Set electron cut.");
	electronCutCmd->SetParameterName("Ecut", false);
	electronCutCmd->SetUnitCategory("Length");
	electronCutCmd->SetRange("Ecut>0.0");
	electronCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	protonCutCmd = new G4UIcmdWithADoubleAndUnit("/physics/setPCut", this);
	protonCutCmd->SetGuidance("Set positron cut.");
	protonCutCmd->SetParameterName("Pcut", false);
	protonCutCmd->SetUnitCategory("Length");
	protonCutCmd->SetRange("Pcut>0.0");
	protonCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	allCutCmd = new G4UIcmdWithADoubleAndUnit("/physics/setCuts", this);
	allCutCmd->SetGuidance("Set cut for all.");
	allCutCmd->SetParameterName("cut", false);
	allCutCmd->SetUnitCategory("Length");
	allCutCmd->SetRange("cut>0.0");
	allCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	pListCmd = new G4UIcmdWithAString("/physics/addPhysics", this);
	pListCmd->SetGuidance("Add physics list.");
	pListCmd->SetParameterName("PList", false);
	pListCmd->AvailableForStates(G4State_PreInit);

	packageListCmd = new G4UIcmdWithAString("/physics/addPackage", this);
	packageListCmd->SetGuidance("Add physics package.");
	packageListCmd->SetParameterName("package", false);
	packageListCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPhysicsListMessenger::~GammaRayTelPhysicsListMessenger() {
	delete gammaCutCmd;
	delete electronCutCmd;
	delete protonCutCmd;
	delete allCutCmd;
	delete pListCmd;
	delete physDir;
	delete packageListCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPhysicsListMessenger::SetNewValue(G4UIcommand *command, G4String newValue) {
	if (command == gammaCutCmd) {
		pPhysicsList->SetCutForGamma(gammaCutCmd->GetNewDoubleValue(newValue));
	} else if (command == electronCutCmd) {
		pPhysicsList->SetCutForElectron(electronCutCmd->GetNewDoubleValue(newValue));
	} else if (command == protonCutCmd) {
		pPhysicsList->SetCutForPositron(protonCutCmd->GetNewDoubleValue(newValue));
	} else if (command == allCutCmd) {
		G4double cut = allCutCmd->GetNewDoubleValue(newValue);
		pPhysicsList->SetCutForGamma(cut);
		pPhysicsList->SetCutForElectron(cut);
		pPhysicsList->SetCutForPositron(cut);
	} else if (command == pListCmd) {
		pPhysicsList->AddPhysicsList(newValue);
	} else if (command == packageListCmd) {
		pPhysicsList->AddPackage(newValue);
	}
}
