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
// $Id: PhysicsListMessenger.cc,v 1.7 2010-01-13 15:53:44 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
/////////////////////////////////////////////////////////////////////////
//
// PhysicsListMessenger
//
// Created: 31.01.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
//
// 

#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
:pPhysicsList(pPhys)
{   
  gammaCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutGamma",this);  
  gammaCutCmd->SetGuidance("Set gamma cut.");
  gammaCutCmd->SetParameterName("Gcut",false);
  gammaCutCmd->SetUnitCategory("Length");
  gammaCutCmd->SetRange("Gcut>=0.0");
  gammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  electCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutEl",this);  
  electCutCmd->SetGuidance("Set electron cut.");
  electCutCmd->SetParameterName("Ecut",false);
  electCutCmd->SetUnitCategory("Length");
  electCutCmd->SetRange("Ecut>=0.0");
  electCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  posCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutPos",this);
  posCutCmd->SetGuidance("Set positron cut.");
  posCutCmd->SetParameterName("Pcut",false);
  posCutCmd->SetUnitCategory("Length");
  posCutCmd->SetRange("Pcut>=0.0");
  posCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  pCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutProt",this);
  pCutCmd->SetGuidance("Set proton cut.");
  pCutCmd->SetParameterName("ProtCut",false);
  pCutCmd->SetUnitCategory("Length");
  pCutCmd->SetRange("ProtCut>=0.0");
  pCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  allCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutsAll",this);
  allCutCmd->SetGuidance("Set cut for all.");
  allCutCmd->SetParameterName("cut",false);
  allCutCmd->SetUnitCategory("Length");
  allCutCmd->SetRange("cut>=0.0");
  allCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  pListCmd = new G4UIcmdWithAString("/testhadr/Physics",this);
  pListCmd->SetGuidance("Add modula physics list.");
  pListCmd->SetParameterName("PList",false);
  pListCmd->AvailableForStates(G4State_PreInit);

  listCmd = new G4UIcmdWithoutParameter("/testhadr/ListPhysics",this);
  listCmd->SetGuidance("Available Physics Lists");
  listCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete gammaCutCmd;
  delete electCutCmd;
  delete posCutCmd;
  delete pCutCmd;
  delete allCutCmd;
  delete pListCmd;
  delete listCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if( command == gammaCutCmd ) {
    if(pPhysicsList) {
      pPhysicsList->SetCutForGamma(gammaCutCmd->GetNewDoubleValue(newValue));
    } else {
      UI->ApplyCommand("/run/setCutForAGivenParticle gamma " + newValue);
    }

  } else if( command == electCutCmd ) {
    if(pPhysicsList) {
      pPhysicsList->SetCutForElectron(electCutCmd->GetNewDoubleValue(newValue));
    } else {
      UI->ApplyCommand("/run/setCutForAGivenParticle e- " + newValue);
    }

  } else if( command == posCutCmd ) {
    if(pPhysicsList) {
      pPhysicsList->SetCutForPositron(posCutCmd->GetNewDoubleValue(newValue));
    } else {
      UI->ApplyCommand("/run/setCutForAGivenParticle e+ " + newValue);
    }

  } else if( command == pCutCmd ) {
    if(pPhysicsList) {
      pPhysicsList->SetCutForProton(pCutCmd->GetNewDoubleValue(newValue));
    } else {
      UI->ApplyCommand("/run/setCutForAGivenParticle proton " + newValue);
    }

  } else if( command == allCutCmd ) {

    if(pPhysicsList) {
      G4double cut = allCutCmd->GetNewDoubleValue(newValue);
      pPhysicsList->SetCutForGamma(cut);
      pPhysicsList->SetCutForElectron(cut);
      pPhysicsList->SetCutForPositron(cut);
      pPhysicsList->SetCutForProton(cut);
    } else {
      UI->ApplyCommand("/run/setCut " + newValue);
    }

  } else if( command == pListCmd ) {
    if(pPhysicsList) {
      G4String name = newValue;
      if(name == "PHYSLIST") {
	char* path = getenv(name);
	if (path) name = G4String(path);
	else {
	  G4cout << "### PhysicsListMessenger WARNING: "
		 << " environment variable PHYSLIST is not defined"
		 << G4endl;
	  return; 
	}
      }
      pPhysicsList->AddPhysicsList(name);
    } else {
      G4cout << "### PhysicsListMessenger WARNING: "
	     << " /testhadr/Physics UI command is not available "
	     << "for reference Physics List" << G4endl;
    }

  } else if( command == listCmd ) {
    if(pPhysicsList) {
      pPhysicsList->List();
    } else { 
      G4cout << "### PhysicsListMessenger WARNING: "
	     << " /testhadr/ListPhysics UI command is not available "
	     << "for reference Physics List" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
