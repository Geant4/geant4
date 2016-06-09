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
// $Id: PhysicsListMessenger.cc,v 1.4 2006/08/11 14:38:11 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
:pPhysicsList(pPhys)
{   
  gammaCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutGamma",this);  
  gammaCutCmd->SetGuidance("Set gamma cut.");
  gammaCutCmd->SetParameterName("Gcut",false);
  gammaCutCmd->SetUnitCategory("Length");
  gammaCutCmd->SetRange("Gcut>0.0");
  gammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  electCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutEl",this);  
  electCutCmd->SetGuidance("Set electron cut.");
  electCutCmd->SetParameterName("Ecut",false);
  electCutCmd->SetUnitCategory("Length");
  electCutCmd->SetRange("Ecut>0.0");
  electCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  posCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutPos",this);
  posCutCmd->SetGuidance("Set positron cut.");
  posCutCmd->SetParameterName("Pcut",false);
  posCutCmd->SetUnitCategory("Length");
  posCutCmd->SetRange("Pcut>0.0");
  posCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  allCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutsAll",this);
  allCutCmd->SetGuidance("Set cut for all.");
  allCutCmd->SetParameterName("cut",false);
  allCutCmd->SetUnitCategory("Length");
  allCutCmd->SetRange("cut>0.0");
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
  delete allCutCmd;
  delete pListCmd;
  delete listCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if( command == gammaCutCmd )
    pPhysicsList->SetCutForGamma(gammaCutCmd->GetNewDoubleValue(newValue));

  if( command == electCutCmd )
    pPhysicsList->SetCutForElectron(electCutCmd->GetNewDoubleValue(newValue));

  if( command == posCutCmd )
   pPhysicsList->SetCutForPositron(posCutCmd->GetNewDoubleValue(newValue));

  if( command == allCutCmd )
    {
      G4double cut = allCutCmd->GetNewDoubleValue(newValue);
      pPhysicsList->SetCutForGamma(cut);
      pPhysicsList->SetCutForElectron(cut);
      pPhysicsList->SetCutForPositron(cut);
    }

  if( command == pListCmd ) {
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
  }

  if( command == listCmd )
    pPhysicsList->List();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
