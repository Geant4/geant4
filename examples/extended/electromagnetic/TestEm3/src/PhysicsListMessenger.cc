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
/// \file electromagnetic/TestEm5/src/PhysicsListMessenger.cc
/// \brief Implementation of the PhysicsListMessenger class
//
// $Id: PhysicsListMessenger.cc 109096 2018-03-26 14:46:51Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
  :G4UImessenger(),fPhysicsList(pPhys),fMaxChargedStep(DBL_MAX)
{
  fPhysDir = new G4UIdirectory("/testem/phys/");
  fPhysDir->SetGuidance("physics list commands");

  fListCmd = new G4UIcmdWithAString("/testem/phys/addPhysics",this);  
  fListCmd->SetGuidance("Add modula physics list.");
  fListCmd->SetParameterName("PList",false);
  fListCmd->AvailableForStates(G4State_PreInit);
  fListCmd->SetToBeBroadcasted(false);      

  fStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/testem/stepMax",this);
  fStepMaxCmd->SetGuidance("Set max allowed step length");
  fStepMaxCmd->SetParameterName("mxStep",false);
  fStepMaxCmd->SetRange("mxStep>0.");
  fStepMaxCmd->SetUnitCategory("Length");
  fStepMaxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete fListCmd;
  delete fPhysDir;
  delete fStepMaxCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == fListCmd )
    { fPhysicsList->AddPhysicsList(newValue); }
  if (command == fStepMaxCmd)
    { fMaxChargedStep = fStepMaxCmd->GetNewDoubleValue(newValue); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
