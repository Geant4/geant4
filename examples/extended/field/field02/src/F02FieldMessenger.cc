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
/// \file field/field02/src/F02FieldMessenger.cc
/// \brief Implementation of the F02FieldMessenger class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F02FieldMessenger.hh"

#include "F02ElectricFieldSetup.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F02FieldMessenger::F02FieldMessenger(F02ElectricFieldSetup* fieldSetup)
 : fElFieldSetup(fieldSetup)
{
  fFieldDir = new G4UIdirectory("/field/");
  fFieldDir->SetGuidance("F02 field tracking control.");

  fStepperCmd = new G4UIcmdWithAnInteger("/field/setStepperType",this);
  fStepperCmd->SetGuidance("Select stepper type for electric field");
  fStepperCmd->SetParameterName("choice",true);
  fStepperCmd->SetDefaultValue(4);
  fStepperCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fUpdateCmd = new G4UIcmdWithoutParameter("/field/update",this);
  fUpdateCmd->SetGuidance("Update calorimeter geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);

  fElFieldZCmd = new G4UIcmdWithADoubleAndUnit("/field/setFieldZ",this);
  fElFieldZCmd->SetGuidance("Define uniform Electric field.");
  fElFieldZCmd->SetGuidance("Electric field will be in Z direction.");
  fElFieldZCmd->SetGuidance("Value of Electric field has to be given in volt/m");
  fElFieldZCmd->SetParameterName("Ez",false,false);
  fElFieldZCmd->SetDefaultUnit("megavolt/m");
  fElFieldZCmd->AvailableForStates(G4State_Idle);

  fElFieldCmd = new G4UIcmdWith3VectorAndUnit("/field/setField",this);
  fElFieldCmd->SetGuidance("Define uniform Electric field.");
  fElFieldCmd->SetGuidance("Value of Electric field has to be given in volt/m");
  fElFieldCmd->SetParameterName("Ex","Ey","Ez",false,false);
  fElFieldCmd->SetDefaultUnit("megavolt/m");
  fElFieldCmd->AvailableForStates(G4State_Idle);

  fMinStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setMinStep",this);
  fMinStepCmd->SetGuidance("Define minimal step");
  fMinStepCmd->SetParameterName("min step",false,false);
  fMinStepCmd->SetDefaultUnit("mm");
  fMinStepCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F02FieldMessenger::~F02FieldMessenger()
{
  delete fStepperCmd;
  delete fElFieldZCmd;
  delete fElFieldCmd;
  delete fMinStepCmd;
  delete fFieldDir;
  delete fUpdateCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F02FieldMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fStepperCmd )
    fElFieldSetup->SetStepperType(fStepperCmd->GetNewIntValue(newValue));
  if( command == fUpdateCmd )
    fElFieldSetup->UpdateIntegrator();
  if( command == fElFieldZCmd )
    fElFieldSetup->SetFieldZValue(fElFieldZCmd->GetNewDoubleValue(newValue));
  if( command == fElFieldCmd )
    fElFieldSetup->SetFieldValue(fElFieldCmd->GetNew3VectorValue(newValue));
  if( command == fMinStepCmd )
    fElFieldSetup->SetMinStep(fMinStepCmd->GetNewDoubleValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
