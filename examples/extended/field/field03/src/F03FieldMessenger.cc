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
/// \file field/field03/src/F03FieldMessenger.cc
/// \brief Implementation of the F03FieldMessenger class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F03FieldMessenger.hh"

#include "F03FieldSetup.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F03FieldMessenger::F03FieldMessenger(F03FieldSetup* fieldSetup)
 : G4UImessenger(),
   fEMfieldSetup(fieldSetup),
   fFieldDir(0),
   fStepperCmd(0),
   fMagFieldZCmd(0),
   fMagFieldCmd(0),
   fLocalMagFieldCmd(0),
   fMinStepCmd(0),
   fUpdateCmd(0)
{
  fFieldDir = new G4UIdirectory("/field/");
  fFieldDir->SetGuidance("F03 field tracking control.");

  fStepperCmd = new G4UIcmdWithAnInteger("/field/setStepperType",this);
  fStepperCmd->SetGuidance("Select stepper type for magnetic field");
  fStepperCmd->SetParameterName("choice",true);
  fStepperCmd->SetDefaultValue(4);
  fStepperCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fUpdateCmd = new G4UIcmdWithoutParameter("/field/update",this);
  fUpdateCmd->SetGuidance("Update calorimeter geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);
 
  fMagFieldZCmd = new G4UIcmdWithADoubleAndUnit("/field/setFieldZ",this);
  fMagFieldZCmd->SetGuidance("Define global magnetic field.");
  fMagFieldZCmd->SetGuidance("Global magnetic field will be in Z direction.");
  fMagFieldZCmd->SetParameterName("Bz",false,false);
  fMagFieldZCmd->SetDefaultUnit("tesla");
  fMagFieldZCmd->AvailableForStates(G4State_Idle);
 
  fMagFieldCmd = new G4UIcmdWith3VectorAndUnit("/field/setField",this);
  fMagFieldCmd->SetGuidance("Define global magnetic field.");
  fMagFieldCmd->SetParameterName("Bx","By","Bz",false,false);
  fMagFieldCmd->SetDefaultUnit("tesla");
  fMagFieldCmd->AvailableForStates(G4State_Idle);
 
  fLocalMagFieldCmd = new G4UIcmdWith3VectorAndUnit("/field/setLocalField",this);
  fLocalMagFieldCmd->SetGuidance("Define local magnetic field.");
  fLocalMagFieldCmd->SetParameterName("Blx","Bly","Blz",false,false);
  fLocalMagFieldCmd->SetDefaultUnit("tesla");
  fLocalMagFieldCmd->AvailableForStates(G4State_Idle);
 
  fMinStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setMinStep",this);
  fMinStepCmd->SetGuidance("Define minimal step");
  fMinStepCmd->SetGuidance("Magnetic field will be in Z direction.");
  fMinStepCmd->SetParameterName("min step",false,false);
  fMinStepCmd->SetDefaultUnit("mm");
  fMinStepCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F03FieldMessenger::~F03FieldMessenger()
{
  delete fStepperCmd;
  delete fMagFieldZCmd;
  delete fMagFieldCmd;
  delete fLocalMagFieldCmd;
  delete fMinStepCmd;
  delete fFieldDir;
  delete fUpdateCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F03FieldMessenger::SetNewValue( G4UIcommand* command, G4String newValue)
{
  if( command == fStepperCmd )
    fEMfieldSetup->SetStepperType(fStepperCmd->GetNewIntValue(newValue));
  if( command == fUpdateCmd )
    fEMfieldSetup->UpdateField();
  if( command == fMagFieldZCmd )
  {
    fEMfieldSetup->SetFieldZValue(fMagFieldZCmd->GetNewDoubleValue(newValue));
    // Check the value
    G4cout << "Set global field value to " <<
      fEMfieldSetup->GetGlobalFieldValue() / gauss << " Gauss " << G4endl;
  }
  if( command == fMagFieldCmd )
  {
    fEMfieldSetup->SetFieldValue(fMagFieldCmd->GetNew3VectorValue(newValue));
    // Check the value
    G4cout << "Set global field value to " <<
      fEMfieldSetup->GetGlobalFieldValue() / gauss << " Gauss " << G4endl;
  }
  if( command == fLocalMagFieldCmd )
  {
    fEMfieldSetup->SetLocalFieldValue(fLocalMagFieldCmd->GetNew3VectorValue(newValue));
    fEMfieldSetup->UpdateField();
    // Check the value
    G4cout << "Set local field value to " <<
      fEMfieldSetup->GetLocalFieldValue() / gauss << " Gauss " << G4endl;
  }
  if( command == fMinStepCmd )
    fEMfieldSetup->SetMinStep(fMinStepCmd->GetNewDoubleValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
