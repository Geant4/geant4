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
/// \file field/field04/src/F04FieldMessenger.cc
/// \brief Implementation of the F04FieldMessenger class
//

#include "F04FieldMessenger.hh"

#include "F04GlobalField.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "F04DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04FieldMessenger::F04FieldMessenger(F04GlobalField* pEMfield,
                                     F04DetectorConstruction* detector)
  : fGlobalField(pEMfield)
{
  fDetector = detector;

  fDetDir = new G4UIdirectory("/field/");
  fDetDir->SetGuidance(" Field tracking control ");

  fCaptureB1Cmd = new G4UIcmdWithADoubleAndUnit("/field/SetCaptureB1",this);
  fCaptureB1Cmd->SetGuidance("Set B1 of the Capture Magnet");
  fCaptureB1Cmd->SetParameterName("CSizeB1",false,false);
  fCaptureB1Cmd->SetDefaultUnit("tesla");
  fCaptureB1Cmd->SetRange("CSizeB1>0.");
  fCaptureB1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fCaptureB2Cmd = new G4UIcmdWithADoubleAndUnit("/field/SetCaptureB2",this);
  fCaptureB2Cmd->SetGuidance("Set B2 of the Capture Magnet");
  fCaptureB2Cmd->SetParameterName("CSizeB2",false,false);
  fCaptureB2Cmd->SetDefaultUnit("tesla");
  fCaptureB2Cmd->SetRange("CSizeB2>0.");
  fCaptureB2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fTransferBCmd = new G4UIcmdWithADoubleAndUnit("/field/SetTransferB",this);
  fTransferBCmd->SetGuidance("Set B of the Transfer Magnet");
  fTransferBCmd->SetParameterName("TSizeB",false,false);
  fTransferBCmd->SetDefaultUnit("tesla");
  fTransferBCmd->SetRange("TSizeB>0.");
  fTransferBCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fStepperCMD = new G4UIcmdWithAnInteger("/field/setStepperType",this);
  fStepperCMD->SetGuidance("Select stepper type for field");
  fStepperCMD->SetParameterName("choice",true);
  fStepperCMD->SetDefaultValue(4);
  fStepperCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  fMinStepCMD = new G4UIcmdWithADoubleAndUnit("/field/setMinStep",this);
  fMinStepCMD->SetGuidance("Define minimal step");
  fMinStepCMD->SetParameterName("min step",false,false);
  fMinStepCMD->SetDefaultUnit("mm");
  fMinStepCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDeltaChordCMD = new G4UIcmdWithADoubleAndUnit("/field/setDeltaChord",this);
  fDeltaChordCMD->SetGuidance("Define delta chord");
  fDeltaChordCMD->SetParameterName("delta chord",false,false);
  fDeltaChordCMD->SetDefaultUnit("mm");
  fDeltaChordCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDeltaOneStepCMD =
                 new G4UIcmdWithADoubleAndUnit("/field/setDeltaOneStep",this);
  fDeltaOneStepCMD->SetGuidance("Define delta one step");
  fDeltaOneStepCMD->SetParameterName("delta one step",false,false);
  fDeltaOneStepCMD->SetDefaultUnit("mm");
  fDeltaOneStepCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDeltaIntersectionCMD =
            new G4UIcmdWithADoubleAndUnit("/field/setDeltaIntersection",this);
  fDeltaIntersectionCMD->SetGuidance("Define delta intersection");
  fDeltaIntersectionCMD->SetParameterName("delta intersection",false,false);
  fDeltaIntersectionCMD->SetDefaultUnit("mm");
  fDeltaIntersectionCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

  fEpsMinCMD = new G4UIcmdWithADoubleAndUnit("/field/setEpsMin",this);
  fEpsMinCMD->SetGuidance("Define eps min");
  fEpsMinCMD->SetParameterName("eps min",false,false);
  fEpsMinCMD->SetDefaultUnit("mm");
  fEpsMinCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

  fEpsMaxCMD = new G4UIcmdWithADoubleAndUnit("/field/setEpsMax",this);
  fEpsMaxCMD->SetGuidance("Define eps max");
  fEpsMaxCMD->SetParameterName("eps max",false,false);
  fEpsMaxCMD->SetDefaultUnit("mm");
  fEpsMaxCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04FieldMessenger::~F04FieldMessenger()
{
  delete fDetDir;

  delete fCaptureB1Cmd;
  delete fCaptureB2Cmd;
  delete fTransferBCmd;

  delete fStepperCMD;
  delete fMinStepCMD;
  delete fDeltaChordCMD;
  delete fDeltaOneStepCMD;
  delete fDeltaIntersectionCMD;
  delete fEpsMinCMD;
  delete fEpsMaxCMD;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04FieldMessenger::SetNewValue( G4UIcommand* command, G4String newValue)
{

  if( command == fCaptureB1Cmd )
    fDetector->SetCaptureMgntB1(fCaptureB1Cmd->GetNewDoubleValue(newValue));

  if( command == fCaptureB2Cmd )
    fDetector->SetCaptureMgntB2(fCaptureB2Cmd->GetNewDoubleValue(newValue));

  if( command == fTransferBCmd )
    fDetector->SetTransferMgntB(fTransferBCmd->GetNewDoubleValue(newValue));

  if( command == fStepperCMD )
  {
    fGlobalField->SetStepperType(fStepperCMD->GetNewIntValue(newValue));
  }
  if( command == fMinStepCMD )
  {
    fGlobalField->SetMinStep(fMinStepCMD->GetNewDoubleValue(newValue));
  }
  if( command == fDeltaChordCMD )
  {
    fGlobalField->SetDeltaChord(fDeltaChordCMD->GetNewDoubleValue(newValue));
  }
  if( command == fDeltaOneStepCMD )
  {
    fGlobalField->
                SetDeltaOneStep(fDeltaOneStepCMD->GetNewDoubleValue(newValue));
  }
  if( command == fDeltaIntersectionCMD )
  {
    fGlobalField->
      SetDeltaIntersection(fDeltaIntersectionCMD->GetNewDoubleValue(newValue));
  }
  if( command == fEpsMinCMD )
  {
    fGlobalField->SetEpsMin(fEpsMinCMD->GetNewDoubleValue(newValue));
  }
  if( command == fEpsMaxCMD )
  {
    fGlobalField->SetEpsMax(fEpsMaxCMD->GetNewDoubleValue(newValue));
  }
}
