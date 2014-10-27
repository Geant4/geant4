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
// $Id: G4ITSteppingMessenger.cc 60427 2012-07-11 16:34:35Z matkara $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// History:
// -----------
//
// -------------------------------------------------------------------
#include <G4ITScheduler.hh>
#include <G4ITSchedulerMessenger.hh>

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

G4ITSchedulerMessenger::G4ITSchedulerMessenger(G4ITScheduler * stepMgr) :
    fITStepManager(stepMgr)
{
  fITDirectory = new G4UIdirectory("/IT/");
  fITDirectory->SetGuidance("IT control commands.");

  // Set end time
  fEndTime = new G4UIcmdWithADoubleAndUnit("/IT/endTime", this);
  fEndTime->SetGuidance("Set end time");
  fEndTime->AvailableForStates(G4State_PreInit, G4State_Idle);
  fEndTime->SetUnitCategory("Time");
  fEndTime->SetDefaultUnit("picosecond");
  fEndTime->SetDefaultValue(1);

  // Set time tolerance
  fTimeTolerance = new G4UIcmdWithADoubleAndUnit("/IT/timeTolerance", this);
  fTimeTolerance->SetGuidance("Set time tolerance");
  fTimeTolerance->AvailableForStates(G4State_PreInit, G4State_Idle);
  fTimeTolerance->SetUnitCategory("Time");
  fTimeTolerance->SetDefaultUnit("picosecond");
  fTimeTolerance->SetDefaultValue(1);

  // Initialize
  fInitCmd = new G4UIcmdWithoutParameter("/IT/initialize", this);
  fInitCmd->SetGuidance("Initialize G4ITStepManager.");
  fInitCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Set Max Null time Step
  fMaxNULLTimeSteps = new G4UIcmdWithAnInteger("/IT/maxNullTimeSteps", this);
  fMaxNULLTimeSteps->SetGuidance("Set Max Null time Step");
  fMaxNULLTimeSteps->SetParameterName("numberOfNullTimeSteps", true);
  fMaxNULLTimeSteps->SetDefaultValue(10);
  fMaxNULLTimeSteps->SetRange("numberOfNullTimeSteps >=0 ");

  fMaxStepNumber = new G4UIcmdWithAnInteger("/IT/maxStepsNumber", this);
  fMaxStepNumber->SetParameterName("maximumNumberOfSteps", true);
  fMaxStepNumber->SetDefaultValue(-1);

  // Beam On
  fProcessCmd = new G4UIcmdWithoutParameter("/IT/process", this);
  fProcessCmd->SetGuidance("Process track stacked in G4ITStepManager");
  fProcessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Verbose
  fVerboseCmd = new G4UIcmdWithAnInteger("/IT/verbose", this);
  fVerboseCmd->SetGuidance("Set the Verbose level of G4ITStepManager.");
  fVerboseCmd->SetGuidance(" 0 : Silent (default)");
  fVerboseCmd->SetGuidance(" 1 : Display reactions");
  fVerboseCmd->SetGuidance(" 2 ");
  fVerboseCmd->SetParameterName("level", true);
  fVerboseCmd->SetDefaultValue(1);
  //  fVerboseCmd->SetRange("level >=0 && level <=4");
}

G4ITSchedulerMessenger::~G4ITSchedulerMessenger()
{
  delete fTimeTolerance;
  delete fITDirectory;
  delete fInitCmd;
  delete fEndTime;
  delete fMaxNULLTimeSteps;
  delete fMaxStepNumber;
  delete fProcessCmd;
  delete fVerboseCmd;
}

void G4ITSchedulerMessenger::SetNewValue(G4UIcommand * command,
                                         G4String newValue)
{
  if (command == fProcessCmd)
  {
    fITStepManager->Process();
  }
  else if (command == fEndTime)
  {
    fITStepManager->SetEndTime(fEndTime->GetNewDoubleValue(newValue));
  }
  else if (command == fTimeTolerance)
  {
    fITStepManager->SetTimeTolerance(
        fTimeTolerance->GetNewDoubleValue(newValue));
  }
  else if (command == fVerboseCmd)
  {
    fITStepManager->SetVerbose(fVerboseCmd->GetNewIntValue(newValue));
  }
  else if (command == fInitCmd)
  {
    fITStepManager->Initialize();
  }
  else if (command == fMaxNULLTimeSteps)
  {
    fITStepManager->SetMaxZeroTimeAllowed(
        fMaxNULLTimeSteps->GetNewIntValue(newValue));
  }
  else if (command == fMaxStepNumber)
  {
    fITStepManager->SetMaxNbSteps(fMaxStepNumber->GetNewIntValue(newValue));
  }
}

G4String G4ITSchedulerMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;

  if (command == fVerboseCmd)
  {
    cv = fVerboseCmd->ConvertToString(fITStepManager->GetVerbose());
  }
  else if (command == fEndTime)
  {
    cv = fEndTime->ConvertToString(fITStepManager->GetEndTime());
  }
  else if (command == fTimeTolerance)
  {
    cv = fTimeTolerance->ConvertToString(fITStepManager->GetTimeTolerance());
  }
  else if (command == fInitCmd)
  {
    cv = fInitCmd->ConvertToString(fITStepManager->IsInitialized());
  }
  else if (command == fMaxNULLTimeSteps)
  {
    cv = fMaxNULLTimeSteps->ConvertToString(
        fITStepManager->GetMaxZeroTimeAllowed());
  }
  else if (command == fMaxStepNumber)
  {
    cv = fMaxStepNumber->ConvertToString(fITStepManager->GetMaxNbSteps());
  }

  return cv;
}

