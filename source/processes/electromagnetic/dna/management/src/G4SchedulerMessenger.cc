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
#include <G4Scheduler.hh>
#include <G4SchedulerMessenger.hh>

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

G4SchedulerMessenger::G4SchedulerMessenger(G4Scheduler * stepMgr) :
    fScheduler(stepMgr)
{
  fITDirectory = new G4UIdirectory("/scheduler/");
  fITDirectory->SetGuidance("Control commands for the time scheduler "
      "(dna chemistry applications).");

  // Set end time
  fEndTime = new G4UIcmdWithADoubleAndUnit("/scheduler/endTime", this);
  fEndTime->SetGuidance("Set time at which the simulation must stop.");
  fEndTime->AvailableForStates(G4State_PreInit, G4State_Idle);
  fEndTime->SetUnitCategory("Time");
  fEndTime->SetDefaultUnit("picosecond");
  fEndTime->SetDefaultValue(1);

  // Set time tolerance
  fTimeTolerance = new G4UIcmdWithADoubleAndUnit("/scheduler/timeTolerance",
                                                 this);
  fTimeTolerance->SetGuidance("This command aims at resolving issues related to"
      " floating points. If two time events are separated by less than the "
      "selected tolerance, they are assumed to belong to the same time step.");
  fTimeTolerance->AvailableForStates(G4State_PreInit, G4State_Idle);
  fTimeTolerance->SetUnitCategory("Time");
  fTimeTolerance->SetDefaultUnit("picosecond");
  fTimeTolerance->SetDefaultValue(1);

  // Initialize
  fInitCmd = new G4UIcmdWithoutParameter("/scheduler/initialize", this);
  fInitCmd->SetGuidance("Initialize G4Scheduler. This is done "
      "for standalone application only (no physics).");
  fInitCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Set Max Null time Step
  fMaxNULLTimeSteps = new G4UIcmdWithAnInteger("/scheduler/maxNullTimeSteps",
                                               this);
  fMaxNULLTimeSteps->SetGuidance("Set maximum allowed zero time steps. After this "
      "threshold, the simulation is stopped.");
  fMaxNULLTimeSteps->SetParameterName("numberOfNullTimeSteps", true);
  fMaxNULLTimeSteps->SetDefaultValue(10);
  fMaxNULLTimeSteps->SetRange("numberOfNullTimeSteps >=0 ");

  fMaxStepNumber = new G4UIcmdWithAnInteger("/scheduler/maxStepNumber", this);
  fMaxStepNumber->SetGuidance("Set the maximum number of time steps. After this "
      "threshold, the simulation is stopped.");
  fMaxStepNumber->SetParameterName("maximumNumberOfSteps", true);
  fMaxStepNumber->SetDefaultValue(-1);

  // Beam On
  fProcessCmd = new G4UIcmdWithoutParameter("/scheduler/process", this);
  fProcessCmd->SetGuidance("Process stacked tracks in G4Scheduler. This is done "
      "for standalone application only (no physics).");
  fProcessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Verbose
  fVerboseCmd = new G4UIcmdWithAnInteger("/scheduler/verbose", this);
  fVerboseCmd->SetGuidance("Set the Verbose level of G4Scheduler.");
  fVerboseCmd->SetGuidance(" 0 : Silent (default)");
  fVerboseCmd->SetGuidance(" 1 : Display reactions");
  fVerboseCmd->SetGuidance(" 2 ");
  fVerboseCmd->SetParameterName("level", true);
  fVerboseCmd->SetDefaultValue(1);
  //  fVerboseCmd->SetRange("level >=0 && level <=4");

  fWhyDoYouStop = new G4UIcmdWithoutParameter("/scheduler/whyDoYouStop",this);
  fWhyDoYouStop->SetGuidance("Will print information on why the scheduler is "
      "stopping the process");

  fUseDefaultTimeSteps = new G4UIcmdWithABool("/scheduler/useDefaultTimeSteps",
                                              this);
  fUseDefaultTimeSteps->SetGuidance("Let the G4 processes decided for the next "
      "time step interval. This command would be interesting if no reaction has "
      "been set and if one will want to track down Brownian objects. "
      "NB: This command gets in conflicts with the declaration of time steps.");
}

G4SchedulerMessenger::~G4SchedulerMessenger()
{
  delete fTimeTolerance;
  delete fITDirectory;
  delete fInitCmd;
  delete fEndTime;
  delete fMaxNULLTimeSteps;
  delete fMaxStepNumber;
  delete fProcessCmd;
  delete fVerboseCmd;
  delete fWhyDoYouStop;
  delete fUseDefaultTimeSteps;
}

void G4SchedulerMessenger::SetNewValue(G4UIcommand * command,
                                         G4String newValue)
{
  if (command == fProcessCmd)
  {
    fScheduler->Process();
  }
  else if (command == fEndTime)
  {
    fScheduler->SetEndTime(fEndTime->GetNewDoubleValue(newValue));
  }
  else if (command == fTimeTolerance)
  {
    fScheduler->SetTimeTolerance(
        fTimeTolerance->GetNewDoubleValue(newValue));
  }
  else if (command == fVerboseCmd)
  {
    fScheduler->SetVerbose(fVerboseCmd->GetNewIntValue(newValue));
  }
  else if (command == fInitCmd)
  {
    fScheduler->Initialize();
  }
  else if (command == fMaxNULLTimeSteps)
  {
    fScheduler->SetMaxZeroTimeAllowed(
        fMaxNULLTimeSteps->GetNewIntValue(newValue));
  }
  else if (command == fMaxStepNumber)
  {
    fScheduler->SetMaxNbSteps(fMaxStepNumber->GetNewIntValue(newValue));
  }
  else if (command == fWhyDoYouStop)
  {
    fScheduler->WhyDoYouStop();
  }
  else if (command == fUseDefaultTimeSteps)
  {
    fScheduler->UseDefaultTimeSteps(fUseDefaultTimeSteps->GetNewBoolValue(newValue));
  }
}

G4String G4SchedulerMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;

  if (command == fVerboseCmd)
  {
    cv = fVerboseCmd->ConvertToString(fScheduler->GetVerbose());
  }
  else if (command == fEndTime)
  {
    cv = fEndTime->ConvertToString(fScheduler->GetEndTime());
  }
  else if (command == fTimeTolerance)
  {
    cv = fTimeTolerance->ConvertToString(fScheduler->GetTimeTolerance());
  }
  else if (command == fInitCmd)
  {
    cv = fInitCmd->ConvertToString(fScheduler->IsInitialized());
  }
  else if (command == fMaxNULLTimeSteps)
  {
    cv = fMaxNULLTimeSteps->ConvertToString(
        fScheduler->GetMaxZeroTimeAllowed());
  }
  else if (command == fMaxStepNumber)
  {
    cv = fMaxStepNumber->ConvertToString(fScheduler->GetMaxNbSteps());
  }
  else if (command == fUseDefaultTimeSteps)
  {
    cv = fUseDefaultTimeSteps->ConvertToString(fScheduler->AreDefaultTimeStepsUsed());
  }

  return cv;
}

