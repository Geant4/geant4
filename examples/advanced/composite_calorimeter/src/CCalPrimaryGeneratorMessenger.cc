///////////////////////////////////////////////////////////////////////////////
// File: CCalPrimaryGeneratorMessenger.cc
// Description: CCalPrimaryGeneratorMessenger adds new commands for
//              primary generator action
///////////////////////////////////////////////////////////////////////////////
#include "CCalPrimaryGeneratorMessenger.hh"
#include "CCalPrimaryGeneratorAction.hh"

#include "G4UImanager.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"

#include "PhysicalConstants.h"

CCalPrimaryGeneratorMessenger::CCalPrimaryGeneratorMessenger(CCalPrimaryGeneratorAction* myGun) : myAction(myGun) {

  verboseCmd = new G4UIcmdWithAnInteger("/CCal/generator/verbose",this);
  verboseCmd->SetGuidance("set Verbosity level ");
  verboseCmd->SetParameterName("value",true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->AvailableForStates(PreInit,Idle);

  rndmCmd = new G4UIcmdWithAString("/CCal/generator/random",this);
  rndmCmd->SetGuidance("Choose randomly energy and direction of the incident particle.");
  rndmCmd->SetGuidance("  Choice : on,off(default)");
  rndmCmd->SetParameterName("choice",true);
  rndmCmd->SetDefaultValue("off");
  rndmCmd->SetCandidates("on off ON OFF");
  rndmCmd->AvailableForStates(PreInit,Idle);

  scanCmd = new G4UIcmdWithAString("/CCal/generator/scan",this);
  scanCmd->SetGuidance("Scan eta and phi ranges with single incident particle");
  scanCmd->SetGuidance("  Choice : on,off(default)");
  scanCmd->SetGuidance("  Ranges : etamin/max, phimin/max are set by other commands ");
  scanCmd->SetParameterName("choice",true);
  scanCmd->SetDefaultValue("off");
  scanCmd->SetCandidates("on off ON OFF");
  scanCmd->AvailableForStates(PreInit,Idle);

  minEnergyCmd = new G4UIcmdWithADoubleAndUnit("/CCal/generator/minEnergy",this);
  minEnergyCmd->SetGuidance("Set minimum Energy for the incident particle.");
  minEnergyCmd->SetParameterName("value",true);
  minEnergyCmd->SetDefaultValue(1.);
  minEnergyCmd->SetDefaultUnit("GeV");
  minEnergyCmd->AvailableForStates(PreInit,Idle);
  
  maxEnergyCmd = new G4UIcmdWithADoubleAndUnit("/CCal/generator/maxEnergy",this);
  maxEnergyCmd->SetGuidance("Set maximum Energy for the incident particle.");
  maxEnergyCmd->SetParameterName("value",true);
  maxEnergyCmd->SetDefaultValue(1.);
  maxEnergyCmd->SetDefaultUnit("TeV");
  maxEnergyCmd->AvailableForStates(PreInit,Idle);

  minPhiCmd = new G4UIcmdWithADoubleAndUnit("/CCal/generator/minPhi",this);
  minPhiCmd->SetGuidance("Set minimum Phi angle for the incident particle direction");
  minPhiCmd->SetGuidance("  Choice : from 0 to 2*pi ");
  minPhiCmd->SetParameterName("value",true);
  minPhiCmd->SetDefaultValue(0);
  minPhiCmd->SetDefaultUnit("radian");
  minPhiCmd->AvailableForStates(PreInit,Idle);

  maxPhiCmd = new G4UIcmdWithADoubleAndUnit("/CCal/generator/maxPhi",this);
  maxPhiCmd->SetGuidance("Set maximum Phi angle for the incident particle direction");
  maxPhiCmd->SetGuidance("  Choice : from 0 to 2*pi ");
  maxPhiCmd->SetParameterName("value",true);
  maxPhiCmd->SetDefaultValue(2.*pi);
  maxPhiCmd->SetDefaultUnit("radian");
  maxPhiCmd->AvailableForStates(PreInit,Idle);

  stepsPhiCmd = new G4UIcmdWithAnInteger("/CCal/generator/stepsPhi",this);
  stepsPhiCmd->SetGuidance("number of steps along Phi for scan ");
  stepsPhiCmd->SetParameterName("value",true);
  stepsPhiCmd->SetDefaultValue(1);
  stepsPhiCmd->AvailableForStates(PreInit,Idle);

  minEtaCmd = new G4UIcmdWithADouble("/CCal/generator/minEta",this);
  minEtaCmd->SetGuidance("Set minimum Eta angle for the incident particle direction");
  minEtaCmd->SetGuidance("  Choice : from 0 to infinity");
  minEtaCmd->SetParameterName("value",true);
  minEtaCmd->SetDefaultValue(0);
  minEtaCmd->AvailableForStates(PreInit,Idle);

  maxEtaCmd = new G4UIcmdWithADouble("/CCal/generator/maxEta",this);
  maxEtaCmd->SetGuidance("Set maximum Eta angle for the incident particle direction");
  maxEtaCmd->SetGuidance("  Choice : from 0 to infinity");
  maxEtaCmd->SetParameterName("value",true);
  maxEtaCmd->SetDefaultValue(3.5);
  maxEtaCmd->AvailableForStates(PreInit,Idle);

  stepsEtaCmd = new G4UIcmdWithAnInteger("/CCal/generator/stepsEta",this);
  stepsEtaCmd->SetGuidance("number of steps along Eta for scan ");
  stepsEtaCmd->SetParameterName("value",true);
  stepsEtaCmd->SetDefaultValue(1);
  stepsEtaCmd->AvailableForStates(PreInit,Idle);

  runNoCmd = new G4UIcmdWithAnInteger("/CCal/generator/runNo",this);
  runNoCmd->SetGuidance("set the run number ");
  runNoCmd->SetParameterName("value",true);
  runNoCmd->SetDefaultValue(0);
  runNoCmd->AvailableForStates(PreInit,Idle);

}

CCalPrimaryGeneratorMessenger::~CCalPrimaryGeneratorMessenger() {

  if (verboseCmd)
    delete verboseCmd;
  if (scanCmd)
    delete rndmCmd;
  if (scanCmd)
    delete scanCmd;
  if (minEnergyCmd)
    delete minEnergyCmd;
  if (maxEnergyCmd)
    delete maxEnergyCmd;
  if (minPhiCmd)
    delete minPhiCmd;
  if (maxPhiCmd)
    delete maxPhiCmd;
  if (stepsPhiCmd)
    delete stepsPhiCmd;
  if (minEtaCmd)
    delete minEtaCmd;
  if (maxEtaCmd)
    delete maxEtaCmd;
  if (stepsEtaCmd)
    delete stepsEtaCmd;
  if (runNoCmd)
    delete runNoCmd;

}

void CCalPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,
						G4String newValues)    { 
  if (command == verboseCmd)
    myAction->SetVerboseLevel(verboseCmd->GetNewIntValue(newValues));
  else if (command == rndmCmd)
    myAction->SetRandom(newValues);
  else if (command == scanCmd)
    myAction->SetScan(newValues);
  else if (command == minEnergyCmd)
    myAction->SetMinimumEnergy(minEnergyCmd->GetNewDoubleValue(newValues));
  else if (command == maxEnergyCmd)
    myAction->SetMaximumEnergy(maxEnergyCmd->GetNewDoubleValue(newValues));
  else if (command == minPhiCmd)
    myAction->SetMinimumPhi(minPhiCmd->GetNewDoubleValue(newValues));
  else if (command == maxPhiCmd)
    myAction->SetMaximumPhi(maxPhiCmd->GetNewDoubleValue(newValues));
  else if (command == stepsPhiCmd)
    myAction->SetStepsPhi(stepsPhiCmd->GetNewIntValue(newValues));
  else if (command == minEtaCmd)
    myAction->SetMinimumEta(minEtaCmd->GetNewDoubleValue(newValues));
  else if (command == maxEtaCmd)
    myAction->SetMaximumEta(maxEtaCmd->GetNewDoubleValue(newValues));
  else if (command == stepsEtaCmd)
    myAction->SetStepsEta(stepsEtaCmd->GetNewIntValue(newValues));
  else if (command == runNoCmd)
    myAction->SetRunNo(runNoCmd->GetNewIntValue(newValues));
 
}
