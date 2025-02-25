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

/// \file G4FieldParametersMessenger.cc
/// \brief Implementation of the G4FieldParametersMessenger class
///
/// This code was initially developed in Geant4 VMC package
/// (https://github.com/vmc-project)
/// and adapted to Geant4.
///
/// \author I. Hrivnacova; IJCLab, Orsay

#include "G4FieldParametersMessenger.hh"
#include "G4FieldParameters.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"

//_____________________________________________________________________________
G4FieldParametersMessenger::G4FieldParametersMessenger(
  G4FieldParameters* fieldParameters)
  : fFieldParameters(fieldParameters)
{
  // Standard constructor

  G4String directoryName = "/field/";
  if (fieldParameters->GetVolumeName() != "") {
    directoryName.append(fieldParameters->GetVolumeName());
    directoryName.append("/");
    fDirectory = new G4UIdirectory(directoryName);
    fDirectory->SetGuidance("Magnetic field control commands.");
  }

  G4String commandName = directoryName;
  commandName.append("fieldType");
  fFieldTypeCmd = new G4UIcmdWithAString(commandName, this);
  G4String guidance = "Select type of the field";
  fFieldTypeCmd->SetGuidance(guidance);
  fFieldTypeCmd->SetParameterName("FieldType", false);
  G4String candidates;
  for (G4int i = kMagnetic; i <= kGravity; i++) {
    G4FieldType ft = (G4FieldType)i;
    candidates += G4FieldParameters::FieldTypeName(ft);
    candidates += " ";
  }
  fFieldTypeCmd->SetCandidates(candidates);
  fFieldTypeCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);

  commandName = directoryName;
  commandName.append("equationType");
  fEquationTypeCmd = new G4UIcmdWithAString(commandName, this);
  guidance = "Select type of the equation of motion of a particle in a field";
  fEquationTypeCmd->SetGuidance(guidance);
  fEquationTypeCmd->SetParameterName("EquationType", false);
  candidates = "";
  for (G4int i = kEqMagnetic; i <= kEqEMfieldWithEDM; i++) {
    G4EquationType et = (G4EquationType)i;
    candidates += G4FieldParameters::EquationTypeName(et);
    candidates += " ";
  }
  fEquationTypeCmd->SetCandidates(candidates);
  fEquationTypeCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);

  commandName = directoryName;
  commandName.append("stepperType");
  fStepperTypeCmd = new G4UIcmdWithAString(commandName, this);
  guidance =
    "Select type of the the integrator of particle's equation of motion in a "
    "field";
  fStepperTypeCmd->SetGuidance(guidance);
  fStepperTypeCmd->SetParameterName("StepperType", false);
  candidates = "";
  for (G4int i = kCashKarpRKF45; i <= kRK547FEq3; i++) {    
    G4StepperType st = (G4StepperType)i;
    if (st == kUserStepper) continue;
    candidates += G4FieldParameters::StepperTypeName(st);
    candidates += " ";
  }
  fStepperTypeCmd->SetCandidates(candidates);
  fStepperTypeCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);

  commandName = directoryName;
  commandName.append("setMinimumStep");
  fSetMinimumStepCmd = new G4UIcmdWithADoubleAndUnit(commandName, this);
  fSetMinimumStepCmd->SetGuidance("Set minimum step in G4ChordFinder");
  fSetMinimumStepCmd->SetParameterName("StepMinimum", false);
  fSetMinimumStepCmd->SetDefaultUnit("mm");
  fSetMinimumStepCmd->SetUnitCategory("Length");
  fSetMinimumStepCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);

  commandName = directoryName;
  commandName.append("setDeltaChord");
  fSetDeltaChordCmd = new G4UIcmdWithADoubleAndUnit(commandName, this);
  fSetDeltaChordCmd->SetGuidance("Set delta chord in G4ChordFinder");
  fSetDeltaChordCmd->SetParameterName("DeltaChord", false);
  fSetDeltaChordCmd->SetDefaultUnit("mm");
  fSetDeltaChordCmd->SetUnitCategory("Length");
  fSetDeltaChordCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);

  commandName = directoryName;
  commandName.append("setDeltaOneStep");
  fSetDeltaOneStepCmd = new G4UIcmdWithADoubleAndUnit(commandName, this);
  fSetDeltaOneStepCmd->SetGuidance(
    "Set delta one step in global field manager");
  fSetDeltaOneStepCmd->SetParameterName("DeltaOneStep", false);
  fSetDeltaOneStepCmd->SetDefaultUnit("mm");
  fSetDeltaOneStepCmd->SetUnitCategory("Length");
  fSetDeltaOneStepCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);

  commandName = directoryName;
  commandName.append("setDeltaIntersection");
  fSetDeltaIntersectionCmd = new G4UIcmdWithADoubleAndUnit(commandName, this);
  fSetDeltaIntersectionCmd->SetGuidance(
    "Set delta intersection in global field manager");
  fSetDeltaIntersectionCmd->SetParameterName("DeltaIntersection", false);
  fSetDeltaIntersectionCmd->SetDefaultUnit("mm");
  fSetDeltaIntersectionCmd->SetUnitCategory("Length");
  fSetDeltaIntersectionCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);

  commandName = directoryName;
  commandName.append("setMinimumEpsilonStep");
  fSetMinimumEpsilonStepCmd = new G4UIcmdWithADouble(commandName, this);
  fSetMinimumEpsilonStepCmd->SetGuidance(
    "Set minimum epsilon step in global field manager");
  fSetMinimumEpsilonStepCmd->SetParameterName("MinimumEpsilonStep", false);
  fSetMinimumEpsilonStepCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);

  commandName = directoryName;
  commandName.append("setMaximumEpsilonStep");
  fSetMaximumEpsilonStepCmd = new G4UIcmdWithADouble(commandName, this);
  fSetMaximumEpsilonStepCmd->SetGuidance(
    "Set maximum epsilon step in global field manager");
  fSetMaximumEpsilonStepCmd->SetParameterName("MaximumEpsilonStep", false);
  fSetMaximumEpsilonStepCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);

  commandName = directoryName;
  commandName.append("setConstDistance");
  fSetConstDistanceCmd = new G4UIcmdWithADoubleAndUnit(commandName, this);
  fSetConstDistanceCmd->SetGuidance(
    "Set the distance within which the field is considered constant.");
  fSetConstDistanceCmd->SetGuidance(
    "Non-zero value will trigger creating a cached magnetic field.");
  fSetConstDistanceCmd->SetParameterName("ConstDistance", false);
  fSetConstDistanceCmd->SetDefaultUnit("mm");
  fSetConstDistanceCmd->SetUnitCategory("Length");
  fSetConstDistanceCmd->SetRange("ConstDistance >= 0");
  fSetConstDistanceCmd->AvailableForStates(G4State_PreInit);

  commandName = std::move(directoryName);
  commandName.append("printParameters");
  fPrintParametersCmd = new G4UIcmdWithoutParameter(commandName, this);
  fPrintParametersCmd->SetGuidance("Prints all accuracy parameters.");
  fPrintParametersCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);
}

//_____________________________________________________________________________
G4FieldParametersMessenger::~G4FieldParametersMessenger()
{
  // Destructor

  delete fDirectory;
  delete fFieldTypeCmd;
  delete fEquationTypeCmd;
  delete fStepperTypeCmd;
  delete fSetMinimumStepCmd;
  delete fSetDeltaChordCmd;
  delete fSetDeltaOneStepCmd;
  delete fSetDeltaIntersectionCmd;
  delete fSetMinimumEpsilonStepCmd;
  delete fSetMaximumEpsilonStepCmd;
  delete fSetConstDistanceCmd;
}

//
// public methods
//

//_____________________________________________________________________________
void G4FieldParametersMessenger::SetNewValue(
  G4UIcommand* command, G4String newValues)
{
  // Apply command to the associated object.

  if (command == fFieldTypeCmd) {
    for (G4int i = kMagnetic; i <= kGravity; i++) {
      G4FieldType ft = (G4FieldType)i;
      if (newValues == G4FieldParameters::FieldTypeName(ft)) {
        fFieldParameters->SetFieldType(ft);
        break;
      }
    }
    return;
  }
  
  if (command == fEquationTypeCmd) {
    for (G4int i = kEqMagnetic; i <= kEqEMfieldWithEDM; i++) {
      G4EquationType et = (G4EquationType)i;
      if (newValues == G4FieldParameters::EquationTypeName(et)) {
        fFieldParameters->SetEquationType(et);
        break;
      }
    }
    return;
  }

  if (command == fStepperTypeCmd) {
    for (G4int i = kCashKarpRKF45; i <= kRK547FEq3; i++) {
      G4StepperType st = (G4StepperType)i;
      if (newValues == G4FieldParameters::StepperTypeName(st)) {
        fFieldParameters->SetStepperType(st);
        break;
      }
    }
    return;
  }

  if (command == fSetMinimumStepCmd) {
    fFieldParameters->SetMinimumStep(
      fSetMinimumStepCmd->GetNewDoubleValue(newValues));
    return;
  }

  if (command == fSetDeltaChordCmd) {
    fFieldParameters->SetDeltaChord(
      fSetDeltaChordCmd->GetNewDoubleValue(newValues));
    return;
  }

  if (command == fSetDeltaOneStepCmd) {
    fFieldParameters->SetDeltaOneStep(
      fSetDeltaOneStepCmd->GetNewDoubleValue(newValues));
    return;
  }

  if (command == fSetDeltaIntersectionCmd) {
    fFieldParameters->SetDeltaIntersection(
      fSetDeltaIntersectionCmd->GetNewDoubleValue(newValues));
    return;
  }

  if (command == fSetMinimumEpsilonStepCmd) {
    fFieldParameters->SetMinimumEpsilonStep(
      fSetMinimumEpsilonStepCmd->GetNewDoubleValue(newValues));
    return;
  }

  if (command == fSetMaximumEpsilonStepCmd) {
    fFieldParameters->SetMaximumEpsilonStep(
      fSetMaximumEpsilonStepCmd->GetNewDoubleValue(newValues));
    return;
  }

  if (command == fSetConstDistanceCmd) {
    fFieldParameters->SetConstDistance(
      fSetConstDistanceCmd->GetNewDoubleValue(newValues));
    return;
  }

  if (command == fPrintParametersCmd) {
    fFieldParameters->PrintParameters();
    return;
  }
}
