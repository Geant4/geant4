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
// G4QSSMessenger implementation
//
// Author: Leandro Gomez Vidal - October 2021.
// --------------------------------------------------------------------

#include "G4QSSMessenger.hh"

G4QSSMessenger::G4QSSMessenger()
{
  commandsShouldBeInMaster = true;

  qssCmdDir = new G4UIdirectory("/QSS/", false);
  qssCmdDir->SetGuidance("G4QSStepper configuration.");

  dQMinCmd = new G4UIcmdWithADoubleAndUnit("/QSS/dQMin",this);
  dQMinCmd->SetDefaultUnit("mm");
  dQMinCmd->SetParameterName("dQMinCmd",false);
  dQMinCmd->SetUnitCategory("Length");

  dQRelCmd = new G4UIcmdWithADouble("/QSS/dQRel",this);
  dQRelCmd->SetGuidance("Default is 1e-5");
  dQRelCmd->SetParameterName("dQRelCmd",false);

  trialProposedStepModifierCmd = new G4UIcmdWithADouble("/QSS/trialProposedStepModifier",this);
  trialProposedStepModifierCmd->SetGuidance("Default is 1");
  trialProposedStepModifierCmd->SetParameterName("trialProposedStepModifier", false);

  stepperSelectorCmd = new G4UIcmdWithAString("/QSS/selectStepper",this);
  stepperSelectorCmd->SetGuidance("Select stepper.");
  stepperSelectorCmd->SetParameterName("choice", false);
  stepperSelectorCmd->SetCandidates("TemplatedDoPri OldRK45 G4QSS2");

}

G4QSSMessenger::~G4QSSMessenger()
{
  delete qssCmdDir;
  delete dQMinCmd;
  delete dQRelCmd;
  delete stepperSelectorCmd;
  delete trialProposedStepModifierCmd;
  //qssStats.print();
}

G4QSSMessenger* G4QSSMessenger::instance()
{
  static G4QSSMessenger theSingleMessengerInstance;
  return &theSingleMessengerInstance;
}

void G4QSSMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
  if ( command == dQMinCmd ) {
    dQMin = dQMinCmd->GetNewDoubleValue(newValue);
  }

  if ( command == dQRelCmd ) {
    dQRel = dQRelCmd->GetNewDoubleValue(newValue);
  }

  if ( command == trialProposedStepModifierCmd ) {
    trialProposedStepModifier = trialProposedStepModifierCmd->GetNewDoubleValue(newValue);
  }

  if(command == stepperSelectorCmd){
    this->selectStepper(newValue);
  }

}

void G4QSSMessenger::selectStepper(const std::string &newValue)
{
  const std::map<std::string, StepperSelection> stepperMapping =
    {{"TemplatedDoPri", TemplatedDoPri}, {"OldRK45", OldRK45}, {"G4QSS2", G4QSS2}};
  _selectedStepper = stepperMapping.at(newValue);
  G4cout << "G4QSSMessenger: Selecting stepper " << newValue << G4endl;
}

G4QSSMessenger::StepperSelection G4QSSMessenger::selectedStepper()
{
  return _selectedStepper;
}
