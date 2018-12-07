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
/// \file electromagnetic/TestEm12/src/StepMaxMessenger.cc
/// \brief Implementation of the StepMaxMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StepMaxMessenger.hh"

#include "StepMax.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepMaxMessenger::StepMaxMessenger(StepMax* stepM)
:G4UImessenger(),fStepMax(stepM),
 fStepMax1Cmd(0),    
 fStepMax2Cmd(0)
{ 
  fStepMax1Cmd = new G4UIcmdWithADoubleAndUnit("/testem/stepMax",this);
  fStepMax1Cmd->SetGuidance("Set max allowed step length");
  fStepMax1Cmd->SetParameterName("mxStep1",false);
  fStepMax1Cmd->SetRange("mxStep1>0.");
  fStepMax1Cmd->SetUnitCategory("Length");
   
  fStepMax2Cmd = new G4UIcmdWithABool("/testem/applyAutomaticStepMax",this);
  fStepMax2Cmd->SetGuidance("apply StepMax computed from histograms");
  fStepMax2Cmd->SetParameterName("mxStep2",true);
  fStepMax2Cmd->SetDefaultValue(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepMaxMessenger::~StepMaxMessenger()
{
  delete fStepMax1Cmd;
  delete fStepMax2Cmd;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StepMaxMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fStepMax1Cmd)
    { fStepMax->SetMaxStep1(fStepMax1Cmd->GetNewDoubleValue(newValue));}
    
  if (command == fStepMax2Cmd)
    { fStepMax->ApplyMaxStep2(fStepMax2Cmd->GetNewBoolValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
