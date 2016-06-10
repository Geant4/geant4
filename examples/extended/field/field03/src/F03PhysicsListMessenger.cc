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
/// \file field/field03/src/F03PhysicsListMessenger.cc
/// \brief Implementation of the F03PhysicsListMessenger class
//
//
// $Id: F03PhysicsListMessenger.cc 76602 2013-11-13 08:33:35Z gcosmo $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F03PhysicsListMessenger.hh"

#include "F03PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F03PhysicsListMessenger::F03PhysicsListMessenger(F03PhysicsList * List)
  : fF03List(List)
{
  fCutGCmd = new G4UIcmdWithADoubleAndUnit("/calor/cutG",this);
  fCutGCmd->SetGuidance("Set cut values by RANGE for Gamma.");
  fCutGCmd->SetParameterName("range",true);
  fCutGCmd->SetDefaultValue(1.);
  fCutGCmd->SetDefaultUnit("mm");
  fCutGCmd->AvailableForStates(G4State_Idle);

  fCutECmd = new G4UIcmdWithADoubleAndUnit("/calor/cutE",this);
  fCutECmd->SetGuidance("Set cut values by RANGE for e- e+.");
  fCutECmd->SetParameterName("range",true);
  fCutECmd->SetDefaultValue(1.);
  fCutECmd->SetDefaultUnit("mm");
  fCutECmd->AvailableForStates(G4State_Idle);

  fSetMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/step/setMaxStep",this);
  fSetMaxStepCmd->SetGuidance("Set max. step length in the detector");
  fSetMaxStepCmd->SetParameterName("mxStep",true);
  fSetMaxStepCmd->SetDefaultUnit("mm");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F03PhysicsListMessenger::~F03PhysicsListMessenger()
{
  delete fSetMaxStepCmd;
  delete fCutGCmd;
  delete fCutECmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F03PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if(command == fCutGCmd)
    {fF03List->SetGammaCut(fCutGCmd->GetNewDoubleValue(newValue));}
  if(command == fCutECmd)
    {fF03List->SetElectronCut(fCutECmd->GetNewDoubleValue(newValue));}
  if(command == fSetMaxStepCmd)
    {fF03List->SetMaxStep(fSetMaxStepCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
