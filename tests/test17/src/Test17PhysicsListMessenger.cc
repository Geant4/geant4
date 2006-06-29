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
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17PhysicsListMessenger.hh"

#include "Test17PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17PhysicsListMessenger::Test17PhysicsListMessenger(Test17PhysicsList * List)
:Test17List(List)
{
  cutGCmd = new G4UIcmdWithADoubleAndUnit("/test17/cutG",this);
  cutGCmd->SetGuidance("Set cut values by RANGE for Gamma.");
  cutGCmd->SetParameterName("cutG",false);
  cutGCmd->SetRange("cutG>0.");
  cutGCmd->SetUnitCategory("Length");
  cutGCmd->AvailableForStates(G4State_Idle);

  cutECmd = new G4UIcmdWithADoubleAndUnit("/test17/cutE",this);
  cutECmd->SetGuidance("Set cut values by RANGE for e- e+.");
  cutECmd->SetParameterName("cutE",false);
  cutECmd->SetRange("cutE>0.");
  cutECmd->SetUnitCategory("Length");
  cutECmd->AvailableForStates(G4State_Idle);

  eCmd = new G4UIcmdWithADoubleAndUnit("/test17/cutGammaEnergy",this);
  eCmd->SetGuidance("Set cut values by ENERGY for secondary gamma.");
  eCmd->SetParameterName("cutenergy",false);
  eCmd->SetRange("cutenergy>0.");
  eCmd->SetUnitCategory("Energy");
  eCmd->AvailableForStates(G4State_Idle);

  eaCmd = new G4UIcmdWithADoubleAndUnit("/test17/cutAugerEnergy",this);
  eaCmd->SetGuidance("Set cut values by ENERGY for Auger electrons.");
  eaCmd->SetParameterName("cutAenergy",false);
  eaCmd->SetRange("cutAenergy>0.");
  eaCmd->SetUnitCategory("Energy");
  eaCmd->AvailableForStates(G4State_Idle);

  setMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/step/setMaxStep",this);
  setMaxStepCmd->SetGuidance("Set max. step length in the detector");
  setMaxStepCmd->SetParameterName("mxStep",false);
  setMaxStepCmd->SetRange("mxStep>0.");
  setMaxStepCmd->SetUnitCategory("Length");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17PhysicsListMessenger::~Test17PhysicsListMessenger()
{
  delete cutGCmd;
  delete cutECmd;
  delete eCmd;
  delete eaCmd;
  delete setMaxStepCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == cutGCmd)
    { Test17List->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}
  if(command == cutECmd)
    { Test17List->SetElectronCut(cutECmd->GetNewDoubleValue(newValue));}
  if(command == eCmd)
    {Test17List->SetCutForSecondaryPhotons(eCmd->GetNewDoubleValue(newValue));}
  if(command == eaCmd)
    {Test17List->SetCutForAugerElectrons(eaCmd->GetNewDoubleValue(newValue));}
  if(command == setMaxStepCmd)
    { Test17List->SetMaxStep(setMaxStepCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

