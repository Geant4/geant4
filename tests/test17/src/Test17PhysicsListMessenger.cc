//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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

  cutPCmd = new G4UIcmdWithADoubleAndUnit("/test17/cutP",this);
  cutPCmd->SetGuidance("Set cut values by RANGE for proton and others.");
  cutPCmd->SetParameterName("cutP",false);
  cutPCmd->SetRange("cutP>0.");
  cutPCmd->SetUnitCategory("Length");    
  cutPCmd->AvailableForStates(G4State_Idle);

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

  rCmd = new G4UIcmdWithADoubleAndUnit("/test17/range",this);
  rCmd->SetGuidance("Display the RANGE of Electron for the current material.");
  rCmd->SetParameterName("range",false);
  rCmd->SetRange("range>0.");
  rCmd->SetUnitCategory("Length");     
  rCmd->AvailableForStates(G4State_Idle);

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
  delete cutPCmd;
  delete rCmd;
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
  if(command == cutPCmd)
    { Test17List->SetProtonCut(cutPCmd->GetNewDoubleValue(newValue));}
  if(command == eCmd)
    {Test17List->SetCutForSecondaryPhotons(eCmd->GetNewDoubleValue(newValue));}
  if(command == eaCmd)
    {Test17List->SetCutForAugerElectrons(eaCmd->GetNewDoubleValue(newValue));}
  if(command == rCmd)
    { Test17List->GetRange(rCmd->GetNewDoubleValue(newValue));}
  if(command == setMaxStepCmd)
    { Test17List->SetMaxStep(setMaxStepCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

